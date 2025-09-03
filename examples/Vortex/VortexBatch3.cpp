#include "base/Type.h"
#include "Mesh/Mesh.h"
#include "Matrix/Matrix.h"
#include "DG/DG_Basis/DG_Basis.h"
// #include "DG/DG_Schemes/ExplicitConvection.h"
// #include "DG/DG_Schemes/ExplicitDiffusion.h"
#include "DG/DG_Schemes/ImplicitConvection.h"
// #include "DG/DG_Schemes/ImplicitDiffusion.h"
#include "EigenSolver/EigenSparseSolver.h"
// #include "unsupported/Eigen/IterativeSolvers"

string sol_filename(uInt Order, uInt Time, uInt Nmesh){
    return "./Order_" + std::to_string(Order) 
            + "/solution/T_" + std::to_string(Time) 
            + "_N_" + std::to_string(Nmesh) + ".txt";
}

#include "problem.h"
#include "submat.h"
ComputingMesh create_mesh(uInt N){
    GeneralMesh mesh = OrthHexMesh({0,0,0},{10,10,10.0/N},{N,N,1});
    mesh.split_hex5_scan();                                   
    mesh.rebuild_cell_topology();                             
    mesh.validate_mesh();                                     
    ComputingMesh cmesh(mesh);                                
    cmesh.m_boundaryTypes.resize(cmesh.m_faces.size());                   
    for(uInt faceId=0;faceId<cmesh.m_faces.size();faceId++){           
        if(cmesh.m_faces[faceId].m_neighbor_cells[1]==uInt(-1)){ 
            const auto& face = cmesh.m_faces[faceId];           
            if(std::abs(face.m_normal[2])>0.5 ) 
                cmesh.m_boundaryTypes[faceId] = BoundaryType::Pseudo3DZ;
            // else if(std::abs(face.m_normal[1])>0.5)          
            //     cmesh.m_boundaryTypes[faceId] = BoundaryType::Pseudo3DY;
            else
                cmesh.m_boundaryTypes[faceId] = BoundaryType::Dirichlet;
        }
    }
    return cmesh;
}



template<uInt Order>
void Run(uInt N){
    auto chrono_start = std::chrono::steady_clock::now();
    auto chrone_clock = [&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-chrono_start).count();};
    auto logging = [&](std::string ss){debug("Time  " + std::to_string(chrone_clock()) + "  \tsec      " + ss);};
    debug("Start   " + std::to_string(chrone_clock()));

    const auto& cmesh = create_mesh(N);
    logging("Split Hex Mesh to Tet");

    using Basis = DGBasisEvaluator<Order>;
    using QuadC = typename AutoQuadSelector<Basis::OrderBasis, GaussLegendreTet::Auto>::type;
    constexpr uInt DoFs = 5*Basis::NumBasis;

    LongVector<DoFs> U_n(cmesh.m_cells.size());

    #pragma omp parallel for schedule(dynamic)
    for(uInt cellId=0;cellId<cmesh.m_cells.size();cellId++){
        const auto& cell = cmesh.m_cells[cellId];
        const auto& rhoU_coef = Basis::func2coef([&](vector3f Xi)->DenseMatrix<5,1>{
            return {rho_Xi(cell,Xi),rhou_Xi(cell,Xi),rhov_Xi(cell,Xi),rhow_Xi(cell,Xi),rhoe_Xi(cell,Xi)};
        });
        for(uInt k=0;k<Basis::NumBasis;k++){
            MatrixView<DoFs,1,5,1>(U_n[cellId],5*k,0) = rhoU_coef[k];
        }
    }
    
    logging("Set Init Value");

    ImplicitConvection<Basis::OrderBasis> convection;
    
    LongVector<DoFs> mass(U_n.size());
    for(uInt cellId=0;cellId<cmesh.m_cells.size();cellId++){
        DenseMatrix<DoFs,DoFs> mass_matrix;
        for(uInt xgId=0; xgId<QuadC::num_points; ++xgId) {
            const auto& p = QuadC::points[xgId];
            auto phi = Basis::eval_all(p[0], p[1], p[2]);
            for(uInt i=0; i<Basis::NumBasis; ++i) {
                for(uInt k=0; k<5; ++k) {
                    mass[cellId](5*i + k, 0) += phi[i] * phi[i] * 
                        QuadC::weights[xgId] * cmesh.m_cells[cellId].compute_jacobian_det();
                }
            }
        }
    }
    
    print(std::array<std::string,8>{"#       time", "rel.err  rho",
                    "rel.err  u", "rel.err  v", "rel.err  w", 
                    "rel.err  e", "rel.err coef", "cpu time"});
    for(uInt kkkk=0;kkkk<1;kkkk++){
        
        Scalar Dt = 0.5;
        Scalar max_dt = Dt * std::pow((5.0/N),(Order+1));
        Scalar dt = max_dt;
        Scalar sub_t = 0;
        do{ 
            dt = std::min(max_dt,Dt-sub_t);
            sub_t += dt;

            
            LongVector<DoFs> U_k = U_n;

            for(uInt picard_iter = 0; picard_iter < 100; picard_iter++){

                LongVector<DoFs> rhs(U_n.size());
                BlockSparseMatrix<DoFs,DoFs> sparse_mat;
                
                convection.assemble(cmesh, U_k, Dt * (kkkk) + sub_t + 0.5 * dt,
                                    sparse_mat, rhs);
                
                for(uInt cellId = 0;cellId<cmesh.m_cells.size();cellId++){
                    const auto& diag = DenseMatrix<DoFs,DoFs>::Diag(mass[cellId]/dt);
                    sparse_mat.add_block(cellId, cellId, diag);
                    rhs[cellId] += diag.multiply(U_n[cellId]);
                }

                sparse_mat.finalize();



                const auto& Arr = get_sub_sparse<Basis::NumBasis,1,1,0,0>(sparse_mat);
                // {
                //     auto& block_val = Arr.storage.ell_blocks[0][0];
                //     auto& block_col = Arr.storage.ell_cols[0][0];
                //     debug("矩阵 Arr 的 0 行 第一个非零块的列号：" + std::to_string(block_col));
                //     debug(block_val);
                // }
                const auto& Arm = get_sub_sparse<Basis::NumBasis,1,3,0,1>(sparse_mat);
                const auto& Are = get_sub_sparse<Basis::NumBasis,1,1,0,4>(sparse_mat);
                const auto& br = get_sub_vector<Basis::NumBasis,1,0>(rhs);
                auto Ur = get_sub_vector<Basis::NumBasis,1,0>(U_k);
                // {
                //     debug(Ur[0]);
                // }
                           
             
                const auto& Amr = get_sub_sparse<Basis::NumBasis,3,1,1,0>(sparse_mat);
                // {
                //     auto& block_val = Amr.storage.ell_blocks[0][0];
                //     auto& block_col = Amr.storage.ell_cols[0][0];
                //     debug("矩阵 Arr 的 0 行 第一个非零块的列号：" + std::to_string(block_col));
                //     debug(block_val);
                // }
                const auto& Amm = get_sub_sparse<Basis::NumBasis,3,3,1,1>(sparse_mat);
                // {
                //     auto& block_val = Amm.storage.ell_blocks[0][0];
                //     auto& block_col = Amm.storage.ell_cols[0][0];
                //     debug("矩阵 Arr 的 0 行 第一个非零块的列号：" + std::to_string(block_col));
                //     debug(block_val);
                // }
                const auto& Ame = get_sub_sparse<Basis::NumBasis,3,1,1,4>(sparse_mat);                
                const auto& bm = get_sub_vector<Basis::NumBasis,3,1>(rhs);
                auto Um = get_sub_vector<Basis::NumBasis,3,1>(U_k);
                // {
                //     debug(Um[0]);
                // }
        
             
                const auto& Aer = get_sub_sparse<Basis::NumBasis,1,1,4,0>(sparse_mat);
                const auto& Aem = get_sub_sparse<Basis::NumBasis,1,3,4,1>(sparse_mat);
                const auto& Aee = get_sub_sparse<Basis::NumBasis,1,1,4,4>(sparse_mat);        
                       
                const auto& be = get_sub_vector<Basis::NumBasis,1,4>(rhs);
                auto Ue = get_sub_vector<Basis::NumBasis,1,4>(U_k);
                // {
                //     debug(Ue[0]);
                // }
                
                // debug(U_k);
                // debug(Ur);
                // debug(Um);
                // debug(Ue);
                
                // debug(rhs);
                // debug(br);
                // debug(bm);
                // debug(be);

                // debug(rhs - sparse_mat.multiply(U_k));
                // debug(be - (Aer.multiply(Ur) + Aem.multiply(Um) + Aee.multiply(Ue)));

                // debug(get_sub_vector<Basis::NumBasis,1,4>(sparse_mat.multiply(U_k)));
                // const auto err_r = get_sub_vector<Basis::NumBasis,1,0>(sparse_mat.multiply(U_k))-(Arr.multiply(Ur) + Arm.multiply(Um) + Are.multiply(Ue));
                // const auto err_m = get_sub_vector<Basis::NumBasis,3,1>(sparse_mat.multiply(U_k))-(Amr.multiply(Ur) + Amm.multiply(Um) + Ame.multiply(Ue));
                // const auto err_e = get_sub_vector<Basis::NumBasis,1,4>(sparse_mat.multiply(U_k))-(Aer.multiply(Ur) + Aem.multiply(Um) + Aee.multiply(Ue));
                // debug(vector3f{{err_r.dot(err_r),err_m.dot(err_m),err_e.dot(err_e)}});


                // debug("000");
                {   
                    const auto& submat_residual = rhs - sparse_mat.multiply(U_k);
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t before = " + oss.str());
                }

                auto sol_Arr = [&](LongVector<DoFs>& Urme_tmp, const auto& Ur, const auto& Um, const auto& Ue){
                    EigenSparseSolver<1*Basis::NumBasis,1*Basis::NumBasis> Arr_solver(Arr,br - Arm.multiply(Um) - Are.multiply(Ue));
                    LongVector<1*Basis::NumBasis> Ur_tmp = Arr_solver.SparseLU(Ur);
                    set_sub_vector<Basis::NumBasis,1,0>(Ur_tmp,Urme_tmp);
                    return Ur_tmp;
                };
                auto sol_Amm = [&](LongVector<DoFs>& Urme_tmp, const auto& Ur, const auto& Um, const auto& Ue){
                    EigenSparseSolver<3*Basis::NumBasis,3*Basis::NumBasis> Amm_solver(Amm,bm - Amr.multiply(Ur) - Ame.multiply(Ue));
                    LongVector<3*Basis::NumBasis> Um_tmp = Amm_solver.SparseLU(Um);
                    set_sub_vector<Basis::NumBasis,3,1>(Um_tmp,Urme_tmp);
                    return Um_tmp;
                };
                auto sol_Aee = [&](LongVector<DoFs>& Urme_tmp, const auto& Ur, const auto& Um, const auto& Ue){
                    EigenSparseSolver<1*Basis::NumBasis,1*Basis::NumBasis> Aee_solver(Aee,be - Aer.multiply(Ur) - Aem.multiply(Um));
                    LongVector<1*Basis::NumBasis> Ue_tmp = Aee_solver.SparseLU(Ue);
                    set_sub_vector<Basis::NumBasis,1,4>(Ue_tmp,Urme_tmp);
                    return Ue_tmp;
                };



                {   
                    LongVector<DoFs> Urme_tmp = U_k;
                    const auto& submat_residual = rhs - sparse_mat.multiply(Urme_tmp);
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t init : " + oss.str());
                }
                {   
                    LongVector<DoFs> Urme_tmp = U_k;
                    sol_Arr(Urme_tmp,Ur,Um,Ue);
                    const auto& submat_residual = rhs - sparse_mat.multiply(Urme_tmp);
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t only Arr : " + oss.str());
                }
                {   
                    LongVector<DoFs> Urme_tmp = U_k;
                    sol_Amm(Urme_tmp,Ur,Um,Ue);
                    const auto& submat_residual = rhs - sparse_mat.multiply(Urme_tmp);
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t only Amm = " + oss.str());
                }
                {   
                    LongVector<DoFs> Urme_tmp = U_k;
                    sol_Aee(Urme_tmp,Ur,Um,Ue);
                    const auto& submat_residual = rhs - sparse_mat.multiply(Urme_tmp);
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t only Aee = " + oss.str());
                }
                {   
                    LongVector<DoFs> Urme_tmp = U_k;
                    const auto& Ur_tmp = sol_Arr(Urme_tmp,Ur,Um,Ue);
                    sol_Aee(Urme_tmp,Ur_tmp,Um,Ue);
                    const auto& submat_residual = rhs - sparse_mat.multiply(Urme_tmp);
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Arr + Aee = " + oss.str());
                }

                
                {   
                 // ========== 新加调试代码：打印各项量级 ==========
                    // logging("  ||br||          = " + std::to_string(std::sqrt(br.dot(br))));
                    // logging("  ||Arr*Ur||      = " + std::to_string(std::sqrt((Arr.multiply(Ur)).dot(Arr.multiply(Ur)))));  //cg 06.17
                    // logging("  ||Arm*Um||      = " + std::to_string(std::sqrt((Arm.multiply(Um)).dot(Arm.multiply(Um)))));  //cg 06.17
                    // logging("  ||Are*Ue||      = " + std::to_string(std::sqrt((Are.multiply(Ue)).dot(Are.multiply(Ue)))));  //cg 06.17

                    // const auto& submat_residual = br - Arr.multiply(Ur) - Arm.multiply(Um) - Are.multiply(Ur);
                    const auto& submat_residual = br - Arr.multiply(Ur) - Arm.multiply(Um) - Are.multiply(Ue);  //cg 06.17
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Submat residual = " + oss.str());
                }
                {   
                    // logging("  ||bm||          = " + std::to_string(std::sqrt(bm.dot(bm))));
                    // logging("  ||Amr*Ur||      = " + std::to_string(std::sqrt((Amr.multiply(Ur)).dot(Amr.multiply(Ur)))));
                    // logging("  ||Amm*Um||      = " + std::to_string(std::sqrt((Amm.multiply(Um)).dot(Amm.multiply(Um)))));
                    // logging("  ||Ame*Ue||      = " + std::to_string(std::sqrt((Ame.multiply(Ue)).dot(Ame.multiply(Ue)))));

                    // const auto& submat_residual = bm - Amr.multiply(Ur) - Amm.multiply(Um) - Ame.multiply(Ur);
                    const auto& submat_residual = bm - Amr.multiply(Ur) - Amm.multiply(Um) - Ame.multiply(Ue);  //cg 06.17
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Submat residual = " + oss.str());
                }
                {   
                    // const auto& submat_residual = be - Aer.multiply(Ur) - Aem.multiply(Um) - Aee.multiply(Ur);
                    // const auto& AU = Aer.multiply(Ur) + Aem.multiply(Um) + Aee.multiply(Ue);
                    // const auto& AU_res = get_sub_vector<Basis::NumBasis,1,4>(sparse_mat.multiply(U_k));
                    // debug((AU-AU_res).dot(AU-AU_res));
                    // const auto& be_res = get_sub_vector<Basis::NumBasis,1,4>(rhs);
                    // debug((be-be_res).dot(be-be_res));
                    
                    // debug((be-AU_res).dot(be-AU_res));
                    // debug((be_res-AU).dot(be_res-AU));
                    // debug((be-AU).dot(be-AU));
                    // debug((be_res-AU_res).dot(be_res-AU_res));
                    const auto& submat_residual = be - Aer.multiply(Ur) - Aem.multiply(Um) - Aee.multiply(Ue);  //cg 06.17
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Submat residual = " + oss.str());
                }


                logging("Start linear solver");


                EigenSparseSolver<DoFs,DoFs> solver(sparse_mat,rhs);
                LongVector<DoFs> U_k_tmp = solver.BiCGSTAB(U_k);

                const auto& delta_x = U_k_tmp - U_k;
                U_k = U_k_tmp;
                Scalar delta = delta_x.dot(delta_x);
                delta = std::sqrt(delta);

                std::ostringstream oss;
                oss << delta;
                {
                    LongVector<DoFs> rhs(U_n.size());
                    BlockSparseMatrix<DoFs,DoFs> sparse_mat;
                    

                    // 用 保存的 U_k 代入、离散
                    convection.assemble(cmesh, U_k, Dt * (kkkk) + sub_t + 0.5 * dt,
                                        sparse_mat, rhs);
                    sparse_mat.finalize();
                    const auto& residual = (mass/dt*U_k + sparse_mat.multiply(U_k)) - (rhs + mass/dt * U_n);
                    std::ostringstream oss;
                    oss << std::sqrt(residual.dot(residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "  " + oss.str());

                }
                logging("Picard iter " + std::to_string(picard_iter) + "  " + oss.str());
                if(delta < 1e-10) break;
            }
            // Picard 结束，赋值，x 为 u^{n+1}
            U_n = U_k;

            logging("Iter  " + std::to_string(kkkk+1) + " \tSub TimeStep \t" + std::to_string(sub_t));

        }while(sub_t < Dt);


        Scalar curr_time = Dt * (kkkk+1);

        std::ofstream fp;
        fp.open(sol_filename(Basis::OrderBasis, kkkk+1, N));
        auto err_integral = [&](LongVector<5*QuadC::num_points> U_h,LongVector<5*QuadC::num_points> U_s){
            DenseMatrix<5,1> err_per_cells = DenseMatrix<5,1>::Zeros();
            DenseMatrix<5,1> sol_per_cells = 1e-47 * DenseMatrix<5,1>::Ones();
            for(uInt cellId=0;cellId<cmesh.m_cells.size();cellId++){
                const auto& cell = cmesh.m_cells[cellId];
                for(uInt xgId=0; xgId<QuadC::num_points; ++xgId) {
                    const DenseMatrix<5,1>& bUh = MatrixView<5*QuadC::num_points,1,5,1>(U_h[cellId],5*xgId,0);
                    const DenseMatrix<5,1>& bUs = MatrixView<5*QuadC::num_points,1,5,1>(U_s[cellId],5*xgId,0);
                    const DenseMatrix<5,1>& bUe = bUh - bUs;
                    const auto& weight = QuadC::weights[xgId] * cell.compute_jacobian_det();
                    // Scalar error_cell = uh[cellId][xgId] - us[cellId][xgId];

                    err_per_cells += pow(bUe,2) * weight;
                    sol_per_cells += pow(bUs,2) * weight;
                }
            }
            return pow(err_per_cells/sol_per_cells,0.5);
        };
        LongVector<5*QuadC::num_points> U_h(cmesh.m_cells.size());
        LongVector<5*QuadC::num_points> U_s(cmesh.m_cells.size());
        LongVector<DoFs> error(cmesh.m_cells.size());
        #pragma omp parallel for schedule(dynamic)
        for(uInt cellId=0;cellId<cmesh.m_cells.size();cellId++){
            const auto& cell = cmesh.m_cells[cellId];
            // rho[cellId] = 0.0;
            const auto& U_func = [&](vector3f Xi)->DenseMatrix<5,1>{
                return {
                    rho_Xi(cell,Xi,curr_time),
                    rhou_Xi(cell,Xi,curr_time),
                    rhov_Xi(cell,Xi,curr_time),
                    rhow_Xi(cell,Xi,curr_time),
                    rhoe_Xi(cell,Xi,curr_time)
                    };
            };
            const auto& U_coef = Basis::func2coef(U_func);
            for(uInt k=0;k<Basis::NumBasis;k++){
                const DenseMatrix<5,1>& block_coef = MatrixView<DoFs,1,5,1>(U_n[cellId],5*k,0) - U_coef[k];
                MatrixView<DoFs,1,5,1>(error[cellId],5*k,0) = block_coef;
            }
            for(uInt xgId=0; xgId<QuadC::num_points; ++xgId) {
                const auto& p = QuadC::points[xgId];
                const auto& pos = cell.transform_to_physical(p);
                const auto& value = Basis::eval_all(p[0],p[1],p[2]);
                const auto& U = Basis::template coef2filed<5,Scalar>(U_n[cellId],p);
                MatrixView<5*QuadC::num_points,1,5,1> block_U_h(U_h[cellId],5*xgId,0);
                MatrixView<5*QuadC::num_points,1,5,1> block_U_s(U_s[cellId],5*xgId,0);
                
                block_U_h = U;
                block_U_h[0] = U[0];
                block_U_s = U_func(p);
            }
        }

        for(uInt cellId=0;cellId<cmesh.m_cells.size();cellId++){
            const auto& cell = cmesh.m_cells[cellId];
            for(uInt xgId=0; xgId<QuadC::num_points; ++xgId) {
                const auto& p = QuadC::points[xgId];
                const auto& pos = cell.transform_to_physical(p);
                const auto& bUh = MatrixView<5*QuadC::num_points,1,5,1>(U_h[cellId],5*xgId,0);
                const auto& bUs = MatrixView<5*QuadC::num_points,1,5,1>(U_s[cellId],5*xgId,0);
                fp <<std::setprecision(16)<< pos[0] << "  " <<std::setprecision(16)<< pos[1] << "  " <<std::setprecision(16)<< pos[2]
                 << "  " <<std::setprecision(16)<<  bUh[0] << "  " <<std::setprecision(16)<<  bUs[0] 
                 << "  " <<std::setprecision(16)<<  bUh[1] << "  " <<std::setprecision(16)<<  bUs[1] 
                 << "  " <<std::setprecision(16)<<  bUh[2] << "  " <<std::setprecision(16)<<  bUs[2]
                 << "  " <<std::setprecision(16)<<  bUh[3] << "  " <<std::setprecision(16)<<  bUs[3] 
                 << "  " <<std::setprecision(16)<<  bUh[4] << "  " <<std::setprecision(16)<<  bUs[4] << std::endl;
            }
        }
        fp.close();

        const auto& U_err = err_integral(U_h,U_s);
        print(std::array<Scalar,8>{curr_time, U_err[0], U_err[1],U_err[2],U_err[3],U_err[4], 
                std::sqrt(error.dot(error)/U_n.dot(U_n)), chrone_clock()});
        // print(vector3f{curr_time, std::sqrt(err.dot(err)/rho.dot(rho)),chrone_clock()});
    }
}

// template void Run<0>(uInt);
// template void Run<1>(uInt);
// template void Run<2>(uInt);
// template void Run<3>(uInt);
// template void Run<4>(uInt);
// template void Run<5>(uInt);
// template void Run<6>(uInt);
// template void Run<7>(uInt);
// template void Run<8>(uInt);
// template void Run<9>(uInt);



int main(int argc, char** argv){
    
    omp_set_num_threads(get_phy_cpu());
    Eigen::setNbThreads(get_phy_cpu());


    int order = std::stoi(argv[1]);
    int meshN = std::stoi(argv[2]);
                             
    // if(order == 0) Run<0>(meshN);
    if(order == 1) Run<1>(meshN);
    if(order == 2) Run<2>(meshN);
    if(order == 3) Run<3>(meshN);
    if(order == 4) Run<4>(meshN);
    // if(order == 5) Run<5>(meshN);
    // if(order == 6) Run<6>(meshN);
    // if(order == 7) Run<7>(meshN);
    // if(order == 8) Run<8>(meshN);
    // if(order == 9) Run<9>(meshN);
}