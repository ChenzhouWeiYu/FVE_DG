#include "Type.h"
#include "OptimizedMesh/OptimizedMesh.h"
#include "Matrix/Matrix.h"
#include "DG_Schemes/DG_Basis.h"
#include "DG_Schemes/Convection.h"
#include "DG_Schemes/ImplicitConvection.h"
#include "Eigen/Sparse"
#include "unsupported/Eigen/IterativeSolvers"

#include "exact.h"

typedef Eigen::SparseMatrix<double,Eigen::RowMajor> EigenSpMat;
typedef Eigen::Triplet<double> Triplet;



constexpr Scalar param_gamma = 1.4;

template<typename Type>
inline Type rho_xyz(Type x, Type y, Type z, Type t){
    return 1.0 + 0.99*std::sin(x+y - 2*t);
}
template<typename Type>
inline Type u_xyz(Type x, Type y, Type z, Type t){
    return 1.0;
}
template<typename Type>
inline Type v_xyz(Type x, Type y, Type z, Type t){
    return 1.0;
}
template<typename Type>
inline Type w_xyz(Type x, Type y, Type z, Type t){
    return 0.0;
}
template<typename Type>
inline Type p_xyz(Type x, Type y, Type z, Type t){
    // return (param_gamma-1)*rho_xyz(x,y,z,t)*e_xyz(x,y,z);
    return 1.0;
}
template<typename Type>
inline Type e_xyz(Type x, Type y, Type z, Type t){
    // Scalar r2 = x*x + y*y;
    // constexpr Scalar r_ds2 = 1.0Q/(2.0Q * 0.1Q * 0.1Q);
    // return 1e-12 + 0.979264*M_1_PI*r_ds2* std::exp(-r2*r_ds2);
    Type p = p_xyz(x,y,z,t);
    Type u = u_xyz(x,y,z,t);
    Type v = v_xyz(x,y,z,t);
    return p/rho_xyz(x,y,z,t)/(param_gamma-1) + 0.5*(u*u+v*v);
    
}

#define Filed_Func(filedname) \
inline Scalar filedname##_xyz(const vector3f& xyz, Scalar t){\
    Scalar x = xyz[0], y = xyz[1], z = xyz[2];\
    return filedname##_xyz(x,y,z,t);\
}\
inline Scalar filedname##_Xi(const CompTetrahedron& cell, const vector3f& Xi, Scalar t){\
    const vector3f& xyz = cell.transform_to_physical(Xi);\
    Scalar x = xyz[0], y = xyz[1], z = xyz[2];\
    return filedname##_xyz(x,y,z,t);\
}\
inline Scalar rho##filedname##_Xi(const CompTetrahedron& cell, const vector3f& Xi, Scalar t){\
    const vector3f& xyz = cell.transform_to_physical(Xi);\
    Scalar x = xyz[0], y = xyz[1], z = xyz[2];\
    return rho_xyz(x,y,z,t)*filedname##_xyz(x,y,z,t);\
}
Filed_Func(rho);
Filed_Func(u);
Filed_Func(v);
Filed_Func(w);
Filed_Func(p);
Filed_Func(e);

#undef Filed_Func



template<uInt Order>
void Run(uInt N){
    auto chrono_start = std::chrono::steady_clock::now();
    auto chrone_clock = [&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-chrono_start).count();};
    auto logging = [&](std::string ss){debug("Time  " + std::to_string(chrone_clock()) + "  \tsec      " + ss);};
    debug("Start   " + std::to_string(chrone_clock()));
    // GeneralMesh mesh = OrthHexMesh({-1.1, -1.1, -1.1/N},{1.1, 1.1, 1.1/N},{N,N,1});
    GeneralMesh mesh = OrthHexMesh({0.0, 0.0, 0.0},{2*M_PI, 2*M_PI, 2*M_PI/N},{N,N,1});
    mesh.split_hex5_scan();                                   
    mesh.rebuild_cell_topology();                             
    mesh.validate_mesh();                                     
    ComputingMesh cmesh(mesh);                                
    cmesh.m_boundary.resize(cmesh.m_faces.size());            
    for(uInt fid=0;fid<cmesh.m_faces.size();fid++){           
        if(cmesh.m_faces[fid].m_neighbor_cells[1]==uInt(-1)){ 
            const auto& face = cmesh.m_faces[fid];            
            if(std::abs(face.m_normal[2])>0.5)                
                cmesh.m_boundary[fid] = 2;
            else
                cmesh.m_boundary[fid] = 1;
        }
    }
    logging("Split Hex Mesh to Tet");

    using Basis = DGBasisEvaluator<Order>;
    using QuadC = typename AutoQuadSelector<Basis::OrderBasis, GaussLegendreTet::Auto>::type;

    LongVector<5*Basis::NumBasis> x(mesh.m_cells.size());

    #pragma omp parallel for schedule(dynamic)
    for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
        const auto& cell = cmesh.m_cells[cid];
        const auto& rho_coef  = Basis::func2coef([&](vector3f Xi){return rho_Xi(cell,Xi); });
        const auto& rhou_coef = Basis::func2coef([&](vector3f Xi){return rhou_Xi(cell,Xi);});
        const auto& rhov_coef = Basis::func2coef([&](vector3f Xi){return rhov_Xi(cell,Xi);});
        const auto& rhow_coef = Basis::func2coef([&](vector3f Xi){return rhow_Xi(cell,Xi);});
        const auto& rhoe_coef = Basis::func2coef([&](vector3f Xi){return rhoe_Xi(cell,Xi);});
        #pragma omp simd
        for(uInt k=0;k<Basis::NumBasis;k++){
            x[cid][5*k+0] = rho_coef[k];
            x[cid][5*k+1] = rhou_coef[k];
            x[cid][5*k+2] = rhov_coef[k];
            x[cid][5*k+3] = rhow_coef[k];
            x[cid][5*k+4] = rhoe_coef[k];
        }
    }
    
    logging("Set Init Value");


    ExplicitConvection<Basis::OrderBasis> explicitconvection;
    ImplicitConvection<Basis::OrderBasis> convection;
    // LongVector<5*Basis::NumBasis> rhs(x.size());
    // BlockSparseMatrix<5*Basis::NumBasis,5*Basis::NumBasis> sparse_mat;
    
    // convection.assemble(cmesh,x,0.0,sparse_mat,rhs);
    // sparse_mat.finalize();


    // const auto& ex_conv = explicitconvection.eval(cmesh, x, 0.0);
    // const auto& im_conv = sparse_mat.multiply(x) - rhs;
    // // print(im_conv);
    // // print(ex_conv);
    // print(max(abs(im_conv - ex_conv)));  //  1e-14

    
    // logging("1111111111");
    
    LongVector<5*Basis::NumBasis> mass(x.size());
    for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
        DenseMatrix<5*Basis::NumBasis,5*Basis::NumBasis> mass_matrix;
        for(uInt g=0; g<QuadC::num_points; ++g) {
            const auto& p = QuadC::points[g];
            auto phi = Basis::eval_all(p[0], p[1], p[2]);
            for(uInt i=0; i<Basis::NumBasis; ++i) {
                for(uInt k=0; k<5; ++k) {
                    mass[cid](5*i + k, 0) += phi[i] * phi[i] * 
                        QuadC::weights[g] * cmesh.m_cells[cid].compute_jacobian_det();
                }
            }
        }
    }
    
    // logging("22222222222");
    // print(sparse_mat);
    // sparse_mat.output_as_scalar("./12345.txt");

    // // auto chrono_start = std::chrono::steady_clock::now();
    // // auto chrone_clock = [&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-chrono_start).count();};
    print(std::array<std::string,8>{"#       time", "rel.err  rho",
                    "rel.err  u", "rel.err  v", "rel.err  w", 
                    "rel.err  e", "rel.err coef", "cpu time"});
    for(uInt kkkk=0;kkkk<10;kkkk++){    
        // Scalar dt = 0.001;
        
        // logging("Iter    " + std::to_string(kkkk+1));
        
        // Scalar h = 10.0/N * 0.5773502691896257645091487805019574556476 / 3.0,   vel = 3.12;
        // Scalar max_dt = 0.9 * h * 1.0/(2*Order+1) * 1.0 / vel;
        
        Scalar Dt = 0.1;


        Scalar max_dt = Dt * std::pow((2*M_PI/N),(Order+1));
        // max_dt = Dt;

        Scalar dt = max_dt;
        uInt kk = 0;
        
        // debug(dt);
        Scalar sub_t = 0;
        do{ 


            // Scalar h = 0.0704416 * 1.0/N,   vel = 0.0;
            // for(uInt fid=0; fid<cmesh.m_faces.size(); ++fid) {
            //     // std::array<Scalar,5> cumsum = {0,0,0,0,0};
            //     const auto& face = cmesh.m_faces[fid];
            //     const auto& cells = face.m_neighbor_cells;
            //     if(cells[1] == uInt(-1)) continue;
            //     const auto& coef_L = x[cells[0]];
            //     const auto& coef_R = x[cells[1]];
            //     auto xi_L = convection.transform_to_cell(face, {1.0/3,1.0/3}, 0);
            //     auto xi_R = convection.transform_to_cell(face, {1.0/3,1.0/3}, 1);
            //     auto basis_L = Basis::eval_all(xi_L[0], xi_L[1], xi_L[2]);
            //     auto basis_R = Basis::eval_all(xi_R[0], xi_R[1], xi_R[2]);

            //     DenseMatrix<5,1> U_L, U_R;
            //     for(uInt bid=0; bid<Basis::NumBasis; ++bid) {
            //         for(uInt k=0; k<5; ++k) {
            //             U_L[k] += basis_L[bid] * coef_L[5*bid + k];
            //             U_R[k] += basis_R[bid] * coef_R[5*bid + k];
            //         }
            //     }
            //     const Scalar lambda = convection.compute_max_wave_speed(U_L, U_R);
            //     vel = std::max(lambda,vel);
            // }


            // Scalar max_dt = std::min(sub_dt * 0.25,   0.25 * h * 1.0/(2*Order+1) * 1.0 / (vel+(2*Order+1)*param_mu/h));
            // if(max_dt<1e-6)break;
            // Scalar dt = std::min(max_dt,sub_dt-sub_t);
            // sub_t += dt;

            dt = 1.0 * max_dt;
            // if(kkkk>1 || kk>3){
            //     dt = 1.0 * max_dt;
            // }else if(kk==0){
            //     dt = 0.1 * max_dt;
            // }else if(kk==1){
            //     dt = 0.2 * max_dt;
            // }else if(kk==2){
            //     dt = 0.3 * max_dt;
            // }else if(kk==3){
            //     dt = 0.4 * max_dt;
            // }
            // kk++;

            dt = std::min(dt,Dt-sub_t);
            sub_t += dt;
            

            
            LongVector<5*Basis::NumBasis> dx = x;

            for(uInt picard_iter = 0; picard_iter < 100; picard_iter++){

                LongVector<5*Basis::NumBasis> rhs(x.size());
                BlockSparseMatrix<5*Basis::NumBasis,5*Basis::NumBasis> sparse_mat;
                

                // 用 保存的 dx 代入、离散
                convection.assemble(cmesh, dx, Dt * (kkkk) + sub_t + 0.5 * dt,
                                    sparse_mat, rhs);
                sparse_mat.finalize();
                // sparse_mat.output_as_scalar("./12345.txt");

                std::vector<Triplet> tripletList;
                // debug(sparse_mat.storage.ell_blocks.size());
                // debug(sparse_mat.storage.ell_cols.size());
                // debug(sparse_mat.storage.ell_max_per_row);
                // debug(sparse_mat.storage.csr_blocks.size());
                // debug(sparse_mat.storage.csr_cols.size());
                // debug(sparse_mat.storage.csr_row_ptr.size());
                // debug(sparse_mat.num_block_rows);
                // break;
                uInt ELL_row_length = sparse_mat.storage.ell_blocks.size();
                uInt ELL_max_row = sparse_mat.storage.ell_max_per_row;
                uInt CSR_nnz = sparse_mat.storage.csr_blocks.size();
                uInt NNZ = (ELL_row_length * ELL_max_row + CSR_nnz) * 5*Basis::NumBasis * 5*Basis::NumBasis;
                uInt ROWs = sparse_mat.num_block_rows * 5*Basis::NumBasis;
                tripletList.reserve(NNZ);

                for(uInt brow=0; brow<sparse_mat.num_block_rows; ++brow){
                    for(uInt i=0; i<sparse_mat.storage.ell_max_per_row; ++i){
                        const uInt bcol = sparse_mat.storage.ell_cols[brow][i];
                        if(bcol == sparse_mat.invalid_index) continue;
                        const auto& block = sparse_mat.storage.ell_blocks[brow][i];
                        for(uInt row=0;row<block.rows();row++){
                            for(uInt col=0;col<block.cols();col++){
                                // if(std::abs(block(row,col))>1e-15)
                                tripletList.push_back(Triplet(brow*block.rows()+row,bcol*block.cols()+col,block(row,col)));
                            }
                        }
                    }

                    const uInt start = sparse_mat.storage.csr_row_ptr[brow];
                    const uInt end = sparse_mat.storage.csr_row_ptr[brow+1];
                    for(uInt idx = start; idx < end; ++idx) {
                        const uInt bcol = sparse_mat.storage.csr_cols[idx];
                        const auto& block = sparse_mat.storage.csr_blocks[idx];

                        for(uInt row=0;row<block.rows();row++){
                            for(uInt col=0;col<block.cols();col++){
                                // if(std::abs(block(row,col))>1e-15)
                                tripletList.push_back(Triplet(brow*block.rows()+row,bcol*block.cols()+col,block(row,col)));
                            }
                        }
                    }
                }
                // for(...)
                // {
                // // ...
                // tripletList.push_back(T(i,j,v_ij));
                // }
                EigenSpMat Eigen_A(ROWs,ROWs);
                // debug(Eigen_A.nonZeros());
                // debug(Eigen_A.rows());
                // debug(Eigen_A.cols());
                // debug(Eigen_A.size());
                Eigen::VectorXd Eigen_b(ROWs), Eigen_x(ROWs);


                // Scalar rdt = 1/dt;
                // if(picard_iter < 8){
                //     if(picard_iter < 2){
                //         rdt = 1e-2 / dt;
                //     }else if(picard_iter < 4){
                //         rdt = 2e-2 / dt;
                //     }else if(picard_iter < 6){
                //         rdt = 3e-2 / dt;
                //     }else if(picard_iter < 8){
                //         rdt = 4e-2 / dt;
                //     }
                // }else if(picard_iter < 10){
                //     rdt = 2e-1 / dt;
                // }else if(picard_iter < 12){
                //     rdt = 3e-1 / dt;
                // }else if(picard_iter < 14){
                //     rdt = 4e-1 / dt;
                // }
                // 但是 时间 u^n 还是旧的 x
                for(uInt r = 0; r < mass.size(); r++){
                    for(uInt rr = 0; rr < mass[r].size(); rr++){
                        tripletList.push_back(Triplet(r*mass[r].size()+rr, r*mass[r].size()+rr, mass[r][rr]/dt));
                    }
                }
                rhs += mass/dt * x;



                Eigen_A.setFromTriplets(tripletList.begin(), tripletList.end());
                for(uInt r=0;r<rhs.size();r++){
                    auto& block = rhs[r];
                    for(uInt rr=0;rr<block.size();rr++){
                        Eigen_b[r * block.size() + rr] = block[rr];
                    }
                }

                Eigen::DGMRES<EigenSpMat> Eigen_solver;
                // Eigen::SparseLU<EigenSpMat> Eigen_solver;
                // Eigen_solver.setMaxIterations(200);
                Eigen_solver.setTolerance(1e-12);
                Eigen_solver.compute(Eigen_A);
                if(Eigen_solver.info()!=Eigen::Success) {
                    sparse_mat.output_as_scalar("12345.txt");
                    print(dx);
                    throw std::runtime_error("Matrix decomposition failed");
                }
                Eigen_x = Eigen_solver.solve(Eigen_b);
                // if(Eigen_solver.info()!=Eigen::Success) {
                //     throw std::runtime_error("Matrix solving failed");
                // }

                // 解出来一个新的 dx
                Scalar delta = 0.0;
                for(uInt r=0;r<dx.size();r++){
                    auto& block = dx[r];
                    for(uInt rr=0;rr<block.size();rr++){
                        Scalar xxx = block[rr] - Eigen_x[r * block.size() + rr];
                        delta += xxx*xxx;
                        block[rr] = Eigen_x[r * block.size() + rr];
                    }
                }
                

                // 保正限制器
                // 对密度 rho，其平均值为 rho_avg = c_0 * phi_0，而 phi_0 = 1，所以 rho_avg = c_0
                // 考虑密度 rho 在每一个积分点上，都 > 0，或者 > eps = 1e-12, 1e-16 这种
                //   具体来说  对每一个积分点处，计算每一个 rho(x_g)，要求 rho(x_g) > eps
                //   实现方法：找 rho(x_g) 的最小值 rho_min
                //   rho_new = theta * (rho - rho_avg) + rho_avg = theta*rho + (1-theta)*rho_avg
                //       theta = min { 1, (rho_avg-eps)/(rho_avg-rho_min) }
                //   守恒型：rho_new_avg = theta*rho_avg + (1-theta)*rho_avg = rho_avg
                #pragma omp parallel for schedule(dynamic)
                for(uInt cid = 0; cid < cmesh.m_cells.size(); cid++){
                    Scalar eps = 1e-16;
                    std::array<Scalar, QuadC::num_points> rho_qp;
                    rho_qp.fill(0.0);
                    for (uInt xgi = 0; xgi < QuadC::num_points; xgi++) {
                        const auto& xg = QuadC::points[xgi];
                        const auto& basis = Basis::eval_all(xg[0],xg[1],xg[2]);
                        for (uInt k = 0; k < Basis::NumBasis; k++) {
                            rho_qp[xgi] += basis[k] * dx[cid][0 + 5*k];  // Density
                        }
                    }
                    Scalar rho_min = *std::min_element(rho_qp.begin(), rho_qp.end());
                    Scalar rho_avg = dx[cid][0];
                    if (rho_min < eps) {
                        Scalar effective_rho_avg = std::max(rho_avg, eps);
                        Scalar numerator = effective_rho_avg - eps;
                        Scalar denominator = effective_rho_avg - rho_min;

                        Scalar theta = (denominator < 1e-32) ? 0.0 : std::min(1.0, numerator / denominator);

                        // Apply to high-order modes of density
                        for (uInt k = 1; k < Basis::NumBasis; k++) {
                            dx[cid][0 + 5*k] *= theta;
                        }
                    }
                    
                    // print(dx[cid][0]);
                    // dx[cid][0] = std::max(dx[cid][0],1e-16);
                    // std::array<Scalar,QuadC::num_points> rho;
                    // rho.fill(0.0);
                    // for(uInt xgi=0;xgi<QuadC::num_points;xgi++){
                    //     const auto& basis = Basis::eval_all(xg[0],xg[1],xg[2]);
                    //     for(uInt k=0;k<Basis::NumBasis;k++){
                    //         rho[xgi] += basis[k]*dx[cid][0+5*k];
                    //     }
                    // }
                    // Scalar rho_min = *std::min_element(rho.begin(), rho.end());
                    // Scalar rho_avg = std::max(dx[cid][0],eps);
                    // Scalar theta = std::min(1.0,std::abs(rho_avg-eps)/(std::abs(rho_avg-rho_min)+1e-32));
                    // 为了实现 rho_new = theta*rho + (1-theta)*rho_avg
                    // 实质上是修改系数：c_k_new = theta*c_k,  c_0_new = c_0
                    // for(uInt k=1;k<Basis::NumBasis;k++){
                    //     dx[cid][0+5*k] *= theta;
                    // }
                    
                    // print(dx[cid][0]);

                    // // 下面开始处理压强
                    // std::array<Scalar, QuadC::num_points> p_qp;
                    // p_qp.fill(0.0);
                    // for (uInt xgi = 0; xgi < QuadC::num_points; xgi++) {
                    //     const auto& basis = Basis::eval_all(QuadC::points[xgi][0], QuadC::points[xgi][1], QuadC::points[xgi][2]);

                    //     Scalar rho = 0.0, rhou = 0.0, rhov = 0.0, rhow = 0.0, rhoE = 0.0;
                    //     for (uInt k = 0; k < Basis::NumBasis; k++) {
                    //         rho  += basis[k] * dx[cid][0 + 5*k];
                    //         rhou += basis[k] * dx[cid][1 + 5*k];
                    //         rhov += basis[k] * dx[cid][2 + 5*k];
                    //         rhow += basis[k] * dx[cid][3 + 5*k];
                    //         rhoE += basis[k] * dx[cid][4 + 5*k];
                    //     }

                    //     Scalar u = rhou / rho;
                    //     Scalar v = rhov / rho;
                    //     Scalar w = rhow / rho;
                    //     Scalar ke = 0.5 * (rhou*u + rhov*v + rhow*w);
                    //     Scalar p = (param_gamma - 1.0) * (rhoE - ke);
                    //     p_qp[xgi] = p;
                    // }
                    // Scalar p_min = *std::min_element(p_qp.begin(), p_qp.end());
                    // Scalar p_avg = (param_gamma - 1.0) * (dx[cid][4] - 0.5 * (dx[cid][1]*dx[cid][1] + dx[cid][2]*dx[cid][2] + dx[cid][3]*dx[cid][3]) / dx[cid][0]);

                    // if (p_min < eps) {
                    //     Scalar effective_p_avg = std::max(p_avg, eps);
                    //     Scalar numerator = effective_p_avg - eps;
                    //     Scalar denominator = effective_p_avg - p_min;

                    //     Scalar theta = (denominator < 1e-32) ? 0.0 : std::min(1.0, numerator / denominator);

                    //     // Apply to high-order modes of energy
                    //     for (uInt k = 1; k < Basis::NumBasis; k++) {
                    //         dx[cid][4 + 5*k] *= theta;
                    //     }
                    // }
                    // std::array<DenseMatrix<5,1>,Basis::NumBasis> U;
                    // // std::array<Scalar,Basis::NumBasis> t_xg;
                    // Scalar t_xg;
                    // U.fill(DenseMatrix<5,1>::Zeros());
                    // Scalar theta_2 = 0.0;
                    // for(uInt xgi=0;xgi<QuadC::num_points;xgi++){
                    //     const auto& xg = QuadC::points[xgi];
                    //     const auto& basis = Basis::eval_all(xg[0],xg[1],xg[2]);
                    //     for(uInt k=0;k<Basis::NumBasis;k++){
                    //         U[xgi][0] += basis[k]*dx[cid][0+5*k];
                    //         U[xgi][1] += basis[k]*dx[cid][1+5*k];
                    //         U[xgi][2] += basis[k]*dx[cid][2+5*k];
                    //         U[xgi][3] += basis[k]*dx[cid][3+5*k];
                    //         U[xgi][4] += basis[k]*dx[cid][4+5*k];
                    //     }
                    //     Scalar  t_xg_i = 0.0;
                    //     std::array<Scalar,11> pp;
                    //     for(uInt pi=0; pi<11; pi++){
                    //         Scalar lambda = 0.1*pi;
                    //         Scalar rho = std::max(lambda*U[xgi][0] + (1-lambda)*U[0][0],1e-16);
                    //         Scalar u = (lambda*U[xgi][1] + (1-lambda)*U[0][1])/rho;
                    //         Scalar v = (lambda*U[xgi][2] + (1-lambda)*U[0][2])/rho;
                    //         Scalar w = (lambda*U[xgi][3] + (1-lambda)*U[0][3])/rho;
                    //         Scalar e = (lambda*U[xgi][4] + (1-lambda)*U[0][4])/rho;
                    //         Scalar p = (param_gamma-1)*rho*(e-0.5*(u*u+v*v+w*w));
                    //         pp[pi] = p;
                    //     }
                    //     // print(pp);
                    //     for(uInt pi=0; pi<10; pi++){
                    //         if((pp[pi]-1e-16)*(pp[pi+1]-1e-16)<0){
                    //             // 二分法，对于存在 零点 的区间，找零点位置
                    //             Scalar xl = 0.1 * pi;
                    //             Scalar xr = 0.1 * (pi + 1);
                    //             Scalar xm = 0.5 * (xl + xr);
                    //             Scalar pl = pp[pi];
                    //             Scalar pr = pp[pi+1];
                    //             for(uInt ii = 0; ii<10; ii++){
                    //                 Scalar lambda = xm;
                    //                 Scalar rho = std::max(lambda*U[xgi][0] + (1-lambda)*U[0][0],1e-16);
                    //                 Scalar u = (lambda*U[xgi][1] + (1-lambda)*U[0][1])/rho;
                    //                 Scalar v = (lambda*U[xgi][2] + (1-lambda)*U[0][2])/rho;
                    //                 Scalar w = (lambda*U[xgi][3] + (1-lambda)*U[0][3])/rho;
                    //                 Scalar e = (lambda*U[xgi][4] + (1-lambda)*U[0][4])/rho;
                    //                 Scalar p = (param_gamma-1)*rho*(e-0.5*(u*u+v*v+w*w));
                    //                 if((pl-1e-16)*(p-1e-16)<0){
                    //                     xr = xm;
                    //                 }
                    //                 else{
                    //                     xl = xm;
                    //                 }
                    //                 // print(vector3f{xl,xm,xr});
                    //                 // 或者 xr-xl 充分小？
                    //             }
                    //             // 二分法结束，最终的 xm 保存，退出
                    //             t_xg_i = xm;
                    //             break;
                    //         }
                    //     }
                    //     // t_xg_i 就是积分点 xg 处，满足要求 p = eps 的 t_xg (theta2)
                    //     if (t_xg_i < theta_2)  theta_2 = t_xg_i;
                    // }
                    // dx[cid] *= theta_2;

                }






                // const auto& xxxx = sparse_mat.multiply(dx);
                // Scalar delta = xxxx.dot(xxxx);
                delta = std::sqrt(delta);
                std::ostringstream oss;
                oss << delta;
                // std::ostringstream oss1;
                // oss1 << Eigen_solver.error();
                // logging("BiCGSTAB it " + std::to_string(Eigen_solver.iterations()) + "  " + oss1.str());
                {
                    LongVector<5*Basis::NumBasis> rhs(x.size());
                    BlockSparseMatrix<5*Basis::NumBasis,5*Basis::NumBasis> sparse_mat;
                    

                    // 用 保存的 dx 代入、离散
                    convection.assemble(cmesh, dx, Dt * (kkkk) + sub_t + 0.5 * dt,
                                        sparse_mat, rhs);
                    sparse_mat.finalize();
                    const auto& residual = (mass/dt*dx + sparse_mat.multiply(dx)) - (rhs + mass/dt * x);
                    std::ostringstream oss;
                    oss << std::sqrt(residual.dot(residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "  " + oss.str());

                }
                logging("Picard iter " + std::to_string(picard_iter) + "  " + oss.str());
                if(delta < 1e-10) break;
            }
            // Picard 结束，赋值，x 为 u^{n+1}
            x = dx;// - dt / mass * dx;

            // const auto& Eigen_r = Eigen_A*Eigen_x-Eigen_b;
            // LongVector<5*Basis::NumBasis> dx(x.size());
            // for(uInt r=0;r<dx.size();r++){
            //     auto& block = dx[r];
            //     for(uInt rr=0;rr<block.size();rr++){
            //         block[rr] = Eigen_b[r * block.size() + rr];
            //     }
            // }
            // print(sparse_mat.multiply(dx)-rhs);
            



            // debug(vector2f{dt,sub_t});
            logging("Iter  " + std::to_string(kkkk+1) + " \tSub TimeStep \t" + std::to_string(sub_t));

            // const auto& k1 = - 1/b * convection.eval(cmesh, x, 
            //                                                 0.1 * (kkkk) + sub_t);
            // const auto& x1 = x + 1.0/2.0 * dt * k1;
            // // for(uInt cid=0;cid<x1.size();cid++){
            // //     x1[cid][5*0+0] = std::max(1e-4,x1[cid][5*0+0]);
            // // }
            // const auto& k2 = - 1/b * convection.eval(cmesh, x1,  
            //                                                 0.1 * (kkkk) + sub_t + 1.0/2.0 * dt); 
            // x = x + dt * k2;
            // for(uInt cid=0;cid<x1.size();cid++){
            //     x[cid][5*0+0] = std::max(1e-4,x[cid][5*0+0]);
            // }

            // const auto& h1 = - 1/b * convection.eval(cmesh, x);
            // const auto& h2 = - 1/b * convection.eval(cmesh, x + 0.5 * dt * h1);
            // const auto& h3 = - 1/b * convection.eval(cmesh, x + 0.5 * dt * h2);
            // const auto& h4 = - 1/b * convection.eval(cmesh, x + 1.0 * dt * h3);
            // x = x + dt * (h1 + 2*h2 + 2*h3 + h4) / 6;

        }while(sub_t < Dt);


        Scalar curr_time = Dt * (kkkk+1);

        std::ofstream fp;
        fp.open("./Order_" + std::to_string(Basis::OrderBasis) + "/rho_T_" + std::to_string(kkkk+1) 
                + "_N_" + std::to_string(N) + ".txt");
        auto err_integral = [&](LongVector<QuadC::num_points> uh,LongVector<QuadC::num_points> us){
            Scalar err_per_cells = 0.0, sol_per_cells = 1e-47;
            for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
                const auto& cell = cmesh.m_cells[cid];
                for(uInt g=0; g<QuadC::num_points; ++g) {
                    const auto& weight = QuadC::weights[g] * cell.compute_jacobian_det();
                    Scalar error_cell = uh[cid][g] - us[cid][g];
                    err_per_cells += std::pow(error_cell,2) * weight;
                    sol_per_cells += std::pow(us[cid][g],2) * weight;
                }
            }
            return std::pow(err_per_cells/sol_per_cells,0.5);
        };

        LongVector<QuadC::num_points> rho_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> rho_s(cmesh.m_cells.size());
        LongVector<QuadC::num_points> u_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> u_s(cmesh.m_cells.size());
        LongVector<QuadC::num_points> v_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> v_s(cmesh.m_cells.size());
        LongVector<QuadC::num_points> w_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> w_s(cmesh.m_cells.size());
        LongVector<QuadC::num_points> e_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> e_s(cmesh.m_cells.size());
        LongVector<5 * Basis::NumBasis> error(cmesh.m_cells.size());

        #pragma omp parallel for schedule(dynamic)
        for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
            const auto& cell = cmesh.m_cells[cid];
            // rho[cid] = 0.0;
            for(uInt g=0; g<QuadC::num_points; ++g) {
                const auto& p = QuadC::points[g];
                const auto& pos = cell.transform_to_physical(p);
                const auto& value = Basis::eval_all(p[0],p[1],p[2]);
                for(uInt k=0;k<Basis::NumBasis;k++){
                    rho_h[cid][g] += value[k] * x[cid][5*k+0];
                    u_h[cid][g] += value[k] * x[cid][5*k+1];
                    v_h[cid][g] += value[k] * x[cid][5*k+2];
                    w_h[cid][g] += value[k] * x[cid][5*k+3];
                    e_h[cid][g] += value[k] * x[cid][5*k+4];
                }
                u_h[cid][g] /= rho_h[cid][g];
                v_h[cid][g] /= rho_h[cid][g];
                w_h[cid][g] /= rho_h[cid][g];
                e_h[cid][g] /= rho_h[cid][g];
                rho_s[cid][g] = rho_xyz(pos,curr_time);
                u_s[cid][g] = u_xyz(pos,curr_time);
                v_s[cid][g] = v_xyz(pos,curr_time);
                w_s[cid][g] = w_xyz(pos,curr_time);
                e_s[cid][g] = e_xyz(pos,curr_time);
            }
            const auto& rho_coef  = Basis::func2coef([&](vector3f Xi){return rho_Xi( cell,Xi,curr_time);});
            const auto& rhou_coef = Basis::func2coef([&](vector3f Xi){return rhou_Xi(cell,Xi,curr_time);});
            const auto& rhov_coef = Basis::func2coef([&](vector3f Xi){return rhov_Xi(cell,Xi,curr_time);});
            const auto& rhow_coef = Basis::func2coef([&](vector3f Xi){return rhow_Xi(cell,Xi,curr_time);});
            const auto& rhoe_coef = Basis::func2coef([&](vector3f Xi){return rhoe_Xi(cell,Xi,curr_time);});
            #pragma omp simd
            for(uInt k=0;k<Basis::NumBasis;k++){
                error[cid][5*k+0] = x[cid][5*k+0] - rho_coef[k];
                error[cid][5*k+1] = x[cid][5*k+1] - rhou_coef[k];
                error[cid][5*k+2] = x[cid][5*k+2] - rhov_coef[k];
                error[cid][5*k+3] = x[cid][5*k+3] - rhow_coef[k];
                error[cid][5*k+4] = x[cid][5*k+4] - rhoe_coef[k];
            }
        }

        for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
            const auto& cell = cmesh.m_cells[cid];
            for(uInt g=0; g<QuadC::num_points; ++g) {
                const auto& p = QuadC::points[g];
                const auto& pos = cell.transform_to_physical(p);
                fp <<std::setprecision(16)<< pos[0] << "  " <<std::setprecision(16)<< pos[1] << "  " <<std::setprecision(16)<< pos[2]
                 << "  " <<std::setprecision(16)<<  rho_h[cid][g] << "  " <<std::setprecision(16)<<  rho_s[cid][g] 
                 << "  " <<std::setprecision(16)<<  u_h[cid][g] << "  " <<std::setprecision(16)<<  u_s[cid][g] 
                 << "  " <<std::setprecision(16)<<  v_h[cid][g] << "  " <<std::setprecision(16)<<  v_s[cid][g] 
                 << "  " <<std::setprecision(16)<<  w_h[cid][g] << "  " <<std::setprecision(16)<<  w_s[cid][g] 
                 << "  " <<std::setprecision(16)<<  e_h[cid][g] << "  " <<std::setprecision(16)<<  e_s[cid][g] << std::endl;
            }
        }
        fp.close();

        print(std::array<Scalar,8>{curr_time, err_integral(rho_h,rho_s), 
                        err_integral(u_h,u_s), err_integral(v_h,v_s), err_integral(w_h,w_s), 
                        err_integral(e_h,e_s), std::sqrt(error.dot(error)/x.dot(x)), chrone_clock()});
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
    if(order == 5) Run<5>(meshN);
    // if(order == 6) Run<6>(meshN);
    // if(order == 7) Run<7>(meshN);
    // if(order == 8) Run<8>(meshN);
    // if(order == 9) Run<9>(meshN);
}