#include <fstream>
#include <sstream>
#include <chrono>
#include "Type.h"
#include "OptimizedMesh/OptimizedMesh.h"
#include "Matrix/Matrix.h"
#include "DG_Schemes/DG_Basis.h"
#include "DG_Schemes/Convection.h"

#include "exact.h"

constexpr Scalar param_gamma = 1.4;

template<typename Type>
inline Type rho_xyz(Type x, Type y, Type z, Type t){
    return 1.0;
}
template<typename Type>
inline Type u_xyz(Type x, Type y, Type z, Type t){
    return 0.0;
}
template<typename Type>
inline Type v_xyz(Type x, Type y, Type z, Type t){
    return 0.0;
}
template<typename Type>
inline Type w_xyz(Type x, Type y, Type z, Type t){
    return 0.0;
}
template<typename Type>
inline Type p_xyz(Type x, Type y, Type z, Type t){
    return (param_gamma-1)*rho_xyz(x,y,z,t)*e_xyz(x,y,z);
}
template<typename Type>
inline Type e_xyz(Type x, Type y, Type z, Type t){
    Scalar r2 = x*x + y*y;
    constexpr Scalar r_ds2 = 1.0Q/(2.0Q * 0.1Q * 0.1Q);
    return 1e-12 + 0.979264*M_1_PI*r_ds2* std::exp(-r2*r_ds2);
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
    GeneralMesh mesh = OrthHexMesh({-1.1, -1.1, -1.1/N},{1.1, 1.1, 1.1/N},{N,N,1});
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
    std::ofstream fp;
    fp.open("./Sedov/Order_" + std::to_string(Basis::OrderBasis) + "/rho_T_" + std::to_string(0) 
            + "_N_" + std::to_string(N) + ".txt");
    fp <<std::setprecision(16)<< ("#       x") << "  " 
        <<std::setprecision(16)<< (" y") << "  " 
        <<std::setprecision(16)<< (" z")
                << "  " <<std::setprecision(16)<<  (" rho")
                << "  " <<std::setprecision(16)<<  (" u")
                << "  " <<std::setprecision(16)<<  (" v")
                << "  " <<std::setprecision(16)<<  (" w")
                << "  " <<std::setprecision(16)<<  (" e") << std::endl;
    for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
        const auto& cell = cmesh.m_cells[cid];
        for(uInt g=0; g<QuadC::num_points; ++g) {
            const auto& p = QuadC::points[g];
            const auto& pos = cell.transform_to_physical(p);
            fp <<std::setprecision(16)<< pos[0] << "  " <<std::setprecision(16)<< pos[1] << "  " <<std::setprecision(16)<< pos[2]
                << "  " <<std::setprecision(16)<<  rho_xyz(pos)
                << "  " <<std::setprecision(16)<<  u_xyz(pos)
                << "  " <<std::setprecision(16)<<  v_xyz(pos)
                << "  " <<std::setprecision(16)<<  w_xyz(pos)
                << "  " <<std::setprecision(16)<<  e_xyz(pos) << std::endl;
        }
    }


    ExplicitConvection<Basis::OrderBasis> convection;
    
    LongVector<5*Basis::NumBasis> b(x.size());
    for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
        DenseMatrix<5*Basis::NumBasis,5*Basis::NumBasis> mass_matrix;
        for(uInt g=0; g<QuadC::num_points; ++g) {
            const auto& p = QuadC::points[g];
            auto phi = Basis::eval_all(p[0], p[1], p[2]);
            for(uInt i=0; i<Basis::NumBasis; ++i) {
                for(uInt k=0; k<5; ++k) {
                    b[cid](5*i + k, 0) += phi[i] * phi[i] * 
                        QuadC::weights[g] * cmesh.m_cells[cid].compute_jacobian_det();
                }
            }
        }
    }

    // auto chrono_start = std::chrono::steady_clock::now();
    // auto chrone_clock = [&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-chrono_start).count();};
    // print(std::array<std::string,8>{"#       time", "rel.err  rho",
    //                 "rel.err  u", "rel.err  v", "rel.err  w", 
    //                 "rel.err  e", "rel.err coef", "cpu time"});
    for(uInt kkkk=0;kkkk<100;kkkk++){    
        // Scalar dt = 0.001;
        
        // logging("Iter    " + std::to_string(kkkk+1));
        
        
        Scalar sub_dt = 0.01;
        // debug(dt);
        Scalar sub_t = 0;
        do{

            Scalar h = 2.2/N * 0.5773502691896257645091487805019574556476 / 3.0,   vel = 0.0;
            for(uInt fid=0; fid<cmesh.m_faces.size(); ++fid) {
                // std::array<Scalar,5> cumsum = {0,0,0,0,0};
                const auto& face = cmesh.m_faces[fid];
                const auto& cells = face.m_neighbor_cells;
                if(cells[1] == uInt(-1)) continue;
                const auto& coef_L = x[cells[0]];
                const auto& coef_R = x[cells[1]];
                auto xi_L = convection.transform_to_cell(face, {1.0/3,1.0/3}, 0);
                auto xi_R = convection.transform_to_cell(face, {1.0/3,1.0/3}, 1);
                auto basis_L = Basis::eval_all(xi_L[0], xi_L[1], xi_L[2]);
                auto basis_R = Basis::eval_all(xi_R[0], xi_R[1], xi_R[2]);

                DenseMatrix<5,1> U_L, U_R;
                for(uInt bid=0; bid<Basis::NumBasis; ++bid) {
                    for(uInt k=0; k<5; ++k) {
                        U_L[k] += basis_L[bid] * coef_L[5*bid + k];
                        U_R[k] += basis_R[bid] * coef_R[5*bid + k];
                    }
                }
                const Scalar lambda = convection.compute_max_wave_speed(U_L, U_R);
                vel = std::max(lambda,vel);
            }


            Scalar max_dt = std::min(sub_dt * 0.25,   0.5 * h * 1.0/(2*Order+1) * 1.0 / vel);
            if(max_dt<1e-6)break;
            Scalar dt = std::min(max_dt,sub_dt-sub_t);
            sub_t += dt;

            // debug(vector2f{dt,sub_t});
            logging("Iter  " + std::to_string(kkkk+1) + " \tSub TimeStep \t" + std::to_string(sub_t));
            const auto& k1 = - 1/b * convection.eval(cmesh, x, 
                                                            0.1 * (kkkk) + sub_t);
            auto x1 = x + 1.0/2.0 * dt * k1;
            for(uInt cid=0;cid<x1.size();cid++){
                x1[cid][5*0+0] = std::max(1e-4,x1[cid][5*0+0]);
            }
            const auto& k2 = - 1/b * convection.eval(cmesh, x1,  
                                                            0.1 * (kkkk) + sub_t + 1.0/2.0 * dt); 
            x = x + dt * k2;
            for(uInt cid=0;cid<x1.size();cid++){
                x[cid][5*0+0] = std::max(1e-4,x[cid][5*0+0]);
            }

            // const auto& h1 = - 1/b * convection.eval(cmesh, x);
            // const auto& h2 = - 1/b * convection.eval(cmesh, x + 0.5 * dt * h1);
            // const auto& h3 = - 1/b * convection.eval(cmesh, x + 0.5 * dt * h2);
            // const auto& h4 = - 1/b * convection.eval(cmesh, x + 1.0 * dt * h3);
            // x = x + dt * (h1 + 2*h2 + 2*h3 + h4) / 6;

        }while(sub_t < sub_dt);


        Scalar curr_time = sub_dt * (kkkk+1);

        std::ofstream fp;
        fp.open("./Sedov/Order_" + std::to_string(Basis::OrderBasis) + "/rho_T_" + std::to_string(kkkk+1) 
                + "_N_" + std::to_string(N) + ".txt");
        fp <<std::setprecision(16)<< ("#       x") << "  " 
            <<std::setprecision(16)<< (" y") << "  " 
            <<std::setprecision(16)<< (" z")
                 << "  " <<std::setprecision(16)<<  (" rho")
                 << "  " <<std::setprecision(16)<<  (" u")
                 << "  " <<std::setprecision(16)<<  (" v")
                 << "  " <<std::setprecision(16)<<  (" w")
                 << "  " <<std::setprecision(16)<<  (" e") << std::endl;
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
        LongVector<QuadC::num_points> u_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> v_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> w_h(cmesh.m_cells.size());
        LongVector<QuadC::num_points> e_h(cmesh.m_cells.size());
        // LongVector<5 * Basis::NumBasis> error(cmesh.m_cells.size());

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
            }
        }

        for(uInt cid=0;cid<cmesh.m_cells.size();cid++){
            const auto& cell = cmesh.m_cells[cid];
            for(uInt g=0; g<QuadC::num_points; ++g) {
                const auto& p = QuadC::points[g];
                const auto& pos = cell.transform_to_physical(p);
                fp <<std::setprecision(16)<< pos[0] << "  " <<std::setprecision(16)<< pos[1] << "  " <<std::setprecision(16)<< pos[2]
                 << "  " <<std::setprecision(16)<<  rho_h[cid][g]
                 << "  " <<std::setprecision(16)<<  u_h[cid][g] 
                 << "  " <<std::setprecision(16)<<  v_h[cid][g] 
                 << "  " <<std::setprecision(16)<<  w_h[cid][g]
                 << "  " <<std::setprecision(16)<<  e_h[cid][g] << std::endl;
            }
        }
        fp.close();

        print(std::array<Scalar,2>{curr_time, chrone_clock()});
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