#include "DG/DG_Schemes/InCompressibleImplicitConvection.h"

#define Explicit_For_Flux(Order) \
template class InCompressibleImplicitConvection<Order+1,Order,AirFluxC>;\
template class InCompressibleImplicitConvection<Order+1,Order,MonatomicFluxC>;\
template class InCompressibleImplicitConvection<Order+2,Order+1,AirFluxC>;\
template class InCompressibleImplicitConvection<Order+2,Order+1,MonatomicFluxC>;\
template class InCompressibleImplicitConvection<Order+3,Order+2,AirFluxC>;\
template class InCompressibleImplicitConvection<Order+3,Order+2,MonatomicFluxC>;

Explicit_For_Flux(0)
Explicit_For_Flux(3)

#undef Explicit_For_Flux


template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_boundary_faces(const ComputingMesh& mesh, 
                      const LongVector<DoFs>& old_solution,
                      const Scalar curr_time,
                      BlockSparseMatrix<DoFs,DoFs>& sparse_mat,
                      LongVector<DoFs>& sparse_rhs) 
{
    /* 边界条件处理通用公式：
       Flux = 0.5[(F(u_in) + F(u_bc))·n + λ(u_in - u_bc)]
       λ = max(|u_in·n|, |u_bc·n|) (Lax-Friedrichs参数)
    */

    #pragma omp parallel for schedule(dynamic)
    for(uInt fid = 0; fid < mesh.m_faces.size(); ++fid) {
        const auto& face = mesh.m_faces[fid];
        if(face.m_neighbor_cells[1] != uInt(-1)) continue;

        const uInt cid = face.m_neighbor_cells[0];
        const auto& coef = old_solution[cid];
        const auto bc_type = mesh.m_boundaryTypes[fid];
        
        DenseMatrix<DoFs,DoFs> face_matrix;
        DenseMatrix<DoFs,1> face_rhs;

        for(uInt g = 0; g < QuadF::num_points; ++g) {
            const auto xi = transform_to_cell(face, QuadF::points[g], 0);
            const Scalar jac_weight = QuadF::weights[g] * face.compute_jacobian_det();
            
            // 重建内部状态 -------------------------------------------------
            /* u_in = Σ (c_u^j, c_v^j, c_w^j)ϕ_j(xi)
               p_in = Σ c_p^jψ_j(xi) */
            const auto u_basis = BasisU::eval_all(xi[0], xi[1], xi[2]);
            const auto p_basis = BasisP::eval_all(xi[0], xi[1], xi[2]);
            
            DenseMatrix<3,1> uvw_inner = DenseMatrix<3,1>::Zeros();
            Scalar p_inner = 0.0;
            for(uInt bid = 0; bid < NU; ++bid) {
                uvw_inner[0] += u_basis[bid] * coef[3*bid];
                uvw_inner[1] += u_basis[bid] * coef[3*bid+1];
                uvw_inner[2] += u_basis[bid] * coef[3*bid+2];
            }
            for(uInt bid = 0; bid < NP; ++bid) {
                p_inner += p_basis[bid] * coef[3*NU + bid];
            }

            // 边界条件处理 -------------------------------------------------
            const DenseMatrix<3,1> normal = face.m_normal;
            DenseMatrix<3,1> uvw_bc;
            Scalar p_bc = 0.0;

            switch(bc_type) {
            case BoundaryType::Dirichlet: {
                /* Dirichlet边界条件处理：
                   u_bc = u_analytic(xi,t)
                   p_bc = p_analytic(xi,t)
                   强加形式：u_ghost = 2u_bc - u_in */
                uvw_bc = uvw_Xi(mesh.m_cells[cid], xi, curr_time); 
                p_bc = p_Xi(mesh.m_cells[cid], xi, curr_time);
                // uvw_bc = 2*uvw_bc - uvw_inner;  // 镜像反射构造ghost cell
                // p_bc = 2*p_bc - p_inner;
                break;
            }
            case BoundaryType::Wall: {
                /* 无滑移壁面边界条件：
                   u_bc = u_wall (静止壁面时为0)
                   p_bc = p_in (Neumann条件) */
                uvw_bc = uvw_Xi(mesh.m_cells[cid], xi, curr_time);
                p_bc = p_inner;  // ∂p/∂n = 0
                break;
            }
            case BoundaryType::Pseudo3DZ: {
                /* 伪三维Z方向处理：
                   u_ghost = (u_in, v_in, -w_in)
                   p_ghost = p_in */
                uvw_bc = {uvw_inner[0], uvw_inner[1], -uvw_inner[2]};
                p_bc = p_inner;
                // uvw_bc = 2*uvw_inner - uvw_bc;  // 构造镜像值
                // p_bc = 2*p_inner - p_bc;
                break;
            }
            default:
                throw std::runtime_error("Unsupported boundary type");
            }

            // 通量计算 -----------------------------------------------------
            /* Lax-Friedrichs参数：λ = max(|u_in·n|, |u_bc·n|)
               平均状态：u_avg = 0.5(u_in + u_bc)
                        p_avg = 0.5(p_in + p_bc) */
            const Scalar lambda = std::max(uvw_inner.norm(), uvw_bc.norm());
            const DenseMatrix<3,1> uvw_avg = 0.5 * (uvw_inner + uvw_bc);
            const Scalar p_avg = 0.5 * (p_inner + p_bc);

            // 动量方程通量项 -------------------------------------------------
            for(uInt j = 0; j < NU; ++j) { // 测试函数ϕ_j
                /* 对流项贡献矩阵：A_uu[j,i] += 0.5[(u_avg·n)ϕ_iϕ_j + λϕ_iϕ_j]I_3
                   压力项贡献矩阵：A_up[j,i] += 0.5nψ_iϕ_j */
                const Scalar conv_coeff = 0.5 * (uvw_avg.dot(normal) + lambda);
                for(uInt i = 0; i < NU; ++i) { // 试探函数ϕ_i
                    assemble_Auu(face_matrix, 
                               conv_coeff * u_basis[i] * u_basis[j] * jac_weight,
                               j, i);
                }

                for(uInt i = 0; i < NP; ++i) { // 试探函数ψ_i
                    assemble_Aup(face_matrix,
                               0.5 * normal * p_basis[i] * u_basis[j] * jac_weight,
                               j, i);
                }

                /* 右端项贡献：rhs += 0.5[(u_bc·n - λ)u_bc + p_bc n]ϕ_j */
                const DenseMatrix<3,1> rhs_term = 
                    (0.5 * (uvw_avg.dot(normal) - lambda) * uvw_bc + 
                    0.5 * normal * p_bc);
                assemble_bu(face_rhs, rhs_term * u_basis[j] * jac_weight, j);
            }

            // 连续性方程通量项 -----------------------------------------------
            for(uInt j = 0; j < NP; ++j) { // 测试函数ψ_j
                /* 散度项贡献矩阵：A_pu[j,i] += 0.5n^Tϕ_iψ_j */
                for(uInt i = 0; i < NU; ++i) { // 试探函数ϕ_i
                    assemble_Apu(face_matrix,
                               0.5 * normal.transpose() * u_basis[i] * p_basis[j] * jac_weight,
                               j, i);
                }

                /* 右端项贡献：rhs += 0.5(u_bc·n)ψ_j */
                assemble_bp(face_rhs, 
                          0.5 * uvw_bc.dot(normal) * p_basis[j] * jac_weight, 
                          j);
            }
        }

        #pragma omp critical
        {
            sparse_mat.add_block(cid, cid, face_matrix);
            sparse_rhs[cid] -= face_rhs;
        }
    }
}
// template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
// void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
// assemble_boundarys(const ComputingMesh& mesh, const LongVector<DoFs>& old_solution,
//             const Scalar curr_time, BlockSparseMatrix<DoFs,DoFs>& sparse_mat,
//             LongVector<DoFs>& sparse_rhs){
//     // 只计算边界的
//     #pragma omp parallel for schedule(dynamic)
//     for(uInt fid=0; fid<mesh.m_faces.size(); ++fid) {
        
//         const auto& face = mesh.m_faces[fid];
//         const auto& cells = face.m_neighbor_cells;

//         // 跳过内部单元
//         if(cells[1] != uInt(-1)) continue;
//         // 只要边界单元

//         const auto bc = mesh.m_boundaryTypes[fid];

//         // 内部单元
//         const auto& coef = old_solution[cells[0]];
//         DenseMatrix<DoFs,DoFs> face_matrix;
//         DenseMatrix<DoFs,1> face_rhs;

//         for(uInt g=0; g<QuadF::num_points; ++g) {
//             // 转换到单元自然坐标
//             const auto& xi = transform_to_cell(face, QuadF::points[g], 0);
//             const auto& u_basis = BasisU::eval_all(xi[0], xi[1], xi[2]);
//             const auto& p_basis = BasisP::eval_all(xi[0], xi[1], xi[2]);
//             const Scalar jac_weight = QuadF::weights[g] * face.compute_jacobian_det();
            
//             // 重建内部状态
//             DenseMatrix<3,1> uvw_inner;
//             for(uInt bid=0; bid<NU; ++bid) {
//                 for(uInt k=0; k<3; ++k)
//                     uvw_inner[k] += u_basis[bid] * coef[3*bid + k];
//             }
//             Scalar p_inner = 0.0;
//             for(uInt bid=0; bid<NP; ++bid) {
//                 p_inner += p_basis[bid] * coef[3*NU + bid];
//             }
//             const DenseMatrix<3,1>& normal = face.m_normal;
//             // 根据边界类型计算通量贡献
//             DenseMatrix<3UL, 3UL> dF;
//             DenseMatrix<3UL, 1UL> drhs;
//             DenseMatrix<3UL, 1UL> dFp;
//             DenseMatrix<3UL, 1UL> drhsp;
//             DenseMatrix<1UL, 3UL> dFu;
//             DenseMatrix<1UL, 1UL> drhsu;
//             if (bc==BoundaryType::Dirichlet) {
//                 // 使用与内部面完全相同
//                 const auto& uvw_D = uvw_Xi(mesh.m_cells[cells[0]],xi,curr_time);
//                 const auto& p_D = p_Xi(mesh.m_cells[cells[0]],xi,curr_time);

//                 auto uvw_ghost = uvw_D + (uvw_D - uvw_inner);
//                 auto p_ghost = p_D + (p_D - p_inner);
//                 // F(uvw) = rho uvw_i delta_ij uvw_j
//                 // auto F_inner = Flux::computeJacobian(U_inner);
//                 // auto F_ghost = Flux::computeJacobian(U_ghost);
//                 // F_inner = F_ghost;
                
//                 // F(uvw) dot n = (u dot n) rho I
//                 Scalar Fn_inner = uvw_inner.dot(normal);
//                 Scalar Fn_ghost = uvw_ghost.dot(normal);
//                 // DenseMatrix<5,5> Fn_inner, Fn_ghost;
//                 // for(int d=0; d<3; ++d){
//                 //     Fn_inner += F_inner[d] * face.m_normal[d];
//                 //     Fn_ghost += F_ghost[d] * face.m_normal[d];
//                 // }
                
//                 const Scalar lambda = std::max(uvw_inner.length(),uvw_ghost.length());
//                 uvw_ghost = 0.5*(uvw_inner+uvw_ghost);
//                 p_ghost = 0.5*(p_inner+p_ghost);
//                 dF = 0.5*(Fn_inner + lambda)*DenseMatrix<3,3>::Identity();
//                 drhs = 0.5*(Fn_inner + lambda)*uvw_ghost;
//                 dFp = normal;
//                 drhsp = normal * p_ghost;
//                 dFu = normal.transpose();
//                 drhsu = dFu.multiply(uvw_ghost);
//             }
//             else if (bc==BoundaryType::Wall || bc==BoundaryType::WallTD
//                      || bc==BoundaryType::WallTN || bc==BoundaryType::WallTR) {
//                 // 使用与内部面完全相同
//                 const auto& uvw_D = uvw_Xi(mesh.m_cells[cells[0]],xi,curr_time);

//                 auto uvw_ghost = uvw_D + (uvw_D - uvw_inner);
//                 auto p_ghost = p_inner;
//                 uvw_ghost = 0.5*(uvw_inner+uvw_ghost);
//                 // F(uvw) = rho uvw_i delta_ij uvw_j
//                 // auto F_inner = Flux::computeJacobian(U_inner);
//                 // auto F_ghost = Flux::computeJacobian(U_ghost);
//                 // F_inner = F_ghost;
                
//                 // F(uvw) dot n = (u dot n) rho I
//                 Scalar Fn_inner = uvw_inner.dot(normal);
//                 Scalar Fn_ghost = uvw_ghost.dot(normal);
//                 // DenseMatrix<5,5> Fn_inner, Fn_ghost;
//                 // for(int d=0; d<3; ++d){
//                 //     Fn_inner += F_inner[d] * face.m_normal[d];
//                 //     Fn_ghost += F_ghost[d] * face.m_normal[d];
//                 // }
                
//                 const Scalar lambda = std::max(uvw_inner.length(),uvw_ghost.length());
//                 dF = 0.5*(Fn_inner + lambda)*DenseMatrix<3,3>::Identity();
//                 drhs = 0.5*(Fn_inner + lambda)*uvw_ghost;
//                 dFp = normal;
//                 drhsp = normal * p_ghost;
//                 dFu = normal.transpose();
//                 drhsu = dFu.multiply(uvw_ghost);
//             }
//             else if (bc==BoundaryType::Neumann) {
//                 // 使用与内部面完全相同
//                 auto uvw_ghost = uvw_inner;
//                 auto p_ghost = p_inner;
//                 // F(uvw) = rho uvw_i delta_ij uvw_j
//                 // auto F_inner = Flux::computeJacobian(U_inner);
//                 // auto F_ghost = Flux::computeJacobian(U_ghost);
//                 // F_inner = F_ghost;
                
//                 // F(uvw) dot n = (u dot n) rho I
//                 Scalar Fn_inner = uvw_inner.dot(normal);
//                 Scalar Fn_ghost = uvw_ghost.dot(normal);
//                 // DenseMatrix<5,5> Fn_inner, Fn_ghost;
//                 // for(int d=0; d<3; ++d){
//                 //     Fn_inner += F_inner[d] * face.m_normal[d];
//                 //     Fn_ghost += F_ghost[d] * face.m_normal[d];
//                 // }
                
//                 const Scalar lambda = std::max(uvw_inner.length(),uvw_ghost.length());
//                 dF = 0.5*(Fn_inner + lambda)*DenseMatrix<3,3>::Identity();
//                 drhs = 0.5*(Fn_inner + lambda)*uvw_ghost;
//                 dFp = normal;
//                 drhsp = normal * p_ghost;
//                 dFu = normal.transpose();
//                 drhsu = dFu.multiply(uvw_ghost);
//             }
//             else if(bc==BoundaryType::Pseudo3DZ){
//                 auto uvw_ghost = uvw_inner * DenseMatrix<3,1>({1,1,-1});
//                 auto p_ghost = p_inner;
//                 uvw_ghost = 0.5*(uvw_inner+uvw_ghost);
//                 p_ghost = 0.5*(p_inner+p_ghost);
//                 Scalar Fn_inner = uvw_inner.dot(normal);
//                 Scalar Fn_ghost = uvw_ghost.dot(normal);
//                 const Scalar lambda = std::max(uvw_inner.length(),uvw_ghost.length());
//                 dF = 0.5*(Fn_inner + lambda)*DenseMatrix<3,3>::Identity();
//                 drhs = 0.5*(Fn_inner + lambda)*uvw_ghost;
//                 dFp = normal;
//                 drhsp = normal * p_ghost;
//                 dFu = normal.transpose();
//                 drhsu = dFu.multiply(uvw_ghost);
//             }
//             else if(bc==BoundaryType::Pseudo3DY){
//                 auto uvw_ghost = uvw_inner * DenseMatrix<3,1>({1,-1,1});
//                 auto p_ghost = p_inner;
//                 uvw_ghost = 0.5*(uvw_inner+uvw_ghost);
//                 p_ghost = 0.5*(p_inner+p_ghost);
//                 Scalar Fn_inner = uvw_inner.dot(normal);
//                 Scalar Fn_ghost = uvw_ghost.dot(normal);
//                 const Scalar lambda = std::max(uvw_inner.length(),uvw_ghost.length());
//                 dF = 0.5*(Fn_inner + lambda)*DenseMatrix<3,3>::Identity();
//                 drhs = 0.5*(Fn_inner + lambda)*uvw_ghost;
//                 dFp = normal;
//                 drhsp = normal * p_ghost;
//                 dFu = normal.transpose();
//                 drhsu = dFu.multiply(uvw_ghost);
//             }
//             else if(bc==BoundaryType::Pseudo3DX){
//                 auto uvw_ghost = uvw_inner * DenseMatrix<3,1>({-1,1,1});
//                 auto p_ghost = p_inner;
//                 uvw_ghost = 0.5*(uvw_inner+uvw_ghost);
//                 p_ghost = 0.5*(p_inner+p_ghost);
//                 Scalar Fn_inner = uvw_inner.dot(normal);
//                 Scalar Fn_ghost = uvw_ghost.dot(normal);
//                 const Scalar lambda = std::max(uvw_inner.length(),uvw_ghost.length());
//                 dF = 0.5*(Fn_inner + lambda)*DenseMatrix<3,3>::Identity();
//                 drhs = 0.5*(Fn_inner + lambda)*uvw_ghost;
//                 dFp = normal;
//                 drhsp = normal * p_ghost;
//                 dFu = normal.transpose();
//                 drhsu = dFu.multiply(uvw_ghost);
//             }
            

//             // 组装到单元矩阵
//             for(uInt j=0; j<NU; ++j){
//                 {
//                     for(uInt i=0; i<NU; ++i){
//                         MatrixView<DoFs,DoFs,3,3> block(face_matrix,3*j,3*i);
//                         block += dF * (u_basis[i] * u_basis[j]) * jac_weight;
//                     }
//                     MatrixView<DoFs,1,3,1> block(face_rhs,3*j,0);
//                     block += drhs * (u_basis[j]) * jac_weight;
//                     }
//                     {
//                     for(uInt i=0; i<NP; ++i) {
//                         MatrixView<DoFs,DoFs,3,1> block(face_matrix, 3*j, 3*NU+i);
//                         block += dFp*(p_basis[i] * u_basis[j]) * jac_weight;
//                     }
//                     MatrixView<DoFs,1,3,1> block(face_rhs,3*j,0);
//                     block += drhsp * (u_basis[j]) * jac_weight;
//                     }
//             }
//             for(uInt j=0; j<NP; ++j) {
//                 for(uInt i=0; i<NU; ++i) {
//                     MatrixView<DoFs,DoFs,1,3> block(face_matrix, 3*NU+j, 3*i);
//                     block += dFu*(u_basis[i] * p_basis[j]) * jac_weight;
//                 }
//                 MatrixView<DoFs,1,1,1> block(face_rhs,3*NU+j,0);
//                 block += drhsu * (p_basis[j]) * jac_weight;
//             }
//         }
//         // debug("22222");
//         // debug(face_matrix);

//         #pragma omp critical
//         {
//         sparse_mat.add_block(cells[0], cells[0], face_matrix);
//         sparse_rhs[cells[0]] -= face_rhs;
//         // debug(vector2u{cells[0], cells[0]});
//         }
//     }

// }

// template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
// void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
// assemble_cells(const ComputingMesh& mesh, const LongVector<DoFs>& old_solution,
//             const Scalar curr_time, BlockSparseMatrix<DoFs,DoFs>& sparse_mat,
//             LongVector<DoFs>& sparse_rhs){
//     // 预计算体积积分点数据
//     constexpr uInt num_vol_points = QuadC::num_points;
//     std::array<std::array<Scalar,3>, num_vol_points> vol_points;
//     std::array<Scalar, num_vol_points> vol_weights;
//     std::array<std::array<Scalar, NU>, num_vol_points> u_basis;
//     std::array<std::array<vector3f, NU>, num_vol_points> u_grads;
//     std::array<std::array<Scalar, NP>, num_vol_points> p_basis;
//     std::array<std::array<vector3f, NP>, num_vol_points> p_grads;
//     for(uInt g=0; g<num_vol_points; ++g) {
//         vol_points[g] = QuadC::points[g];
//         vol_weights[g] = QuadC::weights[g];
//         u_basis[g] = BasisU::eval_all(vol_points[g][0],vol_points[g][1],vol_points[g][2]);
//         u_grads[g] = BasisU::grad_all(vol_points[g][0],vol_points[g][1],vol_points[g][2]);
//         p_basis[g] = BasisP::eval_all(vol_points[g][0],vol_points[g][1],vol_points[g][2]);
//         p_grads[g] = BasisP::grad_all(vol_points[g][0],vol_points[g][1],vol_points[g][2]);
//     }

//     #pragma omp parallel for schedule(dynamic)
//     for(uInt cid=0;cid<mesh.m_cells.size();cid++){
//         const auto& cell = mesh.m_cells[cid];
//         const auto& coef = old_solution[cid];

//         DenseMatrix<DoFs,DoFs> cell_matrix;

//         for(uInt g=0; g<num_vol_points; ++g) {
//             // Auu 部分
//             DenseMatrix<3,1> uvw = DenseMatrix<3,1>::Zeros();
//             for(uInt bid=0; bid<NU; ++bid) {
//                 const Scalar phi = u_basis[g][bid];
//                 for(uInt k=0; k<3; ++k) 
//                     uvw[k] += phi * coef[3*bid + k];
//             }
//             // 原本是需要计算 Jac(F,U) 的，但 F_i(U) = rho u_i u_j = (rho u_i delta_ij ) u_j
//             const Scalar jac_weight = vol_weights[g] * cell.compute_jacobian_det();

//             for(uInt j=0; j<NU; ++j) {
//                 const auto& J = cell.get_invJacMat();
//                 const auto& grad_u_j = J.multiply(u_grads[g][j]);
//                 Scalar mat = uvw[0]*grad_u_j[0] + uvw[1]*grad_u_j[1] + uvw[2]*grad_u_j[2];
//                 mat *= jac_weight;
//                 // 这次改了 for j  for i 顺序
//                 for(uInt i=0; i<NU; ++i) {
//                     MatrixView<DoFs,DoFs,3,3> block(cell_matrix,3*j,3*i);
//                     //  nabla(F(U))  ->  -F(U)*grad V
//                     block -= (mat * u_basis[g][i])*DenseMatrix<3,3>::Identity();
//                 }
//                 for(uInt i=0; i<NP; ++i) {
//                     MatrixView<DoFs,DoFs,3,1> block(cell_matrix,3*j,3*NU+i);
//                     // 压强在等号右边好像是 - grad(p) ？那挪过去就是+号
//                     block -= grad_u_j*jac_weight*p_basis[g][i];
//                 }
//             }
//             for(uInt j=0; j<NP; ++j) {
//                 const auto& J = cell.get_invJacMat();
//                 const auto& grad_p_j = J.multiply(p_grads[g][j]);
//                 for(uInt i=0; i<NU; ++i) {
//                     MatrixView<DoFs,DoFs,1,3> block(cell_matrix,3*NU+j,3*i);
//                     // 速度散度在左边取 + 
//                     block -= grad_p_j.transpose()*jac_weight*u_basis[g][i];
//                 }
//             }
//         }
//         #pragma omp critical
//         sparse_mat.add_block(cid, cid, cell_matrix);
//     }
// }

// template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
// void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
// assemble_internals(const ComputingMesh& mesh, const LongVector<DoFs>& old_solution,
//             const Scalar curr_time, BlockSparseMatrix<DoFs,DoFs>& sparse_mat,
//             LongVector<DoFs>& sparse_rhs){
//     // 只计算内部的
//     #pragma omp parallel for schedule(dynamic)
//     for(uInt fid=0; fid<mesh.m_faces.size(); ++fid) {
//         const auto& face = mesh.m_faces[fid];
//         const auto& cells = face.m_neighbor_cells;
        
//         // 跳过边界单元
//         if(cells[1] == uInt(-1)) continue;
//         // 后面就是只有内部面

//         // 获取左右单元系数
//         const auto& coef_L = old_solution[cells[0]];
//         const auto& coef_R = old_solution[cells[1]];

//         // 左右两个面的分块矩阵
//         DenseMatrix<DoFs,DoFs> face_matrix_LL, face_matrix_LR;
//         DenseMatrix<DoFs,DoFs> face_matrix_RL, face_matrix_RR;

//         for(uInt g=0; g<QuadF::num_points; ++g) {
//             const auto& uv = QuadF::points[g];
//             const Scalar jac_weight = QuadF::weights[g] * face.compute_jacobian_det();

//             // 转换到左右单元自然坐标
//             auto xi_L = transform_to_cell(face, uv, 0);
//             auto xi_R = transform_to_cell(face, uv, 1);

//             // 计算左右单元基函数值
//             auto u_basis_L = BasisU::eval_all(xi_L[0], xi_L[1], xi_L[2]);
//             auto u_basis_R = BasisU::eval_all(xi_R[0], xi_R[1], xi_R[2]);
//             auto p_basis_L = BasisP::eval_all(xi_L[0], xi_L[1], xi_L[2]);
//             auto p_basis_R = BasisP::eval_all(xi_R[0], xi_R[1], xi_R[2]);

//             // 重建左右状态
//             DenseMatrix<3,1> uvw_L, uvw_R;
//             for(uInt bid=0; bid<NU; ++bid) {
//                 for(uInt k=0; k<3; ++k) {
//                     uvw_L[k] += u_basis_L[bid] * coef_L[3*bid + k];
//                     uvw_R[k] += u_basis_R[bid] * coef_R[3*bid + k];
//                 }
//             }
//             Scalar p_L, p_R;
//             for(uInt bid=0; bid<NP; ++bid) {
//                 p_L += p_basis_L[bid] * coef_L[3*NU + bid];
//                 p_R += p_basis_R[bid] * coef_R[3*NU + bid];
//             }

//             // 计算通量雅可比
//             // auto F_L = Flux::computeJacobian(U_L);  // u_i * DenseMatrix<3,3>::Identity();
//             // auto F_R = Flux::computeJacobian(U_R);

//             // 计算法向通量雅可比
//             // DenseMatrix<3,3> Fn_L, Fn_R;
//             // for(int d=0; d<3; ++d) {
//             //     Fn_L += F_L[d] * face.m_normal[d];
//             //     Fn_R += F_R[d] * face.m_normal[d];
//             // }
//             // debug(Fn_L);
            

//             // ( u dot n ) * u，改写为 ((u dot n) * I) * u ，
//             const DenseMatrix<3,1>& normal = face.m_normal;
//             const auto& Fn_L = uvw_L.dot(normal) * DenseMatrix<3,3>::Identity();
//             const auto& Fn_R = uvw_R.dot(normal) * DenseMatrix<3,3>::Identity();

//             // Lax-Friedrichs参数
//             const Scalar lambda = std::max(uvw_L.length(),uvw_R.length());//Flux::computeWaveSpeed(U_L, U_R);
//             //  debug(lambda);
//             // 通量雅可比贡献
//             // 同理，是对角矩阵（数量矩阵）
//             const auto& dF_L = 0.5*(Fn_L + lambda * DenseMatrix<3,3>::Identity());
//             const auto& dF_R = 0.5*(Fn_R - lambda * DenseMatrix<3,3>::Identity());
//             // debug(dF_L);
//             // 组装到左右单元矩阵
//             for(uInt j=0; j<NU; ++j) {
//                 // u_basis[j] 为 测试函数
//                 for(uInt i=0; i<NU; ++i) {
//                     // u_basis[i] 为 u,v,w 的试探函数，也就是 u,v,w 基函数
//                     MatrixView<DoFs,DoFs,3,3> block_LL(face_matrix_LL, 3*j, 3*i);
//                     MatrixView<DoFs,DoFs,3,3> block_LR(face_matrix_LR, 3*j, 3*i);
//                     MatrixView<DoFs,DoFs,3,3> block_RL(face_matrix_RL, 3*j, 3*i);
//                     MatrixView<DoFs,DoFs,3,3> block_RR(face_matrix_RR, 3*j, 3*i);
//                     // ((u dot n + lambda) * u_basis[j] ) 为 u,v,w 共同的系数，所以乘以单位矩阵
//                     block_LL += (dF_L * (u_basis_L[i] * u_basis_L[j]) * jac_weight);
//                     block_LR += (dF_R * (u_basis_R[i] * u_basis_L[j]) * jac_weight);
//                     block_RL += (dF_L * (u_basis_L[i] * u_basis_L[j]) * jac_weight);
//                     block_RR += (dF_R * (u_basis_R[i] * u_basis_L[j]) * jac_weight);
//                 }
//                 for(uInt i=0; i<NP; ++i) {
//                     // p_basis[i] 为 p 的试探函数，也就是 p 基函数
//                     MatrixView<DoFs,DoFs,3,1> block_LL(face_matrix_LL, 3*j, 3*NU+i);
//                     MatrixView<DoFs,DoFs,3,1> block_LR(face_matrix_LR, 3*j, 3*NU+i);
//                     MatrixView<DoFs,DoFs,3,1> block_RL(face_matrix_RL, 3*j, 3*NU+i);
//                     MatrixView<DoFs,DoFs,3,1> block_RR(face_matrix_RR, 3*j, 3*NU+i);
//                     block_LL += ((p_basis_L[i] * u_basis_L[j]) * jac_weight) * normal;
//                     block_LR += ((p_basis_R[i] * u_basis_L[j]) * jac_weight) * normal;
//                     block_RL += ((p_basis_L[i] * u_basis_L[j]) * jac_weight) * normal;
//                     block_RR += ((p_basis_R[i] * u_basis_L[j]) * jac_weight) * normal;
//                 }
//             }
//             for(uInt j=0; j<NP; ++j) {
//                 // 测试函数为 p_basis[j] ，也就是 p 的基函数空间
//                 for(uInt i=0; i<NU; ++i) {
//                     // 对 div(u) * p_basis[j] 进行积分，高斯公式翻转为面积分
//                     MatrixView<DoFs,DoFs,1,3> block_LL(face_matrix_LL, 3*NU+j, 3*i);
//                     MatrixView<DoFs,DoFs,1,3> block_LR(face_matrix_LR, 3*NU+j, 3*i);
//                     MatrixView<DoFs,DoFs,1,3> block_RL(face_matrix_RL, 3*NU+j, 3*i);
//                     MatrixView<DoFs,DoFs,1,3> block_RR(face_matrix_RR, 3*NU+j, 3*i);
//                     block_LL += ((u_basis_L[i] * p_basis_L[j]) * jac_weight) * normal.transpose();
//                     block_LR += ((u_basis_R[i] * p_basis_L[j]) * jac_weight) * normal.transpose();
//                     block_RL += ((u_basis_L[i] * p_basis_L[j]) * jac_weight) * normal.transpose();
//                     block_RR += ((u_basis_R[i] * p_basis_L[j]) * jac_weight) * normal.transpose();
//                 }
//             }

//         }
//         // 这四个矩阵。。。。。终于。。。。。。
//         #pragma omp critical
//         {   
//             sparse_mat.add_block(cells[0], cells[0], face_matrix_LL);
//             sparse_mat.add_block(cells[0], cells[1], face_matrix_LR);
//             sparse_mat.add_block(cells[1], cells[0], (-1.)*face_matrix_RL); 
//             sparse_mat.add_block(cells[1], cells[1], (-1.)*face_matrix_RR);
//             // debug(vector2u{cells[0], cells[1]});
//         }
//     }
// }

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_cells(const ComputingMesh& mesh, 
              const LongVector<DoFs>& old_solution,
              BlockSparseMatrix<DoFs,DoFs>& sparse_mat) 
{
    /* 单元局部矩阵计算公式：
       A_uu[j,i] = -∫_K (u·∇ϕ_j)ϕ_i I_3 dx       (3×3块)
       A_up[j,i] = -∫_K ∇ϕ_j ψ_i dx             (3×1块) 
       A_pu[j,i] = -∫_K ϕ_i ∇ψ_j dx             (1×3块) */
    
    constexpr uInt num_vol_points = QuadC::num_points;
    // 预计算积分点数据
    std::array<vector3f, num_vol_points> vol_points;
    std::array<Scalar, num_vol_points> vol_weights;
    std::array<std::array<Scalar, NU>, num_vol_points> u_basis;  // ϕ值
    std::array<std::array<vector3f, NU>, num_vol_points> u_grads;// ∇ϕ
    std::array<std::array<Scalar, NP>, num_vol_points> p_basis;  // ψ值
    std::array<std::array<vector3f, NP>, num_vol_points> p_grads;// ∇ψ

    // 初始化积分点数据
    for(uInt g = 0; g < num_vol_points; ++g) {
        // 自然坐标到物理坐标转换
        vol_points[g] = QuadC::points[g];
        vol_weights[g] = QuadC::weights[g];
        u_basis[g] = BasisU::eval_all(vol_points[g][0], vol_points[g][1], vol_points[g][2]);
        u_grads[g] = BasisU::grad_all(vol_points[g][0], vol_points[g][1], vol_points[g][2]);
        p_basis[g] = BasisP::eval_all(vol_points[g][0], vol_points[g][1], vol_points[g][2]);
        p_grads[g] = BasisP::grad_all(vol_points[g][0], vol_points[g][1], vol_points[g][2]);
    }

    #pragma omp parallel for schedule(dynamic)
    for(uInt cid = 0; cid < mesh.m_cells.size(); ++cid) {
        const auto& cell = mesh.m_cells[cid];
        const auto& coef = old_solution[cid];
        DenseMatrix<DoFs, DoFs> cell_matrix;
        const Scalar jac_det = cell.compute_jacobian_det();  // |J|
        const auto& invJ = cell.get_invJacMat();             // J^{-T} 已经内蕴了一个转置

        // 重建当前速度和压力
        DenseMatrix<3,1> uvw = DenseMatrix<3,1>::Zeros();
        // DenseMatrix<3,1> grad_u = DenseMatrix<3,1>::Zeros();
        // DenseMatrix<3,1> grad_v = DenseMatrix<3,1>::Zeros();
        // DenseMatrix<3,1> grad_w = DenseMatrix<3,1>::Zeros();
        // DenseMatrix<3,1> grad_p = DenseMatrix<3,1>::Zeros();
        for(uInt g = 0; g < num_vol_points; ++g) {
            const Scalar jac_weight = vol_weights[g] * jac_det; // w_g |J|
            
            // 计算当前积分点物理量
            uvw = DenseMatrix<3,1>::Zeros();
            // grad_u = DenseMatrix<3,1>::Zeros();
            // grad_v = DenseMatrix<3,1>::Zeros();
            // grad_w = DenseMatrix<3,1>::Zeros();
            // grad_p = DenseMatrix<3,1>::Zeros();
            Scalar p = 0.0;
            for(uInt bid = 0; bid < NU; ++bid) {
                const Scalar phi = u_basis[g][bid];
                const auto grad_phi_j = invJ.multiply(u_grads[g][bid]); // ∇ϕ_j = J^{-T}∇̂ϕ_j
                uvw[0] += phi * coef[3*bid];     // u = Σ c_u^i ϕ_i
                uvw[1] += phi * coef[3*bid+1];    // v = Σ c_v^i ϕ_i
                uvw[2] += phi * coef[3*bid+2];    // w = Σ c_w^i ϕ_i
                // grad_u += grad_phi_j * coef[3*bid]; // ∇u = Σ c_u^i ∇ϕ_i
                // grad_v += grad_phi_j * coef[3*bid+1]; // ∇v = Σ c_v^i ∇ϕ_i
                // grad_w += grad_phi_j * coef[3*bid+2]; // ∇w = Σ c_w^i ∇ϕ_i
            }
            for(uInt bid = 0; bid < NP; ++bid) {
                p += p_basis[g][bid] * coef[3*NU + bid]; // p = Σ c_p^i ψ_i
                // grad_p += p_grads[g][bid] * coef[3*NU + bid]; // ∇p = Σ c_p^i ∇ψ_i
            }
            // debug(cell.m_centroid);
            // debug(grad_u.transpose());
            // debug(grad_v.transpose());
            // debug(grad_w.transpose());
            // debug(grad_p.transpose());

            // debug(uvw.transpose());
            // debug(p);

            // 组装动量方程
            for(uInt j = 0; j < NU; ++j) { // 试探函数ϕ_j
                const auto grad_phi_j = invJ.multiply(u_grads[g][j]); // ∇ϕ_j = J^{-T}∇̂ϕ_j
                
                /* 对流项：A_uu[j,i] += -∫(u·∇ϕ_j)ϕ_i I_3 dx
                   = -(u·∇ϕ_j)ϕ_i |J| w_g I_3 */
                const Scalar u_grad_phi = uvw.dot(grad_phi_j);
                for(uInt i = 0; i < NU; ++i) { // 测试函数ϕ_i
                    assemble_Auu(cell_matrix, 
                               -u_basis[g][i] * u_grad_phi * jac_weight,
                               j, i); // 注意i,j顺序
                }

                /* 压力梯度项：A_up[j,i] += -∫∇ϕ_j ψ_i dx 
                   = -∇ϕ_j ψ_i |J| w_g */
                for(uInt i = 0; i < NP; ++i) { // 测试函数ψ_i
                    assemble_Aup(cell_matrix,
                               -p_basis[g][i] * grad_phi_j * jac_weight,
                               j, i);
                }
            }

            // 组装连续性方程
            for(uInt j = 0; j < NP; ++j) { // 试探函数ψ_j
                /* 散度项：A_pu[j,i] += -∫ϕ_i ∇ψ_j dx 
                   = -ϕ_i ∇ψ_j |J| w_g */
                const auto grad_phi_j = invJ.multiply(p_grads[g][j]); // ∇ϕ_j = J^{-T}∇ψ_j
                for(uInt i = 0; i < NU; ++i) { // 测试函数ϕ_i
                    assemble_Apu(cell_matrix,
                               -u_basis[g][i] * grad_phi_j.transpose() * jac_weight,
                               j, i);
                }
            }
            
        
        }
        #pragma omp critical
        sparse_mat.add_block(cid, cid, cell_matrix);
    }
}

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_internal_faces(const ComputingMesh& mesh,
                       const LongVector<DoFs>& old_solution,
                       BlockSparseMatrix<DoFs,DoFs>& sparse_mat) 
{
    /* 内部面通量项公式：
       A_uu贡献：0.5[(u·n ± λ)ϕ_i^± ϕ_j^∓]I_3 
       A_up贡献：0.5ψ_i^± ϕ_j^∓ n
       A_pu贡献：-0.5ϕ_i^± ψ_j^∓ n^T */
    
    #pragma omp parallel for schedule(dynamic)
    for(uInt fid = 0; fid < mesh.m_faces.size(); ++fid) {
        const auto& face = mesh.m_faces[fid];
        if(face.m_neighbor_cells[1] == uInt(-1)) continue;

        const uInt L = face.m_neighbor_cells[0]; // 左单元
        const uInt R = face.m_neighbor_cells[1]; // 右单元
        const auto& coef_L = old_solution[L];
        const auto& coef_R = old_solution[R];
        const DenseMatrix<3,1>& normal = face.m_normal;      // 单位法向量n (L→R)

        DenseMatrix<DoFs,DoFs> mat_LL, mat_LR, mat_RL, mat_RR;

        for(uInt g = 0; g < QuadF::num_points; ++g) {
            const auto& uv = QuadF::points[g];
            const Scalar jac_weight = QuadF::weights[g] * face.compute_jacobian_det();

            // 转换到单元局部坐标
            auto xi_L = transform_to_cell(face, uv, 0);
            auto xi_R = transform_to_cell(face, uv, 1);

            // 计算基函数值
            auto phi_L = BasisU::eval_all(xi_L[0], xi_L[1], xi_L[2]);
            auto phi_R = BasisU::eval_all(xi_R[0], xi_R[1], xi_R[2]);
            auto psi_L = BasisP::eval_all(xi_L[0], xi_L[1], xi_L[2]);
            auto psi_R = BasisP::eval_all(xi_R[0], xi_R[1], xi_R[2]);

            // 重建面两侧状态
            DenseMatrix<3,1> uvw_L, uvw_R;
            uvw_L = DenseMatrix<3,1>::Zeros();
            uvw_R = DenseMatrix<3,1>::Zeros();
            for(uInt bid = 0; bid < NU; ++bid) {
                uvw_L[0] += phi_L[bid] * coef_L[3*bid];
                uvw_L[1] += phi_L[bid] * coef_L[3*bid+1];
                uvw_L[2] += phi_L[bid] * coef_L[3*bid+2];
                
                uvw_R[0] += phi_R[bid] * coef_R[3*bid];
                uvw_R[1] += phi_R[bid] * coef_R[3*bid+1];
                uvw_R[2] += phi_R[bid] * coef_R[3*bid+2];
            }
            Scalar p_L = 0.0, p_R = 0.0;
            for(uInt bid = 0; bid < NP; ++bid) {
                p_L += psi_L[bid] * coef_L[3*NU + bid];
                p_R += psi_R[bid] * coef_R[3*NU + bid];
            }

            // debug(uvw_L.transpose());
            // debug(uvw_R.transpose());
            // debug(p_L);
            // debug(p_R);


            // 计算Lax-Friedrichs参数
            const Scalar lambda = std::max(uvw_L.norm(), uvw_R.norm());

            // 组装通量项
            for(uInt j = 0; j < NU; ++j) { // 测试函数ϕ_j
                for(uInt i = 0; i < NU; ++i) { // 试探函数ϕ_i
                    /* 对流通量项：
                       LL块：0.5(u_L·n + λ)ϕ_i^L ϕ_j^L I_3
                       LR块：0.5(u_R·n - λ)ϕ_i^R ϕ_j^L I_3 */
                    const Scalar flux_LL = 0.5 * (uvw_L.dot(normal) + lambda) 
                                         * phi_L[i] * phi_L[j] * jac_weight;
                    const Scalar flux_LR = 0.5 * (uvw_R.dot(normal) - lambda) 
                                         * phi_R[i] * phi_L[j] * jac_weight;
                    const Scalar flux_RL = 0.5 * (-uvw_L.dot(normal) - lambda) 
                                         * phi_L[i] * phi_R[j] * jac_weight;
                    const Scalar flux_RR = 0.5 * (-uvw_R.dot(normal) + lambda) 
                                         * phi_R[i] * phi_R[j] * jac_weight;

                    assemble_Auu(mat_LL, flux_LL, j, i);
                    assemble_Auu(mat_LR, flux_LR, j, i);
                    assemble_Auu(mat_RL, flux_RL, j, i);
                    assemble_Auu(mat_RR, flux_RR, j, i);
                }

                // 压力通量项
                for(uInt i = 0; i < NP; ++i) { // 试探函数ψ_i
                    /* 压力通量项：
                       LL块：0.5ψ_i^L ϕ_j^L n
                       LR块：0.5ψ_i^R ϕ_j^L n */
                    const DenseMatrix<3,1> pflux_LL = 0.5 * normal * psi_L[i] * phi_L[j] * jac_weight;
                    const DenseMatrix<3,1> pflux_LR = 0.5 * normal * psi_R[i] * phi_L[j] * jac_weight;
                    const DenseMatrix<3,1> pflux_RL = 0.5 * -normal * psi_L[i] * phi_R[j] * jac_weight;
                    const DenseMatrix<3,1> pflux_RR = 0.5 * -normal * psi_R[i] * phi_R[j] * jac_weight;

                    assemble_Aup(mat_LL, pflux_LL, j, i);
                    assemble_Aup(mat_LR, pflux_LR, j, i);
                    assemble_Aup(mat_RL, pflux_RL, j, i);
                    assemble_Aup(mat_RR, pflux_RR, j, i);
                }
            }

            // 连续性方程通量
            for(uInt j = 0; j < NP; ++j) { // 测试函数ψ_j
                for(uInt i = 0; i < NU; ++i) { // 试探函数ϕ_i
                    /* 散度通量项：
                       LL块：-0.5ϕ_i^L ψ_j^L n^T
                       LR块：-0.5ϕ_i^R ψ_j^L n^T */
                    const DenseMatrix<1,3> div_LL = 0.5 * normal.transpose() * phi_L[i] * psi_L[j] * jac_weight;
                    const DenseMatrix<1,3> div_LR = 0.5 * normal.transpose() * phi_R[i] * psi_L[j] * jac_weight;
                    const DenseMatrix<1,3> div_RL = 0.5 * -normal.transpose() * phi_L[i] * psi_R[j] * jac_weight;
                    const DenseMatrix<1,3> div_RR = 0.5 * -normal.transpose() * phi_R[i] * psi_R[j] * jac_weight;

                    assemble_Apu(mat_LL, div_LL, j, i);
                    assemble_Apu(mat_LR, div_LR, j, i);
                    assemble_Apu(mat_RL, div_RL, j, i);
                    assemble_Apu(mat_RR, div_RR, j, i);
                }
            }
        }

        #pragma omp critical
        {
            sparse_mat.add_block(L, L, mat_LL);
            sparse_mat.add_block(L, R, mat_LR);
            sparse_mat.add_block(R, L, mat_RL);
            sparse_mat.add_block(R, R, mat_RR);
        }
    }
}


template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble(const ComputingMesh& mesh, const LongVector<DoFs>& old_solution,
            const Scalar curr_time, BlockSparseMatrix<DoFs,DoFs>& sparse_mat,
            LongVector<DoFs>& sparse_rhs){
    // print("123");
    assemble_cells(mesh,old_solution,sparse_mat);
    // print("234");
    assemble_internal_faces(mesh,old_solution,sparse_mat);
    // print("345");
    assemble_boundary_faces(mesh,old_solution,curr_time,sparse_mat,sparse_rhs);
    // print("456");
}

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
vector3f InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
transform_to_cell(const CompTriangleFace& face, 
                              const vector2f& uv, uInt side) const {
    const auto& nc = face.m_natural_coords[side];
    return nc[0]*(1-uv[0]-uv[1]) + nc[1]*uv[0] + nc[2]*uv[1];
}

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_Auu(DenseMatrix<DoFs,DoFs>& block, 
                const DenseMatrix<3,3>& Auu,
                uInt row, uInt col) {
    // 本质上就是 3*row 到 3*row +2 ，3*col 到 3*col + 2，填入Auu
    block.template View<3,3>(3 * row, 3 * col) += Auu;
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_Auu(DenseMatrix<DoFs,DoFs>& block, 
                const Scalar& Auu,
                uInt row, uInt col) {
    // 数量矩阵 的 特化版本
    block.template View<3,3>(3 * row, 3 * col) += Auu * DenseMatrix<3,3>::Identity();
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_Auu(DenseMatrix<DoFs,DoFs>& block, 
                const DenseMatrix<3,1>& Auu,
                uInt row, uInt col) {
    // 对角矩阵 的 特化版本
    block.template View<3,3>(3 * row, 3 * col) += DenseMatrix<3,3>::Diag(Auu);
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_Aup(DenseMatrix<DoFs,DoFs>& block, 
                const DenseMatrix<3,1>& Aup,
                uInt row, uInt col) {
    block.template View<3,1>(3 * row, 3*NU + col) += Aup; 
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_Apu(DenseMatrix<DoFs,DoFs>& block, 
                const DenseMatrix<1,3>& Apu,
                uInt row, uInt col) {
    block.template View<1,3>(3*NU + row, 3 * col) += Apu;
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_App(DenseMatrix<DoFs,DoFs>& block, 
                const DenseMatrix<1,1>& App,
                uInt row, uInt col) {
    block.template View<1,1>(3*NU + row,3*NU + col) += App;
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_App(DenseMatrix<DoFs,DoFs>& block, 
                const Scalar& App,
                uInt row, uInt col) {
    // 标量特化版本
    block(3*NU + row,3*NU + col) += App;
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_bu(DenseMatrix<DoFs,1>& block_rhs, 
                const DenseMatrix<3,1>& bu, uInt row) {
    // 本质上就是 3*row 到 3*row +2 ，0 列，填入bu
    block_rhs.template View<3,1>(3 * row, 0) += bu;
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_bp(DenseMatrix<DoFs,1>& block_rhs, 
                const DenseMatrix<1,1>& bu, uInt row) {
    block_rhs.template View<1,1>(3*NU + row, 0) += bu;
};

template<uInt OrderU, uInt OrderP, typename Flux, typename GaussQuadCell, typename GaussQuadFace>
void InCompressibleImplicitConvection<OrderU, OrderP, Flux, GaussQuadCell, GaussQuadFace>::
assemble_bp(DenseMatrix<DoFs,1>& block_rhs, 
                const Scalar& bu, uInt row) {
    // 标量特化版本
    block_rhs(3*NU + row, 0) += bu;
};