#pragma once
#include "Matrix/DenseMatrix.h"
#include "Matrix/SparseMatrix.h"
#include "Matrix/LongVector.h"
#include "Matrix/submat.h"
#include "DG/DG_Basis/DG_Basis.h"
#include "Mesh/ComputingMesh.h"
#include "Matrix/block_precondition.h"
#include "EigenSolver/EigenSparseSolver.h"


// 模板类实现 PGMRES
template<uInt BlockDim, uInt Max_dims, bool output_flag>
class PGMRES;


template<uInt BlockDim, uInt InnerMaxIters = 100>
class FVEPreconditioner : public Preconditioner<BlockDim> {
public:
    using Vec = LongVector<BlockDim>;
    static constexpr uInt NumBasis = BlockDim / 5;
    static constexpr uInt Order = (NumBasis < 10 ? (NumBasis == 1 ? 0 : 1) : (NumBasis==10 ? 2 : 3));
// private:
    const BlockSparseMatrix<BlockDim, BlockDim>& A;
    const ComputingMesh& mesh;

    BlockSparseMatrix<5,5*NumBasis> Rmat;
    BlockSparseMatrix<5,5> RAP;
    EigenSparseSolver<5,5> RAP_solver;

public:

    FVEPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs, const ComputingMesh& mesh_)
        : A(A_), mesh(mesh_)
    {
        // Ur_tmp.resize(DoFs);
        // Um_tmp.resize(DoFs);
        // Ue_tmp.resize(DoFs);
        // Arr_DiagVec = Arr.get_inv_diag();
        // Amm_DiagVec = Amm.get_inv_diag();
        // Aee_DiagVec = Aee.get_inv_diag();

        const std::array<std::array<Scalar,NumBasis>,4> phi_N{
            DGBasisEvaluator<Order>::eval_all(static_cast<Scalar>(0),static_cast<Scalar>(0),static_cast<Scalar>(0)),
            DGBasisEvaluator<Order>::eval_all(static_cast<Scalar>(1),static_cast<Scalar>(0),static_cast<Scalar>(0)),
            DGBasisEvaluator<Order>::eval_all(static_cast<Scalar>(0),static_cast<Scalar>(1),static_cast<Scalar>(0)),
            DGBasisEvaluator<Order>::eval_all(static_cast<Scalar>(0),static_cast<Scalar>(0),static_cast<Scalar>(1))
        };

        std::vector<Scalar> weights(mesh.m_points.size());  // 总节点数，记录权重，用于归一化
        for(uInt cellId = 0; cellId < mesh.m_cells.size(); ++cellId){
            const auto& cell = mesh.m_cells[cellId];
            const auto& nodes = cell.m_nodes;  // 四个顶点的 全局节点索引
            // 如果以体积为权重
            #pragma GCC unroll 8
            for(uInt k = 0; k < 4; ++k){
                weights[nodes[k]] += cell.m_volume;
            }
        }
        // BlockSparseMatrix<5,5*NumBasis> Rmat;
        for(uInt cellId = 0; cellId < mesh.m_cells.size(); ++cellId){
            const auto& cell = mesh.m_cells[cellId];
            const auto& nodes = cell.m_nodes;  // 四个顶点的 全局节点索引
            // 如果以体积为权重
            #pragma GCC unroll 8
            for(uInt k = 0; k < 4; ++k){
                Scalar weight = cell.m_volume/weights[nodes[k]];
                uInt Rmat_row = nodes[k];
                uInt Rmat_col = cellId;
                DenseMatrix<5,5*NumBasis> Rmat_blk;
                #pragma GCC unroll 40
                for(uInt t = 0; t < NumBasis; ++t){
                    #pragma GCC unroll 8
                    for(uInt s = 0; s < 5; ++s){
                        Rmat_blk(s,s + 5*t) = weight * phi_N[k][t];
                    }
                }
                Rmat.add_block(Rmat_row,Rmat_col,Rmat_blk);
            }
        }
        Rmat.finalize();
        RAP = compute_XYt(Rmat,compute_XYt(Rmat,A));
        // std::cout << 1111 << std::endl;
        

        Rmat.output_as_scalar("Rmat.txt");
        compute_XYt(Rmat,A).output_as_scalar("RAtmat.txt");
        RAP.output_as_scalar("RAPmat.txt");

        // for(uInt brow=0; brow<RAP.num_block_rows; ++brow){
        //     for(uInt i=0; i<RAP.storage.ell_max_per_row; ++i){
        //         const uInt bcol = RAP.storage.ell_cols[brow][i];
        //         if(bcol == RAP.invalid_index) continue; // 关键跳过
        //         std::cout << "(" << brow << ", " << bcol << ")" << std::endl;
        //     }
        //     const uInt start = RAP.storage.csr_row_ptr[brow];
        //     const uInt end = RAP.storage.csr_row_ptr[brow+1];
            
        //     for(uInt idx = start; idx < end; ++idx) {
        //         // debug(storage.csr_blocks[idx].multiply(vec[idx]));
        //         uInt bcol = RAP.storage.csr_cols[idx];
        //         std::cout << "(" << brow << ", " << bcol << ")" << std::endl;
        //     }
        // }

        RAP_solver = EigenSparseSolver<5,5>(RAP);
        // std::cout << 22222 << std::endl;
    }


    
    Vec apply(const Vec& rhs) override {
        const auto& R_r = Rmat.multiply(rhs);
        // std::cout << 333333 << std::endl;
        // std::cout << R_r.size() <<   "   " <<  R_r.blocks[0].size() << std::endl;
        RAP_solver.set_rhs(R_r);
        // std::cout << 444444444 << std::endl;
        const auto& sol = RAP_solver.SparseLU(R_r); // 无效的初始猜测
        // std::cout << 5555555555555 << std::endl;

        LongVector<5*NumBasis> result(A.num_block_rows); // 预分配空间

        // ELL部分
        for(uInt brow=0; brow<Rmat.num_block_rows; ++brow){
            for(uInt i=0; i<Rmat.storage.ell_max_per_row; ++i){
                const uInt bcol = Rmat.storage.ell_cols[brow][i];
                if(bcol == Rmat.invalid_index) continue; // 关键跳过
                const auto& block = Rmat.storage.ell_blocks[brow][i];
                result[bcol] += block.transpose().multiply(sol[brow]);
            }
        }

        // CSR部分
        for(uInt brow = 0; brow < Rmat.num_block_rows; ++brow) {
            const uInt start = Rmat.storage.csr_row_ptr[brow];
            const uInt end = Rmat.storage.csr_row_ptr[brow+1];
            
            for(uInt idx = start; idx < end; ++idx) {
                uInt bcol = Rmat.storage.csr_cols[idx];
                const auto& block = Rmat.storage.csr_blocks[idx];
                result[bcol] += block.transpose().multiply(sol[brow]);
            }
        }

        return result;

        
        // return sol;
        // return apply(rhs,100,1e-3);
    }

    // Vec apply(const Vec& rhs, uInt max_sweeps = 100, Scalar epsilon = 1e-6) { 

    // }
  















    template<uInt Brow, uInt Bmid, uInt Bcol>
    auto compute_XYt(const BlockSparseMatrix<Brow,Bmid>& X,
                    const BlockSparseMatrix<Bcol,Bmid>& Y) {

        auto gather_row = [](const auto& M, uInt row) {
            using BlockT = typename std::remove_reference_t<decltype(M)>::BlockType;
            std::vector<std::pair<uInt, BlockT>> nz;

            // ELL
            for (uInt i = 0; i < M.storage.ell_max_per_row; ++i) {
                uInt c = M.storage.ell_cols[row][i];
                if (c != M.invalid_index) {
                    nz.emplace_back(c, M.storage.ell_blocks[row][i]);
                }
            }
            // CSR
            uInt s = M.storage.csr_row_ptr[row], e = M.storage.csr_row_ptr[row+1];
            for (uInt idx = s; idx < e; ++idx) {
                nz.emplace_back(M.storage.csr_cols[idx], M.storage.csr_blocks[idx]);
            }
            return nz; // 未排序
        };

        using BlockXY = DenseMatrix<Brow,Bcol>;
        BlockSparseMatrix<Brow,Bcol> result;

        // 可选：缓存 Y 的每一行的哈希表，避免重复构建
        std::vector<std::unique_ptr<std::unordered_map<uInt, DenseMatrix<Bcol,Bmid>>>> Y_maps;
        Y_maps.resize(Y.num_block_rows);

        for (uInt brow = 0; brow < X.num_block_rows; ++brow) {
            // 取 X 的一行（未排序）
            auto Xrow = gather_row(X, brow);

            for (uInt bcol = 0; bcol < Y.num_block_rows; ++bcol) {
                // 构建/取用 Y_row 的哈希
                auto& ptr = Y_maps[bcol];
                if (!ptr) {
                    auto Yrow = gather_row(Y, bcol);
                    auto map_ptr = std::make_unique<
                        std::unordered_map<uInt, DenseMatrix<Bcol,Bmid>>
                    >();
                    map_ptr->reserve(Yrow.size() * 2);
                    for (auto& [c, blk] : Yrow) (*map_ptr)[c] = blk;
                    ptr = std::move(map_ptr);
                }
                const auto& Ymap = *ptr;

                auto acc = BlockXY::Zeros();
                bool hit = false;

                // 遍历 Xrow，哈希查找 Y 的同列
                for (const auto& [cX, BX] : Xrow) {
                    auto it = Ymap.find(cX);
                    if (it != Ymap.end()) {
                        // XY^T：Brow×Bmid 乘 Bcol×Bmid^T（注意你的 multiply 接口）
                        acc += BX.multiply(it->second.transpose());
                        hit = true;
                    }
                }

                if (hit) result.add_block(brow, bcol, acc);
            }
        }

        result.finalize();
        return result;
    }




    // template<uInt Brow, uInt Bmid, uInt Bcol>
    // auto compute_XYt(const BlockSparseMatrix<Brow,Bmid>& X, const BlockSparseMatrix<Bcol,Bmid>& Y){
        
    //     auto get_X = [&](uInt brow, uInt Xcol_ptr) -> std::tuple<uInt,DenseMatrix<Brow,Bmid>> {
    //         if(Xcol_ptr<X.storage.ell_max_per_row) {
    //             uInt i = Xcol_ptr;
    //             const uInt bcol = X.storage.ell_cols[brow][i];
    //             // 这一行的元素数量，没有填满 ell 部分，csr 不用再找了
    //             if(bcol == X.invalid_index) return {uInt(-1), DenseMatrix<Brow,Bmid>::Zeros()};    
    //             const auto& block = X.storage.ell_blocks[brow][i];
    //             return {bcol, block};

    //         }
    //         else{
    //             // ell 部分满了，找 csr 部分
    //             const uInt start = X.storage.csr_row_ptr[brow];
    //             const uInt end = X.storage.csr_row_ptr[brow+1];
    //             uInt idx = (Xcol_ptr - X.storage.ell_max_per_row) + start;
    //             if(idx < end){
    //                 const uInt bcol = X.storage.csr_cols[idx];
    //                 const auto& block = X.storage.csr_blocks[idx]; 
    //                 return {bcol, block};
    //             }
    //             else{
    //                 return {uInt(-1), DenseMatrix<Brow,Bmid>::Zeros()};
    //             }
    //         }
    //     };

    //     auto get_Y = [&](uInt brow, uInt Ycol_ptr) -> std::tuple<uInt,DenseMatrix<Bcol,Bmid>> {
    //         if(Ycol_ptr<Y.storage.ell_max_per_row) {
    //             uInt i = Ycol_ptr;
    //             const uInt bcol = Y.storage.ell_cols[brow][i];
    //             // 这一行的元素数量，没有填满 ell 部分，csr 不用再找了
    //             if(bcol == Y.invalid_index) return {uInt(-1), DenseMatrix<Bcol,Bmid>::Zeros()};    
    //             const auto& block = Y.storage.ell_blocks[brow][i];
    //             return {bcol, block};

    //         }
    //         else{
    //             // ell 部分满了，找 csr 部分
    //             const uInt start = Y.storage.csr_row_ptr[brow];
    //             const uInt end = Y.storage.csr_row_ptr[brow+1];
    //             uInt idx = (Ycol_ptr - Y.storage.ell_max_per_row) + start;
    //             if(idx < end){
    //                 const uInt bcol = Y.storage.csr_cols[idx];
    //                 const auto& block = Y.storage.csr_blocks[idx]; 
    //                 return {bcol, block};
    //             }
    //             else{
    //                 return {uInt(-1), DenseMatrix<Bcol,Bmid>::Zeros()};
    //             }
    //         }
    //     };


    //     BlockSparseMatrix<Brow,Bcol> result;
    //     for(uInt brow = 0; brow < X.num_block_rows; ++brow){
    //         for(uInt bcol = 0; bcol < Y.num_block_rows; ++bcol){
    //             DenseMatrix<Brow,Bcol> acc;
    //             bool nonzeros = false;
    //             // 对矩阵 X
    //             for(uInt Xcol_ptr = 0, Ycol_ptr = 0; ;){
    //                 auto [Xcol, Xmat] = get_X(brow,Xcol_ptr);
    //                 auto [Ycol, Ymat] = get_Y(bcol,Ycol_ptr);
    //                 if ((Xcol == uInt(-1)) || (Ycol == uInt(-1))) break;
    //                 if (Xcol != Ycol){
    //                     if(Xcol < Ycol){
    //                         ++Xcol_ptr;
    //                     }
    //                     else{
    //                         ++Ycol_ptr;
    //                     }
    //                 }
    //                 else{
    //                     acc += Xmat.multiply(Ymat.transpose());
    //                     ++Xcol_ptr;
    //                     ++Ycol_ptr;
    //                     nonzeros = true;
    //                 }
    //             }
                
    //             // std::cout << "(" << brow << "," << bcol << ")" << (nonzeros?"true":"false") << std::endl;
    //             if(nonzeros){
    //                 result.add_block(brow,bcol,acc);
    //             }
                
    //         }
    //     }
    //     result.finalize();
    //     return result;
    // }


};

