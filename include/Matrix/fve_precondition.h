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
template<uInt BlockDim, uInt Max_dims, bool output_flag = false>
class PGMRES;


template<uInt BlockDim, uInt InnerMaxIters = 100>
class FVEPreconditioner : public Preconditioner<BlockDim> {
public:
    using Vec = LongVector<BlockDim>;
    static constexpr uInt NumBasis = BlockDim / 5;
    static constexpr uInt Order = (NumBasis < 10 ? (NumBasis == 1 ? 0 : 1) : (NumBasis==10 ? 2 : 3));
private:
    ComputingMesh mesh;
    BlockSparseMatrix<BlockDim, BlockDim> A;
    BlockSparseMatrix<5,5*NumBasis> Rmat;
    BlockSparseMatrix<5,5> RAP;
public:

    // template<uInt BlockRows, uInt BlockCols>
    // BlockSparseMatrix<BlockRows,BlockCols> block_to_eigen(const BlockSparseMatrix<BlockRows,BlockCols>& mat){
    //     uInt ELL_row_length = mat.storage.ell_blocks.size();
    //     uInt ELL_max_row = mat.storage.ell_max_per_row;
    //     uInt CSR_nnz = mat.storage.csr_blocks.size();
    //     uInt NNZ = (ELL_row_length * ELL_max_row + CSR_nnz) * BlockRows * BlockCols;
    //     uInt ROWs = mat.num_block_rows * BlockRows;

    //     std::vector<Triplet> m_tripletList;
    //     m_tripletList.reserve(NNZ);

    //     for(uInt brow=0; brow<mat.num_block_rows; ++brow){
    //         for(uInt i=0; i<mat.storage.ell_max_per_row; ++i){
    //             const uInt bcol = mat.storage.ell_cols[brow][i];
    //             if(bcol == mat.invalid_index) continue;
    //             const auto& block = mat.storage.ell_blocks[brow][i];
    //             for(uInt row=0;row<block.rows();row++){
    //                 for(uInt col=0;col<block.cols();col++){
    //                     // if(std::abs(block(row,col))>1e-15)
    //                     m_tripletList.push_back(Triplet(brow*block.rows()+row,bcol*block.cols()+col,block(row,col)));
    //                 }
    //             }
    //         }

    //         const uInt start = mat.storage.csr_row_ptr[brow];
    //         const uInt end = mat.storage.csr_row_ptr[brow+1];
    //         for(uInt idx = start; idx < end; ++idx) {
    //             const uInt bcol = mat.storage.csr_cols[idx];
    //             const auto& block = mat.storage.csr_blocks[idx];

    //             for(uInt row=0;row<block.rows();row++){
    //                 for(uInt col=0;col<block.cols();col++){
    //                     // if(std::abs(block(row,col))>1e-15)
    //                     m_tripletList.push_back(Triplet(brow*block.rows()+row,bcol*block.cols()+col,block(row,col)));
    //                 }
    //             }
    //         }
    //     }

    //     EigenCSC m_CSCmat(m_eigenRhs.size(), m_eigenRhs.size());
        
    //     m_CSCmat.setFromTriplets(m_tripletList.begin(), m_tripletList.end());

    // }

    template<uInt Brow, uInt Bmid, uInt Bcol>
    auto compute_XYt(const BlockSparseMatrix<Brow,Bmid>& X, const BlockSparseMatrix<Bcol,Bmid>& Y){
        
        auto get_X = [&](uInt brow, uInt Xcol_ptr) -> std::tuple<uInt,DenseMatrix<Brow,Bmid>> {
            if(Xcol_ptr<X.storage.ell_max_per_row) {
                uInt i = Xcol_ptr;
                const uInt bcol = X.storage.ell_cols[brow][i];
                // 这一行的元素数量，没有填满 ell 部分，csr 不用再找了
                if(bcol == X.invalid_index) return {uInt(-1), DenseMatrix<Brow,Bmid>::Zeros()};    
                const auto& block = X.storage.ell_blocks[brow][i];
                return {bcol, block};

            }
            else{
                // ell 部分满了，找 csr 部分
                const uInt start = X.storage.csr_row_ptr[brow];
                const uInt end = X.storage.csr_row_ptr[brow+1];
                uInt idx = (Xcol_ptr - X.storage.ell_max_per_row) + start;
                if(idx < end){
                    const uInt bcol = X.storage.csr_cols[idx];
                    const auto& block = X.storage.csr_blocks[idx]; 
                    return {bcol, block};
                }
                else{
                    return {uInt(-1), DenseMatrix<Brow,Bmid>::Zeros()};
                }
            }
        };

        auto get_Y = [&](uInt brow, uInt Ycol_ptr) -> std::tuple<uInt,DenseMatrix<Bcol,Bmid>> {
            if(Ycol_ptr<Y.storage.ell_max_per_row) {
                uInt i = Ycol_ptr;
                const uInt bcol = Y.storage.ell_cols[brow][i];
                // 这一行的元素数量，没有填满 ell 部分，csr 不用再找了
                if(bcol == Y.invalid_index) return {uInt(-1), DenseMatrix<Bcol,Bmid>::Zeros()};    
                const auto& block = Y.storage.ell_blocks[brow][i];
                return {bcol, block};

            }
            else{
                // ell 部分满了，找 csr 部分
                const uInt start = Y.storage.csr_row_ptr[brow];
                const uInt end = Y.storage.csr_row_ptr[brow+1];
                uInt idx = (Ycol_ptr - Y.storage.ell_max_per_row) + start;
                if(idx < end){
                    const uInt bcol = Y.storage.csr_cols[idx];
                    const auto& block = Y.storage.csr_blocks[idx]; 
                    return {bcol, block};
                }
                else{
                    return {uInt(-1), DenseMatrix<Bcol,Bmid>::Zeros()};
                }
            }
        };


        BlockSparseMatrix<Brow,Bcol> result;
        for(uInt brow = 0; brow < Brow; ++brow){
            for(uInt bcol = 0; bcol < Bcol; ++bcol){
                DenseMatrix<Brow,Bcol> acc;
                // 对矩阵 X
                for(uInt Xcol_ptr = 0, Ycol_ptr = 0; ((Xcol != uInt(-1)) && (Ycol != uInt(-1)));){
                    auto [Xcol, Xmat] = get_X(brow,Xcol_ptr);
                    auto [Ycol, Ymat] = get_Y(brow,Ycol_ptr);
                    if (Xcol != Ycol){
                        if(Xcol_ptr < Ycol_ptr){
                            ++Xcol_ptr;
                        }
                        else{
                            ++Ycol_ptr;
                        }
                    }
                    else{
                        acc += Xmat.multiply(Ymat.transpose());
                        ++Xcol_ptr;
                        ++Ycol_ptr;
                    }
                }
                result.add_block(brow,bcol,acc);
            }
        }
        result.finalize();
        return result;
    }

    FVEPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs, const ComputingMesh& mesh_)
        : A(A_), mesh(mesh_)
          // 提取子块
        //   Arr(get_sub_sparse<NumBasis,1,1,0,0>(A)),
        //   Arm(get_sub_sparse<NumBasis,1,3,0,1>(A)),
        //   Are(get_sub_sparse<NumBasis,1,1,0,4>(A)),

        //   Amr(get_sub_sparse<NumBasis,3,1,1,0>(A)),
        //   Amm(get_sub_sparse<NumBasis,3,3,1,1>(A)),
        //   Ame(get_sub_sparse<NumBasis,3,1,1,4>(A)),

        //   Aer(get_sub_sparse<NumBasis,1,1,4,0>(A)),
        //   Aem(get_sub_sparse<NumBasis,1,3,4,1>(A)),
        //   Aee(get_sub_sparse<NumBasis,1,1,4,4>(A)),

        //   Arr_solver(Arr), Amm_solver(Amm), Aee_solver(Aee)
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
            #pragma unroll
            for(uInt k = 0; k < 4; ++k){
                weights[nodes[k]] += cell.m_volume;
            }
        }
        // BlockSparseMatrix<5,5*NumBasis> Rmat;
        for(uInt cellId = 0; cellId < mesh.m_cells.size(); ++cellId){
            const auto& cell = mesh.m_cells[cellId];
            const auto& nodes = cell.m_nodes;  // 四个顶点的 全局节点索引
            // 如果以体积为权重
            #pragma unroll
            for(uInt k = 0; k < 4; ++k){
                Scalar weight = cell.m_volume/weights[nodes[k]];
                uInt Rmat_row = node[k];
                uInt Rmat_col = cell;
                DenseMatrix<5,5*NumBasis> Rmat_blk;
                for(uInt t = 0; t < NumBasis; ++t){
                    for(uInt s = 0; s < 5; ++s){
                        Rmat_blk[s,s + 5*t] = weight * phi_N[k][t];
                    }
                }
                Rmat.add_block(Rmat_row,Rmat_col,Rmat_blk);
            }
        }
        Rmat.finalize();
        RAP = compute_XYt(Rmat,compute_XYt(Rmat,A));
    }


    
    Vec apply(const Vec& rhs) override {
        return apply(rhs,100,1e-3);
    }

    Vec apply(const Vec& rhs, uInt max_sweeps = 100, Scalar epsilon = 1e-6) { 

    }
  



};

