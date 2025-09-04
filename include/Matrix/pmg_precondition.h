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


template<uInt BlockDim, uInt InnerMaxIters = 100, uInt LowerOrder = 0>
class PMGPreconditioner : public Preconditioner<BlockDim> {
public:
    using Vec = LongVector<BlockDim>;
    static constexpr uInt NumBasis = BlockDim / 5;
    static constexpr uInt Order = (NumBasis < 10 ? (NumBasis == 1 ? 0 : 1) : (NumBasis==10 ? 2 : 3));
private:
    const BlockSparseMatrix<5*NumBasis, 5*NumBasis>& A;

    // 粗网格的多项式阶为 LowerOrder
    static constexpr uInt LowerNumBasis = DGBasisEvaluator<LowerOrder>::NumBasis;
    BlockSparseMatrix<5*LowerNumBasis, 5*LowerNumBasis> RAP; 
    EigenSparseSolver<5*LowerNumBasis, 5*LowerNumBasis> RAP_solver;

public:


    PMGPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs, const ComputingMesh& mesh_)
        : PMGPreconditioner(A_, DoFs) {}
    
    PMGPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs)
        : A(A_)
    {   
        print("PMG Preconditioner Constructor");
        DenseMatrix<5*LowerNumBasis, 5*LowerNumBasis> sub_matrix;
        for(uInt brow=0; brow<A.num_block_rows; ++brow){
            for(uInt i=0; i<A.storage.ell_max_per_row; ++i){
                const uInt bcol = A.storage.ell_cols[brow][i];
                if(bcol == A.invalid_index) continue; // 关键跳过
                const auto& block = A.storage.ell_blocks[brow][i];
                for(uInt j=0; j<5*LowerNumBasis; ++j){
                    for(uInt k=0; k<5*LowerNumBasis; ++k){
                        sub_matrix(j,k) = block(j,k);
                    }
                }
                RAP.add_block(brow, bcol, sub_matrix);
            }
            const uInt start = A.storage.csr_row_ptr[brow];
            const uInt end = A.storage.csr_row_ptr[brow+1];
            for(uInt idx = start; idx < end; ++idx) {
                uInt bcol = A.storage.csr_cols[idx];
                const auto& block = A.storage.csr_blocks[idx];
                for(uInt j=0; j<5*LowerNumBasis; ++j){
                    for(uInt k=0; k<5*LowerNumBasis; ++k){
                        sub_matrix(j,k) = block(j,k);
                    }
                }
                RAP.add_block(brow, bcol, sub_matrix);
            }
        }
        print("PMG Preconditioner RAP Matrix Constructed");
        RAP.finalize();
        print("PMG Preconditioner RAP Matrix Finalized");
        RAP_solver = EigenSparseSolver<5*LowerNumBasis,5*LowerNumBasis>(RAP);
        print("PMG Preconditioner RAP Solver Constructed");
    }


    
    Vec apply(const Vec& rhs) override {
        const auto& R_r = Rmul(rhs);
        print("PMG Preconditioner Restricted RHS to Coarse Grid");
        RAP_solver.set_rhs(R_r);
        print("PMG Preconditioner Set RHS");
        const auto& sol = RAP_solver.SparseLU(R_r); // 无效的初始猜测
        print("PMG Preconditioner Solve Coarse Grid System");
        return Pmul(sol);
    }

  

private:
    // 限制算子 R 的作用，实际截取高阶向量的前一段
    LongVector<5*LowerNumBasis> Rmul(LongVector<5*NumBasis> vec) {
        LongVector<5*LowerNumBasis> result(vec.size()); // 预分配空间
        for(uInt blkId = 0; blkId < vec.size(); ++blkId) {
            for(uInt i=0; i<5*LowerNumBasis; ++i){
                result[blkId][i] = vec[blkId][i];
            }
        }
        return result;
    }

    // 插值算子 P 的作用，实际是后填充 0
    LongVector<5*NumBasis> Pmul(LongVector<5*LowerNumBasis> vec) {
        LongVector<5*NumBasis> result(vec.size()); // 预分配空间
        for(uInt blkId = 0; blkId < vec.size(); ++blkId) {
            for(uInt i=0; i<5*LowerNumBasis; ++i){
                result[blkId][i] = vec[blkId][i];
            }
            // 初始化就是 0，后续不用填充
        }
        return result;
    }
};

