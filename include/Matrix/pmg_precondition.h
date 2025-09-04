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

template<uInt BlockDim, uInt VcycleMaxIters = 1, uInt LowerOrder = 0, 
    typename Smoother=BJacPreconditioner<BlockDim>, uInt PreIters=0, uInt PostIters=5>
class PMGPreconditioner;

template<uInt BlockDim, uInt VcycleMaxIters, uInt LowerOrder, typename Smoother, uInt PreIters, uInt PostIters>
class PMGPreconditioner : public Preconditioner<BlockDim> {
// public:
private:
    using Vec = LongVector<BlockDim>;
    static constexpr uInt NumBasis = BlockDim / 5;
    static constexpr uInt Order = (NumBasis < 10 ? (NumBasis == 1 ? 0 : 1) : (NumBasis==10 ? 2 : 3));

    const BlockSparseMatrix<5*NumBasis, 5*NumBasis>& A;
    Smoother smother;

    // 粗网格的多项式阶为 LowerOrder
    static constexpr uInt LowerNumBasis = DGBasisEvaluator<LowerOrder>::NumBasis;
    BlockSparseMatrix<5*LowerNumBasis, 5*LowerNumBasis> RAP; 

    EigenSparseSolver<5*LowerNumBasis, 5*LowerNumBasis> RAP_solver;
    std::unique_ptr<PMGPreconditioner<5*LowerNumBasis, VcycleMaxIters, LowerOrder==0?0:LowerOrder-1>> RAP_preconditioner_ptr;
    
    // 自定义的 V-Cycle 次数、前后磨光次数（暂时不生效，以后再弄）
    uInt vcycle_iters = VcycleMaxIters;
    uInt pre_iters = PreIters;
    uInt post_iters = PostIters;
public:
    void set_vcycle_iters(uInt vcycle_iters_){
        vcycle_iters = vcycle_iters_;
    }
    void set_pre_iters(uInt pre_iters_){
        pre_iters = pre_iters_;
    }
    void set_post_iters(uInt post_iters_){
        post_iters = post_iters_;
    }

public:


    PMGPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs, const ComputingMesh& mesh_)
        : PMGPreconditioner(A_, DoFs) {}
    
    PMGPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs)
        : A(A_), smother(Smoother(A_, DoFs)) 
    {   
        // print("PMG Preconditioner Constructor");
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
        // print("PMG Preconditioner RAP Matrix Constructed");
        RAP.finalize();
        // print("PMG Preconditioner RAP Matrix Finalized");
        if constexpr (LowerOrder == 0) RAP_solver = EigenSparseSolver<5*LowerNumBasis,5*LowerNumBasis>(RAP);
        else RAP_preconditioner_ptr = std::make_unique<PMGPreconditioner<5*LowerNumBasis, VcycleMaxIters, LowerOrder==0?0:LowerOrder-1>>(RAP, DoFs);
        // print("PMG Preconditioner RAP Solver Constructed");
    }


    
    Vec apply(const Vec& rhs) override {
        // 设置初始
        Vec x(rhs.size());

        #pragma GCC unroll 20
        for(uInt vcycle = 0; vcycle < VcycleMaxIters; ++vcycle) {
            // 前磨光
            #pragma GCC unroll 20
            for(uInt pre_smoothing = 0; pre_smoothing < PreIters; ++pre_smoothing) {
                x = x + smother.apply(rhs - A.multiply(x));
            }
            // print("PMG Preconditioner Applied Smoother");

            // 限制, 粗网格求解, 插值
            x = x + apply_coarse(rhs - A.multiply(x));
            
            // 后磨光
            #pragma GCC unroll 20
            for(uInt post_smoothing = 0; post_smoothing < PostIters; ++post_smoothing) {
                x = x + smother.apply(rhs - A.multiply(x));
            }

            // 这里是 p=2,3, 时的 W-Cycle 迭代次数有减少，时间提升不明显
            if constexpr (LowerOrder > 0){
                // 限制, 粗网格求解, 插值
                x = x + apply_coarse(rhs - A.multiply(x));
                
                // 后磨光
                #pragma GCC unroll 20
                for(uInt post_smoothing = 0; post_smoothing < PostIters; ++post_smoothing) {
                    x = x + smother.apply(rhs - A.multiply(x));
                }
            }
            // print("PMG Preconditioner Applied Smoother");
            // std::cout << "PMG Preconditioner Iteration " << vcycle << "   " 
            // << "residual norm: " << (rhs - A.multiply(x)).norm() << std::endl;
        }

        return x;
    }

    

  

private:
    // 接收残差的 MG 部分
    inline Vec apply_coarse(const Vec& r) {
        if constexpr (LowerOrder == 0) {
            // 限制
            auto R_r = Rmul(r);
            // print("PMG Preconditioner Restricted RHS to Coarse Grid");

            // 粗网格求解
            RAP_solver.set_rhs(R_r);
            // print("PMG Preconditioner Set RHS");
            auto err = RAP_solver.SparseLU(R_r);
            // print("PMG Preconditioner Solve Coarse Grid System");

            // 插值
            return Pmul(err);
        }
        else{
            // 限制
            auto R_r = Rmul(r);
            // print("PMG Preconditioner Restricted RHS to Coarse Grid");

            // 粗网格求解
            auto err = RAP_preconditioner_ptr->apply(R_r);

            // 插值
            return Pmul(err);
        }

    }

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

