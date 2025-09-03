#pragma once
#include "Matrix/DenseMatrix.h"
#include "Matrix/SparseMatrix.h"
#include "Matrix/LongVector.h"
#include "Matrix/submat.h"

// 模板类实现 PGMRES
template<uInt BlockDim, uInt Max_dims, bool output_flag = false>
class PGMRES;

// Max_dims：控制 Hessenberg 最大维数
// Preconditioner：预条件器类型，接口为 apply(const Vec&, Vec&)
template<uInt BlockDim>
class Preconditioner {
public:
    using Vec = LongVector<BlockDim>;
    virtual Vec apply(const Vec& rhs) = 0;
    virtual ~Preconditioner() = default;
};

template<uInt BlockDim>
class IdentityPreconditioner : public Preconditioner<BlockDim>{
public:
    using Vec = LongVector<BlockDim>;
    IdentityPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs) : A(A_) {     
    }
    Vec apply(const Vec& rhs) override {
        return rhs;
    }
private:
    static constexpr uInt NumBasis = BlockDim/5;
    const BlockSparseMatrix<BlockDim, BlockDim>& A;
};
template<uInt BlockDim>
class DiagPreconditioner : public Preconditioner<BlockDim>{
public:
    using Vec = LongVector<BlockDim>;
    DiagPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs) : A(A_) { 
        DiagVec.resize(DoFs);
        for(uInt brow=0; brow<A.num_block_rows; ++brow){
            for(uInt i=0; i<A.storage.ell_max_per_row; ++i){
                const uInt bcol = A.storage.ell_cols[brow][i];
                if(bcol == brow) {
                    const auto& block = A.storage.ell_blocks[brow][i];
                    for(uInt row = 0; row < BlockDim; row++){
                        DiagVec[brow][row] = 1/block(row,row);
                    }
                }
                
            }
            const uInt start = A.storage.csr_row_ptr[brow];
            const uInt end = A.storage.csr_row_ptr[brow+1];
            
            for(uInt idx = start; idx < end; ++idx) {
                const uInt bcol = A.storage.csr_cols[idx];
                if(bcol == brow) {
                    const auto& block = A.storage.csr_blocks[idx];
                    for(uInt row = 0; row < BlockDim; row++){
                        DiagVec[brow][row] = 1/block(row,row);
                    }
                }
            }
        }
    }
    Vec apply(const Vec& rhs) override {
        return rhs * DiagVec;
    }
private:
    static constexpr uInt NumBasis = BlockDim/5;
    const BlockSparseMatrix<BlockDim, BlockDim>& A;
    Vec DiagVec;
};

template<uInt BlockDim>
class BJacPreconditioner : public Preconditioner<BlockDim>{
public:
    using Vec = LongVector<BlockDim>;
    BJacPreconditioner(const BlockSparseMatrix<BlockDim, BlockDim>& A_, uInt DoFs) : A(A_) { 
        DiagVec.resize(DoFs);
        #pragma omp parallel for schedule(dynamic)
        for(uInt brow=0; brow<A.num_block_rows; ++brow){
            for(uInt i=0; i<A.storage.ell_max_per_row; ++i){
                const uInt bcol = A.storage.ell_cols[brow][i];
                if(bcol == brow) {
                    const auto& block = A.storage.ell_blocks[brow][i];
                    DiagVec[brow] = block.lu();
                }
                
            }
            const uInt start = A.storage.csr_row_ptr[brow];
            const uInt end = A.storage.csr_row_ptr[brow+1];
            
            for(uInt idx = start; idx < end; ++idx) {
                const uInt bcol = A.storage.csr_cols[idx];
                if(bcol == brow) {
                    const auto& block = A.storage.csr_blocks[idx];
                    DiagVec[brow] = block.lu();
                }
            }
        }
    }
    Vec apply(const Vec& rhs) override {
        Vec ret(rhs.size());
        #pragma omp parallel for schedule(dynamic)
        for(uInt brow=0; brow<rhs.size(); ++brow){
            ret[brow] = DiagVec[brow].solve(rhs[brow],DiagVec[brow]);
        }
        return ret;
    }
private:
    static constexpr uInt NumBasis = BlockDim/5;
    const BlockSparseMatrix<BlockDim, BlockDim>& A;
    std::vector<DenseMatrix<BlockDim, BlockDim>> DiagVec;
};