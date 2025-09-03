#pragma once
#include "Matrix/SparseMatrix.h"
#include "Matrix/LongVector.h"

template<uInt Nbasis, uInt SubN, uInt SubM, uInt StartN, uInt StartM>
auto get_sub_sparse(const BlockSparseMatrix<5*Nbasis,5*Nbasis>& sparse_mat){
    BlockSparseMatrix<SubN*Nbasis,SubM*Nbasis> ret;
    for(uInt brow = 0; brow < sparse_mat.num_block_rows; ++brow){
        for(uInt i = 0; i < sparse_mat.storage.ell_max_per_row; ++i){
            const uInt bcol = sparse_mat.storage.ell_cols[brow][i];
            if(bcol == sparse_mat.invalid_index) continue;
            const auto& block = sparse_mat.storage.ell_blocks[brow][i];
            DenseMatrix<SubN*Nbasis,SubM*Nbasis> SubA;
            for(uInt basis_i = 0; basis_i<Nbasis; basis_i++){
                for(uInt basis_j = 0; basis_j<Nbasis; basis_j++){
                    for(uInt row = 0; row < SubN; ++row){
                        for(uInt col = 0; col < SubM; ++col){
                            // if(brow==0 && i==0){
                            //     debug(vector4u{{row,col,SubN*basis + row,SubM*basis + col}});
                            //     debug(block(5*basis + row + StartN, 5*basis + col + StartM));
                            // }
                            // assert(5*basis + row + StartN < 5*Nbasis); // cg 06.17
                            // assert(5*basis + col + StartM < 5*Nbasis); // cg 06.17                                           
                            SubA(SubN*basis_i + row,SubM*basis_j + col) = 
                            block(5*basis_i + row + StartN, 5*basis_j + col + StartM);
                        }
                    }
                }
            }
            ret.add_block(brow,bcol,SubA);
        }

        const uInt start = sparse_mat.storage.csr_row_ptr[brow];
        const uInt end = sparse_mat.storage.csr_row_ptr[brow + 1];
        for(uInt idx = start; idx < end; ++idx){
            const uInt bcol = sparse_mat.storage.csr_cols[idx];
            const auto& block = sparse_mat.storage.csr_blocks[idx];
            DenseMatrix<SubN*Nbasis,SubM*Nbasis> SubA;
            for(uInt basis_i = 0; basis_i<Nbasis; basis_i++){
                for(uInt basis_j = 0; basis_j<Nbasis; basis_j++){
                    for(uInt row = 0; row < SubN; ++row){
                        for(uInt col = 0; col < SubM; ++col){
                            // if(brow==0 && i==0){
                            //     debug(vector4u{{row,col,SubN*basis + row,SubM*basis + col}});
                            //     debug(block(5*basis + row + StartN, 5*basis + col + StartM));
                            // }
                            // assert(5*basis + row + StartN < 5*Nbasis); // cg 06.17
                            // assert(5*basis + col + StartM < 5*Nbasis); // cg 06.17                                           
                            SubA(SubN*basis_i + row,SubM*basis_j + col) = 
                            block(5*basis_i + row + StartN, 5*basis_j + col + StartM);
                        }
                    }
                }
            }
            ret.add_block(brow,bcol,SubA);
        }
    }
    ret.finalize();
    return ret;
}

template<uInt Nbasis, uInt SubN, uInt StartN>
auto get_sub_vector(const LongVector<5*Nbasis>& rhs){
    LongVector<SubN*Nbasis> ret(rhs.size());
    for(uInt blockId = 0; blockId < rhs.size(); blockId++){
        const auto& block = rhs[blockId];
        for(uInt basis = 0; basis<Nbasis; basis++){
            for(uInt row = 0; row < SubN; ++row){
                ret[blockId](SubN*basis + row,0) = block(5*basis + row + StartN, 0);
            }
        }
    }
    return ret;
}

template<uInt Nbasis, uInt SubN, uInt StartN>
void set_sub_vector(const LongVector<SubN*Nbasis>& ret, LongVector<5*Nbasis>& rhs){
    for(uInt blockId = 0; blockId < rhs.size(); blockId++){
        auto& block = rhs[blockId];
        for(uInt basis = 0; basis<Nbasis; basis++){
            for(uInt row = 0; row < SubN; ++row){
                block(5*basis + row + StartN, 0) = ret[blockId](SubN*basis + row,0);
            }
        }
    }
}

// template<uInt Nbasis, uInt SubN, uInt StartN>
// LongVector<5*Nbasis> set_sub_vector(const LongVector<SubN*Nbasis>& ret){
//     LongVector<5*Nbasis> rhs;
//     for(uInt blockId = 0; blockId < rhs.size(); blockId++){
//         auto& block = rhs[blockId];
//         for(uInt basis = 0; basis<Nbasis; basis++){
//             for(uInt row = 0; row < SubN; ++row){
//                 block(5*basis + row + StartN, 0) = ret[blockId](SubN*basis + row,0);
//             }
//         }
//     }
// }