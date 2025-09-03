const auto& Arr = get_sub_sparse<Basis::NumBasis,1,1,0,0>(sparse_mat);                         
             
const auto& Amr = get_sub_sparse<Basis::NumBasis,3,1,1,0>(sparse_mat);
                
const auto& Amm = get_sub_sparse<Basis::NumBasis,3,3,1,1>(sparse_mat);
               
const auto& Ame = get_sub_sparse<Basis::NumBasis,3,1,1,4>(sparse_mat);                
const auto& bm = get_sub_vector<Basis::NumBasis,3,1>(rhs);
auto Um = get_sub_vector<Basis::NumBasis,3,1>(U_k);
const auto& Aer = get_sub_sparse<Basis::NumBasis,1,1,4,0>(sparse_mat);
const auto& Aem = get_sub_sparse<Basis::NumBasis,1,3,4,1>(sparse_mat);
const auto& Aee = get_sub_sparse<Basis::NumBasis,1,1,4,4>(sparse_mat);        
                       
const auto& be = get_sub_vector<Basis::NumBasis,1,4>(rhs);
               auto Ue = get_sub_vector<Basis::NumBasis,1,4>(U_k);
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
                  const auto& submat_residual = br - Arr.multiply(Ur) - Arm.multiply(Um) - Are.multiply(Ue);  //cg 06.17
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Submat residual = " + oss.str());
                }
                {   
                    const auto& submat_residual = bm - Amr.multiply(Ur) - Amm.multiply(Um) - Ame.multiply(Ue);  //cg 06.17
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Submat residual = " + oss.str());
                }
                {   
                    const auto& submat_residual = be - Aer.multiply(Ur) - Aem.multiply(Um) - Aee.multiply(Ue);  //cg 06.17
                    std::ostringstream oss;
                    oss << std::sqrt(submat_residual.dot(submat_residual));
                    logging("Picard iter " + std::to_string(picard_iter) + "\t Submat residual = " + oss.str());
                }

                logging("Start linear solver");