#include "base/Type.h"
#include "Mesh/Mesh.h"
#include "Matrix/Matrix.h"
#include "DG/DG_Basis/DG_Basis.h"

ComputingMesh create_mesh(uInt N){
    GeneralMesh mesh = OrthHexMesh({0.0, 0.0, 0.0},{3.2, 1.0, 1.0/N},{(18*N)/5,N,1});
    mesh.split_hex5_scan();                                   
    mesh.rebuild_cell_topology();                             
    mesh.validate_mesh();                                     
    ComputingMesh cmesh(mesh);                                
    cmesh.m_boundaryTypes.resize(cmesh.m_faces.size());                   
    for(uInt faceId=0;faceId<cmesh.m_faces.size();faceId++){           
        if(cmesh.m_faces[faceId].m_neighbor_cells[1]==uInt(-1)){ 
            const auto& face = cmesh.m_faces[faceId];           
            if(std::abs(face.m_normal[2])>0.5 )          
                cmesh.m_boundaryTypes[faceId] = BoundaryType::Pseudo3DZ;
            else{
                const auto& nodes = face.m_nodes;
                const auto& p0 = cmesh.m_points[nodes[0]];
                const auto& p1 = cmesh.m_points[nodes[1]];
                const auto& p2 = cmesh.m_points[nodes[2]];
                const auto& p3 = cmesh.m_points[nodes[3]];
                const auto& centor = (p0+p1+p2+p3)/4;
                if(centor[0]>1.0/6.0 && centor[1] == 0){
                    cmesh.m_boundaryTypes[faceId] = BoundaryType::Pseudo3DY;
                }
                else
                {
                    cmesh.m_boundaryTypes[faceId] = BoundaryType::Dirichlet;
                }
            }
                
        }
    }
    return cmesh;
}

template<typename QuadC, typename Basis>
std::tuple<LongVector<5*QuadC::num_points>, 
            LongVector<5*QuadC::num_points>, 
            LongVector<5*Basis::NumBasis>> 
reconstruct_solution(const ComputingMesh& cmesh, LongVector<5*Basis::NumBasis>& coef, Scalar curr_time){
    LongVector<5*QuadC::num_points> U_h(cmesh.m_cells.size());
    LongVector<5*QuadC::num_points> U_s(cmesh.m_cells.size());
    LongVector<5*Basis::NumBasis> error(cmesh.m_cells.size());

    #pragma omp parallel for schedule(dynamic)
    for(uInt cellId=0;cellId<cmesh.m_cells.size();cellId++){
        const auto& cell = cmesh.m_cells[cellId];
        const auto& U_func = [&](vector3f Xi)->DenseMatrix<5,1>{
            return U_Xi(cell,Xi,curr_time);
        };
        const auto& U_coef = Basis::func2coef(U_func);
        for(uInt k=0;k<Basis::NumBasis;k++){
            const auto& block_coef = coef[cellId].template SubMat<5,1>(5*k,0) - U_coef[k];
            MatrixView<5*Basis::NumBasis,1,5,1>(error[cellId],5*k,0) = block_coef;
        }
        for(uInt xgId=0; xgId<QuadC::num_points; ++xgId) {
            const auto& p = QuadC::points[xgId];
            const auto& pos = cell.transform_to_physical(p);
            const auto& value = Basis::eval_all(p[0],p[1],p[2]);
            const auto& U = Basis::template coef2filed<5,Scalar>(coef[cellId],p);
            MatrixView<5*QuadC::num_points,1,5,1> block_U_h(U_h[cellId],5*xgId,0);
            MatrixView<5*QuadC::num_points,1,5,1> block_U_s(U_s[cellId],5*xgId,0);
            
            block_U_h = U;
            block_U_h[0] = U[0];
            block_U_s = U_func(p);
        }
    }
    return {U_h, U_s, error};

}