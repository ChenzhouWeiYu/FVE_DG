#pragma once
#include "base/Type.h"
#include "DG/DG_Basis/DG_Basis.h"
#include "Mesh/Mesh.h"
#include "Matrix/Matrix.h"
#include "base/exact.h"
#include "DG/DG_Flux/EulerPhysicalFlux.h"

template<uInt Order=3, typename Flux = AirFluxC, 
        typename GaussQuadCell = GaussLegendreTet::Auto, 
        typename GaussQuadFace = GaussLegendreTri::Auto>
class ExplicitConvection {
private:
    using BlockMat = DenseMatrix<5,5>;
    using Basis = DGBasisEvaluator<Order>;

    using QuadC = typename std::conditional_t<
        std::is_same_v<GaussQuadCell, GaussLegendreTet::Auto>,
        typename AutoQuadSelector<Order, GaussLegendreTet::Auto>::type,
        GaussQuadCell
    >;
    using QuadF = typename std::conditional_t<
        std::is_same_v<GaussQuadFace, GaussLegendreTri::Auto>,
        typename AutoQuadSelector<Order, GaussLegendreTri::Auto>::type,
        GaussQuadFace
    >;

    static constexpr uInt N = Basis::NumBasis;
    vector3f transform_to_cell(const CompTriangleFace& face, const vector2f& uv, uInt side) const;

public:
    LongVector<5*N> eval(const ComputingMesh& mesh, 
                 const LongVector<5*N>& old_solution,
                 const Scalar curr_time = 0.0);
    LongVector<5*N> eval_cells(const ComputingMesh& mesh, 
                 const LongVector<5*N>& old_solution,
                 const Scalar curr_time = 0.0);
    LongVector<5*N> eval_internals(const ComputingMesh& mesh, 
                 const LongVector<5*N>& old_solution,
                 const Scalar curr_time = 0.0);
    LongVector<5*N> eval_boundarys(const ComputingMesh& mesh, 
                 const LongVector<5*N>& old_solution,
                 const Scalar curr_time = 0.0);
                
};


// #define Explicit_For_Flux(Order) \
// extern template class ExplicitConvection<Order,AirFlux>;\
// extern template class ExplicitConvection<Order,MonatomicFlux>;\
// extern template class ExplicitConvection<Order+1,AirFlux>;\
// extern template class ExplicitConvection<Order+1,MonatomicFlux>;\
// extern template class ExplicitConvection<Order+2,AirFlux>;\
// extern template class ExplicitConvection<Order+2,MonatomicFlux>;

// Explicit_For_Flux(0)
// Explicit_For_Flux(3)

// #undef Explicit_For_Flux