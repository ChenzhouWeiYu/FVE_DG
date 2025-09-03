// include/DG/Flux/EulerPhysicalFlux.h
#pragma once
#include "base/Type.h"
#include "Matrix/Matrix.h"

template<uInt GammaNumerator = 7, uInt GammaDenominator = 5>
class EulerPhysicalFlux {
public:
    static constexpr Scalar gamma = 
        static_cast<Scalar>(GammaNumerator) / static_cast<Scalar>(GammaDenominator);
    static constexpr Scalar epslion = 1e-16;
    
    // 改为静态成员函数
    static DenseMatrix<5,3> computeFlux(const DenseMatrix<5,1>& U);
    static std::array<DenseMatrix<5,5>,3> computeJacobian(const DenseMatrix<5,1>& U);
    static Scalar computeWaveSpeed(const DenseMatrix<5,1>& U_L, const DenseMatrix<5,1>& U_R);

private:
    static Scalar positive(Scalar val);
    static Scalar enthalpy(const DenseMatrix<5,1>& U);
    static Scalar soundSpeed(const DenseMatrix<5,1>& U);
    static Scalar waveSpeedMagnitude(const DenseMatrix<5,1>& U);
    static Scalar computePressure(const DenseMatrix<5,1>& U);
};

// #include "EulerPhysicalFlux.cpp"

using MonatomicFluxC = EulerPhysicalFlux<5,3>;  // gamma=5/3
using AirFluxC = EulerPhysicalFlux<7,5>;         // gamma=7/5