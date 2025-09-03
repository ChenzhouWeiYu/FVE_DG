#include "DG/DG_Flux/EulerPhysicalFlux.h"

template class EulerPhysicalFlux<5,3>;  // gamma=5/3
template class EulerPhysicalFlux<7,5>;  // gamma=7/5


template<uInt GammaNumerator, uInt GammaDenominator>
DenseMatrix<5,3> EulerPhysicalFlux<GammaNumerator, GammaDenominator>::computeFlux(const DenseMatrix<5,1>& U){
    const Scalar rho = positive(U[0]);
    const Scalar u = U[1]/rho, v = U[2]/rho, w = U[3]/rho, e = U[4]/rho;
    const Scalar u2 = u*u + v*v + w*w;
    const Scalar p = computePressure(U);
    return {rho*u,rho*v,rho*w,
            rho*u*u+p,rho*v*u  ,rho*w*u  ,
            rho*u*v  ,rho*v*v+p,rho*w*v  ,
            rho*u*w  ,rho*v*w  ,rho*w*w+p,
            u*(rho*e+p),v*(rho*e+p),w*(rho*e+p)};
}

template<uInt GammaNumerator, uInt GammaDenominator>
std::array<DenseMatrix<5,5>,3> EulerPhysicalFlux<GammaNumerator, GammaDenominator>::computeJacobian(const DenseMatrix<5,1>& U) {
    const Scalar rho = positive(U[0]);
    const Scalar u = U[1]/U[0], v = U[2]/U[0], w = U[3]/U[0], e = U[4]/U[0];
    const Scalar u2 = u*u + v*v + w*w;
    const Scalar p = computePressure(U);
    const DenseMatrix<5,5> Jac_x = {0,1,0,0,0,  
                    -u*u+0.5*(gamma-1)*u2, 2*u-(gamma-1)*u, -(gamma-1)*v,-(gamma-1)*w,gamma-1,
                    -u*v,v,u,0,0,
                    -u*w,w,0,u,0,
                    u*(-gamma*e+(gamma-1)*u2),gamma*e-0.5*(gamma-1)*(2*u*u+u2),-(gamma-1)*u*v,-(gamma-1)*u*w,gamma*u};
    const DenseMatrix<5,5> Jac_y = {0,0,1,0,0,  
                    -u*v,v,u,0,0,
                    -v*v+0.5*(gamma-1)*u2, -(gamma-1)*u, 2*v-(gamma-1)*v,-(gamma-1)*w,gamma-1,
                    -v*w,0,w,v,0,
                    v*(-gamma*e+(gamma-1)*u2),-(gamma-1)*u*v,gamma*e-0.5*(gamma-1)*(2*v*v+u2),-(gamma-1)*v*w,gamma*v};
    const DenseMatrix<5,5> Jac_z = {0,0,0,1,0,  
                    -u*w,w,0,u,0,
                    -v*w,0,w,v,0,
                    -w*w+0.5*(gamma-1)*u2,-(gamma-1)*u,-(gamma-1)*v,2*w-(gamma-1)*w,gamma-1,
                    w*(-gamma*e+(gamma-1)*u2),-(gamma-1)*u*w,-(gamma-1)*v*w,gamma*e-0.5*(gamma-1)*(2*w*w+u2),gamma*w};
    return {Jac_x, Jac_y, Jac_z};
}



template<uInt GammaNumerator, uInt GammaDenominator>
Scalar EulerPhysicalFlux<GammaNumerator, GammaDenominator>::computeWaveSpeed(
    const DenseMatrix<5,1>& U_L,
    const DenseMatrix<5,1>& U_R) 
{
    return std::max(
        waveSpeedMagnitude(U_L),
        waveSpeedMagnitude(U_R)
    );
}


template<uInt GammaNumerator, uInt GammaDenominator>
Scalar EulerPhysicalFlux<GammaNumerator, GammaDenominator>::positive(Scalar val) { 
    return std::max(val, epslion); 
}


template<uInt GammaNumerator, uInt GammaDenominator>
Scalar EulerPhysicalFlux<GammaNumerator, GammaDenominator>::enthalpy(const DenseMatrix<5,1>& U) {
    return (U[4] + computePressure(U)) / positive(U[0]);
}


template<uInt GammaNumerator, uInt GammaDenominator>
Scalar EulerPhysicalFlux<GammaNumerator, GammaDenominator>::soundSpeed(const DenseMatrix<5,1>& U) {
    return std::sqrt(gamma * positive(computePressure(U) / positive(U[0])));
}


template<uInt GammaNumerator, uInt GammaDenominator>
Scalar EulerPhysicalFlux<GammaNumerator, GammaDenominator>::waveSpeedMagnitude(const DenseMatrix<5,1>& U) {
    const Scalar rho = positive(U[0]);
    const Scalar u = U[1]/rho, v = U[2]/rho, w = U[3]/rho, e = U[4]/rho;
    return std::sqrt(u*u+v*v+w*w) + soundSpeed(U);
}


template<uInt GammaNumerator, uInt GammaDenominator>
Scalar EulerPhysicalFlux<GammaNumerator, GammaDenominator>::computePressure(const DenseMatrix<5,1>& U) {
    const Scalar rho = positive(U[0]);
    const Scalar u = U[1]/rho, v = U[2]/rho, w = U[3]/rho, e = U[4]/rho;
    return (gamma-1)*(rho*e - 0.5*rho*(u*u + v*v + w*w));
}