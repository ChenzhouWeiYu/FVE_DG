#pragma once
#include "base/Type.h"
#include "Mesh/Mesh.h"


template<typename Type>
Type rho_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type u_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type v_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type w_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type p_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type e_xyz(Type x, Type y, Type z, Type t = 0.0);

template<typename Type>
Type frho_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type fu_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type fv_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type fw_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type fp_xyz(Type x, Type y, Type z, Type t = 0.0);
template<typename Type>
Type fe_xyz(Type x, Type y, Type z, Type t = 0.0);

#define Filed_Func(filedname) \
Scalar filedname##_xyz(const vector3f& xyz, Scalar t = 0.0);\
Scalar filedname##_Xi(const CompTetrahedron& cell, const vector3f& Xi, Scalar t = 0.0);\
Scalar rho##filedname##_Xi(const CompTetrahedron& cell, const vector3f& Xi, Scalar t = 0.0);

Filed_Func(rho);
Filed_Func(u);
Filed_Func(v);
Filed_Func(w);
Filed_Func(p);
Filed_Func(e);
Filed_Func(frho);
Filed_Func(fu);
Filed_Func(fv);
Filed_Func(fw);
Filed_Func(fp);
Filed_Func(fe);
#undef Filed_Func

DenseMatrix<5,1> U_Xi(const CompTetrahedron& cell, const vector3f& Xi, Scalar t);
DenseMatrix<3,1> uvw_Xi(const CompTetrahedron& cell, const vector3f& Xi, Scalar t);