#pragma once
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <vector>
#include <map>
#include <span>
#include <list>
#include <array>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <bitset>
#include <iomanip>
#include <variant>
#include <optional>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <omp.h>



// typedef double Scalar;
using Scalar = double;
using Scalar64 = double;
using Scalar128 = long double;
constexpr bool use_AVX2 = 1;
typedef int64_t Int;
typedef uint64_t uInt;

typedef std::array<Scalar,2> vector2f;
typedef std::array<Scalar,3> vector3f;
typedef std::array<Scalar,4> vector4f;
typedef std::array<uInt,2> vector2u;
typedef std::array<uInt,3> vector3u;
typedef std::array<uInt,4> vector4u;
typedef std::array<uInt,5> vector5u;
typedef std::array<uInt,6> vector6u;
typedef std::array<uInt,8> vector8u;
// typedef std::vector<Scalar> vector3;


template<typename T>
inline void debug(T s){
    std::cerr<<s<<std::endl;
}

// template<>
template<typename T,uInt N>
inline void debug(std::array<T,N> s, uInt stride = 12){
    for(auto&p:s)
    std::cerr<<std::setw(stride)<<p<<"\t";
    std::cerr<<std::endl;
}


template<typename T>
inline void print(T s){
    std::cout<<s<<std::endl;
}

// template<>
template<typename T,uInt N>
inline void print(std::array<T,N> s, uInt stride = 12){
    for(auto&p:s)
    std::cout<<std::setw(stride)<<p<<"  ";
    std::cout<<std::endl;
}



#include "vector_op.h"
#include "GaussQuadrature/GaussLegendre.h"
#include "tool.h"
