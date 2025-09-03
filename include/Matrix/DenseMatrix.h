#pragma once
#include "AVX2.h"
#include "base/Type.h"

template <uInt M, uInt N, uInt SubM, uInt SubN>
class MatrixView;

template <uInt M, uInt N>
class DenseMatrix {
private:
    static constexpr uInt Rows = M;
    static constexpr uInt Cols = N;
    static constexpr uInt Size = M*N;
    
    alignas(32) std::array<Scalar, Size> data; // 行优先存储

public:
    // 访问元素
    Scalar& operator()(uInt row, uInt col) { return data[row*N + col]; }
    const Scalar& operator()(uInt row, uInt col) const { return data[row*N + col]; }

    Scalar& operator[](uInt idx) { return data[idx]; }
    const Scalar& operator[](uInt idx) const { return data[idx]; }

    Scalar* data_ptr() { return data.data(); }
    const Scalar* data_ptr() const { return data.data(); }

    DenseMatrix(){for(auto&v:data)v=0.0;};

    const uInt size() const{return Size;}
    const uInt rows() const{return Rows;}
    const uInt cols() const{return Cols;}
    
    template <uInt SubM, uInt SubN>
    MatrixView<M,N,SubM,SubN> View(uInt sr, uInt sc) {
        return MatrixView<M,N,SubM,SubN>(*this, sr, sc);
    }

    template <uInt SubM, uInt SubN>
    DenseMatrix<SubM,SubN> SubMat(uInt sr, uInt sc) {
        return MatrixView<M,N,SubM,SubN>(*this, sr, sc);
    }

    // template <uInt SubM, uInt SubN>
    // DenseMatrix<SubM,SubN> SubMat(uInt sr, uInt sc) const {
    //     return MatrixView<M,N,SubM,SubN>(*this, sr, sc);
    // }

    // ================= 赋值操作 ================= 
    DenseMatrix& operator=(Scalar val) {
        std::fill(data.begin(), data.end(), val);
        return *this;
    }
    operator std::array<Scalar, Size>() const {
        return data; // data 是 std::array，可以直接返回
    }

    // 转置、迹、LU分解等操作
    DenseMatrix<N,M> transpose() const{
        DenseMatrix<N,M> result;
        for(uInt i=0;i<M;i++) for(uInt j=0;j<N;j++) result(j,i)=(*this)(i,j);
        return result;
    }
    Scalar trace() const{
        if constexpr (M==N){
            Scalar result = 0.0;
            for(uInt i=0;i<M;i++) result+=(*this)(i,i);
            return result;
        }
        else return 0.0;
    }

    // LU分解: 无主元，返回一个新矩阵，内部存储 L 和 U
    template<uInt MM = M, uInt NN = N,typename = std::enable_if_t<MM == NN>>
    DenseMatrix<M,N> lu() const {
        DenseMatrix<M,N> LU = *this;
        for (uInt k = 0; k < N; ++k) {
            for (uInt i = k+1; i < N; ++i) {
                LU(i,k) /= LU(k,k);
                for (uInt j = k+1; j < N; ++j) {
                    LU(i,j) -= LU(i,k) * LU(k,j);
                }
            }
        }
        return LU;
    }

    // 自动进行LU分解后求解
    template<uInt MM = M, uInt NN = N,typename = std::enable_if_t<MM == NN>>
    DenseMatrix<M,1> solve(const DenseMatrix<M,1>& b) const {
        DenseMatrix<M,N> LU = lu();
        return solve(b, LU);
    }
    // 利用LU求解方程 Ax = b
    template<uInt MM = M, uInt NN = N,typename = std::enable_if_t<MM == NN>>
    DenseMatrix<M,1> solve(const DenseMatrix<M,1>& b, const DenseMatrix<M,N>& LU) const{
        DenseMatrix<M,1> y, x;

        // 前向代入 Ly = b
        for (uInt i = 0; i < N; ++i) {
            y(i,0) = b(i,0);
            for (uInt j = 0; j < i; ++j) {
                y(i,0) -= LU(i,j) * y(j,0);
            }
        }

        // 后向代入 Ux = y
        for (int i = N-1; i >= 0; --i) {
            x(i,0) = y(i,0);
            for (uInt j = i+1; j < N; ++j) {
                x(i,0) -= LU(i,j) * x(j,0);
            }
            x(i,0) /= LU(i,i);
        }
        return x;
    }
    // 计算逆矩阵: 自动LU分解
    template<uInt MM = M, uInt NN = N,typename = std::enable_if_t<MM == NN>>
    DenseMatrix<M,N> inverse() const{
        DenseMatrix<M,N> LU = lu();
        return inverse(LU);
    }
    // 已有LU分解，计算逆矩阵
    template<uInt MM = M, uInt NN = N,typename = std::enable_if_t<MM == NN>>
    DenseMatrix<M,N> inverse(const DenseMatrix<M,N>& LU) const{
        DenseMatrix<M,N> inv;
        DenseMatrix<M,1> e, col;
        for (uInt i = 0; i < N; ++i) {
            // 生成单位向量 e_i
            for (uInt j = 0; j < N; ++j) e(j,0) = (i == j ? 1 : 0);
            col = solve(e, LU);
            // 把 col 写入逆矩阵第 i 列
            for (uInt j = 0; j < N; ++j) {
                inv(j,i) = col(j,0);
            }
        }
        return inv;
    }

    // DenseMatrix& operator=(const std::array<Scalar, Size>& rhs) {
    //     for(uInt i=0; i<M*N; ++i)
    //         data[i] = rhs[i];
    //     return *this;
    // }
    DenseMatrix(std::initializer_list<Scalar> init) {
        for(auto&v:data)v=0.0;
        std::copy(init.begin(), init.end(), data.begin());
    }
    DenseMatrix(const std::array<Scalar, Size>& init) {
        for(auto&v:data)v=0.0;
        std::copy(init.begin(), init.end(), data.begin());
    }
    DenseMatrix(const std::array<std::array<Scalar,N>, M>& init) {
        for(uInt i = 0; i < M; i++)
        for(uInt j = 0; j < N; j++)
        data[i*N+j] = init[i][j];
        // std::copy(init.begin(), init.end(), data.begin());
    }

    DenseMatrix& operator=(const DenseMatrix& rhs){
            // debug("222");
        for(uInt i=0; i<M*N; ++i)
            data[i] = rhs[i];
            
            // debug("333");
        return *this;
    }


    // ================= 逐元素运算 =================
    // 宏实现 +-*/ 的快速替换 
    #define ELEMENT_WISE_OP(op)\
    DenseMatrix& operator op##=(const DenseMatrix& rhs) {\
        for(uInt i=0; i<M*N; ++i) data[i] op##= rhs.data[i];\
        return *this;\
    }\
    DenseMatrix operator op(const DenseMatrix& rhs) const {\
        DenseMatrix result(*this);\
        return result op##= rhs;\
    }\
    DenseMatrix& operator op##=(Scalar val) {\
        for(auto& d : data) d op##= val;\
        return *this;\
    }\
    friend DenseMatrix operator op(Scalar val, const DenseMatrix& rhs) { \
        DenseMatrix res(rhs); \
        for (auto& v : res.data) v = val op v; \
        return res; \
    }\
    DenseMatrix operator op(Scalar val) const {\
        DenseMatrix result(*this);\
        return result op##= val;\
    }
    
    ELEMENT_WISE_OP(+)
    ELEMENT_WISE_OP(-)
    ELEMENT_WISE_OP(*)
    ELEMENT_WISE_OP(/)
    #undef ELEMENT_WISE_OP
    DenseMatrix operator -() const {
        return Zeros() - (*this);
    }
    

    // ================= 矩阵乘法 =================
    // 矩阵乘法运算符
    template <uInt K>
    DenseMatrix<M, K> multiply(const DenseMatrix<N, K>& rhs) const {
        DenseMatrix<M, K> result;

        if constexpr (std::is_same_v<Scalar,double> && use_AVX2) {
            if constexpr (M == N && N == K && (M & (M - 1)) == 0) {
                // 仅当尺寸为2的幂时使用Strassen
                avx2::strassen_matrix_mult<M>(this->data_ptr(), rhs.data_ptr(), result.data_ptr());
            } else {
                avx2::matrix_mult<M, N, K>(this->data_ptr(), rhs.data_ptr(), result.data_ptr());
            }
        }
        else{
            for (uInt k = 0; k < N; k++) {
                for (uInt i = 0; i < M; i++) {
                    for (uInt j = 0; j < K; j++) {
                        result(i,j) += (*this)(i,k) * rhs(k,j);
                    }
                }
            }
        }
        return result;
    }
    // 矩阵-向量乘法
    DenseMatrix<M, 1> multiply(const DenseMatrix<N, 1>& rhs) const {
        DenseMatrix<M, 1> result;
        if constexpr (std::is_same_v<Scalar,double> && use_AVX2) {
            avx2::matrix_vector_mult<M, N>(this->data_ptr(), rhs.data_ptr(), result.data_ptr());
        }
        else{
            for (uInt i = 0; i < M; i++) {
                for (uInt k = 0; k < N; k++) {
                    result(i,0) += (*this)(i,k) * rhs(k,0);
                }
            }
        }
        return result;
    }
    DenseMatrix<M, 1> multiply(const std::array<Scalar,N>& rhs) const {
        DenseMatrix<M, 1> result;
        if constexpr (std::is_same_v<Scalar,double> && use_AVX2 && N>=4) {
            avx2::matrix_vector_mult<M, N>(this->data_ptr(), rhs.data(), result.data_ptr());
        }
        else{
            for (uInt i = 0; i < M; i++) {
                for (uInt k = 0; k < N; k++) {
                    result(i,0) += (*this)(i,k) * rhs[k];
                }
            }
        }
        return result;
    }
    // 内积运算
    Scalar dot(const DenseMatrix<M, 1>& rhs) const {
        return vec_dot(data, rhs.data);
    }

    Scalar length() const {
        return vec_length(data);
    }

    Scalar norm() const {
        return vec_length(data);
    }

    // 流输出
    friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& mat) {
        for(uInt i=0; i<M; ++i) {
            for(uInt j=0; j<N; ++j)
                os << mat(i,j) << " ";
            os << "\n";
        }
        return os;
    }

private:
    

public:
    // ================= 数学函数 =================
    #define ELEMENT_WISE_FUNC(FUNC)\
    friend DenseMatrix FUNC(const DenseMatrix& mat) {\
        DenseMatrix result;\
        for(uInt i=0; i<M*N; ++i)\
            result.data[i] = std::FUNC(mat.data[i]);\
        return result;\
    }
    ELEMENT_WISE_FUNC(exp)
    ELEMENT_WISE_FUNC(log)
    ELEMENT_WISE_FUNC(sin)
    ELEMENT_WISE_FUNC(cos)
    ELEMENT_WISE_FUNC(tan)
    // ELEMENT_WISE_FUNC(max)
    // ELEMENT_WISE_FUNC(min)
    ELEMENT_WISE_FUNC(abs)
    ELEMENT_WISE_FUNC(tanh)
    #undef ELEMENT_WISE_FUNC

    friend DenseMatrix pow(const DenseMatrix& mat, Scalar beta) {
        DenseMatrix result;
        for(uInt i=0; i<M*N; ++i)
            result.data[i] = std::pow(mat.data[i], beta);
        return result;
    }

    friend Scalar max(const DenseMatrix& mat) {
        Scalar m;
        for(const auto& b : mat.data)
            m = m > b? m : b;
        return m;
    }

    // 类似实现log, sqrt, sin等...

    // ================= 特殊矩阵生成 =================
    static DenseMatrix Identity() {
        static_assert(M == N, "Identity matrix must be square");
        DenseMatrix mat;
        for(uInt i=0; i<M; ++i)
            mat(i,i) = 1.0;
        return mat;
    }
    static DenseMatrix Diag(DenseMatrix<N,1> diag) {
        static_assert(M == N, "Identity matrix must be square");
        DenseMatrix mat;
        for(uInt i=0; i<M; ++i)
            mat(i,i) = diag(i,0);
        return mat;
    }
    static DenseMatrix Diag(std::array<Scalar,N> diag) {
        static_assert(M == N, "Identity matrix must be square");
        DenseMatrix mat;
        for(uInt i=0; i<M; ++i)
            mat(i,i) = diag[i];
        return mat;
    }
    static DenseMatrix Ones() {
        // static_assert(M == N, "Identity matrix must be square");
        DenseMatrix mat;
        for(uInt i=0; i<Size; ++i)
            mat[i] = 1.0;
            // debug("222");
        return mat;
    }
    static DenseMatrix Zeros() {
        // static_assert(M == N, "Identity matrix must be square");
        DenseMatrix mat;
        for(uInt i=0; i<Size; ++i)
            mat[i] = 0.0;
            // debug("444");
        return mat;
    }
    static DenseMatrix Random() {
        // static_assert(M == N, "Random matrix must be square");
        DenseMatrix mat;
        for(uInt i=0; i<Size; ++i)
            mat[i] = 1.0;
        return mat;
    }

};





// ================= 矩阵视图（零拷贝子矩阵） =================
template <uInt M, uInt N, uInt SubM, uInt SubN>
class MatrixView {
    DenseMatrix<M,N>& parent;
    uInt start_row;
    uInt start_col;

public:
    MatrixView(DenseMatrix<M,N>& mat, uInt sr, uInt sc) 
        : parent(mat), start_row(sr), start_col(sc) 
    {
        if(sr + SubM > M || sc + SubN > N)
            throw std::out_of_range("Submatrix out of bounds");
    }

    MatrixView(const DenseMatrix<M,N>& mat, uInt sr, uInt sc) 
        : parent(mat), start_row(sr), start_col(sc) 
    {
        if(sr + SubM > M || sc + SubN > N)
            throw std::out_of_range("Submatrix out of bounds");
    }

    Scalar& operator()(uInt i, uInt j) {
        return parent(start_row + i, start_col + j);
    }

    Scalar& operator[](uInt idx) {
        return (*this)(idx/SubN, idx%SubN);
    }

    const Scalar& operator()(uInt i, uInt j) const {
        return parent(start_row + i, start_col + j);
    }

    const Scalar& operator[](uInt idx) const {
        return (*this)(idx/SubN, idx%SubN);
    }

    template <uInt K>
    auto operator*(const MatrixView<N,K,SubN,K>& rhs) {
        DenseMatrix<SubM,K> result;
        // 优化后的子矩阵乘法...
        return result;  
    }
        // ================= 逐元素运算 =================
    // 宏实现 +-*/ 的快速替换 
    #define ELEMENT_WISE_OP(op)\
    MatrixView& operator op##=(const DenseMatrix<SubM,SubN>& rhs) {\
        for(uInt i=0; i<SubM; ++i){\
            for(uInt j=0; j<SubN; ++j){\
                parent(start_row + i, start_col + j) op##= rhs(i,j);\
            }\
        } \
        return *this;\
    }\
    DenseMatrix<SubM,SubN> operator op(const DenseMatrix<SubM,SubN>& rhs) {\
        DenseMatrix<SubM,SubN> result;\
        for(uInt i=0; i<SubM; ++i){\
            for(uInt j=0; j<SubN; ++j){\
                result(i,j) = parent(start_row + i, start_col + j) op rhs(i,j);\
            }\
        } \
        return result;\
    }\
    DenseMatrix<SubM,SubN> operator op(const MatrixView& rhs) {\
        DenseMatrix<SubM,SubN> result;\
        for(uInt i=0; i<SubM; ++i){\
            for(uInt j=0; j<SubN; ++j){\
                result(i,j) = parent(start_row + i, start_col + j) op rhs.parent(start_row + i, start_col + j);\
            }\
        } \
        return result;\
    }\
    const DenseMatrix<SubM,SubN>& operator op(const DenseMatrix<SubM,SubN>& rhs) const {\
        DenseMatrix<SubM,SubN> result;\
        for(uInt i=0; i<SubM; ++i){\
            for(uInt j=0; j<SubN; ++j){\
                result(i,j) = parent(start_row + i, start_col + j) op rhs(i,j);\
            }\
        } \
        return result;\
    }\
    const DenseMatrix<SubM,SubN>& operator op(const MatrixView& rhs) const {\
        DenseMatrix<SubM,SubN> result;\
        for(uInt i=0; i<SubM; ++i){\
            for(uInt j=0; j<SubN; ++j){\
                result(i,j) = parent(start_row + i, start_col + j) op rhs.parent(start_row + i, start_col + j);\
            }\
        } \
        return result;\
    }

    MatrixView& operator =(const DenseMatrix<SubM,SubN>& rhs) {
        for(uInt i=0; i<SubM; ++i){
            for(uInt j=0; j<SubN; ++j){
                parent(start_row + i, start_col + j) = rhs(i,j);
            }
        } 
        return *this;
    }

    operator DenseMatrix<SubM,SubN>() const {
        DenseMatrix<SubM,SubN> result;
        for(uInt i=0; i<SubM; ++i){
            for(uInt j=0; j<SubN; ++j){
                result(i,j) = parent(start_row + i, start_col + j);
            }
        } 
        return result;
    }
    static DenseMatrix<SubM,SubN> Sub(const DenseMatrix<M,N>& mat, uInt sr, uInt sc) {
        DenseMatrix<SubM,SubN> result;
        for(uInt i=0; i<SubM; ++i){
            for(uInt j=0; j<SubN; ++j){
                result(i,j) = mat(sr + i, sc + j);
            }
        } 
        return result;
    }
    
    ELEMENT_WISE_OP(+)
    ELEMENT_WISE_OP(-)
    ELEMENT_WISE_OP(*)
    ELEMENT_WISE_OP(/)
    #undef ELEMENT_WISE_OP
};

// ================= 使用示例 =================
/*
DenseMatrix<64,64> a = DenseMatrix<64,64>::Identity();
DenseMatrix<64,64> b = DenseMatrix<64,64>::Ones();

// 自动选择Strassen算法
auto c = a * b; 

// 子矩阵操作
MatrixView<64,64,32,32> sub_a(a, 0, 0);
MatrixView<64,64,32,32> sub_b(b, 16, 16);
auto sub_c = sub_a * sub_b;

// 并行向量运算
LongVector<64> v1, v2;
v1 += v2 * 3.14; // 并行执行
Scalar ip = v1.dot(v2); // SIMD+OpenMP优化
*/

template<>      // 对类模板
template<>      // 对成员函数模板
inline DenseMatrix<3,3> DenseMatrix<3,3>::inverse() const {
    DenseMatrix<3,3> ret;
    const DenseMatrix<3,3>& A = *this;
    Scalar det = 0;
    ret(0,0) = A(1,1)*A(2,2) - A(1,2)*A(2,1);
    ret(1,0) = A(1,2)*A(2,0) - A(1,0)*A(2,2);
    ret(2,0) = A(1,0)*A(2,1) - A(1,1)*A(2,0);
    ret(0,1) = A(2,1)*A(0,2) - A(2,2)*A(0,1);
    ret(1,1) = A(2,2)*A(0,0) - A(2,0)*A(0,2);
    ret(2,1) = A(2,0)*A(0,1) - A(2,1)*A(0,0);
    ret(0,2) = A(0,1)*A(1,2) - A(0,2)*A(1,1);
    ret(1,2) = A(0,2)*A(1,0) - A(0,0)*A(1,2);
    ret(2,2) = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    det += A(0,0)*ret(0,0);
    det += A(1,0)*ret(0,1);
    det += A(2,0)*ret(0,2);
    ret /= det;
    return ret;
}