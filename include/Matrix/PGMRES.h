// 文件名：PGMRES.hpp
#pragma once
#include "Matrix/DenseMatrix.h"
#include "Matrix/SparseMatrix.h"
#include "Matrix/LongVector.h"
#include "Matrix/submat.h"
#include "Matrix/block_precondition.h"
#include "Matrix/simple_precondition.h"
#include "Matrix/fve_precondition.h"
#include "Matrix/pmg_precondition.h"

// 模板类实现 PGMRES
template<uInt BlockDim, uInt Max_dims, bool output_flag = false>
class PGMRES;



template<uInt BlockDim, uInt Max_dims, bool output_flag>
class PGMRES {
public:
    // using Vec = DenseMatrix<BlockDim, 1>;
    using Vec = LongVector<BlockDim>;

    PGMRES(const BlockSparseMatrix<BlockDim, BlockDim>& A_, Preconditioner<BlockDim>& M_)
        : A(A_), M(M_) {}

    // 主迭代函数
    std::tuple<uInt,Scalar> solve(Vec& x, const Vec& b, uInt max_iter, Scalar tol) {
        std::array<Vec, Max_dims + 1> V{};       // 正交基 V_0, ..., V_m
        std::array<Vec, Max_dims> Z{};           // 预条件后向量 z_j = M^{-1} v_j
        DenseMatrix<Max_dims + 1, Max_dims> H{}; // Hessenberg 矩阵
        DenseMatrix<Max_dims + 1, 1> g{};         // g = beta * e1

        // print("spmv 前");
        Vec r0 = b - A.multiply(x);
        // print("spmv后，norm前");
        Scalar beta = r0.norm();
        // print("初始残差 r_0 = " + std::to_string(beta));
        if (beta < tol) return {0, beta};

        V[0] = r0 / beta;
        g(0, 0) = beta;

        // Givens 旋转参数
        std::array<Scalar, Max_dims> cs{}, sn{};

        for (uInt j = 0; j < max_iter && j < Max_dims; ++j) {
            // 预条件：z_j = M^{-1} v_j
            // print("预条件前");
            Z[j] = M.apply(V[j]);
            if constexpr (output_flag){
                printf("iter = %ld ,  ||Az-v|| residual = %le\n",j, (A.multiply(Z[j])-V[j]).norm()/V[j].norm());
            }

            // print("预条件后");

            // w = A z_j
            Vec w = A.multiply(Z[j]);
            // print("预条件的残差向量");

            // Arnoldi 正交化
            for (uInt i = 0; i <= j; ++i) {
                H(i, j) = V[i].dot(w);
                w -= H(i, j) * V[i];
            }
            H(j + 1, j) = w.norm();
            // print("H矩阵");

            if (H(j + 1, j) < 1e-14) break; // 说明基生成失败（已收敛）

            V[j + 1] = w / H(j + 1, j);
            
            // print("开始 GIvens 旋转");
            // Apply Givens rotations to new column H(:,j)
            
            for (uInt i = 0; i < j; ++i)
                apply_givens(H(i, j), H(i + 1, j), cs[i], sn[i]);

            // Create new Givens rotation
            generate_givens(H(j, j), H(j + 1, j), cs[j], sn[j]);
            apply_givens(H(j, j), H(j + 1, j), cs[j], sn[j]);
            apply_givens(g(j, 0), g(j + 1, 0), cs[j], sn[j]);
            if constexpr (output_flag){
                printf("iter = %ld ,   residual = %le\n",j, g(j + 1, 0)/beta);
            }
            // print("Givens旋转完");
            // Check residual
            if (std::abs(g(j + 1, 0))/beta < tol) {
                reconstruct_solution(x, Z, H, g, j + 1);
                return {j,std::abs(g(j + 1, 0))/beta};
            }
        }

        reconstruct_solution(x, Z, H, g, max_iter);
        uInt real_iter = std::min(max_iter,Max_dims);
        return {real_iter,std::abs(g(real_iter, 0))/beta};
    }

private:
    const BlockSparseMatrix<BlockDim, BlockDim>& A;
    Preconditioner<BlockDim>& M;

    // 生成 Givens 旋转
    void generate_givens(Scalar a, Scalar b, Scalar& c, Scalar& s) {
        if (std::abs(b) < 1e-14) {
            c = 1.0; s = 0.0;
        } else if (std::abs(b) > std::abs(a)) {
            Scalar t = -a / b;
            s = 1.0 / std::sqrt(1 + t * t);
            c = s * t;
        } else {
            Scalar t = -b / a;
            c = 1.0 / std::sqrt(1 + t * t);
            s = c * t;
        }
    }

    // 应用 Givens 旋转到 (a,b)
    void apply_givens(Scalar& a, Scalar& b, Scalar c, Scalar s) {
        Scalar temp = c * a - s * b;
        b = s * a + c * b;
        a = temp;
    }

    // 解最小二乘问题并更新解
    void reconstruct_solution(Vec& x, const std::array<Vec, Max_dims>& Z,
                              const DenseMatrix<Max_dims + 1, Max_dims>& H,
                              const DenseMatrix<Max_dims + 1, 1>& g,
                              uInt k) {
        DenseMatrix<Max_dims, 1> y{};
        // 回代解 y: H(0:k,0:k) * y = g(0:k)
        for (int i = k - 1; i >= 0; --i) {
            y(i, 0) = g(i, 0);
            for (uInt j = i + 1; j < k; ++j)
                y(i, 0) -= H(i, j) * y(j, 0);
            y(i, 0) /= H(i, i);
        }

        for (uInt i = 0; i < k; ++i)
            x += Z[i] * y(i, 0);
    }
};