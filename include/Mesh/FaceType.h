#pragma once

#include "base/Type.h"

//---------------- 面类型 ----------------
struct TriangleFace {
    vector3u nodes;
    vector2u neighbor_cells = {uInt(-1), uInt(-1)};
    
    void reorder() {
        // 三角形的节点序不需要处理
    }
};

struct QuadFace {
    vector4u nodes;
    vector2u neighbor_cells = {uInt(-1), uInt(-1)};
    vector2u diagonal_nodes = {uInt(-1), uInt(-1)};; // 剖分对角线信息，记录对角节点的全局索引
    
    // 判断四边形顶点顺序是否满足右手法则
    bool is_rhs(const std::vector<vector3f>& points) const {
        const auto& p0 = points[nodes[0]];
        const auto& p1 = points[nodes[1]];
        const auto& p2 = points[nodes[2]];
        const auto& p3 = points[nodes[3]];
        
        // 计算三个边的向量
        const auto v1 = p1-p0;
        const auto v2 = p2-p1;
        const auto v3 = p3-p2;
        
        // 计算两个三角形的法向量
        const auto n1 = vec_cross(v1, v2);
        const auto n2 = vec_cross(v2, v3);
        
        // 法向量应大致同向
        return vec_dot(n1, n2) > 0;
    }

    // 重新排序四边形顶点保证右手法则和对角线正确
    void reorder(const std::vector<vector3f>& points) {
        if (is_rhs(points)) return;
        
        // 不满足右手法则时调整顶点顺序
        // 方案1：交换最后两个顶点
        std::swap(nodes[2], nodes[3]);
        
        // 二次验证
        if (!is_rhs(points)) {
            // 方案2：逆序整个四边形
            std::reverse(nodes.begin()+1, nodes.end());
        }
    }

    // 选择最短的对角线进行剖分
    void split_diagonal(const std::vector<vector3f>& points) {
        const auto& a = points[nodes[0]];
        const auto& b = points[nodes[1]];
        const auto& c = points[nodes[2]];
        const auto& d = points[nodes[3]];
        
        // 比较对角线长度
        const Scalar len_ac = distance(a, c);
        const Scalar len_bd = distance(b, d);
        
        if (len_ac < len_bd) {
            diagonal_nodes = {nodes[0], nodes[2]}; // 对角线 a-c
        } else {
            diagonal_nodes = {nodes[1], nodes[3]}; // 对角线 b-d
        }
    }
};

using GeneralFace = std::variant<TriangleFace, QuadFace>;
