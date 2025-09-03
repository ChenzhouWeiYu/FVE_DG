#pragma once

#include "FaceType.h"

//---------------- 六面体单元 ----------------
struct Hexahedron {
    vector8u nodes;
    vector6u faces;

    void reorder(const std::vector<vector3f>& points, 
                const std::vector<GeneralFace>& all_faces) {
        // Step 1: 提取所有四边形面的节点
        std::array<vector4u, 6> face_nodes;
        for (uInt i = 0; i < 6; ++i) {
            face_nodes[i] = std::get<QuadFace>(all_faces[faces[i]]).nodes;
        }

        // Step 2: 识别三对相对面
        std::array<std::pair<uInt, uInt>, 3> opposite_pairs = find_opposite_pairs(face_nodes);

        // Step 3: 确定底面（平均Z最低的面）
        uInt bottom_pair_idx = find_bottom_pair(opposite_pairs, face_nodes, points);
        auto& bottom_pair = opposite_pairs[bottom_pair_idx];
        
        // Step 4: 排列底面节点（极角排序）
        vector4u base_nodes = sort_quad_by_polar(face_nodes[bottom_pair.first], points);
        vector4u top_nodes = align_top_nodes(base_nodes, face_nodes[bottom_pair.second], points);

        // Step 5: 构建六面体节点顺序
        nodes = {
            base_nodes[0], base_nodes[1], base_nodes[2], base_nodes[3],
            top_nodes[0],  top_nodes[1],  top_nodes[2],  top_nodes[3]
        };

        // Step 6: 重建面顺序 [底面, 顶面, 前, 右, 后, 左]
        rebuild_face_order(base_nodes, top_nodes, face_nodes);
    }

private:
    // 寻找三对相对面（无共享节点）
    std::array<std::pair<uInt, uInt>, 3> find_opposite_pairs(const std::array<vector4u,6>& faces) {
        std::array<std::pair<uInt, uInt>, 3> pairs;
        std::unordered_set<uInt> matched;
        
        uInt pair_idx = 0;
        for (uInt i = 0; i < 6; ++i) {
            if (matched.count(i)) continue;
            for (uInt j = i+1; j < 6; ++j) {
                if (is_opposite(faces[i], faces[j])) {
                    pairs[pair_idx++] = {i, j};
                    matched.insert(i);
                    matched.insert(j);
                    break;
                }
            }
        }
        return pairs;
    }

    // 判断是否为相对面
    bool is_opposite(const vector4u& a, const vector4u& b) const {
        std::unordered_set<uInt> s(a.begin(), a.end());
        for (uInt n : b) if (s.count(n)) return false;
        return true;
    }

    // 通过平均Z坐标找到最下面的面对
    uInt find_bottom_pair(const std::array<std::pair<uInt, uInt>,3>& pairs,
                         const std::array<vector4u,6>& faces,
                         const std::vector<vector3f>& points) const {
        uInt min_idx = 0;
        Scalar min_z = std::numeric_limits<Scalar>::max();
        
        for (uInt i = 0; i < 3; ++i) {
            Scalar z_sum = 0;
            for (uInt n : faces[pairs[i].first]) z_sum += points[n][2];
            if (z_sum < min_z) {
                min_z = z_sum;
                min_idx = i;
            }
        }
        return min_idx;
    }

    // 极角排序四边形节点（绕几何中心逆时针）
    vector4u sort_quad_by_polar(vector4u quad, const std::vector<vector3f>& points) {
        // 计算几何中心
        vector3f center = {0,0,0};
        for (uInt n : quad) center += points[n];
        center /= 4;

        // 按极角排序
        std::sort(quad.begin(), quad.end(), [&](uInt a, uInt b) {
            vector3f va = points[a] - center;
            vector3f vb = points[b] - center;
            return std::atan2(va[1], va[0]) < std::atan2(vb[1], vb[0]);
        });
        // 确保逆时针（通过顶点0到1到2的向量叉积）
        vector3f v01 = points[quad[1]] - points[quad[0]];
        vector3f v02 = points[quad[2]] - points[quad[0]];
        if (v01[0]*v02[1] - v01[1]*v02[0] < 0) { // 仅需二维判断
            std::swap(quad[1], quad[3]);
        }
        return quad;
    }

    // 对齐顶面节点到底面
    vector4u align_top_nodes(const vector4u& base, const vector4u& top_raw,
                            const std::vector<vector3f>& points) {
        vector4u top = top_raw;
        // 建立顶点映射：每个顶面节点对应最近的底面节点
        std::array<uInt,4> mapping;
        for (uInt i = 0; i < 4; ++i) {
            Scalar min_dist = std::numeric_limits<Scalar>::max();
            for (uInt j = 0; j < 4; ++j) {
                Scalar dx = points[base[i]][0] - points[top[j]][0];
                Scalar dy = points[base[i]][1] - points[top[j]][1];
                Scalar dist = dx*dx + dy*dy; // 仅考虑XY平面距离
                if (dist < min_dist) {
                    min_dist = dist;
                    mapping[i] = j;
                }
            }
        }
        // 重新排列顶面节点
        return {top[mapping[0]], top[mapping[1]], top[mapping[2]], top[mapping[3]]};
    }

    // 根据节点顺序重建面索引
    void rebuild_face_order(const vector4u& base, const vector4u& top,
                           std::array<vector4u,6>& face_nodes) {
        // 预期的标准面连接模式
        const std::array<vector4u,6> expected_faces = {{
            {base[0], base[1], base[2], base[3]},  // 底面
            {top[0],  top[1],  top[2],  top[3]},   // 顶面
            {base[0], base[1], top[1],  top[0]},   // 前
            {base[1], base[2], top[2],  top[1]},   // 右
            {base[2], base[3], top[3],  top[2]},   // 后
            {base[3], base[0], top[0],  top[3]}    // 左
        }};

        // 匹配并更新面节点顺序
        for (auto& face : face_nodes) {
            for (const auto& pattern : expected_faces) {
                if (is_cyclic_match(face, pattern)) {
                    face = pattern;
                    break;
                }
            }
        }
    }

    // 循环匹配检查
    bool is_cyclic_match(const vector4u& a, const vector4u& b) const {
        for (uInt offset = 0; offset < 4; ++offset) {
            bool match = true;
            for (uInt i = 0; i < 4; ++i) {
                if (a[i] != b[(i+offset)%4]) {
                    match = false;
                    break;
                }
            }
            if (match) return true;
        }
        return false;
    }
};