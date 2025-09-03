#pragma once

#include "FaceType.h"

//==================== 三棱柱单元 ====================//
struct Prism {
    vector6u nodes; // [底面0,1,2 | 顶面3,4,5]
    vector5u faces; // [底面tri, 顶面tri, 侧面quad1, quad2, quad3]

    void reorder(const std::vector<vector3f>& points, 
                const std::vector<GeneralFace>& all_faces) {
        // Step 1: 提取两个三角形面
        vector2u tri_faces;
        uInt tri_count = 0;
        for (uInt i = 0; i < 5; ++i) {
            if (std::holds_alternative<TriangleFace>(all_faces[faces[i]])) {
                tri_faces[tri_count++] = faces[i];
                if (tri_count == 2) break;
            }
        }

        // Step 2: 确定底面和顶面（通过面索引位置）
        auto& base_face = std::get<TriangleFace>(all_faces[tri_faces[0]]);
        auto& top_face = std::get<TriangleFace>(all_faces[tri_faces[1]]);
        vector3u base_nodes = base_face.nodes;
        vector3u top_nodes = top_face.nodes;

    // std::cout<<"11"<<std::endl;
        // Step 3: 极角排序底面节点（仅需二维投影）
        {
            vector3f center = {0,0,0};
            for (uInt n : base_nodes) center += points[n];
            center /= 3;

            std::sort(base_nodes.begin(), base_nodes.end(), [&](uInt a, uInt b) {
                vector3f va = points[a] - center;
                vector3f vb = points[b] - center;
                return std::atan2(va[1], va[0]) < std::atan2(vb[1], vb[0]);
            });

            // 验证二维法线方向
            vector3f v1 = points[base_nodes[1]] - points[base_nodes[0]];
            vector3f v2 = points[base_nodes[2]] - points[base_nodes[0]];
            if (v1[0]*v2[1] - v1[1]*v2[0] < 0) {
                std::swap(base_nodes[1], base_nodes[2]);
            }
        }

    // std::cout<<"11"<<std::endl;
        // Step 4: 安全处理侧面四边形
        std::array<vector4u,3> quad_nodes;
        uInt quad_count = 0;
        for (uInt i = 0; i < 5; ++i) { // 仅处理侧面面
            if (std::holds_alternative<QuadFace>(all_faces[faces[i]])) { // 关键修复点
                quad_nodes[quad_count++] = std::get<QuadFace>(all_faces[faces[i]]).nodes;
            } else {
                // throw std::runtime_error("Prism has non-quad lateral face");
            }
        }

    // std::cout<<"11"<<std::endl;
        // Step 5: 建立顶面节点映射
        std::array<uInt,3> top_mapping;
        for (uInt i = 0; i < 3; ++i) {
            const auto& quad = quad_nodes[i];
            uInt a = base_nodes[i];
            uInt b = base_nodes[(i+1)%3];
            
            // 查找对应的顶面边
            for (uInt j = 0; j < 3; ++j) {
                uInt c = top_nodes[j];
                uInt d = top_nodes[(j+1)%3];
                if (contains_edge(quad, c, d)) {
                    top_mapping[i] = j;
                    break;
                }
            }
        }

        // Step 6: 对齐顶面节点
        vector3u ordered_top;
        for (uInt i = 0; i < 3; ++i) {
            ordered_top[i] = top_nodes[top_mapping[i]];
        }
        top_nodes = ordered_top;

        // Step 7: 更新节点和面顺序
        nodes = {base_nodes[0], base_nodes[1], base_nodes[2],
                 top_nodes[0],  top_nodes[1],  top_nodes[2]};
        faces = {tri_faces[0], tri_faces[1], faces[2], faces[3], faces[4]};
    }

private:
    // 检查四边形是否包含指定边
    bool contains_edge(const vector4u& quad, uInt a, uInt b) const {
        for (uInt i = 0; i < 4; ++i) {
            if ((quad[i] == a && quad[(i+1)%4] == b) ||
                (quad[i] == b && quad[(i+1)%4] == a)) {
                return true;
            }
        }
        return false;
    }
};