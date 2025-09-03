#pragma once

#include "FaceType.h"

//==================== 四棱锥单元 ====================//
struct Pyramid {
    vector5u nodes; // [底面0-3, 顶点4]
    vector5u faces; // [底面quad, 侧面tri1-4]

    void reorder(const std::vector<vector3f>& points,
                const std::vector<GeneralFace>& all_faces) {
        // Step 1: 识别底面四边形
        uInt base_face = 0;
        for (uInt i = 0; i < 5; ++i) {
            if (std::holds_alternative<QuadFace>(all_faces[faces[i]])) {
                base_face = i;
                break;
            }
        }

        // Step 2: 调整底面节点顺序
        vector4u base_nodes = std::get<QuadFace>(all_faces[faces[base_face]]).nodes;
        sort_quad_ccw(base_nodes, points);

        // Step 3: 确定顶点
        uInt apex = find_apex(base_nodes, all_faces);

        // Step 4: 重建节点顺序
        nodes = {base_nodes[0], base_nodes[1], base_nodes[2], base_nodes[3], apex};
    }

private:
    // 四边形逆时针排序
    void sort_quad_ccw(vector4u& quad, const std::vector<vector3f>& points) {
        vector3f center = {0,0,0};
        for (uInt n : quad) center += points[n];
        center /= 4;

        std::sort(quad.begin(), quad.end(), [&](uInt a, uInt b) {
            vector3f va = points[a] - center;
            vector3f vb = points[b] - center;
            return std::atan2(va[1], va[0]) < std::atan2(vb[1], vb[0]);
        });

        vector3f v1 = points[quad[1]] - points[quad[0]];
        vector3f v2 = points[quad[2]] - points[quad[0]];
        if (vec_dot(vec_cross(v1, v2), {0,0,1}) < 0) {
            std::reverse(quad.begin(), quad.end());
        }
    }

    // 查找顶点（不属于底面的唯一节点）
    uInt find_apex(const vector4u& base, const std::vector<GeneralFace>& all_faces) const {
        std::unordered_set<uInt> base_set(base.begin(), base.end());
        for (uInt face_id : faces) {
            if (auto* tri = std::get_if<TriangleFace>(&all_faces[face_id])) {
                for (uInt n : tri->nodes) {
                    if (!base_set.count(n)) return n;
                }
            }
        }
        return -1; // Should not reach here
    }
};