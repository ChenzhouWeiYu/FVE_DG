#pragma once

#include "FaceType.h"
#include "HexahedronCell.h"
#include "PrismCell.h"
#include "PyramidCell.h"



struct Tetrahedron {
    vector4u nodes;
    vector4u faces;
    
    void reorder(const std::vector<vector3f>& points, 
                const std::vector<GeneralFace>& all_faces)  {
        // 无需重排面
        auto v1 = points[nodes[1]] - points[nodes[0]];
        auto v2 = points[nodes[2]] - points[nodes[0]];
        auto v3 = points[nodes[3]] - points[nodes[0]];
        if(vec_dot(vec_cross(v1,v2),v3)<0){
            nodes =  {nodes[0],nodes[2],nodes[1],nodes[3]};
        }
    }
};

using GeneralCell = std::variant<Hexahedron, Prism, Pyramid, Tetrahedron>;