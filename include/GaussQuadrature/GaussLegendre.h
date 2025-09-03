#include "base/Type.h"


template <uInt Order, typename DefaultQuad>
struct AutoQuadSelector;

namespace GaussLegendreQuad {

    struct Middle {
        public:
            static constexpr int num_points = 1;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {0., 0.}
            }};
            static constexpr std::array<Scalar, num_points> weights = {4.0};
    };

    #include "QuadQuadLowerOrder.h"
    
};

namespace GaussLegendreTri {
    
    struct Auto {
        public:
            static constexpr int num_points = 0;
    };
    
    struct Middle {
        public:
            static constexpr int num_points = 1;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {1.0/3, 1.0/3}
            }};
            static constexpr std::array<Scalar, num_points> weights = {0.5};
    };

    #include "QuadTriLowerOrder.h"
    #include "QuadTriHighOrder.h"

}

namespace GaussLegendreTet {

    struct Auto {
        public:
            static constexpr int num_points = 0;
    };
    struct Middle {
        public:
            static constexpr int num_points = 1;
            static constexpr std::array<std::array<Scalar,3>, num_points> points = {{
                {1.0/4, 1.0/4, 1.0/4}
            }};
            static constexpr std::array<Scalar, num_points> weights = {1.0/6};
    };
    
    #include "QuadTetLowerOrder.h"
    #include "QuadTetHighOrder.h"
}



template <uInt Order, typename DefaultQuad>
struct AutoQuadSelector {
    using type = DefaultQuad;
};

template <uInt Order>
struct AutoQuadSelector<Order, GaussLegendreTet::Auto> {
    using type = GaussLegendreTet::Degree20Points448;
};

template <uInt Order>
struct AutoQuadSelector<Order, GaussLegendreTri::Auto> {
    using type = GaussLegendreTri::Degree23Points100;
};

#include "AutoQuadSelectorTet.h"
#include "AutoQuadSelectorTri.h"


