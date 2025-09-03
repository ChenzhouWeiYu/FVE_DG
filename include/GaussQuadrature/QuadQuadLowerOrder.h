    struct Degree3Points4 {
        private:
            static constexpr Scalar M_SQRT3 = 1.732050807568877293527446341505872366943;
            static constexpr Scalar M_SQRT1_3 = 0.5773502691896257645091487805019574556476;
        public:
            static constexpr int num_points = 4;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {-M_SQRT1_3, -M_SQRT1_3}, 
                { M_SQRT1_3, -M_SQRT1_3},
                { M_SQRT1_3,  M_SQRT1_3},
                {-M_SQRT1_3,  M_SQRT1_3}
            }};
            static constexpr std::array<Scalar, num_points> weights = {1.0, 1.0, 1.0, 1.0};
    };

    struct Degree5Points7 {
        private:
            static constexpr Scalar r = 0.6831300510639732255480692453680701327157;
            static constexpr Scalar s = 0.8906544217808369920783964022053703474661;
            static constexpr Scalar t = 0.3742566422865147407211610155677170660217;
        public:
            static constexpr int num_points = 7;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {  0,  0}, 
                {  r,  r},
                { -r, -r},
                {  s, -t},
                { -s,  t},
                {  t, -s},
                { -t,  s},

            }};
            static constexpr std::array<Scalar, num_points> weights = {
                4 * 2.0/7.0, 
                4 * 25.0/168.0, 
                4 * 25.0/168.0, 
                4 * 5.0/48.0, 
                4 * 5.0/48.0, 
                4 * 5.0/48.0, 
                4 * 5.0/48.0
                };
    };

    struct Degree7Points12 {
        private:
            static constexpr Scalar r = 0.9258200997725514615665667765839995225293;
            static constexpr Scalar s = 0.3805544332083156563791063590863941355001;
            static constexpr Scalar t = 0.8059797829185987437078561813507442463004;
            static constexpr Scalar Bs = 0.1301482291668486142849798580116827915066;
            static constexpr Scalar Bt = 0.05935794367265755855452631482782338133293;
        public:
            static constexpr int num_points = 12;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {  r,  0}, 
                { -r,  0},
                {  0,  r}, 
                {  0, -r},
                {  s,  s}, 
                { -s,  s},
                {  s, -s}, 
                { -s, -s},
                {  t,  t}, 
                { -t,  t},
                {  t, -t}, 
                { -t, -t}
            }};
            static constexpr std::array<Scalar, num_points> weights = {
                4 * 49.0/810.0, 
                4 * 49.0/810.0, 
                4 * 49.0/810.0, 
                4 * 49.0/810.0, 
                4 * Bs, 
                4 * Bs, 
                4 * Bs, 
                4 * Bs, 
                4 * Bt, 
                4 * Bt, 
                4 * Bt, 
                4 * Bt
                };
    };