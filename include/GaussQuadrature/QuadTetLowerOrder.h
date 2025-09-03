    struct Degree5Points15 {
        private:
            static constexpr Scalar w0 = 0.0197530864197530864197530864198;
            static constexpr Scalar x0 = 0.25;

            static constexpr Scalar w1 = 0.0119895139631697700017306424850;
            static constexpr Scalar x1 = 0.0919710780527230327888451353005;
            static constexpr Scalar y1 = 0.724086765841830901633464594098;
            
            static constexpr Scalar w2 = 0.0115113678710453975467702393492;
            static constexpr Scalar x2 = 0.319793627829629908387625452935;
            static constexpr Scalar y2 = 0.0406191165111102748371236411957;
            
            static constexpr Scalar w3 = 0.00881834215167548500881834215168;
            static constexpr Scalar x3 = 0.0563508326896291557410367300109;
            static constexpr Scalar y3 = 0.443649167310370844258963269989;
        public:
            static constexpr int num_points = 15;
            static constexpr std::array<std::array<Scalar,3>, num_points> points = {{
                {x0,x0,x0},
                {x1,x1,x1},{x1,x1,y1},{x1,y1,x1},{y1,x1,x1},
                {x2,x2,x2},{x2,x2,y2},{x2,y2,x2},{y2,x2,x2},
                {x3,x3,y3},{x3,y3,x3},{y3,x3,x3},
                {y3,y3,x3},{y3,x3,y3},{x3,y3,y3},
            }};
            static constexpr std::array<Scalar, num_points> weights = {
                w0,
                w1,w1,w1,w1,
                w2,w2,w2,w2,
                w3,w3,w3,w3,w3,w3
            };
    };