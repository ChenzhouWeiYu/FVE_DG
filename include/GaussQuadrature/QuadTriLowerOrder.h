    struct Degree3Points4 {
        private:
            static constexpr Scalar w1 = 0.09097930912801141530281549896241817511158;
            static constexpr Scalar w2 = 0.1590206908719885846971845010375818248884;
            static constexpr Scalar x1 = 0.07503111022260811817747559832460309763591;
            static constexpr Scalar x2 = 0.1785587282636164231170351333742224289755;
            static constexpr Scalar y1 = 0.2800199154990740720027959942048077631675;
            static constexpr Scalar y2 = 0.6663902460147013867026932740963667102211;
        public:
            static constexpr int num_points = 4;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {x1,y1}, {y1,x1}, {x2,y2}, {y2,x2}
            }};
            static constexpr std::array<Scalar, num_points> weights = {w1,w1,w2,w2};
    };

    struct Degree3Points4_B {
        private:
            static constexpr Scalar w1 = 0.06510416666666666666666666666666666666667;
            static constexpr Scalar w2 = 0.1921911384555808477492202542009676977468;
            static constexpr Scalar w3 = 0.1776005282110858189174464124656989689199;
            static constexpr Scalar x1 = 0.0;
            static constexpr Scalar y1 = 0.8;
            static constexpr Scalar x2 = 0.4339491425357193931250391681262089353371;
            static constexpr Scalar x3 = 0.1755746669880901306844846413976005884725;
        public:
            static constexpr int num_points = 4;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {x1,y1}, {y1,x1}, {x2,x2}, {x3,x3}
            }};
            static constexpr std::array<Scalar, num_points> weights = {w1,w1,w2,w3};
    };

    struct Degree3Points6 {
        private:
            static constexpr Scalar w1 = 0.5 * 0.0333333333333333333333333333333333333333;
            static constexpr Scalar w2 = 0.5 * 0.3;
            static constexpr Scalar x1 = 0.5;
            static constexpr Scalar y1 = 0.0;
            static constexpr Scalar x2 = 0.1666666666666666666666666666666666666666;
            static constexpr Scalar y2 = 0.6666666666666666666666666666666666666666;
        public:
            static constexpr int num_points = 6;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {x1,x1}, {x1,y1}, {y1,x1}, 
                {x2,x2}, {x2,y2}, {y2,x2}
            }};
            static constexpr std::array<Scalar, num_points> weights = {w1,w1,w1,w2,w2,w2};
    };

    struct Degree4Points7 {
        private:
            static constexpr Scalar w0 = 0.5 * 0.225;
            static constexpr Scalar x0 = 1.0/3;
            
            static constexpr Scalar w1 = 0.5 * 0.1259391805448271525956839455001813336576;
            static constexpr Scalar w2 = 0.5 * 0.1323941527885061807376493878331519996757;
            static constexpr Scalar x1 = 0.1012865073234563388009873619151238280556;
            static constexpr Scalar y1 = 0.7974269853530873223980252761697523438888;
            static constexpr Scalar x2 = 0.4701420641051150897704412095134476005159;
            static constexpr Scalar y2 = 0.05971587178976982045911758097310479896829;
        public:
            static constexpr int num_points = 7;
            static constexpr std::array<std::array<Scalar,2>, num_points> points = {{
                {x0,x0},
                {x1,x1}, {x1,y1}, {y1,x1}, 
                {x2,x2}, {x2,y2}, {y2,x2}
            }};
            static constexpr std::array<Scalar, num_points> weights = {w0,w1,w1,w1,w2,w2,w2};
    };