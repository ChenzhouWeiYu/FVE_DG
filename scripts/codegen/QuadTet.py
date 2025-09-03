import os
import re
import numpy as np


sss = ''
for filename in os.listdir('./QuadData'):
    if(re.search(r'NME_6313_cubature_tetra_p\d+_n\d+',filename)):
        # print(filename)
        result = re.findall(r'\d+',filename)
        polynomial_order = int(result[1])
        points_num = int(result[2])
        # print(polynomial_order,points_num)
        # polynomial_order = int(re.findall(r'\d+',l1[0])[0])
        # points_num = int(re.findall(r'\d+',l2[0])[0])
        # orbits = [int(s) for s in re.findall(r'\d+',l3[0])]
        # # print(polynomial_order,points_num,orbits)
        # points = []
        # weight = []
        with open('./QuadData/'+filename,) as fp:
            lines = fp.readlines()
        s1 = """    struct Degree%dPoints%d {
        public:
            static constexpr int num_points = %d;
            static constexpr std::array<std::array<Scalar,3>, num_points> points = {{
"""%(polynomial_order,points_num,points_num)
        for mm,line in enumerate(lines):
            x,y,z,w = line.split(' ')
            s1 += '                {\n                    %s, \n                    %s, \n                    %s\n                }, \n'%(x,y,z)
            # if((mm)%3==0): s1 += '                '
            # s1 += '{%s, %s, %s},  '%(x,y,z)
            # if((mm)%3==2): s1 += '\n'
        s1 += """            }};
            static constexpr std::array<Scalar, num_points> weights = {
"""
        for mm,line in enumerate(lines):
            x,y,z,w = line.split(' ')
            # if((mm)%3==0): s1 += '                '
            s1 += '                %s/6.0,\n'%(w[:-1])
            # if((mm)%3==2): s1 += '\n'
        s1 += """            };
    };
"""
        # print(s1)
        sss += s1

with open('./QuadTet.h','w') as f:
    print(sss,file=f)