#
## A cube with different bc markers on each face
#
algebraic3d


solid p1 = plane (0, 0, 0; -1, 0, 0) -bc=1;  # yz face through (0,0,0)
solid p2 = plane (0, 0, 0; 0, -1, 0) -bc=2;  # xz face through (0,0,0)
solid p3 = plane (0, 0, 0; 0, 0, -1) -bc=3;  # xy face through (0,0,0)
solid p4 = plane (1, 1, 1; 1, 0, 0)  -bc=4;   
solid p5 = plane (1, 1, 1; 0, 1, 0)  -bc=5;  # xz face through (1,1,1)
solid p6 = plane (1, 1, 1; 0, 0, 1)  -bc=6;  # xy face through (1,1,1)

solid cube = p1 and p2 and p3 and p4 and p5 and p6;

tlo cube;


