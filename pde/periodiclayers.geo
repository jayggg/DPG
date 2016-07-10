algebraic3d

solid xneg = plane (0, 0, 0; -1,  0, 0) -bc=1;
solid xpos = plane (1, 1, 1;  1,  0, 0) -bc=1;
solid yneg = plane (0, 0, 0;  0, -1, 0) -bc=1;
solid ypos = plane (1, 1, 1;  0,  1, 0) -bc=1;

solid zneg = plane (0, 0, 0;  0, 0, -1) -bc=2;
solid zpos = plane (1, 1, 1;  0, 0,  1) -bc=3;

solid layer1 = orthobrick(0, 0, 0; 1, 1, 0.5) and xneg and xpos and yneg and ypos and zneg;

solid layer2 = orthobrick(0, 0, 0.5; 1, 1, 1) and xneg and xpos and yneg and ypos and zpos;

tlo layer1 -maxh=0.3; 
tlo layer2 -maxh=0.3;

identify periodic xneg xpos;
identify periodic yneg ypos;
