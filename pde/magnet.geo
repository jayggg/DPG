algebraic3d

solid bot   = plane (0, 0, 0;  0,  0, -1)  -bc=1;
solid top   = plane (1, 1, 1;  0,  0,  1)  -bc=1;
solid left  = plane (0, 0, 0;  0, -1,  0)  -bc=2;
solid right = plane (1, 1, 1;  0,  1,  0)  -bc=2;
solid front = plane (0, 0, 0; -1,  0,  0)  -bc=3;
solid back  = plane (1, 1, 1;  1,  0,  0)  -bc=3;
 
solid cube = bot and top and left and right and front and back;


solid magnet = cylinder(0.5, 0.7, 0.5; 0.7, 0.7, 0.5; 0.2)
               and plane(0.5, 0.7, 0.5; -1,0,0)
	       and plane(0.7, 0.7, 0.5;  1,0,0);

solid matrix = cube and not magnet;

tlo matrix -transparent -material=matrix -maxh=0.2;
tlo magnet -col=[1,0,0] -material=magnet -maxh=0.2;

identify periodic left right;
identify periodic front back;

