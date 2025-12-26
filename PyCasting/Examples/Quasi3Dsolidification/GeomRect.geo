// Gmsh project created on Sun Oct 15 18:55:25 2023
SetFactory("OpenCASCADE");
W = DefineNumber[ 0.372, Name "Parameters/W" ];
R = DefineNumber[ 0.025, Name "Parameters/R" ];
Point(1) = {0, W/2, 0, 1.0};
Point(2) = {W/2-R, W/2, 0, 1.0};
Point(5) = {W/2, W/2-R, 0, 1.0};
Point(6) = {W/2-R, W/2-R, 0, 1.0};
Point(7) = {W/2, 0, 0, 1.0};
Point(8) = {0, 0, 0, 1.0};
Line(1) = {8, 7};
Line(2) = {7, 5};
Circle(3) = {5, 6, 2};
Line(4) = {2, 1};
Line(5) = {1, 8};
Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Physical Curve("Outer", 6) = {4, 3, 2};
