%test2
clear;
clc;
close all;
%%
% comment: 
% (x_kidney, y_kidney) are coordinates of the domain boundary
% (x_tumor, y_tumor) are coordinates of the tumor ROI
% Note, the list shall not contain repeated points.
global x_kidney y_kidney x_tumor y_tumor
angle = linspace(0,2*pi,100);
angle = angle(1:end-1);
R1 = 1;
R2 = 2;
x_tumor= cos(angle)*R1;
y_tumor = sin(angle)*R1;
x_kidney = cos(angle)*R2;
y_kidney = sin(angle)*R2;
output_filename = 'gmsh_test.geo';
util_generateGmshGeo(output_filename,20,10);
!~/gmsh/gmsh gmsh_test.geo -2 -o gmsh_test.msh
!~/gmsh/gmsh gmsh_test.msh
