% 
clear;
clc;
close all;
%%
% comment: 
% (x_skull, y_skull) are coordinates of the domain boundary
% (x_tumor, y_tumor) are coordinates of the tumor ROI
% Note, THE LIST SHALL NOT CONTAIN REPEATED POINTS
load('rat05skull.mat')
load('rat05tumor.mat')
global x_skull y_skull x_tumor y_tumor
x_skull = rat05skull(1:end-1,1);
y_skull = rat05skull(1:end-1,2);
x_tumor = rat05tumor(1:end-1,1);
y_tumor = rat05tumor(1:end-1,2);
output_filename = 'rat05gmsh.geo';
util_generateGmshGeo(output_filename,20,10);
!~/gmsh2/gmsh rat05gmsh.geo -2 -o rat05gmsh.msh
!~/gmsh2/gmsh rat05gmsh.msh
!dolfin-convert rat05gmsh.msh rat05gmsh.xml
