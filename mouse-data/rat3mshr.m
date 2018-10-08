% 
clear;
clc;
close all;
%%
% comment: 
% (x_skull, y_skull) are coordinates of the domain boundary
% (x_tumor, y_tumor) are coordinates of the tumor ROI
% Note, the list shall not contain repeated points.
load('rat3skull.mat')
load('rat3tumor.mat')
global x_skull y_skull x_tumor y_tumor
x_skull = rat3sbound(1:end-1,1);
y_skull = rat3sbound(1:end-1,2);
x_tumor = rat3tbound(1:end-1,1);
y_tumor = rat3tbound(1:end-1,2);
output_filename = 'rat3gmsh.geo';
util_generateGmshGeo(output_filename,20,10);
!~/gmsh/gmsh rat3gmsh.geo -2 -o rat3gmsh.msh
!~/gmsh/gmsh rat3gmsh.msh
