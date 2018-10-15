% comment: 
% (x_skull, y_skull) are coordinates of the domain boundary
% (x_tumor, y_tumor) are coordinates of the tumor ROI
% Note, THE LIST SHALL NOT CONTAIN REPEATED POINTS

clear;
clc;
close all;
addpath ./

filesNums = %{'01','02','05','06','09','12'};
for i = 1:length(filesNums)
    cd(strcat('./rat',filesNums{i}))
    load('skull_out.mat')
    load('tumor_out.mat')
    global x_skull y_skull x_tumor y_tumor
    x_skull = skull_out(1:end-1,1);
    y_skull = skull_out(1:end-1,2);
    x_tumor = tumor_out(1:end-1,1);
    y_tumor = tumor_out(1:end-1,2);
    output_filename = 'gmsh.geo';
    util_generateGmshGeo(output_filename,20,10);
    !~/gmsh2/gmsh gmsh.geo -2 -o gmsh.msh
    !~/gmsh2/gmsh gmsh.msh
    cd ..
end

rmpath ./