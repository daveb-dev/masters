function generateGmshGeo(output_filename,mesh_refinement_skull)
global x_skull y_skull
if(nargin == 1)
    mesh_refinement_skull = 1;
end
fp=fopen(output_filename,'w');
fprintf(fp,['lc1 = ',num2str(mesh_refinement_skull),';\n']);
%print the nodes on the skull boundary
x = x_skull;
y = y_skull;
Np_skull = length(x);
for n=1:Np_skull
    fprintf(fp,['Point(',num2str(n),')={',num2str(x(n)),',',num2str(y(n)),',',num2str(0),',',num2str(mesh_refinement_skull),'};\n']);
end
%print line connections between the nodes (skull)
Line_count = Np_skull;
for n=1:Np_skull-1
    fprintf(fp,['Line(',num2str(n),')={',num2str(n),',',num2str(n+1),'};\n']);
end
fprintf(fp,['Line(',num2str(Np_skull),')={',num2str(Np_skull),',',num2str(1),'};\n']);
%print line loop
fprintf(fp,['Line Loop(',num2str(Np_skull+1),')={']);
for n=1:Np_skull-1
    fprintf(fp,[num2str(n),',']);
end
fprintf(fp,[num2str(Np_skull),'};\n']);

%print the surface
fprintf(fp,['Plane Surface(',num2str(Np_skull+2),')={',num2str(Np_skull+1),', -',num2str(Np_skull+1),'};\n']);
%print the surface
%name of the area (Surface 1)
%name of the skull boundary (Boundary 1)
fprintf(fp,['Physical Surface("Surface 1")={',num2str(Np_skull+2),'};\n']);
fprintf(fp,['Physical Line("Boundary 1")={']);
for n=1:Np_skull-1
    fprintf(fp,[num2str(n),',']);
end
fprintf(fp,[num2str(Np_skull),'};\n']);
fclose all;
end