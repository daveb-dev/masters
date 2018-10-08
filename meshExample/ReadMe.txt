Quick steps:

Download Gmsh. 
Change the line in wheel_msh.m specifying the directory that Gmsh is installed, then
run wheel_mshr.m

All you need is the coordinates of the tumor ROI and boundary. Do not contain repeated
points in the list. In other word, you do not need to form a loop.

The script shalll generates gmsh_test.msh which you can load in FEniCS by:
mesh = Mesh("gmsh_test.xml")

To get a cell function seperating the two regions, first run command "dolfin-convert gmsh_test.msh gmsh_test.xml" in a terminal.
This generates two xml files with suffix "facet_region" and "physical_region". Then you can load the latter file in FEniCS
cf = MeshFunction("size_t",mesh,"gmsh_test_physical_region.xml")

The tumor region is labeled as 2 while the rest is labeled as 1.

Good luck!

