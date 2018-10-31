#!/bin/bash

for VAR in rat01 rat02 rat05
do
    cd $VAR
    dolfin-convert gmsh.msh gmsh.xml
    cd ..
done
