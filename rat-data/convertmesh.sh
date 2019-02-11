#!/bin/bash

for VAR in rat01 rat02 rat05 rat06 rat09 rat12
do
    cd $VAR
    dolfin-convert gmsh.msh gmsh.xml
    cd ..
done
