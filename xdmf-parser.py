#!/usr/bin/python

import sys
import xml.etree.ElementTree as et

if len(sys.argv) < 2:
    raise RuntimeError('No arguments provided to the parser.')

filename = sys.argv[1]

if len(sys.argv) > 2:
    output = sys.argv[2]
else:
    output = filename[:filename.find('.')] + '-out.xdmf'

tree = et.parse(filename)
root = tree.getroot()
domain = root.find('Domain')
grids = domain.findall('Grid')

timeseries = grids[0]

other_timeseries = grids[1:]

attributes = []

for grid in other_timeseries:
    innergrids = grid.findall('Grid')
    attributes_aux = []
    for x in innergrids:
        attribute = x.find('Attribute')
        attributes_aux.append(attribute)
    attributes.append(attributes_aux)

timeseries_grids = timeseries.findall('Grid')

for grid_attributes in attributes:
    for (grid, attribute) in zip(timeseries_grids, grid_attributes):
        grid.append(attribute)

for x in other_timeseries:
    domain.remove(x)

tree.write(output, encoding="utf-8", xml_declaration=True)