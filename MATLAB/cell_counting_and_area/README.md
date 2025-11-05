Workflow:

1. cell_type_counting_and_area.m

Generally, cell_type_counting_and_area.m  is used to quantify the number of each cell type and the area of each cell. To do so, the code loads a transverse optical section (a .svg file), then segments cells, meshes the cell walls, applies a 2x radially thicker outer cell wall (as observed experimentally), and plots this visually. This visualization is useful because if the .svg file needs editing, the visualization will plot the individual cells that need to be edited. 
