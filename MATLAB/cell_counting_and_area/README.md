Workflow:

1. cell_type_counting_and_area.m

Generally, cell_type_counting_and_area.m  is used to quantify the number of each cell type and the area of each cell. To do so, the code loads a transverse optical section (a .svg file), then segments cells, meshes the cell walls, applies a 2x radially thicker outer cell wall (as observed experimentally), and plots this visually. This visualization is useful because if the .svg file needs editing, the visualization will plot the individual cells that need to be edited. 

2. plot_cell_type_area.m

This code the plots the cell area of each cell type across n samples. 

<img width="1385" height="1007" alt="cell_type_by_area" src="https://github.com/user-attachments/assets/8f3b72b4-6ef8-4086-a0d6-01da125e703e" />
