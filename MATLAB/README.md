# root_image_analysis_and_mechanics MATLAB code

cell_type_counting_and_area.m 

Generally, this code loads a transverse optical section .svg, the user made using Abode Illustrator, then it meshes the cell walls and applies a 2x radially thicker outer cell 
wall and plots this visually. This visualation is useful as if the .svg file needs editing it will show internal cells which need to be editted. Next, the code allows the user
to select cell type by clicking each cell. After cell type selcetion, the number of each cell type and area for each cell is calculate.

If the user only wants to calculate cell-type number and cell-type area. This code will calculate these values. 

If the user wants to go on to make FEM models out of the .svg to solve for cell type stresses and related measurements additional code will be added to perform such.
