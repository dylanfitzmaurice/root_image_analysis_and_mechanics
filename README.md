# root_image_analysis_and_mechanics
This repository contains the code and instructions used to measure cell type specifc area of Arabidopsis thaliana root transverse sections. Additionally, it then creates a finite element model from the cross section and solves for the mechanical stresses at cellular resolution throughout the tissue. 

Code is written in the Matlab coding language. 

Cell counting and area measurements from optical sections.
  1. Use confocal microscope to obtain transverse optical sections.
  2. Use Adobe Illustrator to outline each cell and save file as .svg.
  3. Run Matlab script
       1. The script reads the .svg file (user needs to input the .svg directory and image's scale (micron/pixel).)
       2. The user then selects each cell type (Epidermis, Cortex, Middle Cortex, Endodermis, Stele) by clicking each cell and then pressing enter after every cell of a single cell type has been selected.
       3. The code will then count and calculate the area for every single cell and plot the average and standard deviation for each cell type.  
