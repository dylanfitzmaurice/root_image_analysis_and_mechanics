# Mechanics
This repository contains the code and instructions used to create a finite element models from Arabidopsis thaliana root cross sections and solves for the mechanical stresses at cellular resolution throughout the tissue. 

This is performed after cell counting and area quantification.

Code is written in the Matlab coding language. 

Controls:

A control code to demonstrate how stresses are distributed through a solid disk.
solid_disk_freebc.m

| Radial stress | Circumferential (hoop) stress | Von Mises stress |
|---|---|---|
| <img src="https://github.com/user-attachments/assets/4345a3d4-4673-4de0-baf6-0ca9f2239cfa" width="300"> | <img src="https://github.com/user-attachments/assets/6f1ced39-bf44-483d-aaef-29677d969714" width="300"> | <img src="https://github.com/user-attachments/assets/69c8684e-98d4-45a2-b9c4-258120f5dd0e" width="300"> |



A control to demonstrate how stresses are distributed through a tissue of approximately equal sized cells + 2x thick outer rim. 


| Radial stress | Circumferential (hoop) stress | Von Mises stress |
|---|---|---|
| <img src="https://github.com/user-attachments/assets/c65e0d89-924d-4936-8dcf-4bf9d51fdd54" width="300"> | <img src="https://github.com/user-attachments/assets/ecf4c13b-f50f-4fd2-af3b-8fd27909a70b" width="300"> | <img src="https://github.com/user-attachments/assets/1ddb8198-a43e-4ec2-a54b-1ca9fe81286c" width="300"> |


