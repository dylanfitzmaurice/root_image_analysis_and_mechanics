# Mechanics
This repository contains the code and instructions used to create a finite element models from Arabidopsis thaliana root cross sections and solves for the mechanical stresses at cellular resolution throughout the tissue. 

This is performed after cell counting and area quantification.

Code is written in the Matlab coding language. 

Controls:

A control code to demonstrate how stresses are distributed through a solid disk.
solid_disk_freebc.m

| Radial stress | Circumferential (hoop) stress | Von Mises stress |
|---|---|---|
| <img src="https://github.com/user-attachments/assets/a85366db-6063-4bfc-bfdc-68d057581313" width="300"> | <img src="https://github.com/user-attachments/assets/c810e71b-4fef-4c28-9055-0c9bf5f78a2d" width="300"> | <img src="https://github.com/user-attachments/assets/593a5b88-7105-42d0-bab7-f3b39ca9b763" width="300"> |


A control to demonstrate how stresses are distributed through a tissue of approximately equal sized cells + 2x thick outer rim. 

| Radial stress | Circumferential (hoop) stress | Von Mises stress |
|---|---|---|
| <img src="https://github.com/user-attachments/assets/c65e0d89-924d-4936-8dcf-4bf9d51fdd54" width="300"> | <img src="https://github.com/user-attachments/assets/ecf4c13b-f50f-4fd2-af3b-8fd27909a70b" width="300"> | <img src="https://github.com/user-attachments/assets/1ddb8198-a43e-4ec2-a54b-1ca9fe81286c" width="300"> |


