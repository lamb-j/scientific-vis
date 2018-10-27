# Introduction to Scientific Visualization
CIS 510, Winter 2018,
University of Oregon

This repository contains materials used in the Winter 2018 Introduction to Scientific Visualization course (541). The main goal of this course was to visualize scientific data in different ways. We covered fields, meshes, interpolation, and advection. The primary programing projects involved developing an isosurface filter, a volume renderer, and a ray casting engine. The primary tools used were C++, VTK, and VisIt.

## Contents

* ./differencer - Small program used to compare images and verify correctness
* ./p1 - Install software and run example program
* ./p2 - Field interpolation
* ./p3 - Coloring field values
* ./p4 - Particle advection and streamlines
* ./p4G - Performance evaluation of Euler method and RK4 method for particle advection
* ./p5 - Generating isolines
* ./p6 - (Isosurface filter) Marching cubes cases
* ./p7 - (Isosurface filter) Algorithm implementation
* ./p8 - VTK (Visualization Toolkit) exercises
* ./ray_casting - (Final project 1) Implemented scientific visualization algorithm on ray-casted volume rendering.  Code was ~1000 lines in C++.  Algorithm consisted of calculating ray direction for each pixel on the screen, efficiently searching for overlaps along the ray with a large data set (>500MB), and deriving an image from the resulting intersections.
* ./marching_cubes_osx - (Final project 2) OpenCL implementation of marching cubes algorithm
* ./marching_cubes_linux - (Final project 2) OpenCL implementation of marching cubes algorithm
* ./README.md - This file

## Built With

* C++ - Language used to develop 3D rendering engine
* [VTK](https://www.vtk.org/) - The Visualization Toolkit, used for image file IO and other features
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) - Used for displaying image files

## Authors
* **Jacob Lambert** - jlambert.uo@gmail.com
