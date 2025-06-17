
# Physical Design for Nanometer ICs – Programming Assignments

**Course**: Physical Design for Nanometer ICs, Spring 2025  
**Department**: Electrical Engineering, National Taiwan University  
**Instructor**: Prof. Yao-Wen Chang  
**Student**: Yi-En Wu

This repository contains the programming assignments and corresponding reports developed as part of the Physical Design for Nanometer ICs course. Each assignment focuses on solving a core problem in the physical design flow of integrated circuits using C++ on a Linux platform.

---

## Repository Structure

.
├── pa1/ # Programming Assignment 1 – 2-Way Partitioning
│ ├── bin/ # Compiled binary: fm
│ ├── src/ # Source files (.cpp/.h)
│ ├── Makefile
│ ├── readme.txt
│ └── report.pdf
├── pa2/ # Programming Assignment 2 – Fixed-Outline Floorplanning
│ ├── bin/ # Compiled binary: fp
│ ├── src/ # Source files
│ ├── Makefile
│ ├── readme.txt
│ └── report.pdf
├── pa3/ # Programming Assignment 3 – Global Placement
│ ├── bin/ # Compiled binary: place
│ ├── src/ # Source files
│ ├── Makefile
│ ├── readme.txt
│ └── report.pdf

## Assignments Overview

### PA1 – 2-Way Fiduccia–Mattheyses Circuit Partitioning
- **Objective**: Partition a set of cells into two balanced groups to minimize the cutsize.
- **Algorithm**: Fiduccia–Mattheyses heuristic with a gain bucket list for efficient gain updates and cell movement.
- **Key Components**: Custom bucket list, net and cell structures, balance constraint checking.
- **Executable**: `./fm <input_file> <output_file>`
- **Report**: See `pa1/report.pdf`

### PA2 – Fixed-Outline Floorplanning
- **Objective**: Arrange rectangular macros within a fixed outline to minimize area and wirelength.
- **Techniques**: B*-Tree representation for block placement, simulated annealing for optimization, and contour structure for efficient packing.
- **Parameter**: α value (user-defined) balances area vs. wirelength.
- **Executable**: `./fp <alpha> <input.block> <input.net> <output_file>`
- **Report**: See `pa2/report.pdf`

### PA3 – Analytical Global Placement
- **Objective**: Spread cells across a chip layout to minimize total half-perimeter wirelength while satisfying density constraints.
- **Method**: Analytical placement using weighted-average (WA) wirelength model and bell-shaped density functions, optimized by Conjugate Gradient method.
- **Executable**: `./place -aux <input.aux>`
- **Report**: See `pa3/report.pdf`

---
