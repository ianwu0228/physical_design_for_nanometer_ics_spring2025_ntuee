=========================================
README - Programming Assignment #1: FM Partitioning
=========================================

Student ID: B11901188
Name: Yi-En Wu
Course: Physical Design for Nanometer ICs, Spring 2025
Instructor: Yao-Wen Chang

-----------------------------------------
1. Compilation
-----------------------------------------
To compile the project, open a terminal and navigate to the root directory that contains the Makefile.

Then simply run:

    make

This will compile the source code located in the `src/` directory and generate an executable named `fm` inside the `bin/` directory.

To clean up all compiled files before recompiling:

    make clean
    make

-----------------------------------------
2. Execution
-----------------------------------------
After compilation, you can run the program using the following command format:

    ./bin/fm <input file> <output file>

Example:

    ./bin/fm ./input/input_0.dat ./output/output_0.dat

The program reads a circuit description from the input file, applies the Fiduccia-Mattheyses heuristic for 2-way partitioning, and outputs the partitioning result and cut size into the output file.

-----------------------------------------
3. Input Format
-----------------------------------------
The input file should be in the format:

    <Balance Degree>
    NET <Net Name> [<Cell Name>]+ ;

Example:

    0.5
    NET n1 c2 c3 c4 ;
    NET n2 c3 c6 ;
    ...

-----------------------------------------
4. Output Format
-----------------------------------------
The output will follow the format:

    Cutsize = <value>
    G1 <number of cells>
    <cell names> ;
    G2 <number of cells>
    <cell names> ;

-----------------------------------------
5. File Structure
-----------------------------------------
The submission directory should have this structure:

b11901188_pa1/
├── bin/
│   └── fm            <- The compiled binary
├── src/
│   ├── main.cpp
│   ├── partitioner.cpp
│   ├── partitioner.h
│   ├── cell.h
│   ├── net.h
│   ├── bucketlist.h
├── makefile          <- To build the program
├── readme.txt        <- This file
└── report.pdf        <- Report on data structures and findings

-----------------------------------------
6. Notes
-----------------------------------------
- Make sure your program compiles and runs correctly under Linux.
- Compilation optimization flag `-O3` is recommended.

-----------------------------------------
End of README
