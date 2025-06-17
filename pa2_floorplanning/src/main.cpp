#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>    
#include <cstdlib>   
#include <chrono>

#include "floorplanner.h"
#include "tree.h"
using namespace std;

int main(int argc, char** argv)
{
    srand(static_cast<unsigned int>(time(0)));
    fstream input_blk, input_net, output;
    double alpha;


    // start timing 
    auto start_time = chrono::high_resolution_clock::now();


    if (argc == 5) {
        alpha = stod(argv[1]);
        input_blk.open(argv[2], ios::in);
        input_net.open(argv[3], ios::in);
        output.open(argv[4], ios::out);
        if (!input_blk) {
            cerr << "Cannot open the input file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!input_net) {
            cerr << "Cannot open the input file \"" << argv[3]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[4]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./Floorplanner <alpha> <input block file> " <<
                "<input net file> <output file>" << endl;
        exit(1);
    }

    Floorplanner* fp = new Floorplanner(input_blk, input_net, alpha);
    fp->floorplan();



    // End timing
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cout << "Time taken: " << duration.count()*0.001 << " s" << endl;
    fp->computeMaxX();
    fp->computeMaxY();
    fp->_tree->exportToFile(output, fp->getOutlineWidth(), fp->getOutlineHeight(), duration, fp->getOutputCost(),  fp->getOutputWirelength(), fp->getOutputArea(), fp->_maxX, fp->_maxY);
    cout << "alpha: " << fp->_alpha << endl;

    if(fp->_maxX < fp->getOutlineWidth() && fp->_maxY < fp->getOutlineHeight())
    {
        cout << "Floorplan is feasible." << endl;
    }
    else
    {
        cout << "Floorplan is infeasible." << endl;
    }
    cout << "=============== Floorplan completed. ===============" << endl;

//    fstream output_visual;
//    output_visual.open("visualize.txt", ios::out);
//    fp->_tree->exportToFile_visualize(output_visual, fp->getOutlineWidth(), fp->getOutlineHeight(), duration, fp->getOutputCost(),  fp->getOutputWirelength(), fp->getOutputArea(), fp->_maxX, fp->_maxY);


    return 0;
}


