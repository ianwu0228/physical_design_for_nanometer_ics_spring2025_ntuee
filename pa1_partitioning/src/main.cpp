#include <iostream>
#include <fstream>
#include <vector>
#include "partitioner.h"
#include <chrono>
using namespace std;

int main(int argc, char** argv)
{
    // Start the timer
    auto start_time = chrono::high_resolution_clock::now();
    fstream input, output;

    if (argc == 3) {
        input.open(argv[1], ios::in);
        output.open(argv[2], ios::out);
        if (!input) {
            cerr << "Cannot open the input file \"" << argv[1]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./fm <input file> <output file>" << endl;
        exit(1);
    }

    Partitioner* partitioner = new Partitioner(input);
    partitioner->partition();
    partitioner->cutSize();
    partitioner->printSummary();
    partitioner->writeResult(output);

      // Stop the timer and calculate the duration
      auto end_time = chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
      
      // Print the execution time
      cout << "Execution time: " << duration.count() << " milliseconds" << endl;
      

    return 0;
}


// int main()
// {
//     // Start the timer
//     auto start_time = chrono::high_resolution_clock::now();
    

//     fstream input, output;
//     input.open("./input/input_2.dat", ios::in);
//     output.open("./output/output_2.dat", ios::out);
//     Partitioner* partitioner = new Partitioner(input);
    
//     partitioner->partition();
//     partitioner->cutSize();
//     partitioner->printSummary();
//     printf("done all\n");
//     partitioner->writeResult(output);
    

//     // Stop the timer and calculate the duration
//     auto end_time = chrono::high_resolution_clock::now();
//     auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    
//     // Print the execution time
//     cout << "Execution time: " << duration.count() << " milliseconds" << endl;
    

//     return 0;
// }
