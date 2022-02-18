#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <chrono>

using namespace std;

int main(int argc, char *argv[]) {

    std::vector<std::string> args(argv + 1,argv + argc);

    string caseName = args[0];
    string inputDir = args[1];
    int    nRlz     = args[2];

    cout << "*** PREPARING ***" << endl;
    cout << "Case name: " << caseName << endl;

    system(string("rm -rf ../data/" + caseName + " > /dev/null 2>&1").c_str());     // remove old data dir and contents
    system(string("mkdir ../data/" + caseName).c_str());                            // make new case data dirs
    system(string("mkdir ../data/" + caseName + "/data").c_str());                  // ... /data
    system(string("mkdir ../data/" + caseName + "/input").c_str());                 // ... /input
    system(string("mkdir ../data/" + caseName + "/runtime").c_str());               // ... /runtime

    system(string("cp " + inputDir + "/* ../data/" + caseName + "/input/ > /dev/null 2>&1").c_str());               // copy input file to data input dir
    system(string("cp " + inputDir + "/restart* ../data/" + caseName + "/input/ > /dev/null 2>&1").c_str());

    cout << "*** RUNNING ***" << endl;
    cout << "Output will be written to ../data/" << caseName << "/runtime and ../data/" << caseName << "/data" << endl;

    auto start = chrono::system_clock::now();
    time_t start_time = chrono::system_clock::to_time_t(start);
    cout << "Start simulation time: " << ctime(&start_time) << endl;

    for (int i=0; i<nRlz; i++) {
        cout << "-------------------- REALIZATION " << i << " -----------------------" << endl;
        system(string("../bin/odt-run " + caseName + " " + i + " 2>&1").c_str());
    }

    auto end = chrono::system_clock::now();
    time_t end_time = chrono::system_clock::to_time_t(end);
    cout << "\nEnd simulation time: " << ctime(&end_time);

    chrono::duration<double> elapsed_time  = (end-start);

    if (elapsed_time > chrono::duration<double>(3600.0))
        cout << "Elapsed time: " << elapsed_time.count()/3600.0 << " hr" << endl;
    else if (elapsed_time > chrono::duration<double>(60.0))
        cout << "Elapsed time: " << elapsed_time.count()/60.0 << " min" << endl;
    else
        cout << "Elapsed time: " << elapsed_time.count() << " s" << endl;

    return 0;
}