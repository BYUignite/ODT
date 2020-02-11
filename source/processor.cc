/**
 * @file processor.cc
 * Source file for class processor
 */

#include "processor.h"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int processor::nInst;

///////////////////////////////////////////////////////////////////////////////
/** Constructor */

processor::processor() {

#ifdef DOMPI

    //----------- set MPI Stuff (if on)

    if((++nInst)==1) {                 // Only ever call MPI_* once
        int fake_argc = 0;
        char** fake_argv;
        MPI_Init(&fake_argc, &fake_argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    }
    if(nInst > 1)
        cout << endl << "****** WARNING, more than one processor class object" << endl;

    if(myid==0)
        cout << "\n# MPI IS ON" << "; Nprocs = " << nproc << endl;
    MPI_Barrier(MPI_COMM_WORLD);


#else
    myid  = 0;
    nproc = 1;

#endif

}


///////////////////////////////////////////////////////////////////////////////
/** Destructor */

processor::~processor() {

#ifdef DOMPI

    if((--nInst)==0)             // Only ever finalize mpi once
        MPI_Finalize();

#endif
}

///////////////////////////////////////////////////////////////////////////////





