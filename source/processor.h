/**
 * @file processor.h
 * Header file for class processor: holds basic processor info
 */

#pragma once

#ifdef DOMPI
#include "mpi.h"
#endif

#include <string>
#include <ostream>

///////////////////////////////////////////////////////////////////////////////

/** Class defines MPI properties
 *
 *  @author David O. Lignell
 */

class processor{

  ////////////////////// DATA MEMBERS /////////////////////

   private:

        static int nInst;    ///< number of these class objects (should be 1);

    public:

        int myid;            ///< Process ID
        int nproc;           ///< Number of Processes
        int ierr;            ///< Error flag

#ifdef DOMPI
        MPI_Status   status;
        MPI_Request  request;
#endif

  ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

        processor();            // constructor
        ~processor();           // destructor

};

