/**
 * @file main.cc
 * @brief main driver for ODT 
 */

#include "domain.h"
#include "streams.h"
#include "processor.h"
#include "param.h"
#include "micromixer.h"
#include "meshManager.h"
#include "eddy.h"
#include "solver.h"
#include "randomGenerator.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"

#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

using namespace std;
using namespace Cantera;

//////////////////////////////////////////////////////////////

processor proc;

//////////////////////////////////////////////////////////////

int main(int argc, char*argv[]) {


    if(argc<3) {
        cout << endl << "ERROR: code needs caseName and shift arguments" << endl;
        return 1;
    }
    string caseName= argv[1];            // example: temporalJet (../input/temporalJet, without the ../input/)

    int nShiftFileNumbers = 0;
    stringstream ss1;
    string       s1;
    ss1.clear(); ss1 << argv[2];
    ss1 >> nShiftFileNumbers;

    inputoutput io(caseName, nShiftFileNumbers);
    param       pram(&io);
    streams     strm;
    IdealGasPhase gas("../input/gas_mechanisms/"+pram.chemMechFile);
    Transport   *tran = newTransportMgr("mixture-averaged", &gas);
    eddy        ed;
    meshManager mesher;
    solver      *solv;
    micromixer  *mimx;
    solv = new solver();
    mimx = new micromixer();

    domain domn(NULL,  &pram);
    domain eddl(&domn, &pram);

    // we should increment the seed if we are starting MPI multiple times
    if ( pram.seed >= 0 ) pram.seed += nShiftFileNumbers;
    randomGenerator rand(pram.seed);

    domn.init(&io,  &mesher, &strm, &gas, tran, mimx, &ed, &eddl, solv, &rand);
    eddl.init(NULL, NULL,    NULL,  NULL, NULL, NULL,  NULL,NULL,  NULL,  NULL, true);
    //
    //-------------------

    time_t mytimeStart, mytimeEnd;
    mytimeStart = time(0);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
    *io.ostrm << endl << "#################################################################";


    //-------------------

    domn.solv->calculateSolution();

    //domn.io->outputProperties("../data/init.dat", 0.0); //doldb
    //domn.mimx->advanceOdt(0.0, domn.pram->tEnd);        //doldb

    //-------------------

    //delete mimx;
    //delete solv;

    mytimeEnd = time(0);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
    *io.ostrm << endl << "#  End Time   = " << ctime(&mytimeEnd);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl;

    return 0;


}
