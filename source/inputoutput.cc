
#include "inputoutput.h"
#include "domain.h"
#include "processor.h"
#include <sys/stat.h>             // for mkdir
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>               // max_element

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** inputoutput initialization function
 *
 * @param p_domn  \input set domain pointer with.
 */

void inputoutput::init(domain *p_domn) {
    domn    = p_domn;
    LdoDump = false;
}

///////////////////////////////////////////////////////////////////////////////
/** inputoutput constructor function
 *
 * @param p_caseName \input directory with input files
 * @param p_nShift   \input shift the file numbers by this amount (used for multiple sets of parallel runs).
 */

inputoutput::inputoutput(const string p_caseName, const int p_nShift){

    caseName     = p_caseName;
    inputFileDir = "../data/"+caseName+"/input/";

    nShift = p_nShift;

    inputFile   = YAML::LoadFile(inputFileDir+"input.yaml");     ///< use these "nodes" to access parameters as needed
    params      = inputFile["params"];
    streamProps = inputFile["streamProps"];
    initParams  = inputFile["initParams"];
    radParams   = inputFile["radParams"];
    dvParams    = inputFile["dvParams"];
    dTimes      = inputFile["dumpTimes"];
    bcCond      = inputFile["bcCond"];
    scalarSc    = inputFile["scalarSc"];

    for(int i=0; i<dTimes.size(); i++)
        dumpTimes.push_back(dTimes[i].as<double>());
    dumpTimes.push_back(1.0E200);                       ///< add an extra "infinity" for easy handling of the last time
    iNextDumpTime = 0;

    //----------- set the data directory and runtime file

    string fname;
    stringstream ss1;
    string       s1;

    ss1.clear(); ss1 << setfill('0') << setw(5) << proc.myid + nShift;
    s1 = ss1.str();
    dataDir = "../data/"+caseName+"/data/data_" + s1 + "/";   // e.g., "../data_00001", etc.

    int iflag = mkdir(dataDir.c_str(), 0755);
    if(iflag != 0) {
        cout << "\n********** Error, process " << proc.myid << "failed to create "
            << dataDir << ", or it was already there" << endl;
        exit(0);
    }

    fname = "../data/"+caseName+"/runtime/runtime_" + s1;
    ostrm = new ofstream(fname.c_str());

    //----------- set gnuplot file

    fname = dataDir + "plot_odt.gnu";
    gnufile.open(fname.c_str());
    if(!gnufile) {
        cout << endl << "ERROR OPENING FILE: " << dataDir+"plot_odt.gnu" << endl;
        exit(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Destructor function
*/

inputoutput::~inputoutput() {

    delete ostrm;
    gnufile.close();

}

///////////////////////////////////////////////////////////////////////////////
/** Writes a data file of the domain properties (in order)
 *
 * @param fname \input name of the file to write including path
 *  @param time \input time of the output
 */

void inputoutput::outputProperties(const string fname, const double time) {

    string       s1;
    stringstream ss1;

    //--------------------------

    for(int i=0; i<domn->v.size(); i++)
        domn->v.at(i)->setVar();             // make sure the variable is up to date

    //--------------------------

    ofstream ofile(fname.c_str());
    if(!ofile) {
        *ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
        exit(0);
    }

    *ostrm << endl << "# Writing outputfile: " << fname;
    ostrm->flush();

    //--------------------------

    ofile << "# time = "   << time;

    ofile << "\n# Grid points = "   << domn->ngrd;

    if(!domn->pram->LisHips)
        ofile << "\n# Domain Size = " << domn->Ldomain();

    ofile << "\n# Pressure (Pa) = " << domn->pram->pres << endl;

    // HEWSON setting tecplot friendly output
    if (domn->pram->Ltecplot)
        ofile << "VARIABLES =";
    else
        ofile << "#";
    for(int i=0,j=1; i<domn->v.size(); i++){
        if(domn->v.at(i)->L_output){
            if (domn->pram->Ltecplot)
                ofile << setw(14) << "\"" << j++ << "_" << domn->v.at(i)->var_name << "\"";
            else
                ofile << setw(14) << j++ << "_" << domn->v.at(i)->var_name;
        }
    }

    int iploc;
    ofile << scientific;
    ofile << setprecision(10);
    for(int i=0; i<domn->ngrd; i++) {
        iploc = (domn->pram->LisHips) ? domn->solv->pLoc[i] : i;    // HiPS uses an index array
        ofile << endl;
        for(int k=0; k<domn->v.size(); k++)
            if(domn->v.at(k)->L_output)
                ofile << setw(19) << domn->v.at(k)->d.at(iploc);
    }

    ofile.close();

}

///////////////////////////////////////////////////////////////////////////////
/** Set iNextDumpTime from time. Used for restarts.
 * @param time \input time to use to set iNextDumpTime.
 */

void inputoutput::set_iNextDumpTime(double time) {

    for(int i=0; i<dumpTimes.size(); i++)
        if(dumpTimes[i] > time) {   //set this greater-than dump at the start time
            iNextDumpTime = i;
            break;
        }
}

///////////////////////////////////////////////////////////////////////////////
/** Dumps a domain, sets flag, increments next dump
*/

void inputoutput::dumpDomainIfNeeded(){

    if(!LdoDump) return;

    stringstream ss;
    ss << setfill('0') << setw(5) << iNextDumpTime;
    string fname = dataDir + "dmp_" + ss.str() + ".dat";

    outputProperties(fname, dumpTimes.at(iNextDumpTime));

    iNextDumpTime++;
    LdoDump = false;
}

///////////////////////////////////////////////////////////////////////////////
/** Dumps a domain, sets flag, increments next dump
 *  @param fnameRaw \input file name without the path (just the name).
 *  @param time \input time of the output
 */

void inputoutput::writeDataFile(const string fnameRaw, const double time) {

    string fname = dataDir+fnameRaw;
    outputProperties(fname, time);
    gnufile << "plot '" << fnameRaw << "' us 1:3; pause -1;" << endl;

}

///////////////////////////////////////////////////////////////////////////////
/**Output title of properties displayed to screen. */

void inputoutput::outputHeader() {

    *ostrm << endl << "#--------------------------------------------------"
        << "--------------------------------------------------------------------";
    *ostrm << endl;
    *ostrm << setw(5) << "# EE,"
        << setw(12) << "time,"
        << setw(12) << "t-t0,"
        << setw(10) << "nEtry,"
        << setw(6)  << "ngrd,"
        << setw(12) << "edSize,"
        << setw(12) << "edPos,"
        << setw(12) << "edPa,"
        << setw(12) << "nEposs,"
        << setw(12) << "PaAvg,"
        << setw(12) << "invTauEddy"
        ;

    *ostrm << endl << "#--------------------------------------------------"
        << "--------------------------------------------------------------------";
}

///////////////////////////////////////////////////////////////////////////////
/**Outputs the data corresponding to outputHeader.
 * After a given number of accepted eddies, output this info.
 *
 */

void inputoutput::outputProgress() {

    double dmb = 0.5*(domn->ed->leftEdge + domn->ed->rightEdge);
    if(dmb > domn->posf->d.at(domn->ngrd))
        dmb = dmb-domn->Ldomain();

    *ostrm << scientific << setprecision(3) << endl;
    *ostrm << setw(5)  << domn->solv->neddies                    //  1: EE
        << setw(12) << domn->solv->time                       //  2: time
        << setw(12) << domn->solv->time-domn->solv->t0        //  3: t-t0
        << setw(10) << domn->solv->iEtrials                   //  4: nEtry
        << setw(6)  << domn->ngrd                             //  5: ngrd
        << setw(12) << domn->ed->eddySize                     //  6: edSize
        << setw(12) << dmb                                    //  7: edPos
        << setw(12) << domn->ed->Pa                           //  8: edPa
        << setw(12) << domn->solv->nPaSumC                    //  9: nEposs
        << setw(12) << domn->solv->PaSumC/domn->solv->nPaSumC // 10: PaAvg
        << setw(12) << domn->ed->invTauEddy                   // 11: invTauEddy
        ;
    ostrm->flush();
}

///////////////////////////////////////////////////////////////////////////////
/** Restart
 *  The number of columns in the restart file should match the number and order of domainvariables
 *      that are output to a data file. (That is, the order of the domainvariables in domain).
 */

void inputoutput::loadVarsFromRestartFile() {

    string fname;
    stringstream ss1;
    string       s1;

    for(int k=0; k<domn->v.size(); k++) {
        if(domn->v[k]->L_transported && !domn->v[k]->L_output) {
            cout << endl << "ERROR: to restart, all transported variables need to be in the restart file" << endl;
            exit(0);
        }
    }

    if(domn->pram->rstType == "multiple") {
        ss1.clear(); ss1 << setfill('0') << setw(5) << proc.myid;
        fname = inputFileDir + "restart/restart_" + ss1.str() + ".dat";
    }
    else
        fname = inputFileDir + "/restart.dat";

    ifstream ifile(fname.c_str());
    if(!ifile) {
        cout << endl << "ERROR: reading restart file " << fname << endl;
        exit(0);
    }

    //------------- Get file header information

    getline(ifile, s1);                        // read line "# time = 1.1" (this is the restart time
    ss1.clear();
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> domn->pram->trst;

    getline(ifile, s1);                        // read line "# Grid points = 100"
    ss1.clear();
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> s1 >> domn->ngrd;
    domn->ngrdf = domn->ngrd+1;

    getline(ifile, s1);                        // read line "# Domain Size = 2" (don't use)
    getline(ifile, s1);                        // read line "# Pressure (Pa) = 101325
    getline(ifile, s1);                        // read line "# column headers

    //------------- Get file data columns

    for(int k=0; k<domn->v.size(); k++)
        domn->v[k]->d.resize(domn->ngrd);
    domn->posf->d.resize(domn->ngrdf);

    for(int i=0; i<domn->ngrd; i++) {
        for(int k=0; k<domn->v.size(); k++) {
            if(!domn->v[k]->L_output)
                continue;
            ifile >> domn->v[k]->d[i];
        }
    }

    domn->posf->d[domn->ngrd] = domn->posf->d[0] + domn->pram->domainLength; 

    //------------- Set the variables

    for(int k=0; k<domn->v.size(); k++)
        domn->v[k]->setVar();

}

