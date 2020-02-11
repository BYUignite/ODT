/**
 * @file inputoutput.h
 * Header file for class inputoutput
 */

#pragma once

#include <vector>
#include <string>
#include <ostream>
#include <fstream>
#include "yaml-cpp/yaml.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing inputoutput object
 *
 *  @author David O. Lignell
 */

class inputoutput {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain                   *domn;          ///< pointer to domain object

        ostream                  *ostrm;         ///< Runtime: points to cout or to a file

        string                   caseName;       ///< input file directory
        string                   inputFileDir;   ///< input file directory
        string                   dataDir;        ///< data directory (output)

        YAML::Node               inputFile;      ///< yaml input file object base node
        YAML::Node               params;         ///< yaml sub node
        YAML::Node               sootParams;     ///< yaml sub node
        YAML::Node               streamProps;    ///< yaml sub node
        YAML::Node               initParams;     ///< yaml sub node
        YAML::Node               radParams;      ///< yaml sub node
        YAML::Node               dvParams;       ///< yaml sub node
        YAML::Node               dTimes;         ///< yaml sub node
        YAML::Node               bcCond;         ///< yaml sub node
        YAML::Node               scalarSc;       ///< yaml sub node: for hips

        vector<double>           dumpTimes;      ///< vector of dump times
        int                      iNextDumpTime;  ///< index of next dump time
        bool                     LdoDump;        ///< flag for whether we are dumping a file

        ofstream                 gnufile;        ///< gnuplot script file

        int                      nShift;         ///< for file naming


    //////////////////// MEMBER FUNCTIONS /////////////////

    void outputProperties(const string fname,
                          const double time);    ///< actually write the data file
    void dumpDomainIfNeeded();                     ///< calls outputProperties for dumpTimes
    void writeDataFile(const string fnameRaw,
                       const double time);       ///< writes the gnuplot file and calls outputProperties
    void outputHeader();                         ///< output header info during odt solution
    void outputFlmltHeader();                    ///< output header info during flmlt solution
    void outputProgress();                       ///< output data going with the header info
    void outputFlmltProgress();                  ///< output data going with the header info

    void loadVarsFromRestartFile();

    void set_iNextDumpTime(double time);

    void write_h_mixf_flmlt_profile(const vector<double> &hsens); ///< for flamelets: write h(mixf) to file for use in non-adiabatic cases
    void read_h_mixf_flmlt_profile(vector<double> &mixf, vector<double> &ha, vector<double> &hs);  ///< for flamelets: read data file

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        inputoutput(const string p_inputFileDir, const int nShift);
        void init(domain *p_domn);
        ~inputoutput();

};


////////////////////////////////////////////////////////////////////////////////


