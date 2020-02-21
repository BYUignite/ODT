/** Header file for class radiationProperties
 *
 *  @file radiationProperties.h
 *  @author Victoria B. Lansinger
 *
 */

#pragma once

#include <string>
#include <vector>

class domain;

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** Class implementing radiationProperties, no parent object
 *
 *  Options for calculating absoprtion coefficients:
 *
 *      PLANCK: Planck mean absorption coefficients; see TNF workshop radiation
 *              submodel at www.sandia.gov/TNF/radiation.html
 *      WSGG:   Weighted sum of gray gases method; see Bordbar et al. (2014)
 *              C&F 161:2435-2445
 *      RCSLW:  Rank-correlated spectral-line weighted sum of gray gases
 *              method; see Solovjov et al. (2008) J. Quant. Spectirosc.
 *              Radiat. Transfer 197:26-44
 *
 *  @author Victoria B. Lansinger
 *  @author David O. Lignell
 *
 */

class radiationProperties {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        domain                          *domn;      ///< pointer to domain object

        int                              nRadSp;     ///< number of radiating species
        vector<int>                      iRadIndx;   ///< index of radiating species
        int                              nGG;        ///< number of grey gases for WSGG approaches

        vector<vector<double> >          pm_radCoefs;   ///< Planck mean coef values; [spc][coef]

        vector<vector<vector<double> > > c;          ///< WSGG helper values
        vector<vector<double> >          d;          ///< WSGG helper values


    ////////////////////// MEMBER FUNCTIONS /////////////////

    public:

        void init_planck_mean_coefs();
        void get_planck_mean_coefs(const vector<vector<double> > &xMole,
                                   const vector<double>          &T,
                                   const double                  &pressure,
                                   vector<double>                &Kabs);

        void init_WSGG_coefs();
        void get_WSGG_coefs(const vector<vector<double> > &xMole,
                            const vector<double>          &T,
                            const double                  &pressure,
                            vector<vector<double> >       &Kwsgg,
                            vector<vector<double> >       &awsgg);

        void init_RCSLW_coefs();
        void get_RCSLW_coefs();


    ////////////////////// CONSTRUCTOR FUNCTIONS /////////////

    radiationProperties(){}; // constructor
    void init(domain *p_domn);
    ~radiationProperties(){}                // destructor

};

