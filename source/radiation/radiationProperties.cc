/** Source file for class radiationProperties
 *
 *  @file radiationProperties.cc
 *  @author Victoria B. Lansinger
 *
 */

#include "radiationProperties.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////
/** Constructor for radiationProperties class object
 *
 *  @author Victoria B. Lansinger
 *
 *  @param  p_domn  \input pointer to domain object
 */

void radiationProperties::init(domain *p_domn) {

    // --------------- define constants and values

    domn = p_domn;

    // --------------- populate list of gas species indices

    bool fmissing = false;

    nRadSp = 4;                     // CH4, CO2, H2O, CO
    iRadIndx.resize(nRadSp);

    int isp;

    isp = domn->gas->speciesIndex("CH4");
    isp = (isp > 0) ? isp : domn->gas->speciesIndex("ch4");
    iRadIndx[0] = isp;
    if(isp < 0) fmissing = true;

    isp = domn->gas->speciesIndex("CO2");
    isp = (isp > 0) ? isp : domn->gas->speciesIndex("co2");
    iRadIndx[1] = isp;
    if(isp < 0) fmissing = true;

    isp = domn->gas->speciesIndex("H2O");
    isp = (isp > 0) ? isp : domn->gas->speciesIndex("h2o");
    iRadIndx[2] = isp;
    if(isp < 0) fmissing = true;

    isp = domn->gas->speciesIndex("CO");
    isp = (isp > 0) ? isp : domn->gas->speciesIndex("co");
    iRadIndx[3] = isp;
    if(isp < 0) fmissing = true;

    if(fmissing)
        cout << endl << "Warning: one or more radiating species missing from gas mechanism" << endl;

    // --------------- initialize variables needed to get absorption coefficients

    if (domn->pram->radCoefType == "PLANCKMEAN") {
        init_planck_mean_coefs();
    }
    else if (domn->pram->radCoefType == "WSGG") {
        init_WSGG_coefs();
        nGG = 5;                        // 1 clear gas and 4 gray gases
    }
    else if (domn->pram->radCoefType == "RCSLW") {
        init_RCSLW_coefs();
    }
    else {
        *domn->io->ostrm << endl << "ERROR: radCoefType not recognized" << endl;
        exit(0);
    }


    return;

}

///////////////////////////////////////////////////////////////////////////////
/** Initializes constants needed for calculating Planck mean absorption coefficients
 *
 *
 */
void radiationProperties::init_planck_mean_coefs() {

    pm_radCoefs.resize(nRadSp+1,vector<double>(6,0));  // kp (=) 1/atm*m

    pm_radCoefs[0][0] =  6.6334;     // ch4
    pm_radCoefs[0][1] = -0.0035686;
    pm_radCoefs[0][2] =  1.6682E-08;
    pm_radCoefs[0][3] =  2.5611E-10;
    pm_radCoefs[0][4] = -2.6558E-14;
    pm_radCoefs[0][5] =  0;

    pm_radCoefs[1][0] =  18.741;     // co2
    pm_radCoefs[1][1] = -121.310;
    pm_radCoefs[1][2] =  273.500;
    pm_radCoefs[1][3] = -194.050;
    pm_radCoefs[1][4] =  56.310;
    pm_radCoefs[1][5] = -5.8169;

    pm_radCoefs[2][0] = -0.23093;     // h2o
    pm_radCoefs[2][1] = -1.12390;
    pm_radCoefs[2][2] =  9.41530;
    pm_radCoefs[2][3] = -2.99880;
    pm_radCoefs[2][4] =  0.51382;
    pm_radCoefs[2][5] = -1.86840E-5;

    pm_radCoefs[3][0] =  4.7869;       // co < 750 K
    pm_radCoefs[3][1] = -0.06953;
    pm_radCoefs[3][2] =  2.95775E-4;
    pm_radCoefs[3][3] = -4.25732E-7;
    pm_radCoefs[3][4] =  2.02894E-10;
    pm_radCoefs[3][5] =  0.0;

    pm_radCoefs[4][0] =  10.09;         // co > 750 K
    pm_radCoefs[4][1] = -0.01183;
    pm_radCoefs[4][2] = 4.7753E-6;
    pm_radCoefs[4][3] = -5.87209E-10;
    pm_radCoefs[4][4] = -2.5334E-14;
    pm_radCoefs[4][5] =  0.0;

    return;

}

///////////////////////////////////////////////////////////////////////////////
/** function
 *
 */
void radiationProperties::init_WSGG_coefs() {

    c.resize(4, vector<vector<double> >(5, vector<double>(5, 0.0)));
    d.resize(5, vector<double>(5, 0.0));

    double cc[]={ 7.412956e-001, -5.244441e-001,  5.822860e-001, -2.096994e-001,  2.420312e-002,
                 -9.412652e-001,  2.799577e-001, -7.672319e-001,  3.204027e-001, -3.910174e-002,
                  8.531866e-001,  8.230754e-002,  5.289430e-001, -2.468463e-001,  3.109396e-002,
                 -3.342806e-001,  1.474987e-001, -4.160689e-001,  1.697627e-001, -2.040660e-002,
                  4.314362e-002, -6.886217e-002,  1.109773e-001, -4.208608e-002,  4.918817e-003,
                  1.552073e-001, -4.862117e-001,  3.668088e-001, -1.055508e-001,  1.058568e-002,
                  6.755648e-001,  1.409271e+000, -1.383449e+000,  4.575210e-001, -5.019760e-002,
                 -1.125394e+000, -5.913199e-001,  9.085441e-001, -3.334201e-001,  3.842361e-002,
                  6.040543e-001, -5.533854e-002, -1.733014e-001,  7.916083e-002, -9.893357e-003,
                 -1.105453e-001,  4.646634e-002, -1.612982e-003, -3.539835e-003,  6.121277e-004,
                  2.550242e-001,  3.805403e-001, -4.249709e-001,  1.429446e-001, -1.574075e-002,
                 -6.065428e-001,  3.494024e-001,  1.853509e-001, -1.013694e-001,  1.302441e-002,
                  8.123855e-001, -1.102009e+000,  4.046178e-001, -8.118223e-002,  6.298101e-003,
                 -4.532290e-001,  6.784475e-001, -3.432603e-001,  8.830883e-002, -8.415221e-003,
                  8.693093e-002, -1.306996e-001,  7.414464e-002, -2.029294e-002,  2.010969e-003,
                 -3.451994e-002,  2.656726e-001, -1.225365e-001,  3.001508e-002, -2.820525e-003,
                  4.112046e-001, -5.728350e-001,  2.924490e-001, -7.980766e-002,  7.996603e-003,
                 -5.055995e-001,  4.579559e-001, -2.616436e-001,  7.648413e-002, -7.908356e-003,
                  2.317509e-001, -1.656759e-001,  1.052608e-001, -3.219347e-002,  3.386965e-003,
                 -3.754908e-002,  2.295193e-002, -1.600472e-002,  5.046318e-003, -5.364326e-004};

    double dd[]={ 3.404288e-002,  6.523048e-002, -4.636852e-002,  1.386835e-002, -1.444993e-003,
                  3.509457e-001,  7.465138e-001, -5.293090e-001,  1.594423e-001, -1.663261e-002,
                  4.570740e+000,  2.168067e+000, -1.498901e+000,  4.917165e-001, -5.429990e-002,
                  1.098169e+002, -5.092359e+001,  2.343236e+001, -5.163892e+000,  4.393889e-001};

    int ni = 4;
    int nj = 5;
    int nk = 5;

    for(int i=0; i<ni; i++)
        for(int j=0; j<nj; j++)
            for(int k=0; k<nk; k++)
                c[i][j][k] = cc[ i*nj*nk + j*(nk) + k ];

    for(int i=0; i<ni; i++)
        for(int k=0; k<nk; k++)
            d[i][k] = dd[ i*nk + k ];

    return;

}

///////////////////////////////////////////////////////////////////////////////
/** function
 *
 */
void radiationProperties::init_RCSLW_coefs() {

    return;

}
///////////////////////////////////////////////////////////////////////////////
/** function
 *
 *  @param T        \input temperature (K)
 *  @param P        \input pressure (Pa)
 *  @param xMoleSp  \input species mole fractions
 */
void radiationProperties::get_planck_mean_coefs(const vector<vector<double> > &xMole,
                                                const vector<double>          &T,
                                                const double                  &pressure,
                                                vector<double>                &Kabs) {

    double P = pressure/101325.;         // Pa --> atm

    for (int i=0; i<domn->ngrd; i++) {

        double Ktemp = 0.0;     // temporary storage
        Kabs[i] = 0;            // reset value at grid point

        for(int k=0; k<nRadSp; k++) {

            if(iRadIndx[k] < 0)          // this radiation species k is not in the mechanism
                continue;
            if(k==0) {          // ch4
                Ktemp = pm_radCoefs[k][4];
                for(int j=3; j>=0; j--)
                    Ktemp = Ktemp * T[i] + pm_radCoefs[k][j];
            }
            else if(k==1) {     // co2
                Ktemp = pm_radCoefs[k][5];
                for(int j=4; j>=0; j--)
                    Ktemp = Ktemp * 1000/T[i] + pm_radCoefs[k][j];
            }
            else if(k==2) {     // h2o
                Ktemp = pm_radCoefs[k][5];
                for(int j=4; j>=0; j--)
                    Ktemp = Ktemp * 1000/T[i] + pm_radCoefs[k][j];
            }
            else if(k==3){      // co
                int kk = (T[i]<=750 ? k : k+1);
                Ktemp = pm_radCoefs[kk][4];
                for(int j=3; j>=0; j--)
                    Ktemp = Ktemp * T[i] + pm_radCoefs[kk][j];
            }

            Kabs[i] += xMole[i][iRadIndx[k]]*P*Ktemp;
        }
    }

}


///////////////////////////////////////////////////////////////////////////////
/** function
 *
 *  @param T        \input temperature (K)
 *  @param P        \input pressure (Pa)
 *  @param xMoleSp  \input species mole fractions
 *  @param Kwsgg    \ouput Kwsgg[ngrd, ngg]
 *  @param awsgg    \ouput awsgg[ngrd, ngg]
 *
 */
void radiationProperties::get_WSGG_coefs(const vector<vector<double> > &xMoleSp,
                                         const vector<double>          &T,
                                         const double                  &P,
                                         vector<vector<double> >       &Kwsgg,
                                         vector<vector<double> >       &awsgg){

    for(int igrd=1; igrd<domn->ngrd; igrd++) {

        double XCO2 = xMoleSp.at(igrd).at(iRadIndx[1]);
        double XH2O = xMoleSp.at(igrd).at(iRadIndx[2]);

        // get K's and a's ------------------------

        if(abs(XCO2) < 1E-6) XCO2 = 1E-6;   // check for values that break Mr

        double Mr = XH2O/XCO2;
        double Tr = T[igrd]/1200;

        int ni  = 4;
        int nj  = 5;
        int nk  = 5;

        vector<vector<double> > b(ni, vector<double>(nj,0.0));
        // calculate b_ij, see Bordbar (2014) Eq. 10
        for(int i=0; i<ni; i++)
            for(int j=0; j<nj; j++)
                for(int k=0; k<nk; k++)
                    b[i][j] += c[i][j][k]*pow(Mr,k);

        // calculate a_i, see Bordbar (2014) Eq. 9
        awsgg[0][igrd] = 1.0;
        for(int i=1; i<nGG; i++){
            for(int j=0; j<nj; j++)
                awsgg[i][igrd] += b[i-1][j]*pow(Tr,j);
            awsgg[0][igrd] -= awsgg[i][igrd];
        }

        // calculate K_i, see Bordbar (2014) Eq. 11
        for(int i=1; i<nGG; i++){
            for(int k=0; k<nk; k++)
                Kwsgg[i][igrd] += d[i-1][k]*pow(Mr,k);
            Kwsgg[i][igrd] *= (P/101325)*(XCO2 + XH2O);
        }

    }

}

///////////////////////////////////////////////////////////////////////////////
/** function
 *
 */
void radiationProperties::get_RCSLW_coefs() {

    return;

}

