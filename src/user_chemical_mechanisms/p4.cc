
#include <vector>
#include <cmath>
#include <iostream>

using std::vector;
using std::pow;
using std::exp;

// K. Seshadri and N. Peters, "Asymptotic Structure and Extinction of Methane-Air Diffusion Flames,"
//     Combustion and Flame 73:23-44 (1988).
//     https://doi.org/10.1016/0010-2180(88)90051-X
//
// inputs:
//     rho  [=] kg/m3
//     T    [=] K
//     P    [=] Pa
//     y    [=] species mass fractions 
// outputs:
//     rrsp [=] kmol/m3*s
//
//  Rxn I:     CH4 + 2H + H2O -->   CO + 4H2
//  Rxn II:    CO + H2O        =    CO2 + H2
//  Rxn III:   2H + M         -->   H2 + M
//  Rxn IV:    O2 + 3H2        =    2H + 2H2O

void getProblemSpecificRR(double rho, double T, double P, double *y, double *rrsp) {

    int iCH4  = 0;       // species indicies
    int iO2   = 1;
    int iCO2  = 2;
    int iCO   = 3;
    int iH2O  = 4;
    int iH2   = 5;
    int iH    = 6;
    int iN2   = 7;

    //---------- concentrations

    vector<double> Minv{6.233248145608677e-02, 3.125195324707794e-02, 2.272262491763049e-02, 3.570153516601214e-02, // kg/kmol = g/mol
                        5.550929780738274e-02, 4.960317460317460e-01, 9.920634920634921e-01, 3.569643749553795e-02};

    rho /= 1000.0;   // kg/m3 --> g/cm3

    double CH4 = rho*y[iCH4]*Minv[iCH4];     // mol/cm3
    double O2  = rho*y[iO2 ]*Minv[iO2 ];  
    double CO2 = rho*y[iCO2]*Minv[iCO2];  
    double CO  = rho*y[iCO ]*Minv[iCO ];  
    double H2O = rho*y[iH2O]*Minv[iH2O];  
    double H2  = rho*y[iH2 ]*Minv[iH2 ];  
    double H   = rho*y[iH  ]*Minv[iH  ];  
    double N2  = rho*y[iN2 ]*Minv[iN2 ];  

    double CM = CH4*6.54 + O2*0.4 + CO2*1.5 + CO*0.75 + H2O*6.5 + H2*1.0 + H*1.0 + N2*0.4;

    //---------- Rate constants

    double R = 1.98720425864;     // cal/mol*K

    double k11 = 2.20E4  * pow(T,  3.0)    * exp(-8748. /(R*T));
    double k10 = 4.40E6  * pow(T,  1.5)    * exp(  741. /(R*T));
    double k5  = 2.00E18 * pow(T, -0.8);
    double k1  = 1.20E17 * pow(T, -0.91)   * exp(-16515./(R*T));

    double k10K3    = 1.656E7    * pow(T, 1.5247)  * exp(-14342./(R*T));
    double k10K3KII = -4.326E11  * pow(T, 0.7108)  * exp(-24181./(R*T));
    double k1KIV    = -1.0635E16 * pow(T, -0.6615) * exp(-27915./(R*T));

    //---------- Reaction rates

    double H2inv = (H2 > 0) ? 1.0/H2 : 0.0;

    double wI   = k11*CH4*H;
    double wII  = k10K3*H*H2inv*(CO*H2O) + k10K3KII*H*CO2;
    double wIII = k5*O2*H*CM;
    double wIV  = k1*H*O2 + k1KIV*H*H*H*H2O*H2O*H2inv*H2inv*H2inv;

    //---------- final species rates mol/cm3*s --> kmol/m3*s

    rrsp[iCH4] = -wI;
    rrsp[iO2]  = -wIV;
    rrsp[iCO2] =  wII;
    rrsp[iCO]  =  wI - wII;
    rrsp[iH2O] = -wI - wII + 2*wIV;
    rrsp[iH2]  =  4*wI + wII + wIII - 3*wIV;
    rrsp[iH]   =  2*(-wI - wIII + wIV);
    rrsp[iN2]  =  0.0;

    for(int i=0; i<8; i++)
        rrsp[i] *= 1000;          // mol/cm3*s --> kmol/m3*s
}
