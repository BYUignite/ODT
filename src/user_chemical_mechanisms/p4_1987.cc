
#include <vector>
#include <cmath>
#include <iostream>

using std::vector;
using std::pow;
using std::exp;

// N. Peters and R. J. Key, "The computation of stretched laminar methane-air diffusion flames
//     using a reduced four-step mechanism," Combustion and Flame 68:17-29 (1987).
//     https://doi.org/10.1016/0010-2180(87)90062-9
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

    //---------- rate constants, mol, cm3, s, kcal, K

    double R = 1.98720425864E-3;     // kcal/mol*K

    double k1  = 2.2E4  * pow(T, 3.0)  * exp( -8.76/(R*T));
    double k2  = 1.6E6  * pow(T, 2.1)  * exp( -2.45/(R*T));
    double k6  = 2.0E14;
    double k7  = 7.14E14               * exp( -16.8/(R*T));
    double k8  = 3.0E12;
    double k9  = 4.4E6  * pow(T, 1.5)  * exp(  0.74/(R*T));
    double k10 = 1.2E17 * pow(T,-0.91) * exp(-16.52/(R*T));
    double k14 = 2.0E18 * pow(T,-0.8);
    double k15 = 2.15E22* pow(T,-2.0);
    double k16 = 1.5E14                * exp(  -1.0/(R*T));
    double k17 = 2.5E13                * exp( -0.69/(R*T));
    double k18 = 2.0E13;

    double KC12 = 0.2657   * pow(T,-0.0247) * exp(15.083/(R*T));
    double KII  = 3.828E-5 * pow(T,0.8139)  * exp( 9.839/(R*T));
    double KIV  = 11.283   * pow(T,-0.2484) * exp(11.400/(R*T));

    //---------- reaction rates mol/cm3*s

    //double H2inv = (H2 > 0) ? 1.0/H2 : 0.0;
    double H2inv = 1/(H2 + 1.0E-14);

    double w1  = k1*CH4*H;
    double w2  = k2/KC12*H2O*CH4*H*H2inv;
    double w6  = (k6*CH4*H*H)  * (k1 + k2*H2O/KC12*H2inv) / (k6*H + k7*CM + k8*O2);
    double w8  = (k8*O2*CH4*H) * (k1 + k2*H2O/KC12*H2inv) / (k6*H + k7*CM * k8*O2);
    double w9  = k9/KC12*H     * (CO*H2O*H2inv - CO2/KII);
    double w10 = k10*H         * (O2 - H*H*H2O*H2O/KIV*H2inv*H2inv*H2inv);
    double w14 = k14*H*O2*CM;
    double w15 = k15/KC12*H2O*H*H*CM*H2inv;
    double w16 = w14*k16 / (k16+k17+k18*H2O/KC12*H2inv);

    //---------- global reaction rates mol/cm3*s
    
    double wI   = w1 + w2; 
    double wII  = w9;
    double wIII = w6 + w8 + w14 + w15;
    double wIV  = w10 + w16;

    //---------- final species rates mol/cm3*s --> kmol/m3*s

    rrsp[iCH4] = -wI;
    rrsp[iO2]  = -wIV;
    rrsp[iCO2] =  wII;
    rrsp[iCO]  =  wI - wII;
    rrsp[iH2O] = -wI - wII + 2*wIV;
    rrsp[iH2]  =  4*wI + wII + wIII - 3*wIV;
    rrsp[iH]   =  2*(-wI - wIII + wIV);
    rrsp[iN2]  = 0.0;

    for(int i=0; i<8; i++)
        rrsp[i] *= 1000;          // mol/cm3*s --> kmol/m3*s
}
