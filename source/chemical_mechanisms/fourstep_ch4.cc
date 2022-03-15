/**
 * @file fourstep_ch4.cc
 * @brief Source file for one step ethylene mechanism
 */

#include "fourstep_ch4.h"

////////////////////////////////////////////////////////////////////////////////
/*! getProblemSpecificRR function
 *
 */
void fourstep_ch4::getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const {

//    CH4     + 2 H + H2O     <=> CO     + 4 H2
//    CO     + H2O         <=> CO2 + H2
//    3 H2     + O2         <=> 2 H + 2 H20
//    2 H     + M         <=> H2     + M

// STEP 1: solve for stoichiometric coefficients= v_sum[j][i] = v_backward[][] - v_vorward[][]
// STEP 2: rate-of-progress variable = q[i]
// STEP 3: solve for production rates = rrspc[j]= v_sum[j][i] * q[i]


// declare variables
    int                numEqns = 4;
    int                numSpcs = 8;    //M is not included, it has a 0 reaction rate

    vector < vector<int> >    v_forward;
    vector < vector<int> >    v_backward;
    vector < vector<int> >    v_sum(numEqns);    for(int j=0; j<numEqns; j++){    v_sum[j].resize(numSpcs, 0);    }

    vector < double >        q( numEqns, 0);

    double A, B, C, D, a, b, c, n, T_a;
    double A_inf, n_inf, T_a_inf, Fc, Pr, N, k_0, k_inf, F, dummy, numerator, denominator, denominator2;

// copy of global variables
    double molarMass_H = 1.007947;            //kg/kmol
    double molarMass_C = 12.01078;            //kg/kmol
    double molarMass_O = 15.99943;            //kg/kmol
    double molarMass_N = 14.00672;            //kg/kmol


    //add 1.0e-15 to each mass fraction to avoid division by zero.
    double cCH4    = rho*(yi[0]+1.0e-15)/( molarMass_C     + molarMass_H*4    )*1.0e-3;// (kg/m3) * (kmol/kg) ==> ( kmol /m3 ) * 1e-3 ==> mol/cm3
    double cCO    = rho*(yi[1]+1.0e-15)/( molarMass_C     + molarMass_O    )*1.0e-3;
    double cCO2    = rho*(yi[2]+1.0e-15)/( molarMass_C     + molarMass_O*2    )*1.0e-3;
    double cH2    = rho*(yi[3]+1.0e-15)/( molarMass_H*2              )*1.0e-3;
    double cH    = rho*(yi[4]+1.0e-15)/( molarMass_H              )*1.0e-3;
    double cO2    = rho*(yi[5]+1.0e-15)/( molarMass_O*2            )*1.0e-3;
    double cH2O    = rho*(yi[6]+1.0e-15)/( molarMass_H*2     + molarMass_O     )*1.0e-3;
    double cN2    = rho*(yi[7]+1.0e-15)/( molarMass_N*2            )*1.0e-3;



//STEP 1:  solve for stoichiometric coefficients= v_sum[j][i] = v_backward[][] - v_vorward[][]
    //the following initialization is done in this difficult way to avoid doing the following: temp[0]=1; temp[1]=0; temp[2]=0;...

    //species:         CH4,     CO,     CO2,     H2,     H,     O2,     H2O,     N2,     M(not included-would cancel out)
    //forwards reactions
    int tempList1[]={    1,    0,    0,    0,    2,    0,    1,    0};    vector < int >    temp1(tempList1, tempList1+numSpcs ); v_forward.push_back(temp1);
    int tempList2[]={    0,    1,    0,    0,    0,    0,    1,    0};    vector < int >    temp2(tempList2, tempList2+numSpcs ); v_forward.push_back(temp2);
    int tempList3[]={    0,    0,    0,    3,    0,    1,    0,    0};    vector < int >    temp3(tempList3, tempList3+numSpcs ); v_forward.push_back(temp3);
    int tempList4[]={    0,    0,    0,    0,    2,    0,    0,    0};    vector < int >    temp4(tempList4, tempList4+numSpcs ); v_forward.push_back(temp4);

    //backwards reactions
    int tempList5[]={    0,    1,    0,    4,    0,    0,    0,    0};    vector < int >    temp5(tempList5, tempList5+numSpcs ); v_backward.push_back(temp5);
    int tempList6[]={    0,    0,    1,    1,    0,    0,    0,    0};    vector < int >    temp6(tempList6, tempList6+numSpcs ); v_backward.push_back(temp6);
    int tempList7[]={    0,    0,    0,    0,    2,    0,    2,    0};    vector < int >    temp7(tempList7, tempList7+numSpcs ); v_backward.push_back(temp7);
    int tempList8[]={    0,    0,    0,    1,    0,    0,    0,    0};    vector < int >    temp8(tempList8, tempList8+numSpcs ); v_backward.push_back(temp8);

    for(int j=0; j<numEqns; j++){
      for(int i=0; i<numSpcs; i++)
    v_sum[j][i]= v_backward[j][i] - v_forward[j][i];
    }

//STEP 2: rate-of-progress variable = q[i]

    A    = 3.52e16;    n=-0.7;     T_a=8590;    double k_1f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 7.04e13;    n=-0.264;     T_a=72;        double k_1b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 5.06e4;    n=2.67;     T_a=3166;    double k_2f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 3.03e4;    n=2.633;     T_a=2433;    double k_2b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.17e9;    n=1.3;         T_a=1829;    double k_3f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.29e10;    n=1.196;     T_a=9412;    double k_3b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 7.08e13;    n=0.0;         T_a=148;    double k_5f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.66e13;    n=0.0;         T_a=414;    double k_6f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 2.89e13;    n=0.0;         T_a=-250;    double k_7f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 4.0e22;    n=-2.0;     T_a=0.0;    double k_8f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.3e18;    n=-1.0;     T_a=0.0;    double k_9f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 3.04e17;    n=-0.65;     T_a=52189.0;    double k_9b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 7.6;        n=3.84;     T_a=6431;    double k_10f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 0.41;        n=3.91;     T_a=-1884;    double k_10b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 4.4e6;    n=1.5;         T_a=-373;    double k_11f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 2.41e13;    n=0.222;     T_a=12581;    double k_11b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.3e4;    n=3.0;         T_a=4045;    double k_13f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 7.46;        n=3.51;     T_a=2955;    double k_13b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.6e7;    n=1.83;     T_a=1400;    double k_14f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.01e5;    n=2.237;     T_a=7893;    double k_14b    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 8.43e13;    n=0.0;         T_a=0.0;    double k_15f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 5.74e7;    n=1.9;         T_a=1383;    double k_16f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 3.90e10;    n=0.89;     T_a=205;    double k_17f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 5.0e13;    n=0.0;         T_a=0.0;    double k_18f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.860e17;    n=-1.0;     T_a=8555;    double k_19f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 1.1e13;    n=0.0;         T_a=14000;    double k_20f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 2.0e13;    n=0.0;         T_a=0.0;    double k_21f    = A * pow(temp, n) * exp(-T_a /temp);
    A    = 7.78e13;    n=0.0;         T_a=6800;    double k_22f    = A * pow(temp, n) * exp(-T_a /temp);


    /*----------------------------------------------------------------------------------------------------------*/

    double cM_4f= 2.5*cH2 + 16.0*cH2O + 1.2*cCO + 2.4*cCO2         + 1*(cH+cCH4+    cO2+cN2);
    double cM_8f= 2.5*cH2 + 12.0*cH2O + 1.9*cCO + 3.8*cCO2         + 1*(cH+cCH4+    cO2+cN2);
    double cM_9f= 2.5*cH2 + 12.0*cH2O + 1.9*cCO + 3.8*cCO2         + 1*(cH+cCH4+    cO2+cN2);    double cM_9b    = cM_9f;
    double cM_12f=2.0*cH2 + 6.00*cH2O + 1.5*cCO + 2.0*cCO2 + 2.0*cCH4     + 1*(cH+    cO2+cN2);    double cM_12b    = cM_12f;
    double cM_19f=1.2*cH2 + 12.0*cH2O + 2.5*cCO + 2.5*cCO2        + 1*(cH+cCH4+    cO2+cN2);
    double cM_22f=2.0*cH2 + 6.00*cH2O + 1.5*cCO + 2.0*cCO2 + 2.0*cCH4     + 1*(cH+    cO2+cN2);


    /*----------------------------------------------------------------------------------------------------------*/
    /*
     *Rate coeff calculation for rxns that contain Troe falloff coeff Fc
     *source: http://www.me.berkeley.edu/gri-mech/data/k_form.html
     *
     *k= (k_inf / (1+1/Pr) ) * F
     *
     *     Pr= k_0/k_inf *cM
     *   F= log(Fc) / (   1+ (  dummy/ ( N- 0.14 * dummy )  )²   )
     *     N= 0.75-1.27 * log(Fc)
     *       dummy= log(Pr) + C
     *         C= -0.4-0.67 log(Fc)
    */
    /*----------------------------------------------------------------------------------------------------------*/


    //Calc: k_4f
      A        = 5.75e19;    n    =-1.4;     T_a    =0;        k_0    = A     * pow(temp, n)        * exp(-T_a /temp);
      A_inf    = 4.65e12;    n_inf    =0.44;     T_a_inf    =0;        k_inf    = A_inf * pow(temp, n_inf) * exp(-T_a_inf /temp);
      //Fc    = 0.5;
      Fc    = 0.5* (  exp(-temp/1.0e30) + exp(-temp/1.0e-30)  );

      N = 0.75 - 1.27* log10(Fc);    C = -0.4 - 0.67* log10(Fc); Pr = k_0 / k_inf * cM_4f;    dummy= log10(Pr) +C;
      F = log10(Fc) * (  1 + ( pow(  dummy / ( N-0.14*dummy ) , 2 )));
      //F = exp(F);
      F = pow(10.0, F);

    double k_4f    =( k_inf / (1+1/Pr) ) * F;
    /*----------------------------------------------------------------------------------------------------------*/
    //Calc: k_12f
      A        = 2.47e33;    n    =-4.76; T_a    =1228;        k_0    = A     * pow(temp, n)        * exp(-T_a /temp);
      A_inf    = 1.27e16;    n_inf    =-0.63; T_a_inf    =193;        k_inf    = A_inf * pow(temp, n_inf) * exp(-T_a_inf /temp);
      Fc    = 0.217*exp(-temp /74.0) + 0.783*exp( -temp /2941.0 ) + exp( -6964.0/temp );


      N = 0.75 - 1.27* log10(Fc);    C = -0.4 - 0.67* log10(Fc); Pr = k_0 / k_inf * cM_12f;    dummy= log10(Pr) +C;
      F = log10(Fc) * (  1 + ( pow(  dummy / ( N-0.14*dummy ) , 2 )));
      //F = exp(F);
      F = pow(10.0, F);

    double k_12f=( k_inf / (1+1/Pr) ) * F;
    /*----------------------------------------------------------------------------------------------------------*/
    //Calc: k_12b
      A        = 1.0e36;    n    =-4.92; T_a    =54410;        k_0    = A     * pow(temp, n)        * exp(-T_a /temp);
      A_inf    = 5.17e18;    n_inf    =-0.79; T_a_inf    =53375;        k_inf    = A_inf * pow(temp, n_inf) * exp(-T_a_inf /temp);

      //Fc, N, C    = same as for k_12f
      Pr = k_0 / k_inf * cM_12f;    dummy= log10(Pr) +C;        //change log -->log10
      F = log10(Fc) * (  1 + ( pow(  dummy / ( N-0.14*dummy ) , 2 )));
      //F = exp(F);
      F = pow(10.0, F);

    double k_12b=( k_inf / (1+1/Pr) ) * F;
    //k_12b=6.0e5;    //constant avg value

    /*----------------------------------------------------------------------------------------------------------*/


    double d1    =(k_3b *cH2O * cH+ k_1f *cH *cO2);
    double d2    =k_3f *cH2;
    double d3    =k_1b *k_10b *pow(k_3b,2) *pow(cH,2) * cH2O;
    double d4    =k_10f /pow(k_3f,2) /pow(cH2,2) ;

    //to avoid division by 0, check if numerator = 0, if(yes) set eqn=0:
    double cOH;
    numerator    = ( k_3b *cH2O * cH + k_1f *cH *cO2 ) ;
    denominator2= ( k_10f * k_3f * k_3f *cH2 *cH2 );
    denominator    = ( k_3f *cH2 + k_1b *k_10b *k_3b*k_3b*cH*cH*cH2O / denominator2  );

    if( numerator == 0 )              cOH    = 0.0;
    else if (denominator == 0 || denominator2 ==2 && numerator!=0)    {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                          cOH    = numerator / denominator;

      A        = k_1f *cH *cO2 + k_2b *cOH *cH + k_10b *pow(cOH,2);
      B        = k_1b *cOH + k_2f *cH2 + k_10f *cH2O;
      C        = ( k_13f *cH + k_14f *cOH ) *cCH4;
      D        = k_12f *cH + k_13b *cH2 + k_14b *cH2O;
      a        = k_15f*B;
      b     = B*D + k_15f*(C-A);
      c        = -A*D;

    double cO;
    numerator    = (-b + pow(b*b-4*a*c, 0.5) );
    denominator    = (2*a);
    if( numerator == 0 )              cO    = 0.0;
    else if (denominator == 0 && numerator!=0)                {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                          cO    = numerator / denominator;

    double cHO2;
    numerator    = ( k_4f *cO2 *cH );
    denominator    = ( (k_5f+k_6f) *cH + k_7f *cOH );
    if( numerator == 0 )              cHO2    = 0.0;
    else if (denominator == 0 && numerator!=0)                {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                          cHO2    = numerator / denominator;

    double cCH3;
    numerator    = (k_13f *cH + k_14f *cOH) *cCH4;
    denominator    = (k_12f *cH + k_13b *cH2 + k_14b *cH2O + k_15f *cO);
    if( numerator == 0 )            cCH3    = 0.0;
    else if (denominator == 0 && numerator!=0)                {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                    cCH3    = numerator / denominator;

    double cCH3O;
    numerator    = (k_20f *cCH3 *cO2);
    denominator    = (k_21f *cH + k_22f *cM_22f);
    if( numerator == 0 )            cCH3O    = 0.0;
    else if (denominator == 0 && numerator!=0)                {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                    cCH3O    = numerator / denominator;

    double cCH2O;
    numerator    = (k_15f *cCH3 *cO + ( k_21f *cH + k_22f *cM_22f ) *cCH3O );
    denominator    = ( k_16f *cH + k_17f *cOH );
    if( numerator == 0 )            cCH2O    = 0.0;
    else if (denominator == 0 && numerator!=0)                {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                    cCH2O    = numerator / denominator;

    double cHCO;
    numerator    = ( k_16f *cH + k_17f *cOH ) *cCH2O ;
    denominator    = ( k_18f *cH + k_19f * cM_19f );
    if( numerator == 0 )            cHCO    = 0.0;
    else if (denominator == 0 && numerator!=0)                {    cout<<endl<<"error in oneStepRR"<<endl;    exit(0);    }
    else                    cHCO    = numerator / denominator;


    q[0]    = -k_12f *cCH3 *cH + k_12b *cCH4 + k_13f *cCH4 *cH - k_13b *cCH3 *cH2 + k_14f *cCH4 *cOH - k_14b *cCH3 *cH2O;
    q[1]    = k_11f *cCO *cOH - k_11b *cCO2 *cH;
    q[2]    = k_1f *cO2 *cH - k_1b *cOH *cO + k_5f *cH *cHO2 + k_20f *cCH3 *cO2;
    q[3]    = k_4f *cO2 *cH + k_8f *cM_8f *cH *cOH + k_9f *cM_9f *cH *cH + k_12f *cCH3 *cH - k_12b *cCH4 + k_18f *cHCO *cH - k_20f *cCH3 *cO2 + k_21f *cCH3O *cH - k_9b*cM_9b*cH2;


// STEP 3: solve for production rates = rrspc[j]= v_sum[j][i] * q[i]
    double rrSum= 0.0;
    for(int j=0; j<numSpcs; j++){
      dummy    =0.0;
      for(int i=0; i<numEqns; i++)
    dummy     += v_sum[i][j]*q[i];
      rrsp[j]    = dummy*1.0e3;            // convert mol/cm3*s ==> kmol/(m3*s)
      rrSum    +=dummy;
    }

}