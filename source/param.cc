/**
 * @file param.cc
 * @brief Source file for class \ref param
 */

#include "param.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////
/** param initialization function
 * @param p_domn  \input set domain pointer with.
 */

void param::init(domain *p_domn) {
    domn = p_domn;
}

///////////////////////////////////////////////////////////////////////////////
/** param constructor function
 *
 * @param p_io  \input set inputoutput pointer with.
 */

param::param(inputoutput *p_io) {

    io = p_io;

    seed           = io->params["seed"]           ? io->params["seed"].as<int>()             : -1;
    tEnd           = io->params["tEnd"]           ? io->params["tEnd"].as<double>()          : errMsg<double>("tEnd");
    domainLength   = io->params["domainLength"]   ? io->params["domainLength"].as<double>()  : errMsg<double>("domainLength");
    ngrd0          = io->params["ngrd0"]          ? io->params["ngrd0"].as<int>()            : 1000;     //errMsg<int>("ngrd0");
    rho0           = io->params["rho0"]           ? io->params["rho0"].as<double>()          : 1.0;      //errMsg<double>("rho0");
    kvisc0         = io->params["kvisc0"]         ? io->params["kvisc0"].as<double>()        : 0.001694; //errMsg<double>("kvisc0");
    sdiff0         = io->params["sdiff0"]         ? io->params["sdiff0"].as<double>()        : 0.001694; //errMsg<double>("sdiff0");
    dPdx           = io->params["dPdx"]           ? io->params["dPdx"].as<double>()          : 0.0;
    pres           = io->params["pres"]           ? io->params["pres"].as<double>()          : 101325.0;
    chemMech   = io->params["chemMech"]   ? io->params["chemMech"].as<string>()  : errMsg<string>("chemMech");
    probType       = io->params["probType"]       ? io->params["probType"].as<string>()      : errMsg<string>("probType");

    Z_param        = io->params["Z_param"]        ? io->params["Z_param"].as<double>()       : 400.0;    //errMsg<double>("Z_param");
    A_param        = io->params["A_param"]        ? io->params["A_param"].as<double>()       : 0.666667; //errMsg<double>("A_param");
    C_param        = io->params["C_param"]        ? io->params["C_param"].as<double>()       : 5.0;      //errMsg<double>("C_param");
    LES_type       = io->params["LES_type"]       ? io->params["LES_type"].as<string>()      : "NONE";   //errMsg<string>("LES_type");
    Z_LES          = io->params["Z_LES"]          ? io->params["Z_LES"].as<double>()         : 0.0;
    x0virtual      = io->params["x0virtual"]      ? io->params["x0virtual"].as<double>()     : 0.0;
    diffCFL        = io->params["diffCFL"]        ? io->params["diffCFL"].as<double>()       : errMsg<double>("diffCFL");
    cvode_atol     = io->params["cvode_atol"]     ? io->params["cvode_atol"].as<double>()    : 1.0E-10;
    cvode_rtol     = io->params["cvode_rtol"]     ? io->params["cvode_rtol"].as<double>()    : 1.0E-4;

    Lbuoyant       = io->params["Lbuoyant"]       ? io->params["Lbuoyant"].as<bool>()        : false;
    LPeEddy        = io->params["LPeEddy"]        ? io->params["LPeEddy"].as<bool>()         : false;
    LplanarExpCent0= io->params["LplanarExpCent0"]? io->params["LplanarExpCent0"].as<bool>() : false;
    g              = io->params["g"]              ? io->params["g"].as<double>()             : -9.81;
    LdoDL          = io->params["LdoDL"]          ? io->params["LdoDL"].as<bool>()           : false;
    Lsolver        = io->params["Lsolver"]        ? io->params["Lsolver"].as<string>()       : errMsg<string>("Lsolver");
    Lperiodic      = io->params["Lperiodic"]      ? io->params["Lperiodic"].as<bool>()       : false;
    Lspatial       = io->params["Lspatial"]       ? io->params["Lspatial"].as<bool>()        : false;
    LTMA           = io->params["LTMA"]           ? io->params["LTMA"].as<bool>()            : false;
    LplanarTau     = io->params["LplanarTau"]     ? io->params["LplanarTau"].as<bool>()      : false;
    Lignition      = io->params["Lignition"]      ? io->params["Lignition"].as<bool>()       : false;

    bcType         = io->params["bcType"]         ? io->params["bcType"].as<string>()        : errMsg<string>("bcType");
    cCoord         = io->params["cCoord"]         ? io->params["cCoord"].as<int>()           : 1;
    xDomainCenter  = io->params["xDomainCenter"]  ? io->params["xDomainCenter"].as<double>() : 0.0;


    gDens          = io->params["gDens"]          ? io->params["gDens"].as<double>()         : 30;  //errMsg<double>("gDens");
    dxmin          = io->params["dxmin"]          ? io->params["dxmin"].as<double>()         : 0.0; //errMsg<double>("dxmin");
    dxmax          = io->params["dxmax"]          ? io->params["dxmax"].as<double>()         : 0.2;

    Pmax           = io->params["Pmax"]           ? io->params["Pmax"].as<double>()          : 0.4;    //errMsg<double>("Pmax");
    Pav            = io->params["Pav"]            ? io->params["Pav"].as<double>()           : 0.02;   //errMsg<double>("Pav");
    dtfac          = io->params["dtfac"]          ? io->params["dtfac"].as<double>()         : 2.0;    //errMsg<double>("dtfac");
    nDtSmeanWait   = io->params["nDtSmeanWait"]   ? io->params["nDtSmeanWait"].as<int>()     : 100000; //errMsg<int>("nDtSmeanWait");
    eddyMinCells   = io->params["eddyMinCells"]   ? io->params["eddyMinCells"].as<int>()     : 3;      //errMsg<int>("eddyMinCells");
    DAtimeFac      = io->params["DAtimeFac"]      ? io->params["DAtimeFac"].as<double>()     : 10.0;   //errMsg<double>("DAtimeFac");
    tdfac          = io->params["tdfac"]          ? io->params["tdfac"].as<double>()         : 1.0;
    sLastDA        = io->params["sLastDA"]        ? io->params["sLastDA"].as<int>()          : 100;    //errMsg<int>("sLastDA");
    Lp             = io->params["Lp"]             ? io->params["Lp"].as<double>()            : 0.01;   //errMsg<double>("Lp");
    Lmax           = io->params["Lmax"]           ? io->params["Lmax"].as<double>()          : 1.0;    //errMsg<double>("Lmax");
    Lmin           = io->params["Lmin"]           ? io->params["Lmin"].as<double>()          : dxmin*eddyMinCells; //errMsg<double>("Lmin");

    modDump        = io->params["modDump"]        ? io->params["modDump"].as<int>()          : 1000000; //errMsg<int>("modDump");
    modDisp        = io->params["modDisp"]        ? io->params["modDisp"].as<int>()          : 1;
    Ltecplot       = io->params["Ltecplot"]       ? io->params["Ltecplot"].as<bool>()        : false; //errMsg<int>("modDump");

    LmultiPhase    = io->params["LmultiPhase"]    ? io->params["LmultiPhase"].as<bool>()     : false;
    eSurfTens      = io->params["eSurfTens"]      ? io->params["eSurfTens"].as<double>()     : 0.0;

    Lrestart       = io->params["Lrestart"]       ? io->params["Lrestart"].as<bool>()        : false;
    rstType        = io->params["rstType"]        ? io->params["rstType"].as<string>()       : "single";    // "single" or "multiple"
    trst = 0.0; // (dont read this in, it comes from the restart file

    umin_spatial   = io->params["umin_spatial"]   ? io->params["umin_spatial"].as<double>()  : 0.5;

    // Radiation variables ---------------------

    Lrad         = io->params["Lrad"]             ? io->params["Lrad"].as<bool>()               : false;
    radSolveType = io->radParams["radSolveType"]  ? io->radParams["radSolveType"].as<string>()  : "OPTHIN";
    radCoefType  = io->radParams["radCoefType"]   ? io->radParams["radCoefType"].as<string>()   : "PLANCKMEAN";
    ntheta       = io->radParams["ntheta"]        ? io->radParams["ntheta"].as<int>()           : 40;
    npsi         = io->radParams["npsi"]          ? io->radParams["npsi"].as<int>()             : 80;

    // Soot variables ---------------------

    Lsoot             = io->params["Lsoot"]                ? io->params["Lsoot"].as<bool>()                     : false;
    nsvar             = io->sootParams["nsvar"]            ? io->sootParams["nsvar"].as<int>()                  : 0;
    PSD_method        = io->sootParams["PSD_method"]       ? io->sootParams["PSD_method"].as<string>()          : "MONO";
    rhoSoot           = io->sootParams["rhoSoot"]          ? io->sootParams["rhoSoot"].as<double>()             : 1850.0;    ///< solid soot density, kg/m3
    nucleation_mech   = io->sootParams["nucleation_mech"]  ? io->sootParams["nucleation_mech"].as<string>()     : "NONE";
    growth_mech       = io->sootParams["growth_mech"]      ? io->sootParams["growth_mech"].as<string>()         : "NONE";
    oxidation_mech    = io->sootParams["oxidation_mech"]   ? io->sootParams["oxidation_mech"].as<string>()      : "NONE";
    coagulation_mech  = io->sootParams["coagulation_mech"] ? io->sootParams["coagulation_mech"].as<string>()    : "NONE";
    Cmin              = io->sootParams["Cmin"]             ? io->sootParams["Cmin"].as<int>()                   : 100;       ///< minimum number of carbon atoms in a soot particle
    b_coag            = io->sootParams["b_coag"]           ? io->sootParams["b_coag"].as<double>()              : 0.8536;    ///< coagulation constant; or use 1/sqrt(2), or 1, or avg of those
    for(int i=0; i<io->sootParams["PAH_species"].size(); i++)
        PAH_species.push_back(io->sootParams["PAH_species"][i].as<string>());

    // Premix variables ---------------------

    LisPremix      = io->params["LisPremix"]      ? io->params["LisPremix"].as<bool>()  : false;
    uInflow        = io->params["uInflow"]        ? io->params["uInflow"].as<double>()  : 0.0;

    //---------------------

    // Dirichlet velocity BCs (no-slip)
    uBClo       = io->bcCond["uBClo"]       ? io->bcCond["uBClo"].as<double>()       : 0.0;
    uBChi       = io->bcCond["uBChi"]       ? io->bcCond["uBChi"].as<double>()       : 0.0;
    vBClo       = io->bcCond["vBClo"]       ? io->bcCond["vBClo"].as<double>()       : 0.0;
    vBChi       = io->bcCond["vBChi"]       ? io->bcCond["vBChi"].as<double>()       : 0.0;
    wBClo       = io->bcCond["wBClo"]       ? io->bcCond["wBClo"].as<double>()       : 0.0;
    wBChi       = io->bcCond["wBChi"]       ? io->bcCond["wBChi"].as<double>()       : 0.0;
    // Dirichlet scalar BCs
    sBClo       = io->bcCond["sBClo"]       ? io->bcCond["sBClo"].as<double>()       : 0.0; //errMsg<double>("sBClo");
    sBChi       = io->bcCond["sBChi"]       ? io->bcCond["sBChi"].as<double>()       : 0.0; //errMsg<double>("sBChi");

    // Dirichlet temperature BCs
    hWallBCtype = io->bcCond["hWallBCtype"] ? io->bcCond["hWallBCtype"].as<string>() : "ADIABATIC";
    if(hWallBCtype == "ISOTHERMAL") {
        TBClo   = io->bcCond["TBClo"]       ? io->bcCond["TBClo"].as<double>()       : errMsg<double>("TBClo");
        TBChi   = io->bcCond["TBChi"]       ? io->bcCond["TBChi"].as<double>()       : errMsg<double>("TBChi");
    }
    if(Lrad) {
        TBClo   = io->bcCond["TBClo"]       ? io->bcCond["TBClo"].as<double>()       : errMsg<double>("TBClo");
        TBChi   = io->bcCond["TBChi"]       ? io->bcCond["TBChi"].as<double>()       : errMsg<double>("TBChi");
    }

    //--------------------- un-normalize

    dxmin *= domainLength;
    dxmax *= domainLength;
    Lmax  *= domainLength;
    Lmin  *= domainLength;
    Lp    *= domainLength;

    //--------------------- sanity checks

    if( (cCoord == 2) && (xDomainCenter != 0.0) && ( abs(xDomainCenter) < 0.5*domainLength) ) {
        cout << endl << "ERROR: for cylindrical, set xDomainCenter to be 0.0 (e.g, a pipe) or such as gives an annulus";
        cout << endl << "       That is, xDomainCenter should give a domain, that when rotated about the origin, doesnt overlap itself";
        exit(0);
    }
    if( (cCoord == 3) && (xDomainCenter != 0.0) )  {
        cout << endl << "ERROR: spherical case requires xDomainCenter to be zero";
        exit(0);
    }

    if(LTMA && cCoord == 1)
        cout << endl << "ERROR: don't use LTMA=true for cCoord=1. LTMA is for cylindrical/spherical." << endl;
    if(LplanarTau && cCoord == 1)
        cout << endl << "ERROR: don't use LplanarTau=true for cCoord=1. LPlanarTau is for cylindrical/spherical." << endl;

    if(LdoDL && Lsolver=="STRANG")
        cout << endl << "ERROR: STRANG solver is not set up with Darrieus Landau instability LdoDL" << endl;

}