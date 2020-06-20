
#include "eddy.h"
#include "domain.h"
#include <cmath>
#include <numeric>   //accumulate
#include <algorithm> // min_element

///////////////////////////////////////////////////////////////////////////////
/** eddy initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_eddl  \input set domain pointer with.
 */

void eddy::init(domain *p_domn, domain *p_eddl) {
    domn = p_domn;
    eddl = p_eddl;

    esdp1 = -2.0*domn->pram->Lp;
    esdp2 = exp(-2.0*domn->pram->Lp/domn->pram->Lmax) - exp(-2.0*domn->pram->Lp/domn->pram->Lmin);
    esdp3 = exp(-2.0*domn->pram->Lp/domn->pram->Lmin);
    esdp4 = -esdp1 / esdp2;

    cCoef = vector<double> (3,0.0);
    bCoef = vector<double> (3,0.0);

    LperiodicEddy = false;

    cca = vector<double>(7);
    ccb = vector<double>(5);
    ccc = vector<double>(5);
    ccd = vector<double>(5);

    cca[0] = -6.58411371e+02;
    cca[1] = 1.19518310e+03;
    cca[2] = -8.21169900e+02;
    cca[3] = 2.61583310e+02;
    cca[4] = -3.38649746e+01;
    cca[5] = -6.11440937e-01;
    cca[6] = 1.01054351e+00;

    ccb[0] = -0.48987445;
    ccb[1] = 2.39482488;
    ccb[2] = -4.48286508;
    ccb[3] = 3.89924937;
    ccb[4] = -0.40448995;

    ccc[0] = -0.00169718;
    ccc[1] = 0.02173516;
    ccc[2] = -0.10620329;
    ccc[3] = 0.2399716;
    ccc[4] = 0.77720962;

    ccd[0] = -1.98503860e-07;
    ccd[1] = 1.31972798e-05;
    ccd[2] = -3.21487789e-04;
    ccd[3] = 3.44046508e-03;
    ccd[4] = 9.85964172e-01;

}

///////////////////////////////////////////////////////////////////////////////
/** Sample an eddy size from the eddy size distribution.
 * This could be as simple as a tophat, but this function is
 * more accurate, and more efficient.
 */
void eddy::sampleEddySize() {

    eddySize = domn->pram->Lmin * pow(( 1.0 - ( 1.0 - pow((domn->pram->Lmin/domn->pram->Lmax), 5./3.)) * domn->rand->getRand() ), -0.6);

}

///////////////////////////////////////////////////////////////////////////////
/** Uniformly sample the eddy position on the domain.  For periodic domains the
 * position can be anywhere.  For nonperiodic domains, the position
 * is from 0 to the end where end is the domain size - the eddy size
 * (since the eddy has to fit in the domain).  This also means that
 * the eddy position is the left edge of the eddy.
 * For periodic domains, rightEdge is greater than leftEdge (even for eddies that
 * wrap the domain, and may be outside the range of the base domain.)
 */
void eddy::sampleEddyPosition() {

    if(!domn->pram->Lperiodic) {
        leftEdge  = domn->rand->getRand() * (domn->Ldomain() - eddySize) + domn->posf->d.at(0);
        rightEdge = leftEdge + eddySize;
        if(leftEdge  < domn->posf->d.at(0))          leftEdge  = domn->posf->d.at(0);
        if(rightEdge > domn->posf->d.at(domn->ngrd)) rightEdge = domn->posf->d.at(domn->ngrd);
    }
    else {
        LperiodicEddy = false;             // reset the value
        leftEdge  = domn->rand->getRand() * domn->Ldomain() + domn->posf->d.at(0);
        rightEdge = leftEdge + eddySize;
        if(rightEdge > (domn->Ldomain() + domn->posf->d.at(0)))
            LperiodicEddy = true;
    }

}

/////////////////////////////////////////////////////////////////////////////
/** Triplet map the domain
 *  This version will adjust the TM so that the three
 *      segments are spaced evenly in cylindrical and spherical flows.
 *      (planar = no difference).
 *  This can be converted back to the usual equal volume segments by making fracVleft = 1/3 etc.
 *  @param line \inout either the domn or the eddl
 *  @param iS   \input starting cell
 *  @param iE   \input ending cell
 *  (Note, this could be done by the lv objects separately.)
 *  Applies to an eddy line, or the full odt line.
 */
void eddy::tripMap(domain *line, const int iS, int iE, const double C, const bool LsplitAtEddy) {

    int negrd    = iE-iS+1;
    int negrd2   = 2*negrd;
    int negrd3   = 3*negrd;

    double fracVleft;
    double fracVmidl;
    double fracVrght;
    double pm1;          // +1 or -1 factor to change sign
    double invC = 1.0/C;

    int k, i, ip, im, ie, iw;

    //--------- Set fracVlft, fracVmid, fracVrgt

    double r1 = (2*leftEdge+rightEdge)/3;
    double r2 = (leftEdge+2*rightEdge)/3;

    pm1 =leftEdge*r1 < 0 ? -1.0 : 1.0;             // handles cells that split x=0
    double Vleft = abs(pow(abs(leftEdge), C) - pm1*pow(abs(r1), C));

    pm1 =r1*r2 < 0 ? -1.0 : 1.0;                   // handles cells that split x=0
    double Vmidl = abs(pow(abs(r1), C) - pm1*pow(abs(r2), C));

    pm1 =r2*rightEdge < 0 ? -1.0 : 1.0;            // handles cells that split x=0
    double Vrght = abs(pow(abs(r2), C) - pm1*pow(abs(rightEdge), C));

    double Vtot = Vleft + Vmidl + Vrght;
    fracVleft = Vleft / Vtot;
    fracVmidl = Vmidl / Vtot;
    fracVrght = Vrght / Vtot;

    if(domn->pram->LTMA) {            // use three equal volume segments instead of three equal length sements
        fracVleft = 1.0/3;            // no difference for planar, but affects cylindrical or spherical
        fracVmidl = fracVleft;
        fracVrght = fracVleft;
    }
    
    //---------

    pos0 = line->pos->d;            // for filling kernel
    pos0.at(0) = 0.5*(line->posf->d.at(1)+leftEdge);
    pos0.at(line->ngrd-1) = 0.5*(line->posf->d.at(line->ngrd-1)+rightEdge);

    //----------- Grab cell "volume" array (delta(x^c)) in the eddy region

    dxc.resize(negrd,0.0);

    i=0; ip=1;
    pm1 = line->posf->d.at(iS+ip)*leftEdge < 0 ? -1.0 : 1.0;                     // handles cells that split x=0
    dxc.at(i) = abs(pow(abs(line->posf->d.at(iS+ip)), C) - pm1*pow(abs(leftEdge), C));

    for(int i=1, ip=2; i<negrd-1; i++, ip++) {
        pm1 = line->posf->d.at(iS+ip)*line->posf->d.at(iS+i) < 0 ? -1.0 : 1.0;   // handles cells that split x=0
        dxc.at(i) = abs(pow(abs(line->posf->d.at(iS+ip)), C) - pm1*pow(abs(line->posf->d.at(iS+i) ), C));
    }

    i=negrd-1;
    pm1 = line->posf->d.at(iS+i)*rightEdge < 0 ? -1.0 : 1.0;                     // handles cells that split x=0
    dxc.at(i) = abs(pow(abs(rightEdge), C) - pm1*pow(abs(line->posf->d.at(iS+i)), C));

    //---------- insert cells to be filled with values

    for(k=0; k<line->v.size(); k++)
        line->v.at(k)->d.insert(line->v.at(k)->d.begin()+iE+1, negrd2, 0.0);

    line->ngrd += negrd2;
    line->ngrdf = line->ngrd+1;
    iE += negrd2;

    //---------- fill in variable data (except pos, posf)

    for(k=0; k<line->v.size(); k++) {
        if(line->v.at(k)->var_name=="pos" || line->v.at(k)->var_name=="posf") continue;
        for(i=0; i<negrd; i++) {
            line->v.at(k)->d.at(iS+negrd2+i)   = line->v.at(k)->d.at(iS+i);
            line->v.at(k)->d.at(iS+negrd2-1-i) = line->v.at(k)->d.at(iS+i);
        }
    }

    //---------- update face positions

    if(C==1 || leftEdge > 0.0) {      // case: planar, or eddy on the right of the centerline

        // posf.at(iS) is the same
        line->posf->d.at(iS+1) = pow(pow(leftEdge,C) + fracVleft*dxc.at(0), invC);
        for(ie=2, iw=1, i=1; ie<=negrd; ie++, iw++, i++)
            line->posf->d.at(iS+ie) = pow(pow(line->posf->d.at(iS+iw),C) + fracVleft*dxc.at(i), invC);
        for(ie=negrd+1,iw=negrd,i=negrd-1; ie<=negrd2; ie++, iw++, i--)
            line->posf->d.at(iS+ie) = pow(pow(line->posf->d.at(iS+iw),C) + fracVmidl*dxc.at(i), invC);
        for(ie=negrd2+1,iw=negrd2,i=0; ie<negrd3; ie++, iw++, i++)
            line->posf->d.at(iS+ie) = pow(pow(line->posf->d.at(iS+iw),C) + fracVrght*dxc.at(i), invC);
        // line->posf->d.at(iS+negrd3) is the same due to the inserted cells before.

    }
    else if(rightEdge < 0.0) {                       // case: eddy on the left of centerline

        // posf.at(iS) is the same
        line->posf->d.at(iS+1) = -pow(pow(abs(leftEdge),C) - fracVleft*dxc.at(0), invC);
        for(ie=2, iw=1, i=1; ie<=negrd; ie++, iw++, i++)
            line->posf->d.at(iS+ie) = -pow(pow(abs(line->posf->d.at(iS+iw)),C) - fracVleft*dxc.at(i), invC);
        for(ie=negrd+1,iw=negrd,i=negrd-1; ie<=negrd2; ie++, iw++, i--)
            line->posf->d.at(iS+ie) = -pow(pow(abs(line->posf->d.at(iS+iw)),C) - fracVmidl*dxc.at(i), invC);
        for(ie=negrd2+1,iw=negrd2,i=0; ie<negrd3; ie++, iw++, i++)
            line->posf->d.at(iS+ie) = -pow(pow(abs(line->posf->d.at(iS+iw)),C) - fracVrght*dxc.at(i), invC);
        // line->posf->d.at(iS+negrd3) is the same due to the inserted cells before.

    }
    else {                                           // case: eddy splits the centerline
        // Note, this can directly replace the above two cases!

        double dmb;

        // posf.at(iS) is the same
        pm1 = leftEdge < 0 ? -1 : 1;
        dmb = pow(abs(leftEdge),C) + pm1*fracVleft*dxc.at(0);
        line->posf->d.at(iS+1) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        for(ie=2, iw=1, i=1; ie<=negrd; ie++, iw++, i++) {
            pm1 = line->posf->d.at(iS+iw) < 0 ? -1 : 1;
            dmb = pow(abs(line->posf->d.at(iS+iw)),C) + pm1*fracVleft*dxc.at(i);
            line->posf->d.at(iS+ie) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        }
        for(ie=negrd+1,iw=negrd,i=negrd-1; ie<=negrd2; ie++, iw++, i--) {
            pm1 = line->posf->d.at(iS+iw) < 0 ? -1 : 1;
            dmb = pow(abs(line->posf->d.at(iS+iw)),C) + pm1*fracVmidl*dxc.at(i);
            line->posf->d.at(iS+ie) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        }
        for(ie=negrd2+1,iw=negrd2,i=0; ie<negrd3; ie++, iw++, i++) {
            pm1 = line->posf->d.at(iS+iw) < 0 ? -1 : 1;
            dmb = pow(abs(line->posf->d.at(iS+iw)),C) + pm1*fracVrght*dxc.at(i);
            line->posf->d.at(iS+ie) = dmb < 0.0 ? pow(-dmb, invC) : pm1*pow(dmb, invC);
        }
        // line->posf->d.at(iS+negrd3) is the same due to the inserted cells before.
    }

    //---------- update center positions

    line->pos->setVar();

    //---------- split at eddy faces if needed

    if(LsplitAtEddy) {

        vector<double> interPos(3);                // used for cell split

        if(abs(line->posf->d.at(iE+1) - rightEdge) > 1.0e-15) {
            interPos.at(0) = line->posf->d.at(iE);
            interPos.at(1) = rightEdge;
            interPos.at(2) = line->posf->d.at(iE+1);
            domn->mesher->splitCell(iE, 1, interPos);       // don't interpolate here
        }
        if(abs(line->posf->d.at(iS) - leftEdge) > 1.0e-15) {
            interPos.at(0) = line->posf->d.at(iS);
            interPos.at(1) = leftEdge;
            interPos.at(2) = line->posf->d.at(iS+1);
            domn->mesher->splitCell(iS, 1, interPos);       // again, don't interp
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Compute the inverse of the eddy timescale.
 *  Used to define the eddy acceptance probability.
 *
 *  @param Z_value \input large eddy suppression parameter.
 *  @return false for implausible eddy, true for the usual, proper eddy.
 */
bool eddy::eddyTau(const double Z_value, const double C) {

    double         KK=0;                  // equivalent of the 4/27 fac
    double         rhoK=0, rhoJ=0, rhoKK=0, rhoJK=0;
    vector<double> uRhoJ(3,0.0), uRhoK(3,0.0);

    vector<double> P(3);
    vector<double> Q(3);
    double         S;
    double         A;
    double         eKinEddy     = 0.0;
    double         ePeEddy      = 0.0;
    double         eViscPenalty = 0.0;
    double         eDL = 0.0;
    double         rho_ref = 0.0;

    double         rhoEddy=0.0, viscEddy=0.0;

    vector<double> intRhoKi(eddl->ngrd);
    vector<double> intRhoJi(eddl->ngrd);

    int            i;
    double         Etot;

    //----------- set dxc

    domn->mesher->setGridDxc(eddl, dxc, C);

    if(leftEdge >= 0.0 || eddl->posf->d.at(1) <= 0.0)
        dxc[0] = abs( pow(abs(eddl->posf->d.at(1)), C) - pow(abs(leftEdge),C) );
    else
        dxc[0] = abs( pow(abs(eddl->posf->d.at(1)), C) + pow(abs(leftEdge),C) );

    if(eddl->posf->d.at(eddl->ngrd-1) >= 0.0 || rightEdge <= 0.0)
        dxc[eddl->ngrd-1] = abs( pow(abs(eddl->posf->d.at(eddl->ngrd-1)), C) - pow(abs(rightEdge),C) );
    else
        dxc[eddl->ngrd-1] = abs( pow(abs(eddl->posf->d.at(eddl->ngrd-1)), C) + pow(abs(rightEdge),C) );

    double VolE = accumulate(dxc.begin(), dxc.end(), 0.0);

    //-----------

    if(domn->pram->LplanarTau)
        fillKernel_planarAnalytic();
    else
        fillKernel();

    if(domn->pram->LdoDL) {         // Darrieus Landau instability for variable density flows
        vector<double> tempVec(eddl->ngrd);
        transform(eddl->rho->d.begin(), eddl->rho->d.end(), dxc.begin(), tempVec.begin(), multiplies<double>());
        rho_ref = accumulate(tempVec.begin(), tempVec.end(), 0.0) / VolE;
    }

    ///////////// Fill in integral quantities

    //---------- rhoK, rhoJ, UrhoK, UrhoJ

    for(i=0; i<eddl->ngrd; i++) {
        intRhoKi.at(i)  = K.at(i)*eddl->rho->d.at(i)*dxc.at(i);
        intRhoJi.at(i)  = abs(intRhoKi.at(i));
        KK          += K.at(i)*K.at(i)*dxc.at(i);
        if(!domn->pram->Lspatial) {
            rhoK        += intRhoKi.at(i);
            rhoJ        += abs(intRhoKi.at(i));
            rhoKK       += K.at(i)*intRhoKi.at(i);
            rhoJK       += K.at(i)*intRhoJi.at(i);
            uRhoK.at(0)    += intRhoKi.at(i)*eddl->uvel->d.at(i);
            uRhoK.at(1)    += intRhoKi.at(i)*eddl->vvel->d.at(i);
            uRhoK.at(2)    += intRhoKi.at(i)*eddl->wvel->d.at(i);
            uRhoJ.at(0)    += intRhoJi.at(i)*eddl->uvel->d.at(i);
            uRhoJ.at(1)    += intRhoJi.at(i)*eddl->vvel->d.at(i);
            uRhoJ.at(2)    += intRhoJi.at(i)*eddl->wvel->d.at(i);
            if(domn->pram->LdoDL)
                eDL += eddl->aDL->d.at(i)*K[i]*(eddl->rho->d[i] - rho_ref ) * dxc[i];
        }
        else {       // mass flux, not mass
            rhoK        += eddl->uvel->d.at(i) * intRhoKi.at(i);
            rhoJ        += eddl->uvel->d.at(i) * abs(intRhoKi.at(i));
            rhoKK       += eddl->uvel->d.at(i) * K.at(i)*intRhoKi.at(i);
            rhoJK       += eddl->uvel->d.at(i) * K.at(i)*intRhoJi.at(i);
            uRhoK.at(0)    += eddl->uvel->d.at(i) * intRhoKi.at(i)*eddl->uvel->d.at(i);
            uRhoK.at(1)    += eddl->uvel->d.at(i) * intRhoKi.at(i)*eddl->vvel->d.at(i);
            uRhoK.at(2)    += eddl->uvel->d.at(i) * intRhoKi.at(i)*eddl->wvel->d.at(i);
            uRhoJ.at(0)    += eddl->uvel->d.at(i) * intRhoJi.at(i)*eddl->uvel->d.at(i);
            uRhoJ.at(1)    += eddl->uvel->d.at(i) * intRhoJi.at(i)*eddl->vvel->d.at(i);
            uRhoJ.at(2)    += eddl->uvel->d.at(i) * intRhoJi.at(i)*eddl->wvel->d.at(i);
            if(domn->pram->LdoDL)
                eDL += eddl->uvel->d[i] * eddl->aDL->d.at(i)*K[i]*(eddl->rho->d[i] - rho_ref ) * dxc[i];
        }
    }

    //////////// Compute Eddy Energy

    A = rhoK/rhoJ;
    S = 0.5*(A*A+1.0)*rhoKK - A*rhoJK;
    for(i=0; i<3; i++) {
        P.at(i) = uRhoK.at(i) - A * uRhoJ.at(i);
        Q.at(i) = 0.25*P.at(i)*P.at(i)/S;
    }

    eKinEddy = KK/(eddySize*eddySize*VolE) * (Q.at(0)+Q.at(1)+Q.at(2));

    ePeEddy  = (domn->pram->LPeEddy ? domn->pram->g*rhoK : 0.0);

    //////////// Compute Viscous Energy Penalty

    for(i=0; i<eddl->ngrd; i++) {
        rhoEddy  += eddl->rho->d.at(i)  *dxc.at(i);
        viscEddy += eddl->dvisc->d.at(i)*dxc.at(i);
    }
    rhoEddy  /= VolE;
    viscEddy /= VolE;

    eViscPenalty = Z_value*0.5 * VolE/(eddySize*eddySize)*viscEddy*viscEddy/rhoEddy;

    double Ufavre;
    if(domn->pram->Lspatial) {
        Ufavre = eddyFavreAvgVelocity(dxc);
        eViscPenalty *= Ufavre;
    }

    invTauEddy = 0.0;

    Etot = eKinEddy - eViscPenalty + eDL + ePeEddy;

    if(Etot < 0.0) return false;

    invTauEddy = sqrt(2.0*KK/(rhoKK*VolE*eddySize*eddySize) * Etot);

    if(domn->pram->Lspatial)
        invTauEddy = invTauEddy/Ufavre; // 1/s --> 1/m

    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Compute the kernel coefficients.
 *  This used to be done in eddyTau, but to allow more flexibility, especially for
 *  cylindrical and spherical cases it is split off.
 *  Namely, for cylindrical and spherical, we evaluate a "planar" eddyTau so that
 *  a constant shear will give a uniform eddyTau with respect to position.
 *  But we need to evaluate the kernel coefficients using cylindrical or spherical
 *  so that we are conservative of momentum and energy.
 */

void eddy::set_kernel_coefficients() {

    double         rhoK=0, rhoJ=0, rhoKK=0, rhoJK=0;
    vector<double> uRhoJ(3,0.0), uRhoK(3,0.0);

    vector<double> P(3);
    vector<double> Q(3);
    double         S;
    double         A;

    vector<double> intRhoKi(eddl->ngrd);
    vector<double> intRhoJi(eddl->ngrd);

    int            i;

    //----------- set dxc

    double C = domn->pram->cCoord;
    domn->mesher->setGridDxc(eddl, dxc, C);

    if(leftEdge >= 0.0 || eddl->posf->d.at(1) <= 0.0)
        dxc[0] = abs( pow(abs(eddl->posf->d.at(1)), C) - pow(abs(leftEdge),C) );
    else
        dxc[0] = abs( pow(abs(eddl->posf->d.at(1)), C) + pow(abs(leftEdge),C) );

    if(eddl->posf->d.at(eddl->ngrd-1) >= 0.0 || rightEdge <= 0.0)
        dxc[eddl->ngrd-1] = abs( pow(abs(eddl->posf->d.at(eddl->ngrd-1)), C) - pow(abs(rightEdge),C) );
    else
        dxc[eddl->ngrd-1] = abs( pow(abs(eddl->posf->d.at(eddl->ngrd-1)), C) + pow(abs(rightEdge),C) );

    //-----------

    if(domn->pram->LplanarTau)          // todo: not strictly needed if using the same kernel as in eddyTau
        fillKernel_planarAnalytic();    // todo: probably is needed if using testThirds large eddy suppression

    ///////////// Fill in integral quantities

    //---------- rhoK, rhoJ, UrhoK, UrhoJ

    for(i=0; i<eddl->ngrd; i++) {
        intRhoKi.at(i)  = K.at(i)*eddl->rho->d.at(i)*dxc.at(i);
        intRhoJi.at(i)  = abs(intRhoKi.at(i));
        if(!domn->pram->Lspatial) {
            rhoK        += intRhoKi.at(i);
            rhoJ        += abs(intRhoKi.at(i));
            rhoKK       += K.at(i)*intRhoKi.at(i);
            rhoJK       += K.at(i)*intRhoJi.at(i);
            uRhoK.at(0)    += intRhoKi.at(i)*eddl->uvel->d.at(i);
            uRhoK.at(1)    += intRhoKi.at(i)*eddl->vvel->d.at(i);
            uRhoK.at(2)    += intRhoKi.at(i)*eddl->wvel->d.at(i);
            uRhoJ.at(0)    += intRhoJi.at(i)*eddl->uvel->d.at(i);
            uRhoJ.at(1)    += intRhoJi.at(i)*eddl->vvel->d.at(i);
            uRhoJ.at(2)    += intRhoJi.at(i)*eddl->wvel->d.at(i);
        }
        else {       // mass flux, not mass
            rhoK        += eddl->uvel->d.at(i) * intRhoKi.at(i);
            rhoJ        += eddl->uvel->d.at(i) * abs(intRhoKi.at(i));
            rhoKK       += eddl->uvel->d.at(i) * K.at(i)*intRhoKi.at(i);
            rhoJK       += eddl->uvel->d.at(i) * K.at(i)*intRhoJi.at(i);
            uRhoK.at(0)    += eddl->uvel->d.at(i) * intRhoKi.at(i)*eddl->uvel->d.at(i);
            uRhoK.at(1)    += eddl->uvel->d.at(i) * intRhoKi.at(i)*eddl->vvel->d.at(i);
            uRhoK.at(2)    += eddl->uvel->d.at(i) * intRhoKi.at(i)*eddl->wvel->d.at(i);
            uRhoJ.at(0)    += eddl->uvel->d.at(i) * intRhoJi.at(i)*eddl->uvel->d.at(i);
            uRhoJ.at(1)    += eddl->uvel->d.at(i) * intRhoJi.at(i)*eddl->vvel->d.at(i);
            uRhoJ.at(2)    += eddl->uvel->d.at(i) * intRhoJi.at(i)*eddl->wvel->d.at(i);
        }
    }

    ////////////

    A = rhoK/rhoJ;
    S = 0.5*(A*A+1.0)*rhoKK - A*rhoJK;
    for(i=0; i<3; i++) {
        P.at(i) = uRhoK.at(i) - A * uRhoJ.at(i);
        Q.at(i) = 0.25*P.at(i)*P.at(i)/S;
    }

    //////////// Compute Kernel Coefficients

    cCoef.at(0) = 0.5/S * (-P.at(0) + (P.at(0)>0 ? 1.0 : -1.0)
                 * sqrt( (1-domn->pram->A_param)*P.at(0)*P.at(0)
                         + 0.5*domn->pram->A_param*(P.at(1)*P.at(1)+P.at(2)*P.at(2)) ));
    cCoef.at(1) = 0.5/S * (-P.at(1) + (P.at(1)>0 ? 1.0 : -1.0)
                 * sqrt( (1-domn->pram->A_param)*P.at(1)*P.at(1)
                         + 0.5*domn->pram->A_param*(P.at(0)*P.at(0)+P.at(2)*P.at(2)) ));
    cCoef.at(2) = 0.5/S * (-P.at(2) + (P.at(2)>0 ? 1.0 : -1.0)
                 * sqrt( (1-domn->pram->A_param)*P.at(2)*P.at(2)
                         + 0.5*domn->pram->A_param*(P.at(0)*P.at(0)+P.at(1)*P.at(1)) ));

    for(i=0; i<3; i++)
        bCoef.at(i) = -A*cCoef.at(i);

}

///////////////////////////////////////////////////////////////////////////////
/** Eddy acceptance probability is computed after several previous steps.
 * @param dtSample \input eddy sample time step.
 **/
void eddy::computeEddyAcceptanceProb(const double dtSample) {

    double f, g;

    // f = esdp4/(eddySize*eddySize) * exp(esdp1/eddySize);
    // Pa = dtSample * invTauEddy * domn->pram->C_param / (eddySize * eddySize * f * g);
    // note cancellation of eddySize * eddySize

    f = esdp4 * exp(esdp1/eddySize);

    if(!domn->pram->Lperiodic)
        g = 1.0/(domn->Ldomain() - eddySize);
    else
        g = 1.0/domn->Ldomain();

    Pa = dtSample * invTauEddy * domn->pram->C_param / (f * g);

}

///////////////////////////////////////////////////////////////////////////////
/** Apply kernels K and J to the velocity profile.
 * This is called after the kernel coefficients is computed in eddyTau
 */
void eddy::applyVelocityKernels(domain *line, const int iS, const int iE) {

    ////////// update velocity profiles

    int npts = iE-iS+1;

    set_kernel_coefficients();

    if(domn->pram->Lspatial) {

        bool Lflag = false;                  // if kernel results in negative (or < umin) uvel for Lspatial, then don't apply kernels.
        for(int i=0; i<npts; i++)            //     a better approach would be to find the limiting alpha and use that.
            if ( line->uvel->d.at(i+iS) + cCoef.at(0)*K.at(i) + bCoef.at(0)*abs(K.at(i)) < domn->pram->umin_spatial ) {
                Lflag = true;
                break;
            }
        if(!Lflag) {
            for(int i=0; i<npts; i++) {
                line->uvel->d.at(i+iS) += cCoef.at(0)*K.at(i) + bCoef.at(0)*abs(K.at(i));
                line->vvel->d.at(i+iS) += cCoef.at(1)*K.at(i) + bCoef.at(1)*abs(K.at(i));
                line->wvel->d.at(i+iS) += cCoef.at(2)*K.at(i) + bCoef.at(2)*abs(K.at(i));
            }
        }
    }
    else {
        for(int i=0; i<npts; i++) {
            line->uvel->d.at(i+iS) += cCoef.at(0)*K.at(i) + bCoef.at(0)*abs(K.at(i));
            line->vvel->d.at(i+iS) += cCoef.at(1)*K.at(i) + bCoef.at(1)*abs(K.at(i));
            line->wvel->d.at(i+iS) += cCoef.at(2)*K.at(i) + bCoef.at(2)*abs(K.at(i));
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Fill velocity kernel K (used also for \fun{J=|K|})
 *  this applies the planar analytic definition
 */

void eddy::fillKernel_planarAnalytic() {

    K.resize(eddl->ngrd);
    double L  = eddySize;
    double y0 = leftEdge;
    double y;

    for(int i=0; i<eddl->ngrd; i++) {

        y = eddl->pos->d.at(i);

        if(y >= y0 && y <= y0 + L/3)
            K.at(i) = y- (y0 + 3*(y-y0));
        else if(y >= y0 + L/3 && y <= y0 + 2*L/3)
            K.at(i) = y- (y0 + 2*L - 3*(y-y0));
        else if(y >= y0 + 2*L/3 && y <= y0 + L)
            K.at(i) = y- (y0 + 3*(y-y0)-2*L);
        else
            K.at(i) = 0.0;
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Fill velocity kernel K (used also for \fun{J=|K|})
 *  this defines the kernel in terms of eddy displacements
 */

void eddy::fillKernel() {

    int nseg = eddl->ngrd / 3;

    K.resize(eddl->ngrd);

    //---------- 1st Segment

    K.at(0) = 0.5*(eddl->posf->d.at(1) + leftEdge) - pos0.at(0);  // pos is 0.5*(faces), but we want pos=0.5*(face+eddy_edge)
    for(int i=1, j=1; i<nseg; i++, j++)                         //    pos0 was already updated this way in tripmap
        K.at(i) = eddl->pos->d.at(i) - pos0.at(j);

    //---------- Second Segment

    for(int i=nseg*2-1, j=0; i>=nseg; i--, j++)
        K.at(i) = eddl->pos->d.at(i) - pos0.at(j);

    //---------- Third Segment

    for(int i=nseg*2, j=0; i<nseg*3-1; i++, j++)
        K.at(i) = eddl->pos->d.at(i) - pos0.at(j);
    int i=nseg*3-1, j=nseg-1;
    K.at(i) = 0.5*(rightEdge + eddl->posf->d.at(i)) - pos0.at(j); // pos is 0.5*(faces), but we want pos=0.5*(face+eddy_edge)

}

///////////////////////////////////////////////////////////////////////////////
/** Get the Favre average streamwise velocity for spatial formulations.
*  Used to convert from eddy timescale to the eddy spatial scale.
*
*  @param dxc \input cell volumes
*  @return value of the favre avg velocity
*/

double eddy::eddyFavreAvgVelocity(const vector<double> &dxc) {

    double ufavg = 0.0;
    double ravg  = 0.0;
    double dmb;
    for(int i=0; i<eddl->ngrd; i++) {
        dmb = eddl->rho->d.at(i) * dxc.at(i);
        ufavg += dmb*eddl->uvel->d.at(i);
        ravg  += dmb;
    }
    return ufavg/ravg;
}

