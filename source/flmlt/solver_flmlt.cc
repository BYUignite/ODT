
#include "solver_flmlt.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////
/**
 */

void solver_flmlt::init(domain *p_domn) {
    domn = p_domn;
}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 */

void solver_flmlt::calculateSolution() {

    domn->io->writeDataFile("odt_init.dat", 0.0);

    domn->mesher->adaptGrid(0, domn->ngrd-1);

    domn->io->writeDataFile("odt_init_adpt.dat", 0.0);

    domn->io->outputFlmltHeader();

    domn->mimx->advanceOdt(0.0, domn->pram->tEnd);

    if(domn->pram->heatloss == 0.0) {         
        domn->domc->setCaseSpecificVars();
        vector<double> hsens(domn->ngrd);
        domn->enth->set_hsens(hsens);
        domn->io->write_h_mixf_flmlt_profile(hsens);
    }

}
