
#include "solver_premix.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////
/** Constructor
 *  @param p_domn back pointer to domain object
 */

solver_premix::solver_premix(domain *p_domn) {
    domn = p_domn;
    LES_type3 = false;
}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 */

void solver_premix::calculateSolution() {

    domn->io->writeDataFile("premix_init.dat", 0.0);

    domn->mesher->adaptGrid(0, domn->ngrd-1);

    domn->io->writeDataFile("premix_init_adpt.dat", 0.0);

    domn->io->outputPremixHeader();

    domn->mimx->advanceOdt(0.0, domn->pram->tEnd);

}
