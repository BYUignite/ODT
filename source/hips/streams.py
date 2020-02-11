
##############################################################################
# Mixture Fraction Streams
##############################################################################

from __future__ import division
import cantera as ct
import numpy   as np
import yaml

class streams:

    #-----------------------------------------------------------------------------

    def __init__(self) :

        s = self

        with open("../input/hips.yaml",'r') as yfile:            # open the input file for reading
            yfile = yaml.load(yfile)

        mech         = yfile["streams"]["mechanism"]
        s.pres       = yfile["streams"]["pressure"]
        mole_or_mass = yfile["streams"]["mole_or_mass"]
        s.T0         = yfile["streams"]["T0"]
        s.T1         = yfile["streams"]["T1"]

        s.gas = ct.Solution("../input/"+mech)

        if mole_or_mass == 'mole' :
            s.gas.TPX = s.T0, s.pres, yfile["streams"]["comp0"]
            s.x0 = s.gas.X
            s.y0 = s.gas.Y
            s.h0 = s.gas.h
            s.gas.TPX = s.T1, s.pres, yfile["streams"]["comp1"]
            s.x1 = s.gas.X
            s.y1 = s.gas.Y
            s.h1 = s.gas.h
        else :
            s.gas.TPY = s.T0, s.pres, yfile["streams"]["comp0"]
            s.x0 = s.gas.X
            s.y0 = s.gas.Y
            s.h0 = s.gas.h
            s.gas.TPY = s.T1, s.pres, yfile["streams"]["comp1"]
            s.x1 = s.gas.X
            s.y1 = s.gas.Y
            s.h1 = s.gas.h

        s.compute_Zst()

    #-----------------------------------------------------------------------------

    def compute_Zst(self) :

        s = self

        s.gas.TPY = s.T1, s.pres, s.y1
        mix = ct.Mixture([(s.gas, 1.0)])

        nc1 = mix.element_moles("C")
        nh1 = mix.element_moles("H")
        no1 = mix.element_moles("O")
        nn1 = mix.element_moles("N")

        nfsp = mix.species_moles
        nf = sum(nfsp)

        no2_per_nfuel = (nc1 + nh1/4 - no1/2)/nf

        n1 = 1
        n0 = no2_per_nfuel/s.x0[s.gas.species_index("O2")]
        s.gas.TPY = s.T0, s.pres, s.y0
        M0 = s.gas.mean_molecular_weight
        s.gas.TPY = s.T1, s.pres, s.y1
        M1 = s.gas.mean_molecular_weight

        s.Zst = n1*M1/(n0*M0+n1*M1)

    #-----------------------------------------------------------------------------

    def set_gas_mixf_eq(self, mixf, doEq) :

        s = self

        h = s.h1*mixf + s.h0*(1.0-mixf)
        y = s.y1*mixf + s.y0*(1.0-mixf)
        s.gas.HPY = h, s.pres, y

        if doEq : s.gas.equilibrate("HP")






