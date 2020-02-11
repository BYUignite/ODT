from __future__ import division
import numpy as np
#import yaml
from scipy.integrate import odeint
import sys
import streams
sys.path.append('c_rates')
#import getrates_py


###########################################################################################

class hips :

    #----------------------------------------------------------------------------------------------------------

    def __init__(s) :

        #with open("../input/hips.yaml",'r') as yfile:            # open the input file for reading
        #    yfile = yaml.load(yfile)

        s.Nlevels = 9           # yfile["Nlevels"]                             # number of levels in the tree: 0, 1, 2, ... N-1
        s.L0      = 1.0         # yfile["L0"]                                  # tree lengthscale
        s.tau0    = 1.0         # yfile["tau0"]                                # integral timescale
        #s.Da1     = 100         # yfile["Da1"]
        #s.Da2     = 100         # yfile["Da2"]
        s.Cparam  = 1.0         # yfile["Cparam"]                              # eddy rate parameter
        s.tRun    = 40.0        # yfile["tRun"]                                # run time
        s.fMix    = 1.0         # yfile["fMix"]                                # run time
        s.cfl     = 1.0         # yfile["cfl"]                                 # dt = cfl * tMix
        s.Afac    = 0.5         # yfile["Afac"]                                # level lengthscale reduction factor (0.5)
        s.Lrxn    = True        # yfile["Lrxn"]                                # true if reaction is on
        s.myMech  = "cold_rxn"    # yfile["myMech"] if "myMech" in yfile else "canteraRR"   # string: 4step_ch4, oneStep_c2h4, red_c2h4
        s.dtDump  = 0.100       # yfile["dtDump"]

        s.nParcels = int(2**(s.Nlevels-1))
        s.Nm1      = s.Nlevels-1
        s.Nm3      = s.Nlevels-3

        s.strm     = streams.streams()

        if len(sys.argv) < 3:
            print('error in # args')
            exit(0)
        s.myid = sys.argv[1]
        s.Da1 = float(sys.argv[2])
        s.Da2 = float(sys.argv[3])

        #------------------- Set the level lengthscales, timescales, rates

        i = np.arange(s.Nlevels)                          # index array

        s.Llevels = np.ones(s.Nlevels)*s.L0*s.Afac**i     # lengthscale at each level (per node)
        s.taus    = np.ones(s.Nlevels)*s.tau0* \
                    (s.Llevels/s.L0)**(2.0/3.0)/s.Cparam  # tau at each level (per node)
        s.rates   = np.ones(s.Nlevels)/s.taus*2**i        # lambda rate at each level (for all nodes) (note, last two levels not used for eddy events)
        s.tMix    = s.fMix * s.taus[-1]                   # parcel mixing timescale

        #------------------- Set the parcel addresses (index array)

        s.pLoc = np.arange(s.nParcels,dtype='uint64')

        #------------------- set the eddy event times at all levels

        s.setEddyEventTimes()

        #------------------- set output times

        if s.dtDump > 0:
            s.tDump = np.arange(s.dtDump, s.tRun, s.dtDump)
            s.tDump = np.append(s.tDump, s.tRun)
            s.itDumpNext = 0

        #------------------- Initialize the parcel properties

        s.mixf = np.zeros(s.nParcels)
        s.enth = np.zeros(s.nParcels)
        s.ysp  = np.zeros((s.nParcels,s.strm.gas.n_species))

        #------------------- Set the parcel properties

        s.strm.set_gas_mixf_eq(0.0, doEq=False)
        s.enth[:int(s.nParcels/2)] = s.strm.gas.h
        s.ysp[:int(s.nParcels/2),:] = s.strm.gas.Y

        s.strm.set_gas_mixf_eq(1.0, doEq=False)
        s.enth[int(s.nParcels/2):] = s.strm.gas.h
        s.ysp[int(s.nParcels/2):,:] = s.strm.gas.Y

        s.mixf[int(s.nParcels/2):] = 1.0              # simple step profile:  00000001111111

    #----------------------------------------------------------------------------------------------------------

    def setEddyEventTimes(s) :
        """
        Sets arrays s.eTimes and s.eLevels which hold the eddy event times and the corresponding tree base level.
        First make an array of arrays for times at each level.
            These are from a Poisson process with the given rate at each level.
        Then collapse these into a single array.
        Then sort these times, along with the corresponding level array.
        """

        eTimes  = np.empty(0)                                     # each element is an array of times at the given level
        eLevels = np.empty(0, dtype='int64')                      # each element is an array of levels. [000000000],[1111111],[2222], etc. (see below)
        for i,rate in enumerate(s.rates[:-2]) :                   # don't do the last three levels since these can't be the base node for a swap
            nEapprox = int(s.tRun*rate)
            r = np.random.rand(nEapprox * 10)
            times = np.cumsum(-np.log(r)/rate)
            times = times[np.where(times < s.tRun)]
            eTimes = np.append(eTimes, times)
            eLevels = np.append(eLevels, np.ones(np.size(times), dtype='int64') * i)

        eTimes  = np.hstack(eTimes)
        eLevels = np.hstack(eLevels)

        isort = np.argsort(eTimes)
        s.eTimes = eTimes[isort]
        s.eLevels = eLevels[isort]

    #----------------------------------------------------------------------------------------------------------

    def select_and_swap_two_subtrees(s, iLevel) :
        """
        Pass in the level of the tree for the base of the swap: iLevel.
        Randomly select a node on iLevel.
        Go down two levels and select nodes 0q and 1r, where q, r are randomly 0 or 1
        Find the starting index of the Q-tree and R-tree to swap and the number of parcels.
        Then swap the cells.
        Return the starting index for the Q-tree and R-tree, and the number of parcels swaped.

        For a 6 level tree: 0, 1, 2, 3, 4, 5:
        If iLevel = 1, then suppose i=1, 0q = 00 and 1r = 11:
        Then we are swaping 0100** with 0111** or (01|00|**) with (01|11|**)
           or i0qs with i1rs, where i = 01; 0q = 00; 1r = 11; and s = **

        We use bitwise shifts for easy powers of 2.
        The swap is done by adding or subtracting a value (shift),
            which should be equivalent to flipping the swapping the two 0q bits and 1r bits.
                                                                                                                  Level
                                                                                                                ---------
                                                        *                                                           0
                                                     /     \
                                                  /           \
                                               /                 \
                                            /                       \
                                         /                             \
                                      /                                   \
                                   /                                         \
                                /                                               \
                               *                                                (*)  01|0000                        1
                              / \                                               / \
                            /     \                                           /     \
                          /         \                                       /         \
                        /             \                                   /             \
                      /                 \                               /                 \
                    /                     \                           /                     \
                   *                       *                         *                       *                      2
                  / \                     / \                       / \                     / \
                /     \                 /     \                   /     \                 /     \
              /         \             /         \               /         \             /         \
             *           *           *           *            [*] 00|**    *           *          [*] 11|**         3
            / \         / \         / \         / \           / \         / \         / \         / \
           /   \       /   \       /   \       /   \         /   \       /   \       /   \       /   \
          *     *     *     *     *     *     *     *       *     *     *     *     *     *     *     *             4
         / \   / \   / \   / \   / \   / \   / \   / \     / \   / \   / \   / \   / \   / \   / \   / \
        00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31           5
                                                          ^^^^^^^^^^^                         ^^^^^^^^^^^
        """

        iNode = np.random.randint(0, 1<<iLevel)                       # base node for the swap
        if iLevel == s.Nm3 :                                           # if swapping parcels (not nodes),
            zero_q = 0                                                #     it doesn't matter which of the two we choose
            one_r  = 2                                                #     so choose the left to make simple mixing easier
        else :
            zero_q = np.random.randint(0,2)                               # 0q where q is 0 or 1
            one_r  = 2+np.random.randint(0,2)                             # 1r where r is 0 or 1

        Qstart = (zero_q << s.Nm3-iLevel) + (iNode << s.Nm1-iLevel)   # starting index of Q parcels
        Rstart = (one_r  << s.Nm3-iLevel) + (iNode << s.Nm1-iLevel)   # starting index of R parcels
        nPswap = 1 << s.Nm3-iLevel                                    # number of parcels that will be swapped

        Qend = Qstart + nPswap
        Rend = Rstart + nPswap

        aa = s.pLoc[Qstart:Qend].copy()
        s.pLoc[Qstart:Qend] = s.pLoc[Rstart:Rend]
        s.pLoc[Rstart:Rend] = aa

        return Qstart, Rstart, nPswap                                 # aids in applying local mixing processes

    #----------------------------------------------------------------------------------------------------------

    def stepper(s) :
        """
        Perform each eddy event in the list.
        If the eddy event involved single parcels, then perform a simple mixing operation.
        Collect and save the mixture fraction variable to an output file.
        """

        s.mixfAll = s.mixf.copy()
        s.YaAll   = s.ysp[:,0]
        s.YbAll   = s.ysp[:,1]
        s.YrAll   = s.ysp[:,2]
        s.YpAll   = s.ysp[:,3]


        tlast = 0.0
        for i, time in enumerate(s.eTimes) :

            #print("EE at time %g" %(time))
            (QS, RS, nPs) = s.select_and_swap_two_subtrees(s.eLevels[i])

            if s.Lrxn :
                s.mix_react_all_parcels(time-tlast)
            else :
                if nPs == 1 : s.simple_mixing_local(QS, RS)

            if time >= s.tDump[s.itDumpNext] or s.dtDump < 0.0:
                s.itDumpNext += 1
                print("time = ", time)

                s.mixfAll = np.vstack([s.mixfAll, s.mixf[s.pLoc]])
                s.YaAll   = np.vstack([s.YaAll,   s.ysp[s.pLoc,0]])
                s.YbAll   = np.vstack([s.YbAll,   s.ysp[s.pLoc,1]])
                s.YrAll   = np.vstack([s.YrAll,   s.ysp[s.pLoc,2]])
                s.YpAll   = np.vstack([s.YpAll,   s.ysp[s.pLoc,3]])

            tlast = time

        s.mixfAll = np.vstack([s.mixfAll, s.mixf[s.pLoc]])
        s.YaAll   = np.vstack([s.YaAll,   s.ysp[s.pLoc,0]])
        s.YbAll   = np.vstack([s.YbAll,   s.ysp[s.pLoc,1]])
        s.YrAll   = np.vstack([s.YrAll,   s.ysp[s.pLoc,2]])
        s.YpAll   = np.vstack([s.YpAll,   s.ysp[s.pLoc,3]])

        if s.dtDump < 0.0:
            header = ""
        else:
            header = "1_0.0, " + ", ".join([str(i+2)+"_"+str(t) for i,t in enumerate(s.tDump)])
        np.savetxt("../data/mixf_"+str(s.myid) +".out", s.mixfAll.T, fmt='%.5E', header=header)
        np.savetxt("../data/ya_"+str(s.myid)+".out",   s.YaAll.T,   fmt='%.5E', header=header)
        np.savetxt("../data/yb_"+str(s.myid)+".out",   s.YbAll.T,   fmt='%.5E', header=header)
        np.savetxt("../data/yr_"+str(s.myid)+".out",   s.YrAll.T,   fmt='%.5E', header=header)
        np.savetxt("../data/yp_"+str(s.myid)+".out",   s.YpAll.T,   fmt='%.5E', header=header)

        print("done with myid=", s.myid)

    #----------------------------------------------------------------------------------------------------------

    def simple_mixing_local(s, QS, RS) :

        mix = 0.5*(s.mixf[s.pLoc[QS]] + s.mixf[s.pLoc[QS+1]])
        s.mixf[s.pLoc[QS:QS+2]] = mix

        mix = 0.5*(s.mixf[s.pLoc[RS]] + s.mixf[s.pLoc[RS+1]])
        s.mixf[s.pLoc[RS:RS+2]] = mix

        for k in range(4):

            mix = 0.5*(s.ysp[s.pLoc[QS],k] + s.ysp[s.pLoc[QS+1],k])
            s.ysp[s.pLoc[QS:QS+2],k] = mix

            mix = 0.5*(s.ysp[s.pLoc[RS],k] + s.ysp[s.pLoc[RS+1],k])
            s.ysp[s.pLoc[RS:RS+2],k] = mix

    #----------------------------------------------------------------------------------------------------------

    def mix_react_all_parcels(s, tadv) :
        """
        Strang splitting: half step mixing, full step reaction, half step mixing.
        Equations are dp_me/dt = (p_nb - p_me)/tMix + S.
        Here, p is some property, p_nb is the property at the neighboring parcel, and
            p_me is the property at the node in question.
        Solve the mixing term analytically.
        """

        s.ip_me = s.pLoc[np.arange(s.nParcels)]              # parcel indices
        s.ip_nb = np.empty(s.nParcels)               # parcel neighbor indicies
        for i in range(s.nParcels) :
            s.ip_nb[i] = s.pLoc[i+1] if s.pLoc[i]%2 == 0 else s.pLoc[i-1]


        s.ip_nb = s.ip_nb.astype(int)                             # converts ip_nb to integers
        t = 0.0
        dt  = s.tMix * s.cfl
        dt2 = dt*0.5
        keepOnSteppin = True

        while keepOnSteppin :
            if t+dt >= tadv :
                dt = tadv - t
                dt2 = dt*0.5
                keepOnSteppin = False

            #----------  Diffuse a half step

            s.mixf[s.ip_me] = s.mixf[s.ip_me] + (1.0-np.exp(-2.0*dt2/s.tMix)) * 0.5 * \
                                       (s.mixf[s.ip_nb] - s.mixf[s.ip_me])

            s.ysp[s.ip_me,:]  = s.ysp[s.ip_me,:] + (1.0-np.exp(-2.0*dt2/s.tMix)) * 0.5 * \
                                        (s.ysp[s.ip_nb,:] - s.ysp[s.ip_me,:])

            if s.myMech != "cold_rxn":
                s.enth[s.ip_me] = s.enth[s.ip_me] + (1.0-np.exp(-2.0*dt2/s.tMix)) * 0.5 * \
                                           (s.enth[s.ip_nb] - s.enth[s.ip_me])

            #----------  React a full step

            for s.ip_now in range(s.nParcels) :
                y0 = s.ysp[s.ip_now,:]
                times = np.array([0,dt])
                s.ysp[s.ip_now,:] = odeint(s.rhsf_rxr, y0, times, mxstep=10000)[1,:]

            #----------  Diffuse a half step

            s.mixf[s.ip_me] = s.mixf[s.ip_me] + (1.0-np.exp(-2.0*dt2/s.tMix)) * 0.5 * \
                                       (s.mixf[s.ip_nb] - s.mixf[s.ip_me])

            s.ysp[s.ip_me,:]  = s.ysp[s.ip_me,:] + (1.0-np.exp(-2.0*dt2/s.tMix)) * 0.5 * \
                                        (s.ysp[s.ip_nb,:] - s.ysp[s.ip_me,:])

            if s.myMech != "cold_rxn":
                s.enth[s.ip_me] = s.enth[s.ip_me] + (1.0-np.exp(-2.0*dt2/s.tMix)) * 0.5 * \
                                           (s.enth[s.ip_nb] - s.enth[s.ip_me])

            #----------

            t = t+dt

    #----------------------------------------------------------------------------------------------------------

    def rhsf_rxr(s, y, t) :

        s.strm.gas.HPY = s.enth[s.ip_now], s.strm.pres, y                       # changed strm.P to strm.pres
        rr = np.empty(s.strm.gas.n_species)                                     # changed s.gas.n_species to s.strm.gas.species

        #if s.myMech == "4step_ch4" :
        #    getrates_py.getrates_py_4step_ch4( s.strm.gas.density, s.strm.gas.T, s.strm.gas.P, s.strm.gas.U, rr )
        #if s.myMech == "oneStep_c2h4" :
        #    getrates_py.getrates_py_oneStep_c2h4( s.strm.gas.density, s.strm.gas.T, s.strm.gas.P, s.strm.gas.U, rr )
        #if s.myMech == "red_c2h4" :
        #    getrates_py.getrates_py_red_c2h4( s.strm.gas.density, s.strm.gas.T, s.strm.gas.P, s.strm.gas.U, rr )
        if s.myMech == "cold_rxn":
            rr = np.empty(s.strm.gas.n_species)                                     # changed s.gas.n_species to s.strm.gas.species
            iA = 0
            iB = 1
            iR = 2
            iP = 3

            rr[iA] = -s.Da1*y[iA]*y[iB] - 0.5*s.Da2*y[iA]*y[iR]
            rr[iB] = -s.Da1*y[iA]*y[iB]
            rr[iR] =  2.0*s.Da1*y[iA]*y[iB] - s.Da2*y[iA]*y[iR]
            rr[iP] =                      1.5*s.Da2*y[iA]*y[iR]

            return rr

        else :
            rr = s.strm.gas.net_production_rates

        return rr * s.strm.gas.molecular_weights / s.strm.gas.density


###########################################################################################
# RUN THE CODE

h = hips()       # setup a HiPS object

h.stepper()      # step through the eddy events

