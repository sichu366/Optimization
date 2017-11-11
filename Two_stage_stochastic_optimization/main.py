from Two_stage_stochastic_optimization.power_flow_modelling import idx_bus
from os.path import dirname, join
from pypower.loadcase import loadcase
from pypower.ext2int import ext2int
from numpy import zeros, c_, shape, ix_,ones,r_,arange,sum
from numpy import flatnonzero as find
from scipy.sparse.linalg import inv

from pypower.idx_brch import F_BUS, T_BUS, BR_X, TAP, SHIFT, BR_STATUS,RATE_A
from pypower.idx_cost import MODEL, NCOST, PW_LINEAR, COST, POLYNOMIAL
from pypower.idx_bus import BUS_TYPE, REF, VA, VM, PD, GS, VMAX, VMIN
from pypower.idx_gen import GEN_BUS, VG, PG, QG, PMAX, PMIN, QMAX, QMIN

from scipy.sparse import csr_matrix as sparse

from Two_stage_stochastic_optimization.solvers.mix_integer_solvers import miqp_gurobi# The gurobi solver

def optimal_power_flow(*args):
    casedata = args[0] # Target power flow modelling
    mpc = loadcase(casedata) # Import the power flow modelling
    ## convert to internal indexing
    mpc = ext2int(mpc)
    baseMVA, bus, gen, branch = mpc["baseMVA"], mpc["bus"], mpc["gen"], mpc["branch"] #

    nb = shape(mpc['bus'])[0]  ## number of buses
    nl = shape(mpc['branch'])[0]  ## number of branches
    ng = shape(mpc['gen'])[0]  ## number of dispatchable injections

    ## Formualte the
    stat = branch[:, BR_STATUS]  ## ones at in-service branches
    b = stat / branch[:, BR_X]  ## series susceptance
    tap = ones(nl)  ## default tap ratio = 1
    i = find(branch[:, TAP])  ## indices of non-zero tap ratios
    tap[i] = branch[i, TAP]  ## assign non-zero tap ratios

    ## build connection matrix Cft = Cf - Ct for line and from - to buses
    f = branch[:, F_BUS]  ## list of "from" buses
    t = branch[:, T_BUS]  ## list of "to" buses
    i = r_[range(nl), range(nl)]  ## double set of row indices
    ## connection matrix
    Cft = sparse((r_[ones(nl), -ones(nl)], (i, r_[f, t])), (nl, nb))

    ## build Bf such that Bf * Va is the vector of real branch powers injected
    ## at each branch's "from" bus
    Bf = sparse((r_[b, -b], (i, r_[f, t])), shape=(nl, nb))  ## = spdiags(b, 0, nl, nl) * Cft

    ## build Bbus
    Bbus = Cft.T * Bf
    # The distribution factor
    Distribution_factor = Bf*inv(Bbus)

    Cg = sparse((ones(ng), (gen[:, GEN_BUS], arange(ng))), (nb, ng)) # Sparse index generation method is different from the way of matlab

    Pd = sum(bus[:,PD])# Total power demand

    # Formulate the problem


    return mpc

if __name__=="__main__":
    from Two_stage_stochastic_optimization.power_flow_modelling.case30 import case30
    casedata = case30()
    result = optimal_power_flow(casedata)
    # print(result)