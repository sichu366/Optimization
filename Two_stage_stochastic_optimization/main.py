from Two_stage_stochastic_optimization.power_flow_modelling import idx_bus
from os.path import dirname, join
from pypower.loadcase import loadcase
from pypower.ext2int import ext2int
from numpy import zeros, c_, shape, ix_,ones
from numpy import flatnonzero as find
from pypower.idx_brch import F_BUS, T_BUS, BR_X, TAP, SHIFT, BR_STATUS

def optimal_power_flow(*args):
    casedata = args[0] # Target power flow modelling
    mpc = loadcase(casedata) # Import the power flow modelling
    ## convert to internal indexing
    mpc = ext2int(mpc)
    baseMVA, bus, gen, branch = mpc["baseMVA"], mpc["bus"], mpc["gen"], mpc["branch"] #

    nb = shape(mpc['bus'])[0]  ## number of buses
    nl = shape(mpc['branch'])[0]  ## number of branches
    ng = shape(mpc['gen'])[0]  ## number of dispatchable injections

    stat = branch[:, BR_STATUS]  ## ones at in-service branches
    b = stat / branch[:, BR_X]  ## series susceptance
    tap = ones(nl)  ## default tap ratio = 1
    i = find(branch[:, TAP])  ## indices of non-zero tap ratios
    tap[i] = branch[i, TAP]  ## assign non-zero tap ratios
    return mpc

if __name__=="__main__":
    from Two_stage_stochastic_optimization.power_flow_modelling.case30 import case30
    casedata = case30()
    result = optimal_power_flow(casedata)
    # print(result)