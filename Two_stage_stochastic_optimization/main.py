from Two_stage_stochastic_optimization.power_flow_modelling import idx_bus
from os.path import dirname, join
from pypower.loadcase import loadcase

def optimal_power_flow(*args):
    casedata = args[0] # Target power flow modelling
    mpc = loadcase(casedata) # Import the power flow modelling


    return mpc

if __name__=="__main__":
    from Two_stage_stochastic_optimization.power_flow_modelling.case30 import case30
    casedata = case30()
    result = optimal_power_flow(casedata)
    # print(result)