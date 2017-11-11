from Two_stage_stochastic_optimization.power_flow_modelling import idx_bus
from os.path import dirname, join

def optimal_power_flow(*args):
    casedata = args[0] # Target power flow modelling



if __name__=="__main__":
    casedata = join(dirname(__file__), 'power_flow_modelling/case9')
    result = optimal_power_flow(casedata)