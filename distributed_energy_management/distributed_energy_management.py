# The distributed energy management system

from numpy import zeros,vstack

def problem_formulation_local_recovery(*args):
    from distributed_energy_management.configuration import configuration_time_line
    from distributed_energy_management.modelling.power_flow.idx_ed_recovery_format import PG, RG, PUG, RUG, PBIC_AC2DC, PBIC_DC2AC, \
        PESS_C, PESS_DC, RESS, EESS, PMG, PPV, PWP, PL_AC, PL_UAC, PL_DC, PL_UDC, NX

    model = args[0]  # If multiple models are inputed, more local ems models will be formulated
    ## The infeasible optimal problem formulation
    T = configuration_time_line.default_look_ahead_time_step["Look_ahead_time_ed_time_step"]
    nx = T * NX
    lb = zeros(nx)
    ub = zeros(nx)

    for i in range(T):
        ## Update lower boundary
        lb[i * NX + PG] = model["DG"]["PMIN"]
        lb[i * NX + RG] = model["DG"]["PMIN"]
        lb[i * NX + PUG] = model["UG"]["PMIN"]
        lb[i * NX + RUG] = model["UG"]["PMIN"]
        lb[i * NX + PBIC_AC2DC] = 0
        lb[i * NX + PBIC_DC2AC] = 0
        lb[i * NX + PESS_C] = 0
        lb[i * NX + PESS_DC] = 0
        lb[i * NX + RESS] = 0
        lb[i * NX + EESS] = model["ESS"]["SOC_MIN"] * model["ESS"]["CAP"]
        lb[
            i * NX + PMG] = 0  # The line flow limitation, the predefined status is, the transmission line is off-line
        lb[i * NX + PPV] = 0
        lb[i * NX + PWP] = 0
        lb[i * NX + PL_AC] = 0
        lb[i * NX + PL_UAC] = 0
        lb[i * NX + PL_DC] = 0
        lb[i * NX + PL_UDC] = 0
        ## Update lower boundary
        ub[i * NX + PG] = model["DG"]["PMAX"]
        ub[i * NX + RG] = model["DG"]["PMAX"]
        ub[i * NX + PUG] = model["UG"]["PMAX"]
        ub[i * NX + RUG] = model["UG"]["PMAX"]
        ub[i * NX + PBIC_AC2DC] = model["BIC"]["CAP"]
        ub[i * NX + PBIC_DC2AC] = model["BIC"]["CAP"]
        ub[i * NX + PESS_C] = model["ESS"]["PMAX_CH"]
        ub[i * NX + PESS_DC] = model["ESS"]["PMAX_DIS"]
        ub[i * NX + RESS] = model["ESS"]["PMAX_DIS"] + model["ESS"]["PMAX_CH"]
        ub[i * NX + EESS] = model["ESS"]["SOC_MAX"] * model["ESS"]["CAP"]
        ub[
            i * NX + PMG] = 0  # The line flow limitation, the predefined status is, the transmission line is off-line
        ub[i * NX + PPV] = model["PV"]["PG"][i]
        ub[i * NX + PWP] = model["WP"]["PG"][i]
        ub[i * NX + PL_AC] = model["Load_ac"]["PD"][i]
        ub[i * NX + PL_UAC] = model["Load_uac"]["PD"][i]
        ub[i * NX + PL_DC] = model["Load_dc"]["PD"][i]
        ub[i * NX + PL_UDC] = model["Load_udc"]["PD"][i]

    ## Constraints set
    # 1) Power balance equation
    Aeq = zeros((T, nx))
    beq = []
    for i in range(T):
        Aeq[i][i * NX + PG] = 1
        Aeq[i][i * NX + PUG] = 1
        Aeq[i][i * NX + PBIC_AC2DC] = -1
        Aeq[i][i * NX + PBIC_DC2AC] = model["BIC"]["EFF_DC2AC"]
        Aeq[i][i * NX + PL_AC] = -1
        Aeq[i][i * NX + PL_UAC] = -1
        beq.append(0)
    # 2) DC power balance equation
    Aeq_temp = zeros((T, nx))
    for i in range(T):
        Aeq_temp[i][i * NX + PBIC_AC2DC] = model["BIC"]["EFF_AC2DC"]
        Aeq_temp[i][i * NX + PBIC_DC2AC] = -1
        Aeq_temp[i][i * NX + PESS_C] = -1
        Aeq_temp[i][i * NX + PESS_DC] = 1
        Aeq_temp[i][i * NX + PMG] = -1
        Aeq_temp[i][i * NX + PL_DC] = -1
        Aeq_temp[i][i * NX + PL_UDC] = -1
        Aeq_temp[i][i * NX + PPV] = 1
        Aeq_temp[i][i * NX + PWP] = 1
        beq.append(0)
    Aeq = vstack([Aeq, Aeq_temp])

    # 3) Energy storage system
    # 3) Energy storage system
    Aeq_temp = zeros((T, nx))
    for i in range(T):
        if i == 0:
            Aeq_temp[i][i * NX + EESS] = 1
            Aeq_temp[i][i * NX + PESS_C] = -model["ESS"]["EFF_CH"] * configuration_time_line.default_time[
                "Time_step_ed"] / 3600
            Aeq_temp[i][i * NX + PESS_DC] = 1 / model["ESS"]["EFF_DIS"] * configuration_time_line.default_time[
                "Time_step_ed"] / 3600
            beq.append(model["ESS"]["SOC"] * model["ESS"]["CAP"])
        else:
            Aeq_temp[i][(i - 1) * NX + EESS] = -1
            Aeq_temp[i][i * NX + EESS] = 1
            Aeq_temp[i][i * NX + PESS_C] = -model["ESS"]["EFF_CH"] * configuration_time_line.default_time[
                "Time_step_ed"] / 3600
            Aeq_temp[i][i * NX + PESS_DC] = 1 / model["ESS"]["EFF_DIS"] * configuration_time_line.default_time[
                "Time_step_ed"] / 3600
            beq.append(0)
    Aeq = vstack([Aeq, Aeq_temp])

    # Inequality constraints
    # Inequality constraints
    # 1) PG + RG <= PGMAX
    Aineq = zeros((T, nx))
    bineq = []
    for i in range(T):
        Aineq[i][i * NX + PG] = 1
        Aineq[i][i * NX + RG] = 1
        bineq.append(model["DG"]["PMAX"])
    # 2) PG - RG >= PGMIN
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + PG] = -1
        Aineq_temp[i][i * NX + RG] = 1
        bineq.append(-model["DG"]["PMIN"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 3) PUG + RUG <= PUGMAX
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + PUG] = 1
        Aineq_temp[i][i * NX + RUG] = 1
        bineq.append(model["UG"]["PMAX"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 4) PUG - RUG >= PUGMIN
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + PUG] = -1
        Aineq_temp[i][i * NX + RUG] = 1
        bineq.append(-model["UG"]["PMIN"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 5) PESS_DC - PESS_C + RESS <= PESS_DC_MAX
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + PESS_DC] = 1
        Aineq_temp[i][i * NX + PESS_C] = -1
        Aineq_temp[i][i * NX + RESS] = 1
        bineq.append(model["ESS"]["PMAX_DIS"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 6) PESS_DC - PESS_C - RESS >= -PESS_C_MAX
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + PESS_DC] = -1
        Aineq_temp[i][i * NX + PESS_C] = 1
        Aineq_temp[i][i * NX + RESS] = 1
        bineq.append(model["ESS"]["PMAX_CH"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 7) EESS - RESS*delta >= EESSMIN
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + EESS] = -1
        Aineq_temp[i][i * NX + RESS] = configuration_time_line.default_time["Time_step_ed"] / 3600
        bineq.append(-model["ESS"]["SOC_MIN"] * model["ESS"]["CAP"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 8) EESS + RESS*delta <= EESSMAX
    Aineq_temp = zeros((T, nx))
    for i in range(T):
        Aineq_temp[i][i * NX + EESS] = 1
        Aineq_temp[i][i * NX + RESS] = configuration_time_line.default_time["Time_step_ed"] / 3600
        bineq.append(model["ESS"]["SOC_MAX"] * model["ESS"]["CAP"])
    Aineq = vstack([Aineq, Aineq_temp])
    # 9) RG + RUG + RESS >= sum(Load)*beta + sum(PV)*beta_pv + sum(WP)*beta_wp

    # No reserve requirement
    c = zeros(NX)
    if model["DG"]["COST_MODEL"] == 2:
        c[PG] = model["DG"]["COST"][1]
    else:
        c[PG] = model["DG"]["COST"][0]
    c[PUG] = model["UG"]["COST"][0]
    c[PESS_C] = model["ESS"]["COST_CH"][0]
    c[PESS_DC] = model["ESS"]["COST_DIS"][0]
    # The sheding cost
    c[PPV] = -model["PV"]["COST"]
    c[PWP] = -model["WP"]["COST"]
    c[PL_AC] = -model["Load_ac"]["COST"][0]
    c[PL_UAC] = -model["Load_uac"]["COST"][0]
    c[PL_DC] = -model["Load_dc"]["COST"][0]
    c[PL_UDC] = -model["Load_udc"]["COST"][0]

    lb_list = [ ]
    ub_list = [ ]
    c_list = [ ]
    for i in range(NX):
        lb_list.append(lb[i])
        ub_list.append(ub[i])
        c_list.append(c[i])

    C = c_list * T
    # Generate the quadratic parameters
    Q = zeros((nx, nx))
    for i in range(T):
        if model["DG"]["COST_MODEL"] == 2:
            Q[i * NX + PG][i * NX + PG] = model["DG"]["COST"][1]
    mathematical_model = {"Q": Q,
                          "c": C,
                          "Aeq": Aeq,
                          "beq": beq,
                          "A": Aineq,
                          "b": bineq,
                          "lb": lb,
                          "ub": ub}

    return mathematical_model


