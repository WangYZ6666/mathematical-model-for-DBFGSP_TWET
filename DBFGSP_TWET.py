import sys

import gurobipy as gp
import numpy as np
from gurobipy import GRB, max_
import pandas as pd
from pygantt.pygantt import *
import pandas as pd
import xlwt
from xlwt import Workbook

sys.path.append("../")

try:
    # Read Instance
    data_csv = pd.read_csv("../Instances/DBFGSP_TWET_2x6x2.csv", skip_blank_lines=True, index_col=False)

    position_dict = table_find_pos(data_csv, ['Factories', 'Families', 'Machines', 'Total number of jobs',
                                              'Number of Jobs in each Family', 'Jobs in each Family',
                                              'Processing times of jobs', 'On machine', 'Time Window', 'Weight'])

    h = 0x0000ffff
    num_jobs = int(data_csv.iloc[position_dict['Total number of jobs'][0][1] + 1,
                                 data_csv.columns.get_loc(position_dict['Total number of jobs'][0][0])])
    num_machines = int(data_csv.iloc[position_dict['Machines'][0][1],
                                     data_csv.columns.get_loc(position_dict['Machines'][0][0]) + 1])
    num_groups = int(data_csv.iloc[position_dict['Families'][0][1],
                                   data_csv.columns.get_loc(position_dict['Families'][0][0]) + 1])
    num_factories = int(data_csv.iloc[position_dict['Factories'][0][1],
                                      data_csv.columns.get_loc(position_dict['Factories'][0][0]) + 1])

    data_row = position_dict['Jobs in each Family'][0][1] + 1

    jobs_in_each_group = [np.array([0])]
    for l in np.arange(0, num_groups):
        jobs_in_each_group.append(np.insert(np.array(data_csv.loc[data_row, :].dropna().astype(int)) + 1, 0, 0))
        data_row = data_row + 1

    jobs_in_each_group = np.array(jobs_in_each_group, dtype=object)

    data_row = position_dict['Processing times of jobs'][0][1] + 1

    process_time = np.array(data_csv.iloc[data_row: data_row + num_jobs,
                            0:num_machines].astype(int))

    process_time = np.insert(process_time, 0, np.zeros((1, num_machines), dtype=int), 0)
    process_time = np.insert(process_time, 0, np.zeros((1, num_jobs + 1), dtype=int), 1)

    setup_time = np.array([np.zeros((num_groups + 1, num_groups + 1), dtype=int)])
    for setup_time_position in position_dict['On machine']:
        data_row = setup_time_position[1] + 1
        setup_time_on_one_machine = np.array(data_csv.iloc[data_row: data_row + num_groups,
                                             0:num_groups].astype(int))
        setup_time_on_one_machine = np.insert(setup_time_on_one_machine, 0,
                                              setup_time_on_one_machine[np.diag_indices(num_groups)], 0)
        setup_time_on_one_machine = np.insert(setup_time_on_one_machine, 0,
                                              np.zeros((1, num_groups + 1), dtype=int), 1)
        setup_time = np.append(setup_time, [setup_time_on_one_machine], 0)

    data_row = position_dict['Time Window'][0][1] + 1
    TimeWindow = np.array(data_csv.iloc[data_row: data_row + num_groups, 0:2].astype(int))
    TimeWindow = np.insert(TimeWindow, 0, np.zeros((1, 2), dtype=int), 0)

    data_row = position_dict['Weight'][0][1] + 1
    Weight = np.array(data_csv.iloc[data_row: data_row + num_groups, 0:2].astype(int))
    Weight = np.insert(Weight, 0, np.zeros((1, 2), dtype=int), 0)

    group_array = np.arange(0, num_groups + 1)
    machine_array = np.arange(1, num_machines + 1)

    # Create a new model
    model = gp.Model("DBFGSP_TWET")
    model.setParam(GRB.Param.TimeLimit, 60)

    # Create decision variables
    c = {}
    for l in group_array[1:]:
        c[l] = model.addVars(jobs_in_each_group[l][1:], [l], machine_array, vtype=GRB.INTEGER, name="c")

    d = {}
    for l in group_array[1:]:
        d[l] = model.addVars(jobs_in_each_group[l][1:], [l], machine_array, vtype=GRB.INTEGER, name="d")

    # 每个组的完工时间
    c_fam = model.addVars(group_array[1:], vtype=GRB.INTEGER, name="c_fam")

    x = model.addVars(group_array, group_array, vtype=GRB.BINARY, name="x")

    y = {}
    for l in group_array[1:]:
        y[l] = model.addVars(jobs_in_each_group[l], jobs_in_each_group[l], [l], vtype=GRB.BINARY, name="y")

    Earliness = model.addVars(group_array[1:], vtype=GRB.INTEGER, name="Earliness")
    temp_Earliness = model.addVars(group_array[1:], vtype=GRB.INTEGER, name="temp_Earliness")

    Tardiness = model.addVars(group_array[1:], vtype=GRB.INTEGER, name="Tardiness")
    temp_Tardiness = model.addVars(group_array[1:], vtype=GRB.INTEGER, name="temp_Tardiness")

    TWET = model.addVar(vtype=GRB.INTEGER, name="TWET")

    model.setParam(GRB.Param.IntFeasTol, 1e-9)

    # Set objective
    # (1)
    model.setObjective(TWET, GRB.MINIMIZE)

    # Add constraints

    # (2)
    model.addConstrs(gp.quicksum(x[l, l1] for l1 in group_array if l1 != l) == 1
                     for l in group_array[1:])

    # (3)
    model.addConstrs(gp.quicksum(x[l, l1] for l in group_array if l != l1) == 1
                     for l1 in group_array[1:])

    # (4)
    model.addConstr(gp.quicksum(x[0, l1] for l1 in group_array[1:]) <= num_factories)

    # (5)
    model.addConstr(gp.quicksum(x[l, 0] for l in group_array[1:]) <= num_factories)

    # (6)
    model.addConstr(gp.quicksum(x[0, l1] for l1 in group_array[1:]) ==
                    gp.quicksum(x[l, 0] for l in group_array[1:]))

    # (7)
    model.addConstrs(gp.quicksum(y[l][j, j1, l] for j1 in jobs_in_each_group[l] if j1 != j) == 1
                     for l in group_array[1:]
                     for j in jobs_in_each_group[l])

    # (8)
    model.addConstrs(gp.quicksum(y[l][j, j1, l] for j in jobs_in_each_group[l] if j != j1) == 1
                     for l in group_array[1:]
                     for j1 in jobs_in_each_group[l])

    # (9)
    model.addConstrs(d[l][j, l, i] >= c[l][j, l, i]
                     for l in group_array[1:]
                     for j in jobs_in_each_group[l][1:]
                     for i in machine_array)

    # (10)
    model.addConstrs(c[l][j1, l, i] >= d[l][j, l, i] + process_time[j1][i] + (y[l][j, j1, l] - 1) * h
                     for l in group_array[1:]
                     for j in jobs_in_each_group[l][1:]
                     for j1 in jobs_in_each_group[l][1:]
                     if j1 != j
                     for i in machine_array)

    # (11)
    model.addConstrs(
        c[l1][j1, l1, i] >= d[l][j, l, i] + setup_time[i, l, l1] + process_time[j1][i] + (x[l, l1] - 1) * h
        for l in group_array[1:]
        for l1 in group_array[1:]
        if l != l1
        for j in jobs_in_each_group[l][1:]
        for j1 in jobs_in_each_group[l1][1:]
        for i in machine_array)

    # (12)
    model.addConstrs(c[l][j, l, i] >= setup_time[i, 0, l] + process_time[j][i] + (x[0, l] - 1) * h
                     for l in group_array[1:]
                     for j in jobs_in_each_group[l][1:]
                     for i in machine_array)

    # (13)
    model.addConstrs(c[l][j, l, i + 1] == d[l][j, l, i] + process_time[j][i + 1]
                     for l in group_array[1:]
                     for j in jobs_in_each_group[l][1:]
                     for i in machine_array[:-1])

    # (14)
    # model.addConstrs(c_fam[l] >= d[l][j, l, num_machines]
    #                  for l in group_array[1:]
    #                  for j in jobs_in_each_group[l][1:])
    for l in group_array[1:]:
        model.addGenConstrMax(c_fam[l], [d[l][j, l, num_machines]
                                         for j in jobs_in_each_group[l][1:]], 0)

    # (15)
    # for l in group_array[1:]:
    #     model.addConstr(temp_Earliness[l] >= TimeWindow[l][0] - c_fam[l])

    for l in group_array[1:]:
        model.addConstr(Earliness[l] >= TimeWindow[l][0] - c_fam[l])

    # for l in group_array[1:]:
    #     model.addGenConstrMax(Earliness[l], [temp_Earliness[l], 0], 0)

    # (16)
    # for l in group_array[1:]:
    #     model.addConstr(temp_Tardiness[l] >= c_fam[l] - TimeWindow[l][1])

    for l in group_array[1:]:
        model.addConstr(Tardiness[l] >= c_fam[l] - TimeWindow[l][1])

    # for l in group_array[1:]:
    #     model.addGenConstrMax(Tardiness[l], [temp_Tardiness[l], 0], 0)

    # (17)
    model.addConstr(TWET == gp.quicksum((Weight[l][0] * Earliness[l] + Weight[l][1] * Tardiness[l])
                                        for l in group_array[1:]))

    # Optimize model
    model.optimize()
    #
    # model.computeIIS()
    # model.write("model.ilp")


    for v in model.getVars():
        print('%s %g' % (v.Varname, v.X))

    print('Obj: %g' % model.ObjVal)

except gp.GurobiError as e:
    print('Error code ' + str(e.errno) + ': ' + str(e))
except AttributeError as e:
    print('Encountered an attribute error')
