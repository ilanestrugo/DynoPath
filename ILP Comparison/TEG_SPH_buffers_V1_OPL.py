"""
    Tal's take on the TEG shortest-path heuristic with funnel outputs
    Creation Data: April 1, 2023,  (a version with funnels)
    Update April 4, 2023 - add funnels and script mode (v2)
    Update April 6, 2023 - add belts in script mode (v3)
    Update April 7, 2023 - disallow south movements, add block output for statistics (v4)
    Update April 10, 2023 - allow sparse inputs cells, block and std calculation bug fix, internal documentation (v5)
    Update June 9, 2023 - conversion to version w/o funnels and with general grid structure


    TEG is stored in a T x Lx x Ly x 5 array  (where LX and Ly are the dimensions of the blocking rectangle)
    the five elements of the last axis are
    0 - stay arc
    1 - up (north) arc
    2 - right (east) arc
    3 - down (south) arc
    4 - left (west) arc

    There are four types of cells in the grid
    1. Regular cells from which the items can move in all four cardinal directions or stay (except at the edges of the grid)
    2. Input cells like regular cells but to which the items can not move in
    3. Output cells where the items are vanished from the system.
       From these cells the item cannot move out to grid cells and cannot stay
    4. Inactive cells. Items cannot move in or out of these cells.
       If the input and output cells are out of the main rectangle, we have inactive cells between them

"""

import numpy as np
#import pandas as pd
import time
import sys
import pickle
import random
import itertools
import subprocess
import os

res_file_name = sys.argv[0][:-3] + "_res.csv"
block_len = 100
option_argv_position = 6  # position of the option list in argv

if len(sys.argv) < option_argv_position:
    print(f"Usage: python {sys.argv[0]} grid_file_name T1 warmup [OPTIONS]")
    print("The grid file is constructed as follow")
    print("   Line 1: Lx,Ly the dimensions of the blocking rectangle of the grid including input and output cells")
    print("   Line 2: List of input cells (e.g., (1,1),(2,1),(3,1)")
    print("   Line 3: List of output cells")
    print("   Line 4: List of inactive cells (can be empty)")
    print("T1 -  number of time steps including warmup (but not including cool down period)")
    print("warmup - number of time steps in the warmup period")
    print("OPTIONS is an optional string parameter that can contain any subset of SHBD (not case sensitive)")
    print(f"  'S' - a script data structure is created and saved as a pickle file (for the animation)")
    print(f"  'H' - an header row is written to the result file {res_file_name}")
    print(f"  'B' - set block length to the value of the next parameter (default {block_len})")
    print("   'D' - save detailed in/out data, one for each time step in a pickle file")
    exit(1)

# Parameters from the command line
grid_file_name, T1, warmup, seed, time_limit = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
T = T1 + warmup  # full graph including cooldown time

do_script = False
data_save = False
sort_option = 'Stay' #options for sorting of options: Stay,Move,Random,Location

if len(sys.argv) >= option_argv_position + 1:
    if "H" in sys.argv[option_argv_position].upper():
        f = open(res_file_name, "a")
        f.write(
            "Grid File Name, Mode, Lx, Ly, 4WC, Inputs, Outputs, OutputGroups, T1, warmup, seed, count, moves / item, left before cd, exit rate, std, entered after wu, entrance Rate, std, W,L, cpu time, max max buffer, makespan, flow time, flow count, OPL flow time, OPL obj, OPL lb\n")
        f.close()

    if "S" in sys.argv[option_argv_position].upper():
        do_script = True
        script = np.empty(T, dtype=object)
        for t in range(T):
            script[t] = []
        item_destination = {}

    if "B" in sys.argv[option_argv_position].upper():
        if len(sys.argv) >= option_argv_position + 2:
            block_len = int(sys.argv[option_argv_position + 1])

    if T // block_len != T / block_len or warmup // block_len != warmup / block_len:
        print(f"Panic: block len ({block_len}) must be a divider of T ({T1}) and warmup time ({warmup})")
        exit(1)

    if "D" in sys.argv[option_argv_position].upper():
        data_save = True

    #options for sorting of options: Stay,Move,Random
    if "M" in sys.argv[option_argv_position].upper():
        sort_option = 'Move'

# table for translating the move ids to coordinates and vice versa
delta = [(0, 0), (0, 1), (1, 0), (0, -1), (-1, 0)]  # the five directions
delta_inv = {(0, 0): 0, (0, 1): 1, (1, 0): 2, (0, -1): 3, (-1, 0): 4}

# get layout from file
with open(grid_file_name, "r") as file:
    newline_break = ""
    for readline in file:
        line_strip = readline.strip()
        newline_break += line_strip + "|"

arr = newline_break.split('|')
L = arr[0].split(',') # lx,Ly
Is = arr[1][1:-1].split('),(') # list of inputs
Os = arr[2][1:-1].split('),(') # list of final outputs
IAs = []
if ")" in arr[3][1:-1]:
    IAs = arr[3][1:-1].split('),(') # list if inactive cells

Lx = int(L[0])
Ly = int(L[1])

I=[]
O=[]
IA=[]
OT=[] #list of output with types
Z=[] #list of optional z

I0_tuple = '{'
for i in range(len(Is)):
    x=int(Is[i].split(',')[0])
    y=int(Is[i].split(',')[1])
    I.append((x, y))
    node_str = '<' + str(x) +' '+str(y)+'>'
    I0_tuple = I0_tuple + node_str
I0_tuple = I0_tuple +'}'

OT_tuple = '{'
for i in range(len(Os)):
    x=int(Os[i].split(',')[0])
    y=int(Os[i].split(',')[1])
    z=int(Os[i].split(',')[2])
    O.append((x, y))
    OT.append((x, y, z)) #04-07 add Z
    if (z not in Z):
        Z.append(z)

    node_str = '<' + str(x) +' '+str(y)+ ' '+str(z) +'>'
    OT_tuple = OT_tuple + node_str
OT_tuple = OT_tuple +'}'


for i in range(len(IAs)):
    x=int(IAs[i].split(',')[0])
    y=int(IAs[i].split(',')[1])
    IA.append((x, y))

N_out = np.empty((Lx, Ly, 5), dtype=object)  # array of neighbors
N_in = np.empty((Lx, Ly, 5), dtype=object)  # array of cells that I am a neighbor of.

""" Create set of neighbors for the 4WC  """
C_tuple = '{'
for x in range(Lx):
    for y in range(Ly):
        if (x, y) not in IA:
            if (x, y) not in O:
                N_out[x, y, 0] = (x, y)  # stay
                N_in[x, y, 0] = (x, y)
            if y < Ly - 1 and (x, y+1) not in I and (x, y+1) not in IA:
                N_out[x, y, 1] = (x, y + 1)  # north
                N_in[x, y + 1, 1] = (x, y)
            if x < Lx - 1 and (x+1, y) not in I and (x+1, y) not in IA:
                N_out[x, y, 2] = (x + 1, y)  # East
                N_in[x + 1, y, 2] = (x, y)
            if y > 0 and (x, y-1) not in I and (x, y-1) not in IA:
                N_out[x, y, 3] = (x, y - 1)  # South
                N_in[x, y - 1, 3] = (x, y)
            if x > 0 and (x-1, y) not in I and (x-1, y) not in IA:
                N_out[x, y, 4] = (x - 1, y)  # West
                N_in[x - 1, y, 4] = (x, y)
        #build aldo C tuple for OPL
        if (x, y) not in O and (x, y) not in IA:
            node_str = '<' + str(x) +' '+str(y)+'>'
            C_tuple = C_tuple + node_str
C_tuple = C_tuple + '}'

print(f"Layout {Lx} x {Ly} with {len(I)} inputs,  {len(O)} outputs cells with {len(Z)} conveyors")

start = time.time()
# construct and initialized TEG
TEG = np.ones((T, Lx, Ly, 5), dtype=np.int8)
for x in range(Lx):
    for y in range(Ly):
        if (x, y) in IA or (x, y) in O:
            TEG[:, x, y, :] = 0
        else:
            for d in range(1, 5):
                if N_out[x, y, d] is None:
                    TEG[:, x, y, d] = 0

from_node = np.empty((T, Lx, Ly), dtype=object)
print(f"Graph building time {(time.time() - start):.3f} sec.   ({T} layers)")

"""  Find shortest path from input_num to out_out num on the TEG starting at time t
    Assumes TEG, T, Lx, Ly, from_node are global variables """


def find_path(input_num, output_num, t, sort_option, verbal=False):
    C = np.full((T, Lx, Ly), np.iinfo(np.int16).max, dtype=np.int16)
    curr_t = t
    x, y = I[input_num][0], I[input_num][1]
    C[t, x, y] = 0
    from_node[t, x, y] = (-1, -1)
    open_nodes = [I[input_num]]
    not_found_yet = True
    x_out = 0
    y_out = 0
    #find line from I to O
    a = 1
    b = 0
    #18-08-23 several move options
    d_range = [0,1,2,3,4] #default order stay,up,right,down.keft
    d1_range = [1,2,4,0,3] #set order up,right,left,stay,down
    d2_range = [2,4,1,0,3] #set order right,left,stay,up,down
    d_range_up_Side = [1,2] #array with 2 direction only
    if (sort_option == 'Move'):
        d_range = d1_range #set order up,right,left,stay,down

    #till here
    while not_found_yet:
        new_nodes = []
        for (x, y) in open_nodes:
            for d in d_range:
                if TEG[curr_t, x, y, d]:
                    # print(curr_t,d, x,y)
                    new_x, new_y = N_out[x, y, d]
                    # print(new_x, new_y)

                    #V1 - add 1 instead of (d > 0)
                    if C[curr_t + 1, new_x, new_y] > C[curr_t, x, y] + (d > 0):
                        C[curr_t + 1, new_x, new_y] = C[curr_t, x, y] + (d > 0)
                    #if C[curr_t + 1, new_x, new_y] > C[curr_t, x, y] + 1: #(d > 0):
                    #    C[curr_t + 1, new_x, new_y] = C[curr_t, x, y] + 1 #(d > 0)
                        new_nodes.append((new_x, new_y))
                        from_node[curr_t + 1, new_x, new_y] = (x, y)
                        if (new_x, new_y, output_num) in OT: #04-07 changed
                            not_found_yet = False
                            x_out = new_x
                            y_out = new_y
        open_nodes = new_nodes
        curr_t += 1

    #num_of_moves = C[curr_t, O[output_num][0], O[output_num][1]]
    num_of_moves = C[curr_t, x_out, y_out]
    if verbal:
        print(f"Found path of {curr_t - t} time steps and {num_of_moves} moves")
    # recover path

    #x, y = O[output_num][0], O[output_num][1]
    x, y = x_out, y_out
    P = [(x, y)]
    while True:
        x, y = from_node[curr_t, x, y]
        if x == -1:
            break
        P.insert(0, (x, y))
        curr_t -= 1

    if verbal:
        print(f"The path is {P}")

    return P, num_of_moves


""" Update the TEG to reflect a path P starting at time t
    assume P is a legtimate path on TEG
"""


def update_TEG(P, t):
    for i in range(len(P) - 1):
        x, y, next_x, next_y = P[i][0], P[i][1], P[i + 1][0], P[i + 1][1]

        # delete all arcs out of (x,y) at time t+i
        TEG[t + i, x, y, :] = 0

        # delete all arcs into (next_x, next_y) at time t+i
        for dd in range(5):  # scan all neighboring locations
            if N_in[next_x, next_y, dd]:
                neighbor_x, neighbor_y = N_in[next_x, next_y, dd]  # find the neighbor location
                TEG[t + i, neighbor_x, neighbor_y, dd] = 0  # delete  arc

        # delete arcs into (x,y) at time t+i except block movement
        if i > 0:  # no need to delete arcs at the first movement
            #replced 26-03-24
            d = delta_inv[(next_x - x, next_y - y)]  # direction of the current movement
            #if (next_x, next_y) == P[-1]:  # ejection step
            #    d = 1  # north
            #else:
            #    d = delta_inv[(next_x - x, next_y - y)]  # direction of the current movement
            # Eliminate all entrances to the current location  at the current interval
            for dd in range(1, 5):  # scan all neighboring locations (not including current)
                if dd != d and N_in[x, y, dd]:  # we disallow movements from all neighbors expect from the one that are moving at the same direction
                    neighbor_x, neighbor_y = N_in[x, y, dd]  # find the neighbor location
                    TEG[
                        t + i, neighbor_x, neighbor_y, dd] = 0  # delete the arc from the neighbor location to the current
        # delete arcs out of (next_x,next_y) at time t+i except block movement
        if i < len(P) - 2:  # not last move
            d = delta_inv[(next_x - x, next_y - y)]  # direction of the current movement
            for dd in range(1, 5):  # scan all neighboring locations (not including current)
                if dd != d:
                    TEG[t + i, next_x, next_y, dd] = 0  # delete the arc from the destination location to its neighbors

np.random.seed(seed)
count = 0
T1 = T - warmup
items_enter = np.zeros(T1)
items_leave = np.zeros((T, len(Z)))  # detailed by time step and exit buffrt
total_time_after_warmup = 0
total_moves_after_warmup = 0
makespan = 0
flowTime = 0
flowCount = 0
for t in range(T1):
    for i in range(len(I)):
        if TEG[t, I[i][0], I[i][1], 0]:  # Input i is not occupied
            dest = np.random.randint(0, len(Z)) #04-07 get type of ouput from Z
            P, m = find_path(i, dest, t, sort_option)
            update_TEG(P, t)

            # script structure
            # For each time step we have a list with one item for each movement
            # Each movement is described by a 4 tuple (id, orig, dest, last_move flag)
            if do_script:
                for j in range(len(P) - 1):
                    script[t + j].append((count, P[j], P[j + 1], False))
                item_destination[count] = dest
                x, y = P[-1]
                tt = t + len(P) - 1

            count += 1
            items_enter[t] += 1
            items_leave[t + len(P) - 1, dest] += 1
            if t >= warmup:
                total_time_after_warmup += (len(P) - 1)
                total_moves_after_warmup += m

            if makespan < t + len(P) - 1:
                makespan = t + len(P) - 1

            #get flowTime for items exiting after T1-1
            if (t + len(P) - 1 > T1-1):
                leftFlowTime = t + len(P) - 1 - (T1-1)
                flowTime = flowTime + leftFlowTime
                flowCount = flowCount + 1

            # if do_blocks:
            #     block_moves[t//block_len] += m
            #     block_exit[(t+len(P)-1)//block_len] += 1
            #     block_enter[t//block_len] += 1
            #

    #get last entry point - build status
    I_tuple = '{'
    if (t == T1-1):
        print(f"count: {count:.3f}")
        #f = open('Status_for_OPL.csv', "a")
        makespan = makespan - t
        first=0
        scrip_t = script[t]
        for s in range(len(scrip_t)):
            node = scrip_t[s]
            i=node[0]
            xy=node[1]
            x=xy[0]
            y=xy[1]
            dest = item_destination[i]
            node_str = '<' + str(x) +' '+str(y)+ ' '+str(dest) +'>'
            I_tuple = I_tuple + node_str
            #print (node_str)
            #if (first==0):
                #f.write(f"{node_str}")
                #first=1
            #else:
                #f.write(f",{node_str}")
            #f.write("\n")
        #f.write(f",t={t}")
        #f.write(f",Makespan={makespan}")
        #f.write(f",FlowTime={flowTime}")
        #f.write(f",FlowCount={flowCount}")
        #f.write("\n")
        #f.close()
        I_tuple = I_tuple + '}'

cpu_time = (time.time() - start)
print(flowTime)

#run cplex to caluclate same problem - set T as makespan
#build C
file_export = ''
f = open("Sorting_Data_buffers.dat","w")
f.write('file_export = "%s";\n'%file_export)
f.write('time_limit = %d;\n'% time_limit)
f.write('alpha=0.000000;') #parameter for lower weight for arcs from Input points
f.write('beta=1.000000;') #parameter for weight of all arcs
f.write('gamma=0.010000;') #parameter for added weight of movement arcs
f.write('T=%d;\n'%makespan)
f.write('type_num=%d;\n'%len(Z))
f.write('C=%s;\n'%C_tuple)
f.write('I0=%s;\n'%I0_tuple)
f.write('OT=%s;\n'%OT_tuple)
f.write('I=%s;\n'%I_tuple)
f.close()

# Delete the file to make sure that we read the new one
if os.path.exists("sorting_out.txt"):
    os.remove("sorting_out.txt")

try:
    subprocess.run(["oplrun", "Sorting_Model_buffers.mod", "Sorting_Data_buffers.dat"], check=True)

    if os.path.exists("sorting_out.txt"):
        with open("sorting_out.txt", "r") as file:
            newline_break = ""
            for readline in file:
                line_strip = readline.strip()
                newline_break += line_strip + ","
            print(newline_break)

    arr = newline_break.split(',')
    OplFlowTime = int(arr[0].replace('FlowTime=',''))
    OplObj = float(arr[2].replace('obj=',''))
    OplLb = float(arr[3].replace('lb=',''))
except:
    print("Could not solve the model")
    OplFlowTime = ''
    OplObj = ''
    OplLb = ''

enter_after_wu = np.sum(items_enter[warmup:T1])
leave_before_cd = np.sum(items_leave[warmup:T1])

T2 = T1 - warmup  # periods after warmup and before cooldown
L = total_time_after_warmup / T2
W = total_time_after_warmup / enter_after_wu

# exit_rate_blocks = [np.average(items_leave[i*block_len:(i+1)*block_len])/block_len for i in range(warmup// block_len,T1//block_len)]
exit_rate_blocks = [np.average(items_leave[i * block_len:(i + 1) * block_len]) for i in range(T // block_len)]
entrance_rate_blocks = [np.average(items_enter[i * block_len:(i + 1) * block_len]) for i in
                        range((T - warmup) // block_len)] + [0] * (warmup // block_len)

std_exit_rate = np.std(exit_rate_blocks[int(warmup / block_len):int((T - warmup) / block_len)])
std_entrance_rate = np.std(entrance_rate_blocks[int(warmup / block_len):int((T - warmup) / block_len)])

#Grid File Name, Lx, Ly, 4WC, Inputs, Outputs, T1, warmup, count, moves / item, left before cd, exit rate, std, entered after wu, entrance Rate, std, W,L, cpu time\n")
buffer_level = np.zeros(len(Z))
max_buffer_level = np.zeros(len(Z))
for t in range(T):
    buffer_level += items_leave[t]
    buffer_level = np.maximum(0, buffer_level - 1)   # reduce 1 assuming the conveyor rate is 1 (i.e., the same speed as the grid)
    max_buffer_level = np.maximum(buffer_level, max_buffer_level)

sort_option = sort_option + 'V1'
f = open(res_file_name, "a")
f.write(f"{grid_file_name}, {sort_option}, {Lx}, {Ly}, {Lx * Ly}, {len(I)}, {len(O)}, {len(Z)}, {T}, {warmup}, {seed}, {count}, "
        f"{total_moves_after_warmup / enter_after_wu:.3f}, {leave_before_cd}, "
        f"{leave_before_cd / T2:.3f}, {std_exit_rate:.4f},{enter_after_wu},"
        f"{enter_after_wu / T2:.3f}, {std_entrance_rate:.4f},{W:.3f},{L:.3f}, {cpu_time:.2f}")
f.write(f",{np.max(max_buffer_level)}, {makespan},{flowTime},{flowCount},{OplFlowTime},{OplObj},{OplLb}")
#for x in max_buffer_level:
#    f.write(f",{x}")
f.write("\n")
f.close()

print(f"{Lx}x{Ly} with {len(O)} ")
print(f"Throughput: {count / T1:.3f}")
print(f"Throughput after warmup based on items leaving the system: {leave_before_cd / T2:.3f} ({std_exit_rate:.4f})")
print(
    f"Throughput after warmup based on items entering the system: {enter_after_wu / T2:.3f} ({std_entrance_rate:.4f})")
print(f"Mean time in system (W) after warmup: {W:.3f}")
print(f"Mean WIP (L) after warmup: {L:.3f}")
print(f"Entrance blocks {entrance_rate_blocks}")
print(f"Entrance blocks {exit_rate_blocks}")
print(f"Total time {(time.time() - start):.3f} sec.  ({count / (time.time() - start):.1f} items per seconds)")
print(f"Buffer max levels {max_buffer_level}")

if do_script:
    f = open(f"script_{Lx}x{Ly}_Z{len(Z)}_T{T}_{sort_option}_{seed}.p", "wb")
    pickle.dump((Lx, Ly, I, OT, IA, count, script, item_destination), f)
    f.close()
    print("script saved")

if data_save:
    f = open(f"in_out_{Lx}x{Ly}_Z{len(Z)}_T{T}_{sort_option}_{seed}.p", "wb")
    pickle.dump((Lx, Ly, I, OT, IA, items_leave, items_enter), f)
    f.close()
    print("in/out data saved")
