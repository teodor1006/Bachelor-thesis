import re
import random
from collections import defaultdict
import numpy as np
import simple_colors
import networkx as nx
import matplotlib.pyplot as plt


def get_data(str):
    pattern = re.compile("\d+")
    result = re.findall(pattern,str)
    return result

def get_qubit(str):
    pattern = re.compile("q\[(\d+)")
    result = re.findall(pattern,str)
    return result


def converter_circuit_from_qasm(input_file_name):
    gate_list = []
    num_qubits = 0

    with open(input_file_name,'r') as qasm_file:
        for line in (l.strip() for l in qasm_file):
            if line[0:4] == 'creg':
                num_qubits = int(get_data(line)[0])
            if line[0:1] == 'x':
                x = get_qubit(line)
                x_target = int(x[0])
                listSingle = [x_target]
                gate_list.append(listSingle)
            if line[0:1] == 'h':
                h = get_qubit(line)
                h_target = int(h[0])
                listSingle = [h_target]
                gate_list.append(listSingle)
            if line[0:2] == 'rz':
                rz = get_qubit(line)
                rz_target = int(rz[0])
                listSingle = [rz_target]
                gate_list.append(listSingle)
            if line[0:2] == 'ry':
                ry = get_qubit(line)
                ry_target = int(ry[0])
                listSingle = [ry_target]
                gate_list.append(listSingle)
            if line[0:2] == 'rx':
                rx = get_qubit(line)
                rx_target = int(rx[0])
                listSingle = [rx_target]
                gate_list.append(listSingle)
            if line[0:2] == 'cx':
                cnot = get_qubit(line)
                cnot_control = int(cnot[0])
                cnot_target = int(cnot[1])
                listSingle = [cnot_control,cnot_target]
                gate_list.append(listSingle)

    return gate_list, num_qubits



def list_str_to_int(gate_list):
    return [list(map(int, gate)) for gate in gate_list]


def count_global_gates(gate_list,part_qbits):
    global_gates = 0
    for gate in gate_list:
        gate_qbits = set(gate)
        is_global = True
        for qubits in part_qbits:
            if gate_qbits.issubset(qubits):
                is_global = False
                break
        if is_global:
            global_gates += 1
    return global_gates

def find_global_gates(gate_list,part_qbits):
    global_gates = []
    for gate in gate_list:
        gate_qbits = set(gate)
        is_global = True
        for qubits in part_qbits:
            if gate_qbits.issubset(qubits):
                is_global = False
                break
        if is_global:
            global_gates.append(gate)
    return global_gates



def find_local_gates(gate_list,part_qbits):
    local_gates = []
    for gate in gate_list:
        gate_qbits = set(gate)
        is_local = False
        for qubits in part_qbits:
            if gate_qbits.issubset(qubits):
                is_local = True
                break
        if is_local:
            local_gates.append(gate)
    return local_gates


def count_local_gates(gate_list, part_qbits):
    local_gates = 0
    for gate in gate_list:
        gate_qbits = set(gate)
        is_local = False
        for qubits in part_qbits:
            if gate_qbits.issubset(qubits):
                is_local = True
                break
        if is_local:
            local_gates += 1
    return local_gates


def kernighan_lin_partition(gate_list):

    G = nx.Graph()
    # Count the number of times each edge appears in the gate_list
    edge_weights = defaultdict(int)
    for gate in gate_list:
        for i in range(1,len(gate)):
            edge = tuple(sorted([gate[i - 1],gate[i]]))
            if edge in edge_weights:
                edge_weights[edge] += 1
            else:
                edge_weights[edge] = 1

    # Add edges to the graph with their corresponding counts
    for edge,count in edge_weights.items():
        G.add_edge(edge[0],edge[1],weight=count)

    partition = nx.community.kernighan_lin_bisection(G,seed=0)

    return partition


def chromosome_coding(gate_list, part_qbits):
    global_gates = count_global_gates(gate_list, part_qbits)
    code = [(random.randint(0, 1)) for _ in range(global_gates)]
    return code



def initialize_population(gate_list,part_qbits):
    population_list = []
    for i in range(population_size):
        population_list.append(chromosome_coding(gate_list,part_qbits))
    return population_list


def non_execute(gate1, gate2):
    local_gates = find_local_gates(gate_list,part_qbits)
    global_gates = find_global_gates(gate_list,part_qbits)
    chromosome = chromosome_coding(gate_list,part_qbits)

    # 1. Check if gate2 is a local gate whose one of the qubits is the same as the migrated qubit of gate1
    if (gate2 in local_gates) and (gate1 in global_gates):
        if len(gate2) == 2:
            if (gate2[0] == gate1[0]) or (gate2[1] == gate1[1]):
                return True
        elif len(gate2) == 1:
            if (gate2[0] == gate1[0]):
                return True

    # 2. Different labels -> different partition execution
    elif gate1 in global_gates and gate2 in global_gates:
        i = global_gates.index(gate1)
        j = global_gates.index(gate2)
        if chromosome[i] != chromosome[j]:
            return True

    # 3. No common qubits and same partition (label) execution -> requires another teleportation
    elif gate1 in global_gates and gate2 in global_gates:
        i = global_gates.index(gate1)
        j = global_gates.index(gate2)
        if chromosome[i] == chromosome[j] and ((gate1[0] != gate2[0]) and (gate1[1] != gate2[1])):
            return True

    return False


def non_commute(gate1, gate2):

    if len(gate1) == 2 and len(gate2) == 2:
        if gate1[0] == gate2[1] or gate2[0] == gate1[1]:
            return True

    elif len(gate1) == 2 and len(gate2) == 1:
        if gate2[0] == gate1[0] or gate2[0] == gate1[1]:
            return True

    elif len(gate2) == 2 and len(gate1) == 1:
        if gate1[0] == gate2[0] or gate1[0] == gate2[1]:
            return True

    return False


def min_teleport(gate_list, chromosome, nt):

    global_gates = find_global_gates(gate_list, part_qbits)
    local_gates = find_local_gates(gate_list, part_qbits)

    # Base condition
    if not gate_list:
        return nt

    if gate_list[0] in local_gates:
        gate_list.pop(0)
        return min_teleport(gate_list, chromosome, nt)

    temp = gate_list[0]
    nt += 1
    gate_list.remove(gate_list[0])
    global_gates.remove(global_gates[0])

    for i in range(len(gate_list)):
        found = False
        if i < len(gate_list)  and not non_execute(temp, gate_list[i]):
            for k in range(i,-1,-1):
                if non_commute(gate_list[i], gate_list[k]):
                    found = True
                    break
            if found:
                gate_list.remove(gate_list[i])
                if i < len(global_gates):
                    global_gates.remove(global_gates[i])

    nt += 1

    return min_teleport(gate_list, chromosome, nt)



def population_fitness(gate_list, chromosomes,nt):
    fitness_list = []
    for chromosome in chromosomes:
        gate_list_copy = gate_list.copy()
        temp = min_teleport(gate_list_copy, chromosome, nt)
        fitness_list.append(temp)
        #print(f"Chromosome: {chromosome}, Fitness: {temp}")
    #print(f"Fitness List: {fitness_list}")
    return fitness_list


def graph_qubits(gate_list):

    G = nx.Graph()
    # Count the number of times each edge appears in the gate_list
    edge_weights = defaultdict(int)
    for gate in gate_list:
        for i in range(1,len(gate)):
            edge = tuple(sorted([gate[i - 1],gate[i]]))
            if edge in edge_weights:
                edge_weights[edge] += 1
            else:
                edge_weights[edge] = 1

    # Add edges to the graph with their corresponding counts
    for edge,count in edge_weights.items():
        G.add_edge(edge[0],edge[1],weight=count)

    pos = nx.spring_layout(G)

    # Draw graph with edge weights
    nx.draw(G,pos,with_labels=True,font_weight='bold')
    edge_labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)

    plt.show()

def graph_after_partitioning(gate_list):

    G = nx.Graph()
    # Count the number of times each edge appears in the gate_list
    edge_weights = defaultdict(int)
    for gate in gate_list:
        for i in range(1,len(gate)):
            edge = tuple(sorted([gate[i - 1],gate[i]]))
            if edge in edge_weights:
                edge_weights[edge] += 1
            else:
                edge_weights[edge] = 1

    # Add edges to the graph with their corresponding counts
    for edge,count in edge_weights.items():
        G.add_edge(edge[0],edge[1],weight=count)

    partitions = kernighan_lin_partition(gate_list)
    colors = []
    for node in G.nodes():
        if node in partitions[0]:
            colors.append("red")
        elif node in partitions[1]:
            colors.append("green")
    color = colors
    pos = nx.spring_layout(G)
    nx.draw(G,pos,node_color=color,with_labels=True)
    edge_labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels,font_weight='bold')
    nx.draw_networkx_labels(G,pos)
    plt.show()


def roulette_wheel_selection(fitness_list, chromosome_list, elite_count):
    selected_individuals = []

    # Select the top N individuals as elites and add them to the selected individuals list
    elite_indices = sorted(range(len(fitness_list)),key=lambda k:fitness_list[k])[:elite_count]
    for i in elite_indices:
        selected_individuals.append(chromosome_list[i])

    # Calculate the fitness sum of the remaining individuals
    fitness_sum = sum(fitness_list) - sum([fitness_list[i] for i in elite_indices])

    # Select the remaining individuals using roulette wheel selection
    while len(selected_individuals) < len(chromosome_list):
        r = random.uniform(0,fitness_sum)
        partial_sum = 0
        for i in range(len(fitness_list)):
            if i not in elite_indices:
                partial_sum += fitness_list[i]
                if partial_sum > r:
                    selected_individuals.append(chromosome_list[i])
                    break

    return selected_individuals


def two_point_crossover(selected_population_list, population_size, elite_count):
    elite = selected_population_list[:elite_count]
    non_elite = selected_population_list[elite_count:]
    new_population = []

    offspring_count = population_size - len(elite)
    # Generate offspring
    for i in range(offspring_count):
        father = random.choice(non_elite)
        mother = random.choice(non_elite)
        child = father
        if np.random.rand() < Crossover_rate:
            # Set crossover points
            point1 = random.randint(0, len(father) - 1)
            point2 = random.randint(point1, len(father))
            child = father[:point1] + mother[point1:point2] + father[point2:]
        child = bitflip(child)
        new_population.append(child)
    population = elite + new_population
    return population


def bitflip(child):
   for i in range(len(child)):
       if np.random.rand() < Mutation_rate:
           child[i] = 1 - child[i]
   return child



if __name__ == '__main__':

    input_filename = 'C:/Users/Teodor/QC_GA_Partition/qasm/parity_247.qasm'
    gate_list = converter_circuit_from_qasm(input_filename)[0]
    gate_list = list_str_to_int(gate_list)
    print(gate_list)

    print('Number of Gates: ' + simple_colors.green((len(gate_list))))

    num_qubits = converter_circuit_from_qasm(input_filename)[1]
    num_qubits = int(num_qubits)
    print('Number of Qubits: ' + simple_colors.green(num_qubits))

    part_qbits = kernighan_lin_partition(gate_list)
    print('Partitioned qubits: ' + simple_colors.green(part_qbits))

    global_gates_num = count_global_gates(gate_list,part_qbits)
    print('Number of Global Gates: ' + simple_colors.red(global_gates_num))
    global_gates_list = find_global_gates(gate_list,part_qbits)
    print('List of Global Gates: ' + simple_colors.red(str(global_gates_list)))

    local_gates_num = count_local_gates(gate_list,part_qbits)
    print('Number of Local Gates: ' + simple_colors.blue(local_gates_num))
    local_gates_list = find_local_gates(gate_list,part_qbits)
    print('List of Local Gates: ' + simple_colors.blue(str(local_gates_list)))

    Crossover_rate = 0.9
    Mutation_rate = 0.1
    N_generations = 100
    population_size = 30
    elite_count = 2

    #init_graph = graph_qubits(gate_list)
    #graph_after = graph_after_partitioning(gate_list)

    nt = 0
    best_fitness = []

    population_list = initialize_population(gate_list,part_qbits)
    fitness_list = population_fitness(gate_list,population_list,nt)


    for _ in range(N_generations):
        selected_population_list = roulette_wheel_selection(fitness_list,population_list,elite_count)
        population_list = two_point_crossover(selected_population_list,population_size,elite_count)
        fitness_list = population_fitness(gate_list,population_list,nt)
        best_fitness.append(np.min(fitness_list))

    print('Best Teleportation Cost: ',end='')
    print(simple_colors.green(np.min(best_fitness)))
