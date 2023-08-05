import re
import numpy as np
import simple_colors
import random
import math
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



def count_global_gates(gate_list,chromosome):
    global_gates = 0
    for gate in gate_list:
        gate_qbits = set(gate)
        is_global = False  # initialize to False
        partitions = set(chromosome[q] for q in gate_qbits)
        if len(partitions) > 1:  # check if there are qubits in different partitions
            is_global = True
        if is_global:
            global_gates += 1

    return global_gates


def count_local_gates(gate_list, chromosome):
    local_gates = 0
    for gate in gate_list:
        gate_qbits = set(gate)
        is_local = True
        partitions = set(chromosome[q] for q in gate_qbits)
        if len(partitions) > 1:  # check if there are qubits in different partitions
            is_local = False
        if is_local:
            local_gates += 1
    return local_gates


def find_global_gates(gate_list, chromosome):
    global_gates = []
    for gate in gate_list:
        gate_qbits = set(gate)
        is_global = False
        partitions = set(chromosome[q] for q in gate_qbits)
        if len(partitions) > 1:  # check if there are qubits in different partitions
            is_global = True
        if is_global:
            global_gates.append(gate)
    return global_gates


def find_local_gates(gate_list,chromosome):
    local_gates = []
    for gate in gate_list:
        gate_qbits = set(gate)
        is_local = True  # initialize to False
        partitions = set(chromosome[q] for q in gate_qbits)
        if len(partitions) > 1:  # check if there are qubits in different partitions
            is_local = False
        if is_local:
            local_gates.append(gate)
    return local_gates


def chromosome_coding(num_qubits,num_partitions):

    qubits_per_partition = int(math.ceil(num_qubits / num_partitions))
    # Shuffle the qubits to ensure each qubit appears only once
    qubits = list(range(num_qubits))
    random.shuffle(qubits)
    # Assign each qubit to a partition in the shuffled order
    partition_indices = list(range(num_partitions)) * (qubits_per_partition // num_partitions + 1)
    encoding = [None] * num_qubits
    for qubit in qubits:
        partition = partition_indices.pop(0)
        encoding[qubit] = partition + 1
        partition_indices.append(partition)
    return encoding


def initialize_population(num_qubits,num_partitions):
    population_list = []
    for i in range(population_size):
        population_list.append(chromosome_coding(num_qubits,num_partitions))
    return population_list


def check_partition_crossing(gate1, gate2, chromosome):

    partitions1 = []
    partitions2 = []
    for qubit1 in gate1:
        partitions1.append(chromosome[qubit1])
    for qubit2 in gate2:
        partitions2.append(chromosome[qubit2])
    # Check if the partitions are the same for both gates
    return all(partition1 == partition2 for partition1,partition2 in zip(sorted(partitions1), sorted(partitions2)))


def fitness_function(gate_list, chromosome):

    global_gates = find_global_gates(gate_list,chromosome)
    local_gates = find_local_gates(gate_list,chromosome)
    num_teleportations = 0
    executed_gates = []
    #print(chromosome)

    for gate in gate_list:
        if gate in local_gates:
            executed_gates.append(gate)

    i = 0
    while i < len(global_gates):
        j = i + 1
        while j < len(global_gates):

            # Check if the global gates cross the same partitions and share the same qubit
            if check_partition_crossing(global_gates[i],global_gates[j],chromosome) and set(
                    global_gates[i]).intersection(set(global_gates[j])):
                j += 1
            else:
                break
        if j > i + 1:
            executed_gates += global_gates[i:j]
            num_teleportations += 2
            i = j
        else:
            executed_gates.append(global_gates[i])
            num_teleportations += 2
            i += 1

    return num_teleportations


def population_fitness(gate_list, population_list):
    fitness_list = []
    for i in range(len(population_list)):
        fitness_list.append(fitness_function(gate_list,population_list[i]))
    return fitness_list


def tournament_selection(fitness_list, chromosome_list, tournament_size, elite_count):

    selected_individuals = []
    elite_indices = sorted(range(len(fitness_list)),key=lambda k:fitness_list[k])[:elite_count]
    for i in elite_indices:
        selected_individuals.append(chromosome_list[i])

    # Perform tournament selection for the remaining individuals
    while len(selected_individuals) < len(chromosome_list):
        tournament_indices = random.sample(range(len(chromosome_list)),tournament_size)
        tournament_fitness = [fitness_list[i] for i in tournament_indices]
        tournament_chromosomes = [chromosome_list[i] for i in tournament_indices]

        best_index = tournament_fitness.index(min(tournament_fitness))
        best_chromosome = tournament_chromosomes[best_index]
        selected_individuals.append(best_chromosome)

    return selected_individuals


def pmx(parent1, parent2):
    cut1 = random.randint(0, len(parent1)-2)
    cut2 = random.randint(cut1+1, len(parent1)-1)

    offspring = parent1.copy()

    # hash table to map elements in parent1 to parent2
    mapping = {parent1[i]: parent2[i] for i in range(cut1, cut2+1)}

    # Map the segment between the crossover points from parent1 to parent2
    for i in range(cut1, cut2+1):
        if parent1[i] != parent2[i]:
            # Find the corresponding element in parent2 using the hash table
            j = parent2.index(mapping[parent1[i]])
            offspring[i], offspring[j] = offspring[j], offspring[i]

    return offspring


def crossover(selected_population_list, population_size, elite_count):
    elite = selected_population_list[:elite_count]
    non_elite = selected_population_list[elite_count:]
    new_population = []
    offspring_count = population_size - len(elite)
    # Generate the offspring
    for i in range(offspring_count):
        father = random.choice(non_elite)
        mother = random.choice(non_elite)
        if np.random.rand() < Crossover_rate:
            child = pmx(father, mother)
            child = swap_mutation(child)
            new_population.append(child)
    population = elite + new_population
    return population

def swap_mutation(child):
    if np.random.rand() < Mutation_rate:
        pos1, pos2 = random.sample(range(len(child)), 2)
        child[pos1], child[pos2] = child[pos2], child[pos1]
    return child

def visualize_global_gates():

    circuits = ['Sym9_147','Rd73_140','QFT4','QFT8','QFT16','4gt5_76','4mod7','Figure4','Parity_247']
    k2_global_gates = [45,29,8,32,128,13,18,4,8]
    k3_global_gates = [65,40,10,42,170,21,24,5,11]
    k4_global_gates = [76,47,12,48,192,24,30,7,12]

    x_indexes = range(len(circuits))
    width = 0.25
    plt.bar(x_indexes,k2_global_gates,width=width,label='K=2')
    plt.bar([i + width for i in x_indexes],k3_global_gates,width=width,label='K=3')
    plt.bar([i + width * 2 for i in x_indexes],k4_global_gates,width=width,label='K=4')

    for i,v in enumerate(k2_global_gates):
        plt.text(i - width / 2 - 0.15,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k3_global_gates):
        plt.text(i + width / 2 - 0.15,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k4_global_gates):
        plt.text(i + width * 1.5,v + 1,str(v),color='black',fontweight='bold',fontsize='small')

    plt.ylabel('Anzahl an globalen Gatter')
    plt.title('Genetische Algorithmen')

    plt.xticks([i + width for i in x_indexes],circuits)
    plt.xticks(rotation=45)

    plt.legend()
    plt.show()



def visualize_teleportation_cost():

    circuits = ['Sym9_147','Rd73_140','QFT4','QFT8','QFT16','4gt5_76','4mod7','Figure4','Parity_247']
    k2_telep_cost = [12,10,4,8,16,4,4,2,2]
    k3_telep_cost = [34,18,6,14,30,10,10,4,4]
    k4_telep_cost = [52,32,12,24,52,22,32,12,6]

    x_indexes = range(len(circuits))
    width = 0.25
    plt.bar(x_indexes,k2_telep_cost,width=width,label='K=2')
    plt.bar([i + width for i in x_indexes],k3_telep_cost,width=width,label='K=3')
    plt.bar([i + width * 2 for i in x_indexes],k4_telep_cost,width=width,label='K=4')

    for i,v in enumerate(k2_telep_cost):
        plt.text(i - width / 2,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k3_telep_cost):
        plt.text(i + width / 2 - 0.05,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k4_telep_cost):
        plt.text(i + width * 1.5,v + 1,str(v),color='black',fontweight='bold',fontsize='small')

    plt.ylabel('Teleportationskosten')
    plt.title('Genetische Algorithmen Teleportationskosten')
    plt.xticks([i + width for i in x_indexes],circuits)
    plt.xticks(rotation=45)
    plt.legend()
    plt.show()



if __name__ == '__main__':

    input_filename = 'C:/Users/Teodor/QC_GA_Partition/qasm/sym9_147_basic.qasm'
    gate_list = converter_circuit_from_qasm(input_filename)[0]
    gate_list = list_str_to_int(gate_list)
    print('Original Quantum Circuit: ' + str(simple_colors.green(gate_list)))

    print('Number of Gates: ' + simple_colors.green((len(gate_list))))

    num_qubits = converter_circuit_from_qasm(input_filename)[1]
    num_qubits = int(num_qubits)
    print('Number of Qubits: ' + simple_colors.green(num_qubits))

    Crossover_rate = 0.9
    Mutation_rate = 0.1
    N_generations = 100
    population_size = 100
    elite_count = 2
    tournament_size = 2
    num_partitions = 3

    best_fitness = []
    mean_fitness = []
    population_list = initialize_population(num_qubits,num_partitions)
    fitness_list = population_fitness(gate_list,population_list)

    for _ in range(N_generations):
        selected_population_list = tournament_selection(fitness_list,population_list,elite_count,tournament_size)
        population_list = crossover(selected_population_list,population_size,elite_count)
        fitness_list = population_fitness(gate_list,population_list)
        best_fitness.append(np.min(fitness_list))

    plt.plot(range(N_generations),best_fitness,label="Beste Fitness")
    plt.xlabel("Generation")
    plt.ylabel("Fitnesswert")
    plt.title("Fitnesswert über Generationen")
    plt.legend()
    plt.show()

    print(simple_colors.magenta(
        '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'))
    print('Fitness_List: ',end='')
    print(simple_colors.green(best_fitness))
    print('Best Teleportation Cost: ',end='')
    print(simple_colors.green(min(best_fitness)))

    #bar_gg = visualize_global_gates()
    #bar_tc = visualize_teleportation_cost()



