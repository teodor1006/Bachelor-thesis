import re
import simple_colors
import matplotlib.pyplot as plt
import networkx as nx
import heapq
from collections import defaultdict


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



def multi_kerlin(gate_list, num_partitions):

    G = nx.Graph()
    edges = defaultdict(int)
    for gate in gate_list:
        for i in range(1,len(gate)):
            edge = tuple(sorted([gate[i - 1],gate[i]]))
            if edge in edges:
                edges[edge] += 1
            else:
                edges[edge] = 1

    for edge,count in edges.items():
        G.add_edge(edge[0],edge[1],weight=count)

    if num_partitions == 2:
        partition = nx.community.kernighan_lin_bisection(G,seed=0)
    else:
        subgraphs = [G]
        # Keep partitioning the subgraphs until there are num_partitions subgraphs
        while len(subgraphs) < num_partitions:
            # For each subgraph, partition it using Kernighan-Lin bisection
            new_subgraphs = []
            for subgraph in subgraphs:
                for part in nx.community.kernighan_lin_bisection(subgraph,seed=0):
                    new_subgraphs.append(subgraph.subgraph(part))
            subgraphs = new_subgraphs

            # If there are more than num_partitions subgraphs, merge the smallest ones
            if len(subgraphs) > num_partitions:
                # Get the two smallest subgraphs using a heap
                smallest_subgraphs = heapq.nsmallest(2, subgraphs, key=lambda sg: sg.number_of_nodes())
                merged_subgraph = nx.compose_all(smallest_subgraphs)
                # Replace the two smallest subgraphs with the merged subgraph
                subgraphs.remove(smallest_subgraphs[0])
                subgraphs.remove(smallest_subgraphs[1])
                subgraphs.append(merged_subgraph)

        partition = [set(subgraph.nodes) for subgraph in subgraphs]

    return partition

def graph_qubits(gate_list):
    G = nx.Graph()

    # Count the number of times each edge appears in the gate_list
    edge_counts = defaultdict(int)
    for gate in gate_list:
        for i in range(1,len(gate)):
            edge = tuple(sorted([gate[i - 1],gate[i]]))
            if edge in edge_counts:
                edge_counts[edge] += 1
            else:
                edge_counts[edge] = 1

    for edge,count in edge_counts.items():
        G.add_edge(edge[0],edge[1],weight=count)

    pos = nx.spring_layout(G)

    nx.draw(G,pos,with_labels=True,font_weight='bold')
    edge_labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)

    plt.show()

def graph_after_partitioning(gate_list,num_partitions):

    G = nx.Graph()
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

    partitions = multi_kerlin(gate_list,num_partitions)
    colors = []
    for node in G.nodes():
        if node in partitions[0]:
            colors.append("red")
        elif node in partitions[1]:
            colors.append("green")
        elif node in partitions[2]:
            colors.append("gray")
        else:
            colors.append("orange")
    color = colors
    pos = nx.spectral_layout(G)
    nx.draw(G,pos,node_color=color,with_labels=True)
    edge_labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels,font_weight='bold')
    nx.draw_networkx_labels(G,pos)
    plt.show()


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


def count_local_gates(gate_list,part_qbits):
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


def check_partition_crossing(gate1, gate2, part_qbits):

    partitions1 = []
    partitions2 = []
    for qubit1 in gate1:
        for part_index,partition in enumerate(part_qbits):
            if qubit1 in partition:
                partitions1.append(part_index)
                break
    for qubit2 in gate2:
        for part_index,partition in enumerate(part_qbits):
            if qubit2 in partition:
                partitions2.append(part_index)
                break
    # Check if the partitions are the same for both gates
    return all(partition1 == partition2 for partition1,partition2 in zip(sorted(partitions1), sorted(partitions2)))


def min_telep(gate_list,part_qbits):

    global_gates = find_global_gates(gate_list,part_qbits)
    local_gates = find_local_gates(gate_list,part_qbits)
    executed_gates = []
    num_teleportations = 0

    for gate in gate_list:
        if gate in local_gates:
            executed_gates.append(gate)

    i = 0
    while i < len(global_gates):
        j = i + 1
        while j < len(global_gates):

            # Check if the global gates cross the same partitions and share the same qubit
            if check_partition_crossing(global_gates[i],global_gates[j],part_qbits) and set(
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

def visualize_global_gates():

    circuits = ['Sym9_147','Rd73_140','QFT4','QFT8','QFT16','4gt5_76','4mod7','Figure4','Parity_247']
    k2_global_gates = [45,29,8,32,128,13,18,4,8]
    k3_global_gates = [58,47,10,40,160,24,29,6,12]
    k4_global_gates = [76,47,12,48,192,24,30,7,12]

    x_indexes = range(len(circuits))
    width = 0.25
    plt.bar(x_indexes,k2_global_gates,width=width,label='K=2')
    plt.bar([i + width for i in x_indexes],k3_global_gates,width=width,label='K=3')
    plt.bar([i + width * 2 for i in x_indexes],k4_global_gates,width=width,label='K=4')

    # Add exact number of global gates to each bar
    for i,v in enumerate(k2_global_gates):
        plt.text(i - width / 2 - 0.15,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k3_global_gates):
        plt.text(i + width / 2 - 0.15,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k4_global_gates):
        plt.text(i + width * 1.5,v + 1,str(v),color='black',fontweight='bold',fontsize='small')

    plt.ylabel('Anzahl an globalen Gatter')
    plt.title('Kernighan-Lin-Algorithmus')

    plt.xticks([i + width for i in x_indexes],circuits)
    plt.xticks(rotation=45)

    plt.legend()
    plt.show()


def visualize_teleportation_cost():

    circuits = ['Sym9_147','Rd73_140','QFT4','QFT8','QFT16','4gt5_76','4mod7','Figure4','Parity_247']
    k2_telep_cost = [14,10,4,12,32,6,4,2,2]
    k3_telep_cost = [26,36,6,30,108,30,42,10,14]
    k4_telep_cost = [60,36,12,48,158,30,50,12,16]

    x_indexes = range(len(circuits))
    width = 0.25
    plt.bar(x_indexes,k2_telep_cost,width=width,label='K=2')
    plt.bar([i + width for i in x_indexes],k3_telep_cost,width=width,label='K=3')
    plt.bar([i + width * 2 for i in x_indexes],k4_telep_cost,width=width,label='K=4')

    for i,v in enumerate(k2_telep_cost):
        plt.text(i - width / 2 - 0.05,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k3_telep_cost):
        plt.text(i + width / 2 - 0.15,v + 1,str(v),color='black',fontweight='bold',fontsize='small')
    for i,v in enumerate(k4_telep_cost):
        plt.text(i + width * 1.5,v + 1,str(v),color='black',fontweight='bold',fontsize='small')

    plt.ylabel('Teleportationskosten')
    plt.title('Kernighan-Lin-Algorithmus Teleportationskosten')

    plt.xticks([i + width for i in x_indexes],circuits)
    plt.xticks(rotation=45)

    plt.legend()
    plt.show()



if __name__ == '__main__':

    input_filename = 'C:/Users/Teodor/QC_GA_Partition/qasm/4mod7-v0_94_basic.qasm'
    gate_list = converter_circuit_from_qasm(input_filename)[0]
    gate_list = list_str_to_int(gate_list)
    print('Original Quantum Circuit: ' + str(simple_colors.green(gate_list)))

    print('Number of Gates: ' + simple_colors.green((len(gate_list))))

    num_qubits = converter_circuit_from_qasm(input_filename)[1]
    num_qubits = int(num_qubits)
    print('Number of Qubits: ' + simple_colors.green(num_qubits))

    num_partitions = 2

    part_qbits = multi_kerlin(gate_list,num_partitions)
    print('Partitioned qubits: ' + simple_colors.green(part_qbits))


    global_gates_num = count_global_gates(gate_list,part_qbits)
    print('Number of Global Gates: ' + simple_colors.red(global_gates_num))
    global_gates_list = find_global_gates(gate_list,part_qbits)
    print('List of Global Gates: ' + simple_colors.red(str(global_gates_list)))


    local_gates_num = count_local_gates(gate_list,part_qbits)
    print('Number of Local Gates: ' + simple_colors.blue(local_gates_num))
    local_gates_list = find_local_gates(gate_list,part_qbits)
    print('List of Local Gates: ' + simple_colors.blue(str(local_gates_list)))


    tc = min_telep(gate_list,part_qbits)
    print('Best Teleportation Cost: ',end='')
    print(simple_colors.green(tc))

    #graph = graph_qubits(gate_list)
    #graph_part = graph_after_partitioning(gate_list,num_partitions)

    #bar_gg = visualize_global_gates()
    #bar_tc = visualize_teleportation_cost()