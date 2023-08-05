import re
import simple_colors
import random

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



def random_partition(num_qubits, num_partitions):
    for _ in range(200):  # Run for 200 iterations
        qubit_indices = list(range(num_qubits))
        random.shuffle(qubit_indices)
        partition_size = num_qubits // num_partitions
        partitions = [qubit_indices[i*partition_size:(i+1)*partition_size] for i in range(num_partitions)]
        for i in range(num_qubits % num_partitions):
            partitions[i].append(qubit_indices[num_partitions*partition_size + i])
        return partitions


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


if __name__ == '__main__':

    input_filename = 'C:/Users/Teodor/QC_GA_Partition/qasm/sym9_147_basic.qasm'
    gate_list = converter_circuit_from_qasm(input_filename)[0]
    gate_list = list_str_to_int(gate_list)
    print('Original Quantum Circuit: ' + str(simple_colors.green(gate_list)))

    print('Number of Gates: ' + simple_colors.green((len(gate_list))))

    num_qubits = converter_circuit_from_qasm(input_filename)[1]
    num_qubits = int(num_qubits)
    print('Number of Qubits: ' + simple_colors.green(num_qubits))

    num_partitions = 2

    part_qbits = random_partition(num_qubits,num_partitions)
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