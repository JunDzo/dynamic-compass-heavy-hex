import stim
import sys
import numpy as np

def index_to_coord(i: int, d: int):
    # data qubits are at (x = odd, y = even)
    # z ancilla are at (x = even, y = even) for 0 < x < 2d
    # x ancilla are at (x = even, y = odd) for x = 2 mod 4, y = 1 mod 4 or x = 0 mod 4 and y = 3 mod 4 
    coord = np.array([i % (2*d+1), i // (2*d+1)])
    return coord

def coord_to_index(coord: np.ndarray, d: int):
    i = coord[0] + coord[1]*(2*d+1)
    return i

def translate_index(i: int, translation: list[int,int], d: int):
    translation = np.array(translation)
    new_coord = index_to_coord(i,d) + translation
    if new_coord[0] < 0 or new_coord[0] > 2*d:
        return -1
    return coord_to_index(new_coord, d)

def list_qubits(d: int):
    data_qubits = []
    x_ancillas = []
    z_ancillas = []
    for i in range((2*d+1)*(2*d-1)):
        coord = index_to_coord(i,d)
        if (coord[0] % 2 == 1 and coord[1] % 2 == 0):
            data_qubits.append(i)
        elif ((coord[0] % 4 == 2 and coord[1] % 4 == 1) 
          or (coord[0] % 4 == 0 and coord[1] % 4 == 3)):
            x_ancillas.append(i)
        elif (coord[0] % 2 == coord[1] % 2 == 0 and 0 < coord[0] < 2*d):
            z_ancillas.append(i)
    return data_qubits, x_ancillas, z_ancillas

def apply_cnot(q1: int, q2: int, p_2q: float):
    circuit = stim.Circuit()
    if q1 == -1 or q2 == -1:
        return circuit
    circuit.append("CNOT", [q1, q2])
    circuit.append("DEPOLARIZE2", [q1, q2], p_2q)
    return circuit

def idle_errors(i: int, idle_qubits: list[list[int,int]], p_idle: float, d: int):
    circuit = stim.Circuit()
    if index_to_coord(i,d)[0] % 4 == 2:
        idle_qubits = [translate_index(i, translation, d) for translation in idle_qubits]
        idle_qubits = [i for i in idle_qubits if i != -1]
        circuit.append("DEPOLARIZE1", idle_qubits, p_idle)
    return circuit

def make_x_circuit(
    p_1q: float,
    p_2q: float,
    p_idle: float,
    p_reset: float,
    p_meas: float,
    x_ancillas: list[int],
    d: int
):
    circuit = stim.Circuit()
    for i in x_ancillas:
        circuit.append("RX", i)
        circuit.append("Z_ERROR", i, p_reset)
        circuit.append("RZ", translate_index(i,[0,1],d))
        circuit.append("X_ERROR", translate_index(i,[0,1],d), p_reset)
        circuit += idle_errors(i, [[-1,1],[1,1],[-1,-1],[1,-1]], p_idle, d)
    for i in x_ancillas:
        circuit.append("RZ", translate_index(i,[0,-1],d))
        circuit.append("X_ERROR", translate_index(i,[0,-1],d), p_reset)
        circuit += apply_cnot(i, translate_index(i,[0,1],d), p_2q)
        circuit += idle_errors(i, [[-1,1],[1,1],[-1,-1],[1,-1]], p_idle, d)
    circuit.append("TICK")
    for i in x_ancillas:
        circuit += apply_cnot(i, translate_index(i,[0,-1],d), p_2q)
        circuit += apply_cnot(translate_index(i,[0,1],d), translate_index(i, [-1,1],d), p_2q)
        circuit += idle_errors(i, [[1,1],[-1,-1],[1,-1]], p_idle, d)
    circuit.append("TICK")
    for i in x_ancillas:
        circuit += apply_cnot(translate_index(i,[0,-1],d), translate_index(i,[1,-1],d), p_2q)
        circuit += apply_cnot(translate_index(i,[0,1],d), translate_index(i,[1,1],d), p_2q)
        circuit += idle_errors(i, [[-1,1],[-1,-1]], p_idle, d)
        circuit.append("DEPOLARIZE1", i, p_idle)
    circuit.append("TICK")
    for i in x_ancillas:
        circuit += apply_cnot(i, translate_index(i,[0,1],d), p_2q)
        circuit += apply_cnot(translate_index(i,[0,-1],d), translate_index(i,[-1,-1],d), p_2q)
        circuit += idle_errors(i, [[-1,1],[1,1],[1,-1]], p_idle, d)
    circuit.append("TICK")
    for i in x_ancillas:
        circuit.append("X_ERROR", translate_index(i,[0,1],d), p_meas)
        circuit.append("MZ", translate_index(i,[0,1],d))
        circuit.append("DETECTOR", [stim.target_rec(-1)])
        circuit += apply_cnot(i, translate_index(i,[0,-1],d), p_2q)
        circuit += idle_errors(i, [[-1,1],[1,1],[-1,-1],[1,-1]], p_idle, d)
    circuit.append("TICK")
    for i in x_ancillas:
        circuit.append("X_ERROR", translate_index(i,[0,-1],d), p_meas)
        circuit.append("MZ", translate_index(i,[0,-1],d))
        circuit.append("DETECTOR", [stim.target_rec(-1)])
        circuit += idle_errors(i, [[-1,1],[1,1],[-1,-1],[1,-1]], p_idle, d)
    for i in x_ancillas:
        #this is the same timestep as above
        #but i separate them for easier measurement indexing
        circuit.append("Z_ERROR", i, p_meas)
        circuit.append("MX", i)
    circuit.append("TICK")
    return circuit

def make_z_circuit(
    p_1q: float,
    p_2q: float,
    p_idle: float,
    p_reset: float,
    p_meas: float,
    z_ancillas: list[int],
    d: int
):
    circuit = stim.Circuit()
    for i in z_ancillas:
        circuit.append("RZ", i)
        circuit.append("X_ERROR", i, p_reset)
        circuit += idle_errors(i, [[-1,0],[1,0]], p_idle, d)
    circuit.append("TICK")
    for i in z_ancillas:
        circuit += apply_cnot(i, translate_index(i,[-1,0],d), p_2q)
        circuit += idle_errors(i, [[1,0]], p_idle, d)
    circuit.append("TICK")
    for i in z_ancillas:
        circuit += apply_cnot(i, translate_index(i,[1,0],d), p_2q)
        circuit += idle_errors(i, [[-1,0]], p_idle, d)
    circuit.append("TICK")
    for i in z_ancillas:
        circuit.append("X_ERROR", i, p_meas)
        circuit.append("MZ", i)
        circuit += idle_errors(i, [[-1,0],[1,0]], p_idle, d)
    circuit.append("TICK")
    return circuit

def x_readout(
    p_meas: float, 
    data_qubits: list[int], 
    x_ancillas: list[int], 
    n_x: int,
    n_z: int, 
    d: int
):
    circuit = stim.Circuit()
    circuit.append("Z_ERROR", data_qubits, p_meas)
    circuit.append("RX", 0) #this is part of a hack to make all stabilisers weight 4
    for i in range(n_x):
        q = x_ancillas[i]
        support = [[-1,1],[1,1],[1,-1],[-1,-1]]
        support = [translate_index(q,j,d) for j in support]
        support = [j for j in support if j != -1]
        while len(support) < 4:
            #add previously unused qubit 0 to support repeatedly until length is 4
            #so it is easier to count back over the measurement record when building detectors
            support.append(0)
        circuit.append("MX", support)
        circuit.append("DETECTOR", [stim.target_rec(-1), 
                                    stim.target_rec(-2),  
                                    stim.target_rec(-3),  
                                    stim.target_rec(-4),  
                                    stim.target_rec(-(4*(i+1) + n_z + n_x - i))])
    return circuit

def z_readout(
    p_meas: float,
    data_qubits: list[int],
    z_ancillas: list[int],
    n_x: int,
    n_z: int,
    d: int
):
    circuit = stim.Circuit()
    circuit.append("X_ERROR", data_qubits, p_meas)
    for i in range(n_z):
        q = z_ancillas[i]
        circuit.append("MZ", translate_index(q, [-1,0]))
        circuit.append("MZ", translate_index(q, [1,0]))
        circuit.append("DETECTOR", [stim.target_rec(-1),
                                    stim.target_rec(-2),
                                    stim.target_rec(-(2*(i+1) + n_z - i))])
    return circuit

def make_full_circuit(
    p_1q: float,
    p_2q: float,
    p_idle: float,
    p_reset: float,
    p_meas: float,
    d: int,
    n_cycles: int,
    init_basis: str
):
    #setup
    circuit = stim.Circuit()
    data_qubits, x_ancillas, z_ancillas = list_qubits(d)
    n_x = len(x_ancillas)
    n_z = len(z_ancillas)
    x_circuit = make_x_circuit(p_1q, p_2q, p_idle, p_reset, p_meas, x_ancillas, d)
    z_circuit = make_z_circuit(p_1q, p_2q, p_idle, p_reset, p_meas, z_ancillas, d)
    
    #init
    if init_basis == 'z':
        circuit.append("RZ", data_qubits)
        circuit.append("X_ERROR", data_qubits, p_reset)
    elif init_basis == 'x':
        circuit.append("RX", data_qubits)
        circuit.append("Z_ERROR", data_qubits, p_reset)
    else:
        sys.exit("Error: initialisation basis must be x or z")

    #first cycle
    circuit += x_circuit
    if init_basis == 'x':
        for i in range(1,n_x+1):
            circuit.append("DETECTOR", [stim.target_rec(-i)])
    circuit += z_circuit
    if init_basis == 'z':
        for i in range(1,n_z+1):
            circuit.append("DETECTOR", [stim.target_rec(-i)]) 

    #other cycles
    if n_cycles > 1:
        for i in range(1,n_x+1):
            x_circuit.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-(3*n_x + n_z + i))])
        for i in range(1,n_z+1):
            z_circuit.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-(3*n_x + n_z + i))])
        joint_circuit = x_circuit + z_circuit
        repeated_circuit = stim.Circuit(
            """REPEAT """ + str(n_cycles-1) + """ {\n""" + str(joint_circuit) + """\n}"""
        )
        circuit += repeated_circuit

    #readout
    if init_basis == 'x':
        circuit += x_readout(p_meas, data_qubits, x_ancillas, n_x, n_z, d)
        x_observable = [coord_to_index(np.array([1+2*i,0]),d) for i in range(d)]
        circuit.append("MX", x_observable)
        circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-(i+1)) for i in range(d)])
    elif init_basis == 'z':
        circuit += z_readout(p_meas, data_qubits, z_ancillas, n_x, n_z, d)
        z_observable = [coord_to_index(np.array([1,2*i]),d) for i in range(d)]
        circuit.append("MZ", z_observable)
        circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-(i+1)) for i in range(d)])

    return circuit

if __name__ == "__main__":
    d = 3
    n_cycles = 2
    init_basis = 'x'

    p_1q = 0.001
    p_2q = 0.001
    p_idle = 0.001
    p_reset = 0.001
    p_meas = 0.001

    circuit = make_full_circuit(p_1q, p_2q, p_idle, p_reset, p_meas, d, n_cycles, init_basis)  
    with open('circuit.svg', 'w') as f:
        f.write(str(circuit.diagram('timeline-svg')))
