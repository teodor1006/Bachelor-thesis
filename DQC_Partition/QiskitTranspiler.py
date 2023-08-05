import matplotlib.pyplot as plt
from qiskit import *
#import matplotlib.pyplot as plt
#from math import pi

qc = QuantumCircuit(4,4)

qc.cx(0,1)
qc.cx(2,0)
qc.cx(0,3)
qc.h(3)
qc.cx(3,1)
qc.h(1)
qc.cx(1,3)
qc.cx(0,1)
qc.cx(3,0)

#qc.draw(output='mpl')
#plt.show()
#print(qc)


target_basis = ['h','cx','ry','rz','rx','x']
decomposed = transpile(qc,basis_gates=target_basis,optimization_level=3)
decomposed.draw(output='mpl')


qasm_str = decomposed.qasm()
with open("qasm/4gt5_76_basic.qasm","w") as f:
    f.write(qasm_str)
