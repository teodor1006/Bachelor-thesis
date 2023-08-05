from qiskit import *
import matplotlib.pyplot as plt

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

fig, ax = plt.subplots(figsize=(8, 8))
qc.draw(output='mpl', ax=ax, vertical_compression='tight')

line_1 = -1.5
line_2 = -2.5

ax.axhline(line_1, color='red', linestyle='dashed', linewidth=2)
ax.axhline(line_2, color='red', linestyle='dashed', linewidth=2)

ax.text(-2.5, line_1 + 0.8, 'P1', fontsize=14, color='red')
ax.text(-2.5, line_2 + 0.4, 'P2', fontsize=14, color='red')
ax.text(-2.5, -3.1, 'P3', fontsize=14, color='red')

plt.show()



