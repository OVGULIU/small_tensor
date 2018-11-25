import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('strain_stress.txt')

strain = data[:,0]
stress = data[:,1]

plt.plot(strain, stress)
plt.show()


