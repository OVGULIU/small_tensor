import numpy as np 
import matplotlib.pyplot as plt 
# #############################
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
# #############################

data = np.loadtxt('strain_stress.txt')

strain = data[:,0]
stress = data[:,1]

plt.plot(strain, stress)
plt.xlabel('Compressive Strain')
plt.ylabel('Compressive Stress')
plt.savefig("ConcreteCompressCyclic.jpg")
plt.show()



