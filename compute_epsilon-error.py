import numpy as np
import os
import sys
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

##########################################
# Parse command line

if not(os.path.isfile('runtime.inpt')) or not(os.path.isfile('data.inpt')) or len(sys.argv) < 2:
   print('''
Usage: python compute_capacitance.py polarization_filename
''')
   quit()
##########################################
# Constants
k = 1.3806485279e-23 #J.K-1
e = 1.602176620898e-19 #C
bohr2ang = 0.52917721067
ang2m = 1e-10
eps0 = 8.85418782e-12
##########################################

def read_info():
   # Read temperature and potential difference
   with open('runtime.inpt') as run:
     for line in run:
       if 'temperature' in line:
         temp = float(line.split()[1])

   # Read box parameters and num elec
   with open('data.inpt') as run:
     for line in run:
       if 'box' in line:
         line2 = run.readline()
         cellx = float(line2.split()[0])
         celly = float(line2.split()[1])
         cellz = float(line2.split()[2])
   volume = cellx * celly * cellz * (bohr2ang * ang2m)**3 #m^3

   return temp, volume
##########################################

temp, volume = read_info()
print("System parameters:\nTemperature {} K\nVolume {} m^3".format(temp, volume))
print("=========================================")

for polname in sys.argv[1:]:
   if os.path.isfile(polname):
      print("Reading polarization from {}".format(polname))
      total_pol = np.loadtxt(polname, comments = '#', unpack=True) 
   else:
      print("Didn't find file {}".format(polname))
      quit()

nframes = np.shape(total_pol)[0]
print(nframes)
#total_pol = []
#for i in range(0, nframes):
#   total_pol.append(sqrt(pol_x[i]**2 + pol_y[i]**2 + pol_z[i]**2))

mean_pol = np.mean(total_pol)

# Autocorrelation function
print("Computing autocorrelation function")
print("==================================")
Cqq = np.fft.fftn(total_pol - mean_pol)
CC = Cqq[:] * np.conjugate(Cqq[:])
CC[:] = np.fft.ifftn(CC[:])
autocorr = (CC[:len(CC)//2]).real

# Integration of the squared normalized autocorrelation function
Tccum = cumtrapz((autocorr[:]/autocorr[0])**2, dx = 1.0)
plt.plot(np.array(range(nframes//2)), (autocorr[:]/autocorr[0])**2)
plt.grid(True)
plt.ylabel(r'Correlation function $<Q(0) Q(t)>^2/<Q^2>^2$')
plt.xlabel('Time (frames)')
plt.title("What is the correlation time Tc?")
plt.show()
Tc = int(input("Choose the correlation time Tc (in frames):"))
print("==========================================")

# Compute Differential Capacitance
print("Dielectric constant calculation")
epsilon = 1 + (np.mean(total_pol*total_pol) * volume) / (3 * k * temp * eps0)
error_epsilon = np.sqrt(4 * Tccum[Tc] / nframes) * epsilon
print("Epsilon = {} +/- {} ".format(epsilon, error_epsilon))
print("====================================")

#def gaussian(x, mu, var):
#    return 1/np.sqrt(2*math.pi*var) * np.exp(-(x - mu)**2 / (2 * var))
#
#plt.figure(figsize=(20,10))
#plt.hist(equil_charges, bins=200, histtype='step', density=True, linewidth=2, label='Total charge on left electrode')
#x = np.array(np.arange(mean_charge - 20*variance_charge, mean_charge + 20*variance_charge, 0.01))
#plt.plot(x, gaussian(x, mean_charge, variance_charge), label='Gaussian distribution of mean {} and variance {}'.format(mean_charge, variance_charge))
#plt.xlabel(r'Total_charge $Q=\sum_i q_i$')
#plt.legend()
#plt.savefig("charge_distribution.png", bbox_inches='tight')
#plt.show()
#
# Uncomment to obtain logarithmic plot of the charge distribution
# as in Limmer et al. Charge Fluctuations in Nanoscale Capacitors (PRL 2013)
#hist, bin_edge = np.histogram(equil_charges - mean_charge, bins=200, density=True)
#xvalues = (bin_edge[1:] + bin_edge[:-1]) / (2 * np.sqrt(variance_charge))
#yvalues = -((bin_edge[1:] + bin_edge[:-1])/2)**2 / (2 * variance_charge)

#plt.figure(figsize=(10,10))
#plt.semilogy(xvalues, hist)
#plt.semilogy(xvalues, 1/math.sqrt(2* math.pi * variance_charge) * np.exp(yvalues), label='Gaussian function')
#plt.ylabel(r'Total charge distribution $P(Q)$')
#plt.xlabel(r'$\delta Q / \sqrt{<\delta Q^2>}$')
#plt.legend()
#plt.savefig("log_charge_distribution.png", bbox_inches='tight')
#plt.show()

