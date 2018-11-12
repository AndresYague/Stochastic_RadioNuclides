# Python packages
import numpy as np
import glob

# Path results
path = './'

# Read list of tau
with open(path+'tauList.in', 'r') as ff:
    ff.readline() # Skip header
    line = ff.readline() # Actual values
    
    line = (line.replace('D','E')).split()
    nb_tau = len(line)
    tauList = np.array(map(float, line))

# Read list of times
with open(path+'tArray.in', 'r') as ff:
    ff.readline() # Skip header
    line = ff.readline() # Actual values
    
    line = (line.replace('D','E')).split()
    nb_t = len(line)
    tArray = np.array(map(float, line))

# Get output files
files = sorted(glob.glob("%s/Output*" %path))

# Declare the radioactive abundances array
radio_ab = []
for i_tau in range(nb_tau):
    radio_ab.append([])
    
# Read evolution of radioactive abundances
for the_file in files:
    with open(the_file, 'r') as ff:
        not_header = False
        i_tau = -1
        for line in ff:
            if '#' in line:
                i_tau += 1
                not_header = True
            elif not_header:
                line = map(float, line.split()[1:])
                radio_ab[i_tau].append(line)

# Transpose here
for i_tau in range(nb_tau):
    radio_ab[i_tau] = np.transpose(radio_ab[i_tau])

nb_run = len(radio_ab[0][0])

# Get statistics
minn = np.zeros((nb_tau,nb_t))
perc_5  = np.zeros((nb_tau,nb_t))
perc_32 = np.zeros((nb_tau,nb_t))
perc_50 = np.zeros((nb_tau,nb_t))
perc_68 = np.zeros((nb_tau,nb_t))
perc_95 = np.zeros((nb_tau,nb_t))
maxx = np.zeros((nb_tau,nb_t))
mean = np.zeros((nb_tau,nb_t))
for i_tau in range(nb_tau):
    for i_t in range(nb_t):
        minn[i_tau][i_t] = min(radio_ab[i_tau][i_t])
        perc_5[i_tau][i_t]  = np.percentile(radio_ab[i_tau][i_t], 5)
        perc_32[i_tau][i_t] = np.percentile(radio_ab[i_tau][i_t], 32)
        perc_50[i_tau][i_t] = np.percentile(radio_ab[i_tau][i_t], 50)
        perc_68[i_tau][i_t] = np.percentile(radio_ab[i_tau][i_t], 68)
        perc_95[i_tau][i_t] = np.percentile(radio_ab[i_tau][i_t], 95)
        mean[i_tau][i_t] = np.mean(radio_ab[i_tau][i_t])
        maxx[i_tau][i_t] = np.max(radio_ab[i_tau][i_t])
np.save('radioStat', [minn,perc_5,perc_32,perc_50,perc_68,perc_95,maxx,mean])
