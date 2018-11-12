# Python packages
import numpy as np
import glob

# Path results
path = './'

# Read list of tau
with open(path+'tauList.in', 'r') as ff:
    is_first = True
    for line in ff:
        if is_first:
            is_first = False
        line = (line.replace('D','E')).split()
        nb_tau = len(line)
        tauList = np.zeros(nb_tau)
        for i in range(nb_tau):
            tauList[i] = float(line[i])
ff.close()

# Read list of times
with open(path+'tArray.in', 'r') as ff:
    is_first = True
    for line in ff:
        if is_first:
            is_first = False
        line = (line.replace('D','E')).split()
        nb_t = len(line)
        tArray = np.zeros(nb_t)
        for i in range(nb_t):
            tArray[i] = float(line[i])
ff.close()

# Get output files
files = sorted(glob.glob("%s/Output*" %path))
nb_files = len(files)

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
                line = [float(x) for x in line.split()[1:]]
                radio_ab[i_tau].append(np.array(line))
    ff.close()
nb_run = len(radio_ab[0])

# Get statistics
minn = np.zeros((nb_tau,nb_t))
perc_5  = np.zeros((nb_tau,nb_t))
perc_32 = np.zeros((nb_tau,nb_t))
perc_50 = np.zeros((nb_tau,nb_t))
perc_68 = np.zeros((nb_tau,nb_t))
perc_95 = np.zeros((nb_tau,nb_t))
maxx = np.zeros((nb_tau,nb_t))
mean = np.zeros((nb_tau,nb_t))
#list_abundance_theorem = []
#f_ini = 0.6
#i_t_start = int(nb_t*f_ini)
for i_tau in range(nb_tau):
#    list_abundance_theorem.append([])
    for i_t in range(nb_t):
        list_abundance = np.zeros(nb_run)
        for i_run in range(nb_run):
            list_abundance[i_run] = radio_ab[i_tau][i_run][i_t]
#            if i_t >= i_t_start:
#                list_abundance_theorem[i_tau].append(radio_ab[i_tau][i_run][i_t])
        minn[i_tau][i_t] = min(list_abundance)
        perc_5[i_tau][i_t]  = np.percentile(list_abundance, 5)
        perc_32[i_tau][i_t] = np.percentile(list_abundance, 32)
        perc_50[i_tau][i_t] = np.percentile(list_abundance, 50)
        perc_68[i_tau][i_t] = np.percentile(list_abundance, 68)
        perc_95[i_tau][i_t] = np.percentile(list_abundance, 95)
        mean[i_tau][i_t] = np.mean(list_abundance)
        maxx[i_tau][i_t] = max(list_abundance)
np.save('radioStat', [minn,perc_5,perc_32,perc_50,perc_68,perc_95,maxx,mean])

