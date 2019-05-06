# Python packages
import numpy as np
import random

# ======================================================

# Define PDF - Square Box
t_min_sb, t_max_sb = '5e7', '1e10'

# Define the radioactive setup
gamma = '1.00e8'

# Define the run time
maxTime = '15e9'
nb_run = '1e3'

# Output file name
output_name = 'tEvents_box_'+t_min_sb+'_'+t_max_sb+'_gamma_'+gamma+'_tend_'+maxTime+'_'+nb_run+'runs.in'
t_min_sb = float(t_min_sb)
t_max_sb = float(t_max_sb)
kwargs = {"t_min":t_min_sb, "t_max":t_max_sb}
gamma = float(gamma)
maxTime = float(maxTime)
nb_run = int(float(nb_run))

print('Estimated file size:', '%.2E' %((maxTime/gamma)*nb_run*15.0))

# ======================================================

# Keep track of the time
import time as t_module
def get_time(start_time):
    out = 'Run time: ' + \
    str(round((t_module.time() - start_time),10))+"s"
    return out


# Power-Law PDF
def sample_power_law(t_min, t_max, pw_index, norm):
    '''
    Randomly sample a delay time from a power-law distribution.
    
    Arguments
    =========
      - t_min, t_max: Time boundaries of the distribution.
      - pw_index: Power-law index (a) in PDF(t) = t^a.
      - norm: Normalization factor, such that max(PDF) = 1.0.
    '''
    
    # Calculate quantities to limit the number of operations in the loop
    dt_boundary = t_max - t_min
    
    # Loops that searches a valid delay time
    while True:
        
        # Randomly select a coordinate in the t-PDF space
        t_delay = random.random()*dt_boundary + t_min
        yy = random.random()
        
        # Return the delay time if under the PDF curve
        if yy < norm*t_delay**pw_index:
            return t_delay


# Squate-box PDF, where each time value is equally probable
def sample_square_box(t_min, t_max):
    '''
    Randomly sample a delay time from a uniform distribution.
    
    Argument
    =========
      - t_min, t_max: Time boundaries of the distribution.
    '''
    
    # Calculate quantities to limit the number of operations in the loop
    dt_boundary = t_max - t_min
    return random.random()*dt_boundary + t_min


def get_list_t_events(delta, nb_events, pdf, **kwargs):
    '''
    Get list of times where events occurs, assuming a constant
    star formation history.  Events refer to yields ejection
    events following the formation of a progenitor object.
    
    Arguments
    =========
      - delta: Constant recurrence time for the formation of progenitor.
      - nb_events: Total number of events to be calculated.
      - pdf: Probability distribution function for the event times.
      - kwargs: Input parameters for the pdf function.
    '''
    
    # Declare the list of times where events occur
    t_events = np.zeros(nb_events)
    
    # Get fill the list of times
    # "i_e*delta" is the formation time of the i_e progenitor
    for i_e in range(nb_events):
        t_events[i_e] = i_e*delta + pdf(**kwargs)
        
    # Sort the arrays to create a chronology
    t_events = np.sort(t_events)
    t_events -= t_events[0]
    t_events += random.random()*delta

    return t_events

# = = = = = = = = = = = 

# Define PDF - Power Law
#t_min_pw, t_max_pw = 1.0e7, 1.0e10
#pw_index = -1.0
#kwargs = {"t_min":t_min_pw, "t_max":t_max_pw, 
#          "pw_index":pw_index, "norm":(1.0 / t_min_pw**pw_index)}

t_events = []
start_time = t_module.time()
for i in range(nb_run):
#    t_events.append(get_list_t_events(gamma, int(maxTime/gamma), sample_power_law, **kwargs))
    t_events.append(get_list_t_events(gamma, int(maxTime/gamma), sample_square_box, **kwargs))
print(get_time(start_time))
print("nb_events per line:",(int(maxTime/gamma)))
print("nb_events in total:",(int(maxTime/gamma))*nb_run)

# Open the file for writting all lists of t_events
tEvents_file = open(output_name, "w")

# Write the header
# Number of events, Number of runs
tEvents_file.write(str(int(maxTime/gamma))+' '+str(nb_run)+'\n')

# Write each run
for i in range(nb_run):
    for number in t_events[i]:
        tEvents_file.write(str('%.8E' %number)+' ')
    tEvents_file.write('\n')
tEvents_file.close()



