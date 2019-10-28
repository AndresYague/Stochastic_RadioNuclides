# Python packages
import matplotlib.pyplot as plt
import os.path
import numpy as np
import glob

def cutRuns(tArray, runs, minTime = None, maxTime = None):
    '''Cut down the runs and combine them for the statistics'''
    
    if minTime is None:
        minTime = 0
    if maxTime is None:
        maxTime = tArray[-1]*10
    
    # Look for the indices
    indxMin, indxMax = 0, 0
    
    # Up
    ii = 0
    while ii < len(tArray) - 1:
        if tArray[ii] > minTime:
            indxMin = ii
            break
        
        ii += 1
    
    # Down
    ii = len(tArray)
    while ii > 0:
        if tArray[ii - 1] < maxTime:
            indxMax = ii
            break
        
        ii -= 1
    
    # Cut and combine the runs
    tArray = tArray[indxMin:indxMax]
    cutRuns = []
    for run in runs:
        cutRuns.append(run[indxMin:indxMax])
    
    return tArray, cutRuns

def getRuns(filName):
    '''Retrieve all the runs in filName as a list, in order'''
    
    li = []
    with open(filName, "r") as fread:
        tArray = [float(x) for x in fread.readline().split()]
        for line in fread:
            li.append([float(x) for x in line.split()])
    
    # Cut and combine the runs
    tArray, li = cutRuns(tArray, li, 1e10, 1e14)
    
    return np.array(tArray), np.array(li)

def getAllRuns():
    '''Get and save all the cut runs'''
    
    allFiles = sorted(glob.glob("output_data/*.txt"))
    lenAllFiles = len(allFiles)
    while len(allFiles) > 0:
        fil = allFiles.pop()
        
        # Read the file and combine the runs
        tArray, runs = getRuns(fil)
        
        # Save this file
        fileName = os.path.split(fil)[1][:-4]
        np.save("npyFiles/"+fileName, runs)
        
        # Print progress
        print("Done {}%".format(round(100*(lenAllFiles - len(allFiles))/lenAllFiles)))

def plots():
    # Define minimum
    minimum = 1e-100

    # Select tau1 and tau2
    tau1 = '1.00e6'
    tau2 = '1.00e8'

    # Get the list of files with those taus
    allFiles = sorted(glob.glob("output_data/*.txt"))
    while len(allFiles) > 0:
        fil1 = allFiles.pop()
        if not fil1.endswith(tau1 + ".txt"):
            continue
        
        # Filename found, create the other pair
        indx = fil1.index("tau")
        fil2 = fil1[:indx] + "tau_" + tau2 + ".txt"
        if not os.path.isfile(fil2):
            continue
        
        # Now read each file, combine the runs and get the statistics
        tArray, runs1 = getRuns(fil1)
        tArray, runs2 = getRuns(fil2)
        
        print(fil1)
        print(fil2)
        
        # Combine the runs
        combinedRuns = []
        for run1, run2 in zip(runs1, runs2):
            ratio = run1/(run2 + minimum)
            combinedRuns += list(ratio)
            
            plt.plot(tArray, ratio)
        
        plt.show()
        
        plt.hist(combinedRuns, bins = 1000)
        plt.show()

if __name__ == "__main__":
    '''Run plots or get all runs'''
    
    #plots()
    
    # Get and save all runs
    print("Getting all the runs...")
    getAllRuns()
    print("All runs saved!")
