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
    #tArray, li = cutRuns(tArray, li, 1e10, 1e14)
    tArray, li = cutRuns(tArray, li, 0, 1e14)
    
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
        np.save(fil, runs)
        
        # Print progress
        print("Done {}%".format(100*(lenAllFiles - len(allFiles))/lenAllFiles))

def plots():
    # Define minimum
    minimum = 1e-100

    # Select tau
    tau = '1.00e9'

    # Get the list of files with those taus
    allFiles = sorted(glob.glob("output_data/*.txt"))
    while len(allFiles) > 0:
        fil = allFiles.pop()
        #if not "3.16e8" in fil:
            #continue
        if not fil.endswith(tau + ".txt"):
            continue
        
        # Now read each file, combine the runs and get the statistics
        tArray, runs = getRuns(fil)
        
        print(fil)
        
        # Combine the runs
        combinedRuns = []
        parallelRuns = []
        gray = (0.5, 0.5, 0.5, 0.2)
        for run in runs:
            #combinedRuns += list(run)
            parallelRuns.append(run)
            
            plt.plot(tArray, run, color = gray)
        
        parallelRuns = np.transpose(parallelRuns)
        medianOfRuns = []
        for time in parallelRuns:
            medianOfRuns.append(np.median(time))
        
        plt.plot(tArray, medianOfRuns, "k")
        plt.plot(tArray, tArray*0 + 3.16, "r--")
        plt.ylim([2.2, 4])
        plt.show()
        
        #plt.hist(combinedRuns, bins = 1000)
        #plt.show()

if __name__ == "__main__":
    '''Run plots or get all runs'''
    
    plots()
    
    # Get and save all runs
    #print("Getting all the runs...")
    #getAllRuns()
    #print("All runs saved!")
