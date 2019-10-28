# Python packages
import matplotlib.pyplot as plt
import os.path
import numpy as np
import glob

def writeArray(fwrite, array):
    '''Write array to the file fwrite'''
    
    s = ["{:E}".format(x) for x in array]
    s = " ".join(s)
    s += "\n"
    
    fwrite.write(s)

def main():
    '''Plot spreads, sigmas and formula'''
    
    # Lower and upper limits for delta statistics
    lowAge, highAge = 1e10, 1.4e10
    
    # Get list of gammas and taus
    taus = ["1.00e9"]
    
    # Initialize variables for the reading
    oldFilNam = None
    allFiles = sorted(glob.glob("npyFiles/*.npy"))
    for tauStr in taus:
        for fil in allFiles:
            # Get the file for this tau
            if not fil.endswith(tauStr + ".npy"):
                continue
            
            tau = float(tauStr)
            runs = np.load(fil)
            
            # Get gamma
            lnlst = fil.split("_")
            indxGamm = lnlst.index("gamma")
            gammVal = float(lnlst[indxGamm + 1])
            
            # Get the statistics for the deltas but only if it is a new file
            filNam = fil[fil.find("box"):fil.find("runs") + 4] + ".in"
            filNam = "../database/tEvents_" + filNam
            if oldFilNam is None or filNam != oldFilNam:
                allDeltas = []
                with open(filNam, "r") as fread:
                    # Skip first line
                    fread.readline()
                    
                    # Now get all deltas between lowAge and highAge
                    for line in fread:
                        tEvents = [float(x) for x in line.split()]
                        
                        # Get the deltas between lowAge and highAge
                        deltas = []
                        for ii in range(len(tEvents) - 1):
                            if tEvents[ii] >= lowAge and tEvents[ii] <= highAge:
                                deltas.append(tEvents[ii + 1] - tEvents[ii])
                            elif tEvents[ii] > highAge:
                                break
                        
                        allDeltas += deltas
                    
                    allDeltas = np.array(allDeltas)
            
            oldFilNam = filNam
            
            # Get and plot stats
            combinedRuns = []
            for run in runs:
                combinedRuns += list(run)
            
            combinedRuns = np.array(combinedRuns)
            deltMean = np.mean(allDeltas)
            deltStd = np.std(allDeltas)
            mean = np.mean(combinedRuns)
            median = np.median(combinedRuns)
            spread1 = np.percentile(combinedRuns, 50 - 34.1)
            spread2 = np.percentile(combinedRuns, 50 + 34.1)
            spread = spread2 - spread1
            std = np.std(combinedRuns)
            
            print("")
            print(fil[:fil.index("gamma") - 1])
            print("tau = {}, gamm = {:e}, tau/gamma = {}".format(tauStr, gammVal, tau/gammVal))
            print("Median = {}; median/mean = {}".format(median, median/mean))
            print("Mean = {}; spread = {}; std = {}".format(mean, spread/mean*100, std/mean*100))
            print("spread- = {}; spread+ = {}".format(spread1 - mean, spread2 - mean))

if __name__ == "__main__":
    main()
