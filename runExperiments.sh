# Get the filename
filName=$1

# Now copy it to tEvents.in
echo "Copying tEvents"
cp "database/"$filName".in" "tEvents.in"

# Now run the experiment, change the number of CPU if desired
echo "Running experiment"
mpirun -np 10 ./radioCalc

# Now combine outputs
echo "Combining outputs"
python combine_outputs.py

# Removing redundant outputs
echo "Removing the Output*.txt files"
rm Output*.txt

# Renaming the tau files
echo "Renaming files"
rename 's/tau_/'$filName'_tau_/' tau_*.txt
