Amir Bayegan, Peter Clote -- Boston Colleg
************************************************************************
This program simulates kinetics of folding with event-driven Monter Carlo or Gillespie algorithms. The available options are as follows:
************************************************************************
USAGE: ./mc [-s sequence]|[-f FastaFile] [-m method] [-y random seed] [-r number of runs] [-l trajectory length] [-d 0|1|2] [-t temperature(c)] [-e 99|04] [-o output file prefix ] [-v]
try ./mc -h for help

The available options are as follows:

-s : input RNA sequence
-f : input fasta file containing the RNA. This flag and '-s' are mutually exclusive
-m : Define which method should be used to produce trajectories. Use "emc" or "gil" for event-driven monte carlo or gillespie respecively
-y : fix the random seed to the given value. This allows to compare the results from the methods
-r : define the number of trajectories for the given RNA. The default is 1.
-l : define the number of states after which trajectory stops. The default is stopping at the MFE structure
-d : set the dangles to be 0,1 or 2 to compute the energy of each structures
-t : set the temprature (C) to compute the energy of structures
-e : set the energy model to be Turner 99, Turner 2004 or Andronescu 2007  by entering 99|04|07 respectively
-o : set the output folder prefix name to save the mean first passage times. The default prefix is 'rna'. The prefix will be followed by '_firstPassageTimes'
-v : print the trajectory on the output. Save the printed output if you need the trajectories. The output fields are (1)the state,(2)energy of state, (3)time stayed in the state, (4)number of neighbors of the state
