Train Simulation
================
Files: `train.cpp` `Dimensions.h` `plot.py` `run.py` `run_all.sh`
This code simulates the flow around a train. This inmplements periodic boundary conditions along the top and bottom boundaries, and appropriate boundary conditions at the inlet and outlet. The given file `run.py` runs the test case and asks the user to implement parameters. The file `run_all.sh` is a `bash` script file that will iterate over the separations that were used in the paper, and will send an email to a user-specified email address after finishing each velocity. Both `run.py` and `run_all.sh` will check to see if the executable exists before running. If it doesn't, it will make it. If the user wants to run it by themselves, they will have to make the executable first by running the `make` command with the given `Makefile`.

When running `run.py` it asks the user for the following input
1. Velocity: Integer velocity that you want the train to go at
2. Index of the domain grid where you want the first train to start at in X
3. Index of the domain grid where you want both trains to be at in Y
4. Number of units/discretizations you want for train width
5. Number of units/discretizations you want for train height
6. Number of unitz/discretizations for the separation between the two train cars

This code outputs a file based on the following filename format: `train_U_width_height_separation` where each of these were user defined. It also prints out comments (lines starting with `#`) in the file which contain all the parameters that were in the given simulation, including the system start time, system end time, and the total runtime.

After the pressure converges on the domain, the program runs `plot.py` to output the plots for the vorticity and streamfunction, the pressure and streamfunction, the x-velocity, and the y-velocity into the current directory. These are output in pdf format and come from utilizing `streamfunction` and `quiver` plots in `matplotlib`.