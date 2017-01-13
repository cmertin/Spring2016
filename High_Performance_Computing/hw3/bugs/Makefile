CC = g++
MPI = mpic++
CFLAGS = -O2 -fopenmp -std=c++11

bug1: mpi_solved1.c
	$(MPI) $(CFLAGS) -o bug1 mpi_solved1.c

bug2: mpi_solved2.c
	$(MPI) $(CFLAGS) -o bug2 mpi_solved2.c

bug3: mpi_solved3.c
	$(MPI) $(CFLAGS) -o bug3 mpi_solved3.c

bug4: mpi_solved4.c
	$(MPI) $(CFLAGS) -o bug4 mpi_solved4.c

bug5: mpi_solved5.c
	$(MPI) $(CFLAGS) -o bug5 mpi_solved5.c

bug6: mpi_solved6.c
	$(MPI) $(CFLAGS) -o bug6 mpi_solved6.c

bug7: mpi_solved7.c
	$(MPI) $(CFLAGS) -o bug7 mpi_solved7.c

clean:
	rm -f bug* *~ \#* 
