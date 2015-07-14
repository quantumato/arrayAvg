import os

rank = int(os.environ["OMPI_COMM_WORLD_RANK"])

if rank == 0:
	os.system("gdb mpirun arrayAvg6")
else:
	os.system("mpirun -np 16 arrayAvg6")
