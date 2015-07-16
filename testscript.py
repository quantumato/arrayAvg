import os

procs = [4, 9, 25, 64]
params = [10, 50, 100, 1000]

os.system("echo starting arrayAvg4")
os.system("echo ")

for i in procs:
	os.system("echo number of processes: {}".format(i))
	os.system("echo ")
	for j in params:
		os.system("echo args: {}".format(j))
		os.system("echo ")
		for k in range(5):
			os.system("mpirun -np {} arrayAvg4 {} {}".format(i, j, j))

os.system("echo starting arrayAvg5\n")
os.system("echo ")

for i in procs:
	os.system("echo number of processes: {}".format(i))
	os.system("echo ")
	for j in params:
		os.system("echo args: {}".format(j))
		os.system("echo ")
		for k in range(5):
			os.system("mpirun -np {} arrayAvg5 {} {}".format(i, j, j))

os.system("echo starting arrayAvg6\n")
os.system("echo ")
		
for i in procs:
	os.system("echo number of processes: {}".format(i))
	os.system("echo ")
	for j in params:
		os.system("echo args: {}".format(j))
		os.system("echo ")
		for k in range(5):
				os.system("mpirun -np {} arrayAvg6 {} {}".format(i, j, j))
