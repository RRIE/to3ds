def sba_parse(test_txt, results_file):
	ftest = open(test_txt)
	fresults = open(results_file, "w+")
	lines = ftest.readlines()
	derivs = []
	A = []
	U = []
	V = []
	damping = []
	S = []
	system = []
	db_i = []
	fn = []
	for i in range(len(lines)):
		if lines[i].split("\t")[0] == "derivs":
			derivs.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "A":
			A.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "U":
			U.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "V":
			V.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "damping":
			damping.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "S":
			S.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "system":
			system.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "db_i's":
			db_i.append(lines[i].split("\t")[2].split("\n")[0])
		elif lines[i].split("\t")[0] == "fn":
			fn.append(lines[i].split("\t")[2].split("\n")[0])
	
	fresults.write("Computing Derivatives\t")
	for i in range(len(derivs)):
		fresults.write(derivs[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Computing A and B\t")
	for i in range(len(A)):
		fresults.write(A[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Computing U\t")
	for i in range(len(U)):
		fresults.write(U[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Computing V\t")
	for i in range(len(V)):
		fresults.write(V[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Computing damping factor\t")
	for i in range(len(damping)):
		fresults.write(damping[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Computing S and e_j\t")
	for i in range(len(S)):
		fresults.write(S[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Solving Linear System\t")
	for i in range(len(system)):
		fresults.write(system[i]+"\t")
	fresults.write("\n")
	
	fresults.write("Computing db_i\t")
	for i in range(len(db_i)):
		fresults.write(db_i[i]+"\t")
	fresults.write("\n")
	
	fresults.write("fn eval\t")
	for i in range(len(fn)):
		fresults.write(fn[i]+"\t")
	fresults.write("\n")
	
	fresults.close()
	ftest.close()