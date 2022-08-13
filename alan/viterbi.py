from readfasta import read_record
from viterbi_tools import *
import random

# Transition probabilities
# Exon = 0	Intron = 1
# [0][0] -> exon to exon
# [0][1] -> exon to intron
# [1][0] -> intron to exon
# [1][1] -> intron to intron
tpf = sys.argv[1]
tp = read_transition(tpf)
# State probabilities
# Exon = 0 Intron = 1
# kmer represented in int
exspf = sys.argv[2]
inspf = sys.argv[3]
exsp = read_mm(exspf, 'exon')
insp = read_mm(inspf, 'intron')
sp = [exsp, insp]
# Probabilities to logodds
tp = logodds(tp)
sp = logodds(sp)
# Get seq
letter = {}
seq = None
for idn, s in read_record('worm.fa'):
	seq = s
	break
# Process seq and turn into int value
seq = seq[0:1000] # subseq for developing
seq = list(seq)
for i in range(len(seq)):
	if seq[i] == 'A' or seq[i] == 'a' : seq[i] = 0
	if seq[i] == 'C' or seq[i] == 'c' : seq[i] = 1
	if seq[i] == 'G' or seq[i] == 'g' : seq[i] = 2
	if seq[i] == 'T' or seq[i] == 't' : seq[i] = 3
	if seq[i] == 'N' or seq[i] == 'n' : seq[i] = random.randint(0,3)

# Initialize matrices 
pmat = [[0]*(len(seq)+1) for _ in range(len(sp))]
tmat = [[-1]*(len(seq)+1) for _ in range(len(sp))]

# Initialize state probabilities
pmat[0][0] = math.log(0.67) # initial prob for exon
pmat[1][0] = math.log(0.3)  # initial prob for intron

# Fill pmat
for i in range(1, len(seq)+1): # Iterate over each base
	for j in range(len(sp)): # Iterate over each state against current base
		cur_seq = seq[i-1]
		max_score = None
		max_state = None
		for s in range(len(sp)):
			# s=0 -> exon ; s=1 -> intron
			# score = previous state p * transition state p * state p
			score = pmat[s][i-1] + tp[s][j] + sp[j][cur_seq]
			if max_score == None or score > max_score:
				max_score = score
				max_state = s
		pmat[j][i] = max_score
		tmat[j][i] = max_state

# Trace back
max_state = None
max_score = None
for i in range(len(sp)):
	score = pmat[i][len(seq)]
	if max_score == None or score > max_score:
		max_score = score
		max_state = i

tb = []
curi = max_state
curj = len(seq)
while tmat[curi][curj] != -1:
	tb.insert(0, tmat[curi][curj])
	curj-=1
	curi=tmat[curi][curj]

#print_mat(pmat)
#print_mat(pmat, vert=True)
#print_mat(tmat)
#print(tb)
