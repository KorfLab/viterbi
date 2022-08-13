import math
import sys
import random

def read_coordinates(coorfile, genlen):
	coor = [1 for _ in range(genlen)]
	with open(coorfile) as fh:
		for line in fh.readlines():
			f = line.split()
			beg = int(f[1])-1
			end = int(f[2])-1
			for i in range(beg, end+1):
				coor[i] = 0
	return coor


def evgenome(tb, k):
	features = []
	cur_state = tb[0]
	beg = 1
	end = None
	for i, f in enumerate(tb):
		if f != cur_state:
			end = i + k - 1
			features.append((beg, end, cur_state))
			cur_state = f
			beg = end+1
	
	features.append((beg,len(tb),cur_state))
	return features
	
def evfeat(tb, feature):
	count = sum(tb)
	if feature == 'exon':
		if 1 - count/len(tb) > 0.5: return True
		return False
	if feature == 'intron':
		if count/len(tb) > 0.5: return True
		return False
		
def seq2int(seq):
	seq = list(seq)
	for i in range(len(seq)):
		if seq[i] == 'A' or seq[i] == 'a' : seq[i] = 0
		if seq[i] == 'C' or seq[i] == 'c' : seq[i] = 1
		if seq[i] == 'G' or seq[i] == 'g' : seq[i] = 2
		if seq[i] == 'T' or seq[i] == 't' : seq[i] = 3
		if seq[i] == 'N' or seq[i] == 'n' : seq[i] = random.randint(0,3)
	return seq

def ntidx(kmer, k):
	idx = 0
	for i in range(k):
		nt = kmer[i]
		if nt == 'A' or nt == 'a': idx += pow(4, k-i-1) * 0
		if nt == 'C' or nt == 'c': idx += pow(4, k-i-1) * 1
		if nt == 'G' or nt == 'g': idx += pow(4, k-i-1) * 2
		if nt == 'T' or nt == 't': idx += pow(4, k-i-1) * 3
	return idx

def numidx(kmer, k):
	idx = 0
	for i in range(k): idx += pow(4, k-i-1) * kmer[i]
	return idx
	
def read_mm(mmf, feat):
	kmers = None
	k = None
	with open(mmf) as fh:
		for line in fh.readlines():
			line = line.rstrip()
			if line == '': continue
			f = line.split()
			if f[0] == '#' and f[1] != feat:
				sys.exit('Error: exon mm file error!')
			if f[0] != '#' and len(f) > 1:
				kmer = f[0]
				if k == None:
					k = len(kmer)
					kmers = ([0 for _ in range(pow(4,k))])
				prob = float(f[1])
				kmers[ntidx(kmer, k)] = prob
	
	return kmers, k

def read_transition(tpf):
	tp = [[0,0],[0,0]]
	count = 0
	with open(tpf) as fh:
		for line in fh.readlines():
			line = line.rstrip()
			f = line.split(' ')
			if f[0] == '#' and f[1] != 'transition_states':
				sys.exit('Error: transition state file error!')
			if len(f) == 2:
				tp[count][0] = float(f[0])
				tp[count][1] = float(f[1])
				count+=1
	return tp
	
def print_mat(mat, vert = False):
	if vert:
		row_num = len(mat)
		col_num = len(mat[0])
		for i in range(col_num):
			for j in range(row_num):
				val = mat[j][i]
				print(f'{val:.4f}\t', end= '')
			print('')
	else:
		for row in mat:
			for val in row:
				print(f'{val:.3f}\t', end= '')
			print('')


def logodds(mat):
	for row in range(len(mat)):
		for val in range(len(mat[0])):
			mat[row][val] = math.log(mat[row][val])
	return mat
