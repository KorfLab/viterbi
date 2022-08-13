import math
import sys

def get_idx(kmer, k):
	idx = 0
	for i in range(k):
		nt = kmer[i]
		if nt == 'A' or nt == 'a': idx += pow(4, k-i-1) * 0
		if nt == 'C' or nt == 'c': idx += pow(4, k-i-1) * 1
		if nt == 'G' or nt == 'g': idx += pow(4, k-i-1) * 2
		if nt == 'T' or nt == 't': idx += pow(4, k-i-1) * 3
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
				kmers[get_idx(kmer, k)] = prob
	
	return kmers

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
