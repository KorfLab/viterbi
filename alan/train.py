import sys
from readfasta import read_record
from itertools import product

def make_kmer(k):
	kmers = {}
	for kmer in product(list('ACGT'), repeat = k):
		kmer = ''.join(kmer)
		kmers[kmer] = 0
	return kmers

def train(fasta, k):
	numseq = 0
	totlen = 0
	kmers = make_kmer(k)
	for idn, seq in read_record(fasta):
		for i in range(len(seq)-k+1):
			kmer = seq[i:i+k]
			if kmer in kmers: kmers[kmer]+=1
		numseq += 1
		totlen += len(seq)
	if 0 in kmers.values():
		sys.exit('Error: 1 or more kmers contain 0 count!')
	avglen = totlen / numseq
	tot_kmers = sum(kmers.values())
	for kmer in kmers: kmers[kmer] /= tot_kmers
	
	return avglen, kmers
	
efile = sys.argv[1]
ifile = sys.argv[2]
k = int(sys.argv[3])

# Exon
elen, ekmers = train(efile, k)
ei = 1 / elen
ee = 1 - ei

# Intron
ilen, ikmers = train(ifile, k)
ie = 1 / ilen
ii = 1 - ie

with open('data/exon.mm', 'w') as fh:
	fh.write('# exon\n')
	count = 0
	for kmer in ekmers:
		fh.write(f'{kmer} {ekmers[kmer]}\n')
		count+=1
		if count%k == 0: fh.write('\n')

with open('data/intron.mm', 'w') as fh:
	fh.write('# intron\n')
	count = 0
	for kmer in ikmers:
		fh.write(f'{kmer} {ikmers[kmer]}\n')
		count+=1
		if count%k == 0: fh.write('\n')

with open('data/transition.prob', 'w') as fh:
	fh.write('# transition_states exon=0 intron=1\n')
	fh.write(f'{ee} {ei}\n')
	fh.write(f'{ie} {ii}')
