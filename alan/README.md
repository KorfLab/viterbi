# Evaluation

Evalution for my viterbi algorithm

## Preprocess sequence

Lowercase bases are converted to their corresponding uppercase bases

Ns are randomly converted to A, C, G, or T

## Training and Testing

Up to kmers with k=4 are trained and tested. Training data contains 1 or more kmers with 0 count when k=5.

### Features (Distinguishing exon from intron, vice versa)

Training set: first half of exons.fa and introns.fa

Testing set: second half of exons.fa and introns.fa

### Genome (Indentify exon)

Training set: all of exons.fa and introns.fa

Testing data: worms.fa

Sequence preprocessing: random seed = 1

## Results

### Features

A sequence is determined to be predicted successfully if more than 70% of its bases are predicted to be the correct feature.

Exon:

|  k  | success rate |
|:---:|:------------:|
|  1  |    0.9451    |
|  2  |    0.9333    |
|  3  |    0.9333    |
|  4  |    0.9176    |

Intron:

|  k  | success rate |
|:---:|:------------:|
|  1  |    0.7079    |
|  2  |    0.8218    |
|  3  |    0.8465    |
|  4  |    0.8663    |

### Mini worm genome

Success rate is determined from the number of bases in overlapped features divided by total number of bases

|  k  | success rate |
|:---:|:------------:|
|  1  |    0.6583    |
|  2  |    0.7421    |
|  3  |    0.7748    |
|  4  |    0.7798    |
