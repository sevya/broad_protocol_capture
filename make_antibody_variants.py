#! /usr/bin/env python2.7

from Bio import AlignIO, SeqIO
import random, sys

aln_indices = range(385, 392) + range(489, 490) + range(492, 498) + range(575, 581) + range(605, 610) + range(618, 621) + range(648, 652)

aminoacids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

if len(sys.argv) < 3:
	print 'Usage: python2.7 make_antibody_variants.py [mutations per ab] [# total abs]\nEx. python2.7 make_antibody_variants.py 5 500'
	exit()

## Make random antibody variants and output to FASTA file
num_muts = int(sys.argv[1])
num_variants = int(sys.argv[2])
num_digits = len(str(num_variants))

native_seq = 'EWMGWIKPERGAVSYAPQRDLYRDASW'
with open('antibody_variants.fasta', 'w') as out:
	for ii in range(num_variants):
		new_seq = [aa for aa in native_seq]
		for jj in range(num_muts):
			position = random.choice(range(len(new_seq)))
			newaa = random.choice(aminoacids)
			new_seq[position] = newaa
		new_seq = ''.join(new_seq)
		variant_name = 'ab_%s'%(str(ii+1).zfill(num_digits))
		out.write( ''.join(['>', variant_name, '\n', new_seq, '\n']) )

# Convert antibody variants fasta file into resfiles
threading_indices = range(45, 63) + range(71, 75) + range(101, 106)
for handle in SeqIO.parse('antibody_variants.fasta', 'fasta'):
	name = str(handle.id)
	sequence = str(handle.seq)
	with open('ab_resfiles/'+name+'.resfile', 'w') as out:
		out.write('NATAA\nstart\n')
		for thread_seqpos, aa in zip(threading_indices, sequence):
			out.write(' '.join([str(thread_seqpos), 'H', 'PIKAA', aa, '\n']))	
