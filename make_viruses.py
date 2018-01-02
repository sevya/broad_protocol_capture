#! /usr/bin/env python2.7

aln_indices = range(295, 303) + range(384, 391) + range(445, 451) + range(475, 483) + range(498, 502)
aln_indices = range(385, 392) + range(489, 490) + range(492, 498) + range(575, 581) + range(605, 610) + range(618, 621) + range(648, 652)

aminoacids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

from Bio import AlignIO
import random

aln_dict = {}

handle = AlignIO.parse('181_viruses_plus_4J6R.aln', 'clustal')
for record in handle:
	for item in record:
		aln_dict[str(item.id)] = str(item.seq)

threading_indices = range(160, 168) + range(222, 229) + range(276, 282) + range(306, 314) + range(324, 328)

def write_resfile():
	for name, sequence in aln_dict.items():
		with open('resfiles/'+name+'.resfile', 'w') as out:
			out.write('NATAA\nstart\n')
			for aln_seqpos, thread_seqpos in zip(aln_indices, threading_indices):
				aa = sequence[aln_seqpos-1]
				out.write(' '.join([str(thread_seqpos), 'G', 'PIKAA', aa, '\n']))

sequence_dict = {}
for name, sequence in aln_dict.items():
	myseq = ''
	for aln_seqpos, thread_seqpos in zip(aln_indices, threading_indices):
		aa = sequence[aln_seqpos-1]
		if aa == '-':
			aa = aln_dict['4J6R'][aln_seqpos-1]
		myseq += aa
	sequence_dict[name] = myseq
#	print name.ljust(10), myseq.ljust(40)

sequence_dict.pop('4J6R')
keys = sequence_dict.keys()

num_muts = 5
num_variants = 3

# Make 40 variants per virus, each variant with 5 point mutations
def make_viral_variants_fasta():
	with open('viral_variants.fasta', 'w') as out:
		for virus in keys:
			base = sequence_dict[virus]
			virus = virus.replace('.', '-').replace('_', '-')
			out.write( ''.join(['>', virus, '\n', base, '\n']) )
			for ii in range(num_variants):
				new_seq = [aa for aa in base]
				for jj in range(num_muts):
					position = random.choice(range(len(new_seq)))
					newaa = random.choice(aminoacids)
					new_seq[position] = newaa
				new_seq = ''.join(new_seq)
				variant_name = virus+'-v%02d'%ii
				out.write( ''.join(['>', variant_name, '\n', new_seq, '\n']) )	

# Convert viral variants fasta file into resfiles
from Bio import SeqIO
def make_viral_variants_resfile():
	for handle in SeqIO.parse('viral_variants.fasta', 'fasta'):
		name = str(handle.id)
		sequence = str(handle.seq)
		with open('resfiles/'+name+'.resfile', 'w') as out:
			out.write('NATAA\nstart\n')
			for thread_seqpos, aa in zip(threading_indices, sequence):
				out.write(' '.join([str(thread_seqpos), 'G', 'PIKAA', aa, '\n']))
# make_viral_variants_resfile()

def make_antibody_variants_fasta():
	num_muts = 5
	num_variants = 500
	base = 'EWMGWIKPERGAVSYAPQRDLYRDASW'
	with open('antibody_variants.fasta', 'w') as out:
		for ii in range(num_variants):
			new_seq = [aa for aa in base]
			for jj in range(num_muts):
				position = random.choice(range(len(new_seq)))
				newaa = random.choice(aminoacids)
				new_seq[position] = newaa
			new_seq = ''.join(new_seq)
			variant_name = 'ab_%02d'%(ii+1)
			out.write( ''.join(['>', variant_name, '\n', new_seq, '\n']) )

# Convert antibody variants fasta file into resfiles
def make_antibody_variants_resfile():
	threading_indices = range(45, 63) + range(71, 75) + range(101, 106)
	for handle in SeqIO.parse('antibody_variants.fasta', 'fasta'):
		name = str(handle.id)
		sequence = str(handle.seq)
		with open('resfiles/'+name+'.resfile', 'w') as out:
			out.write('NATAA\nstart\n')
			for thread_seqpos, aa in zip(threading_indices, sequence):
				out.write(' '.join([str(thread_seqpos), 'H', 'PIKAA', aa, '\n']))	

make_viral_variants_fasta()
make_viral_variants_resfile()
exit()

threading_dict = {}
for key, value in aln_dict.items():
	threading_dict[key] = ''.join([value[index-1] for index in aln_indices])

## Find my unique viruses
virus_keys = [key for key in threading_dict if 'renum' not in key]
virus_seqs = [threading_dict[key] for key in virus_keys]
unique_viral_seqs = list(set(virus_seqs))
unique_virus_keys = []
for seq in unique_viral_seqs:
	index = virus_seqs.index(seq)
	unique_virus_keys.append(virus_keys[index])

## Find threading positions based on the sequence alignment
threading_indices_dict = {}
for key in aln_dict:
	if 'renum' in key:
		threading_indices = []
		for index in aln_indices:
			dash_sequence = aln_dict[key][:index]
			nodash_sequence = aln_dict[key][:index].replace('-', '')
			num_dashes = len(dash_sequence) - len(nodash_sequence)
			threading_indices.append( index - num_dashes )
		threading_indices_dict[key] = threading_indices
		# print key, ' '.join([str(i) for i in threading_indices])
## Make pymol-style selection
for structure_key in threading_indices_dict:
	print 'cmd.select(\''+structure_key+'_bind_site\', \''+structure_key+' and chain G and resi', '+'.join([str(i) for i in threading_indices_dict[structure_key]])+'\')'
	print 'cmd.color(\'red\', \''+structure_key+'_bind_site\')'
exit()

for structure_key in threading_indices_dict:
	for vir_key in unique_virus_keys:
		if 'renum' in vir_key:
			continue
		# print structure_key, vir_key
		structure_name = structure_key.split('_')[0]
		with open('resfiles/'+structure_name+'_'+vir_key+'.resfile', 'w') as out:
			out.write('NATAA\nstart\n')
			for seqpos, aa in zip(threading_indices_dict[structure_key], threading_dict[vir_key]):
				out.write(' '.join([str(seqpos), 'G', 'PIKAA', aa, '\n']))
