#! /usr/bin/env python2.7

from Bio import SeqIO
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import sys

if len(sys.argv) < 2:
	print 'Usage: python2.7 compile_results.py [filename]'
	exit()

ab_names = dict([(str(item.id), str(item.seq)) for item in SeqIO.parse('antibody_mutants_test.fasta', 'fasta')])
virus_names = [line.strip() for line in open('virus_list.txt')]


## Get ab-virus pair names from the description
def get_key( line ):
	name = line.split()[-1]
	return '_'.join(name.split('_')[1:3])

design_data = pd.read_table( sys.argv[1], sep=r"\s*", header=1, skiprows=0, engine='python' )
pool = mp.Pool( processes=10 )
ab_virus_pairs = pool.map( get_key, design_data['description'] )
design_data['ab_virus_pairs'] = ab_virus_pairs


## Set up a fxn to get the ddG of the lowest scoring model for a given key
def get_ddg_lowest_scoring( key, df, top_n=1 ):
	subset_df = df.loc[df['ab_virus_pairs'] == key]
	subset_df = subset_df.sort_values(by='total_score')
	no_models = len(list(subset_df['ddG']))
	ddg = np.mean(list(subset_df['ddG'])[:top_n])
	ab, virus = key.split('_')
	return (ab, virus, float(ddg), no_models)

## Set up a fxn to get the score of the lowest scoring models for a given key
def get_score_lowest_scoring( key, df, top_n=5 ):
	subset_df = df.loc[df['ab_virus_pairs'] == key]
	subset_df = subset_df.sort_values(by='total_score')
	no_models = len(list(subset_df['total_score']))
	score = np.mean(list(subset_df['total_score'])[:top_n])
	ab, virus = key.split('_')
	return (ab, virus, float(score), no_models)

native_ddg_results = pool.map(
	partial(get_ddg_lowest_scoring, df=data),
	map(lambda item: 'ab-native_'+item, virus_names) 
	)
native_ddg_results.sort(key=lambda arr: arr[1])

native_score_results = pool.map(
	partial(get_score_lowest_scoring, df=data),
	map(lambda item: 'ab-native_'+item, virus_names) 
	)
native_score_results.sort(key=lambda arr: arr[1])

for key in ab_names.keys():
	design_ddg_results = pool.map(
		partial(get_ddg_lowest_scoring, df=data),
		map(lambda item: key+'_'+item, virus_names) 
		)
	design_ddg_results.sort(key=lambda arr: arr[1])

	design_score_results = pool.map(
		partial(get_score_lowest_scoring, df=data),
		map(lambda item: key+'_'+item, virus_names) 
		)
	design_score_results.sort(key=lambda arr: arr[1])

	cutoff = -28.5
	print 'Cutoff of %.1f'%cutoff
	native_breadth 	= np.mean([c<cutoff for a,b,c,d in native_results])*100
	design_breadth	= np.mean([c<cutoff for a,b,c,d in design_results])*100

	print '%0.1f %0.1f '%(native_breadth,design_breadth)
	print '============================'

	plt.plot(
		[c for a,b,c,d in native_results],
		[c for a,b,c,d in design_results],
		'bo')
	plt.xlabel('Native DDG')
	plt.ylabel(key+' DDG')
	plt.title(key+' design results')
	xmin,xmax,ymin,ymax = plt.axis()
	axmin = min(xmin,ymin)
	axmax = max(xmax,ymax)
	plt.axis((axmin,axmax,axmin,axmax))
	plt.plot(np.arange(axmin,axmax+1),
			 np.arange(axmin,axmax+1),
			 'k--', alpha=0.5)
	plt.show()

	plt.plot(
		[c for a,b,c,d in native_score_results],
		[c for a,b,c,d in design_score_results],
		'bo')
	plt.xlabel('Native Score')
	plt.ylabel(key+' Score')
	plt.title(key+' design results')
	xmin,xmax,ymin,ymax = plt.axis()
	axmin = min(xmin,ymin)
	axmax = max(xmax,ymax)
	plt.axis((axmin,axmax,axmin,axmax))
	plt.plot(np.arange(axmin,axmax+1),
			 np.arange(axmin,axmax+1),
			 'k--', alpha=0.5)
	plt.show()

	

	# for nativeDDG, designDDG, nativeScore, designScore in zip(native_results, design_results, native_score_results, design_score_results):
	# 	deltaDDG = designDDG[2] - nativeDDG[2]
	# 	deltaScore = designScore[2] - nativeScore[2]
	# 	outfile.write(str(deltaDDG)+'\t'+str(deltaScore)+'\n')
