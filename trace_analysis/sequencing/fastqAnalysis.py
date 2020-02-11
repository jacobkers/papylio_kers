# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:59:54 2018

@author: Ivo Severins
"""

import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import re
import gc
from tabulate import tabulate
from pathlib import Path

class FastqData:
    def __init__(self, path):
        self.path = Path(path)
        self.write_path = self.path.parent
        self.adapter_sequences = ['AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'ATCTCGTATGCCGTCTTCTGCTTG']
        self._text_file = None

        file = open(path, 'r')

        lineCount = 0
        recordCount = 0
        
        instrument = list()
        run = list()
        flowcell = list()
        lane = list()
        tile = list()
        x = list()
        y = list()
        sample = list()
        
        sequence = list()
        quality = list()
        
        
        for line in file:
            if lineCount % 4 == 0: 
                splitHeader = re.split('[@: -]',line.strip())
                
                instrument.append(splitHeader[1])
                run.append(splitHeader[2])
                flowcell.append(splitHeader[4])
                lane.append(int(splitHeader[5]))
                tile.append(int(splitHeader[6]))
                x.append(int(splitHeader[7]))
                y.append(int(splitHeader[8]))
                sample.append(int(splitHeader[12]))
                
            if lineCount % 4 == 1: 
                sequence.append(line.strip())
            if lineCount % 4 == 3: 
                quality.append(line.strip())
                       
                recordCount +=1
        
            lineCount +=1
        
        
        file.close()
        
                
        self.instrument = np.array(instrument)
        self.run = np.array(run)
        self.flowcell = np.array(flowcell)
        self.lane = np.array(lane)
        self.tile = np.array(tile)
        self.x = np.array(x)
        self.y = np.array(y)
        self.sample = np.array(sample)
        self.sequence = np.array(sequence, dtype = bytes).view('S1').reshape((len(sequence),-1))
        self.quality = np.array(quality, dtype = bytes).view('S1').reshape((len(quality),-1))

        self.write_to_text_file('Total number of sequences: ' + str(self.sequence.shape[0]) + '\n\n')

    @property
    def tile_numbers(self):
        return np.unique(self.tile)

    @property
    def coordinates(self):
        return np.vstack([self.x, self.y]).T

    @property
    def text_file(self):
        if not self._text_file:
            self._text_file = self.write_path.joinpath('Output.txt')
            with self._text_file.open('w') as f:
                f.write('Analysis\n\n')
        return self._text_file

    def write_to_text_file(self, input_text):
        with self.text_file.open('a') as f:
            f.write(input_text + '\n')


    def select(self,indices,copyData = False):
        if copyData:
            import copy
            selection = copy.deepcopy(self)
        else:
            selection = self
        
        selection.instrument = self.instrument[indices]
        selection.run = self.run[indices]
        selection.flowcell = self.flowcell[indices]
        selection.lane = self.lane[indices]
        selection.tile = self.tile[indices]
        selection.x = self.x[indices]
        selection.y = self.y[indices]
        selection.sample = self.sample[indices]
        selection.sequence = self.sequence[indices,:]
        selection.quality = self.quality[indices,:]
        
        if copyData:
            return selection

    def selection(self, **kwargs):
        selection = np.zeros((self.sequence.shape[0], len(kwargs)), dtype=bool)
        for i, (key, value) in enumerate(kwargs.items()):
            if key == 'sequence':
                selection[:,i] = np.all(self.sequence[:, 0:len(value)] == np.array(list(value), dtype = bytes), axis = 1)
            else:
                selection[:,i] = getattr(self, key) == value

        return np.all(selection, axis = 1)

    def number_of_matches(self, sequence):
        # sequence must be a string
        return np.sum(self.sequence[:,0:len(sequence)]==np.array(list(sequence), dtype = bytes),1)


    def matches_per_tile(self, sequence):
        if len(sequence) < 5: return

        header = [''] + list(self.tile_numbers)
        emptyRow = ['' for i in np.arange(9)]

        number_of_matches = self.number_of_matches(sequence)

        allClusters = ['All clusters'] + [np.sum(self.selection(tile = tile)) for tile in self.tile_numbers]
        fullMatchCount = ['Full matches'] + [np.sum(self.tile[number_of_matches == (len(sequence))] == tile) for tile in self.tile_numbers]
        oneMisMatch = ['1 mismatch'] + [np.sum(self.tile[number_of_matches == (len(sequence) - 1)] == tile) for tile in self.tile_numbers]
        twoMisMatches = ['2 mismatches'] + [np.sum(self.tile[number_of_matches == (len(sequence) - 2)] == tile) for tile in self.tile_numbers]
        threeMisMatches = ['3 mismatches'] + [np.sum(self.tile[number_of_matches == (len(sequence) - 3)] == tile) for tile in self.tile_numbers]
        fourMisMatches = ['4 mismatches'] + [np.sum(self.tile[number_of_matches == (len(sequence) - 4)] == tile) for tile in self.tile_numbers]
        fiveMisMatches = ['5 mismatches'] + [np.sum(self.tile[number_of_matches == (len(sequence) - 5)] == tile) for tile in self.tile_numbers]

        # lessThan3mismatches = ['<=2 mismatches'] + [np.sum(data.tile[Nmatch > 47] == tile) for tile in self.tile_numbers]]

        table = tabulate([allClusters,
                          emptyRow,
                          fullMatchCount,
                          oneMisMatch,
                          twoMisMatches,
                          threeMisMatches,
                          fourMisMatches,
                          fiveMisMatches,
                          emptyRow
                          ], header)

        print(table)
        self.write_to_text_file('Matches per tile with sequence: ' + sequence + '\n\n')
        self.write_to_text_file(table)

    def export_positions_per_tile(self):
        # Export clusterpositions
        for tile in self.tile_numbers:
            x = self.x[self.selection(tile=tile)]
            y = self.y[self.selection(tile=tile)]
            np.savetxt(self.write_path.joinpath(str(tile)+'.loc'), np.transpose([x,y]), fmt='%u', delimiter='\t')

    def show_tiles(self):
        # Figure with matching cluster positions
        fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex='all', sharey='all')
        for i, tile in enumerate(self.tile_numbers):
            ax = axes.flatten()[i]
            x = self.x[self.selection(tile = tile)]
            y = self.y[self.selection(tile = tile)]
            ax.scatter(x, y, c='k', marker='.')
            ax.set_title('Tile ' + str(tile))
            ax.set_aspect('equal')

            ax.set_xlim([0, 31000])
            ax.set_ylim([0, 31000])

            if i in [4, 5, 6, 7]: ax.set_xlabel('x (FASTQ)')
            if i in [0, 4]: ax.set_ylabel('y (FASTQ)')

        fig.savefig(self.write_path.joinpath('clusterPositionsWithMatchingSequence.pdf'), bbox_inches='tight')
        fig.savefig(self.write_path.joinpath('clusterPositionsWithMatchingSequence.png'), bbox_inches='tight')



#https://stackoverflow.com/questions/9476797/how-do-i-create-character-arrays-in-numpy

#dataDictionary = {
#		'instrument': instrument,
#		'run': run,
#		'flowcell': flowcell,
#		'lane': lane,
#		'tile': tile,
#		'x': x,
#		'y': y,
#		'sample': sample,
#		'sequence': sequence,
#		'quality': quality
#		}
#
#
#df = pd.DataFrame(columns=['instrument','run','flowcell','lane','tile','x','y','sample','sequence','quality'])
#df = pd.DataFrame(dataDictionary)

# A way to get out specific sequences
#regex = re.compile('(\w){0}AA(\w){49}')
#len(list(filter(regex.match, sequence)))


#
#
#[df['sequence'][:][i] == df['sequence'][1][i] for i in range(len(df['sequence'][1]))]
#
#gc.collect()
#
#
#seqList = list(map(list, sequence))
#seqArray = np.array(list(map(np.array,seqList)))
#
#compSeq = np.array(list('ACTGTTTTTTTTTTTTTTTTACTACCTCTTTTTTTTTTTTTTT'))
#df['mmCountCompSeq'] = np.sum(seqArray[:,0:43]==compSeq, axis=1)
#df['firstMmPosition'] = np.logical_not((seqArray[:,0:43]==compSeq)).argmax(axis=1)
#
#
#df[df['mmCountCompSeq']==42]['firstMmPosition'].plot.hist(bins = 43)


# =============================================================================
# def matchCount(row):
# 	return sum([x==y for x,y in zip(row['sequence'],'CCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTC')])
# 
# 
# df.apply(lambda row: matchCount(row), axis=1)
# =============================================================================

# =============================================================================
# 
# 
# def matchCount(row):
#	return sum([x==y for x,y in zip(row['sequence'],'CCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTC')])
# test = list()
# for seq in sequence:
# 	test.append(sum([x==y for x,y in zip(seq,'GATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTC')]))
# 	
# 	
# 		
# test = list()
# for seq in sequence:
# 	test.append([ord(char) for char in seq])
# 	
# a = np.array(test)
# 
# b = np.sum(a[1,:]==a,axis=1)
# =============================================================================

# =============================================================================
# df.sequence.str[1:3]
# 
# 
# 
# a = np.array([list(x) for x in df.sequence[1:3].str[1:40]])
# np.sum(df.sequence.str[1:3]==['T','T'],axis=1)
# =============================================================================





# =============================================================================
# 
# from Bio import SeqIO
# #record_dict = SeqIO.to_dict(SeqIO.parse(r"C:\Users\Ivo Severins\Desktop\Sequencing data\20180705\Undetermined_S0_L001_I1_001.fastq","fastq"))
# 

# 
# 
# 
# records = list(SeqIO.parse(path+file, "fastq"))
# 
# x = [o.id[33:38] for o in records]
# 
# 
# x = [re.search('(?=:)...',o.id) for o in records]
# 
# re.split(
# 
# 
# test = [o.letter_annotations['phred_quality'] for o in records]
# 
# 
# 
# 
# 
# 
# file = open(path+file, 'r') 
# print file.readline(): 
# =============================================================================






