#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 19:59:54 2017

@author: ruhshan
"""
from Bio import SeqIO
from Bio import pairwise2
import concurrent.futures
import pandas as pd
import math

class Proteome:
    def __init__(self, filepath, name):
        self.name = name
        self.protein_dict = {}
        for s in SeqIO.parse(filepath, 'fasta'):
            self.protein_dict[s.id.split('|')[1]]=str(s.seq)

    def __len__(self):
        return len(self.protein_dict)
    @property
    def ids(self):
        return list(self.protein_dict.keys())
    @property
    def proteins(self):
        return list(self.protein_dict.values())

    def seq(self,id):
        return self.protein_dict[id]

class Operations:
    @classmethod
    def get_score(cls,s):
        sequence=cls.proteome.seq(s)
        alignment = pairwise2.align.globalxx(cls.target,sequence)[0]
        return  100 * alignment[2]/alignment[4], s

    @classmethod
    def concurrent_scoring(cls,target, proteome):
        cls.target = target
        cls.proteome = proteome
        bscore = 0
        seq_id = ''
        ct=0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for score,s in executor.map(cls.get_score, proteome.ids):
                if score>bscore:
                    bscore = score
                    seq_id = s
                ct+=1
                if ct%200==0:
                    percentage = (ct/len(proteome))*100
                    print(proteome.name+':'+str(math.ceil(percentage))+'%',end='\r')
        print()
        return round(bscore,3), seq_id

    @staticmethod
    def get_proteomes():
        # proteomes=[Proteome('proteus_sample.fasta', 'Proteus'),
        #             Proteome('staphyllo_sample.fasta', 'Staphylo'),
        #             Proteome('ecoli_sample.fasta', 'Ecoli')]
        proteomes=[Proteome('proteus_sample', 'Proteus'),
                    Proteome('staphyllo_sample', 'Staphylo'),
                    Proteome('ecoli_sample', 'Ecoli')]
        return proteomes
    @classmethod
    def append_score(cls,protein_id,scores):
        try:
            csv = open('similarity_scores.csv','r')
        except:
            csv_header = cls.largest_proteome.name+","+",".join([x.name for x in cls.other_proteomes])+"\n"
            csv = open('similarity_scores.csv', 'w')
            csv.write(csv_header)
            csv.close()
        with open('similarity_scores.csv','a') as f:
            row = protein_id+','+','.join([str(s) for s in scores])+'\n'
            f.write(row)

    @classmethod
    def get_last_protein(cls):
        try:
            df=pd.read_csv('similarity_scores.csv')
            return df.iloc[-1][0]
        except:
            return 0
    @classmethod
    def perform_calulation(cls, proteomes):
        cls.largest_proteome = proteomes[0]
        for proteome in proteomes:
            if len(proteome) > len(cls.largest_proteome):
                cls.largest_proteome = proteome

        lp_index=proteomes.index(cls.largest_proteome)
        proteomes.pop(lp_index)

        cls.other_proteomes = proteomes
        try:
            done_up_to = cls.largest_proteome.ids.index(cls.get_last_protein())
        except:
            done_up_to = 0
        print("####### {} of {} completed #######".format(done_up_to+1, len(cls.largest_proteome)))
        if done_up_to>0:
            done_up_to+=1
        ct=1
        for protein in cls.largest_proteome.ids[done_up_to:]:
            sequence = cls.largest_proteome.seq(protein)
            scores=[]
            protein_ids=[]
            print("Searchin for {} of {}:".format(protein, cls.largest_proteome.name))
            for proteome in proteomes:
                score, protein_id=cls.concurrent_scoring(sequence, proteome)
                print("     For {} best similarity is {}{} with {}".format(proteome.name, score,'%',protein_id))
                scores.append(score)
                protein_ids.append(protein_id)
            cls.append_score(protein, scores)
            print("####### {} of {} completed #######".format(done_up_to+ct, len(cls.largest_proteome)))
            ct+=1

if __name__=='__main__':
    proteomes = Operations.get_proteomes()
    #Operations.get_last_index()
    Operations.perform_calulation(proteomes)
    #print(Operations.concurrent_scoring(proteomes[1].proteins[3], proteomes[2]))
