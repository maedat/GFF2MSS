#!/usr/bin/env python3
# coding: UTF-8

# 
# ======================================================================
# Project Name    : GFF2MSS
# File Name       : tRNA2gff3.py
# Version       : 1.0.0
# Encoding        : python
# Creation Date   : 2019\08\29
# Author : Taro Maeda 
# license     MIT License (http://opensource.org/licenses/mit-license.php)
# Copyright (c) 2019 Taro Maeda
# ======================================================================
# 

import argparse
import pandas as pd
import numpy as np

def GET_ARGS():
    parser = argparse.ArgumentParser(description="Conveter from tsv-based tRNA data to gff3")
    parser.add_argument('-i','--inf',  help="Input File (tsv file modified from tRNAscan-SE 2.0 structure file)", required=True)
    parser.add_argument('-o','--out',  help="Output gff3-like File", required=True)
    return parser.parse_args()


def DF_Extract(df):
    OUT_CHA = ""
    SOURCE =""
    ID_counter =0
    for index, row in df.iterrows():
#        print(row['intron_start'])
        ID_counter += 1
        SOURCE = row['source']
        Type = row['Type']
        Anticodon = row['Anticodon']
        INTRON = False
        if row['intron_start'] == row['intron_start'] :  #intron NaN
            INTRON = True
        PSEUDO = False
        if row['Possible_pseudogene'] == row['Possible_pseudogene']:  #pseudogene NaN
            PSEUDO = True
        #裏表
        if row['start']-row['end']<0: #strand[+]
            tRNA_START     = str(row['start'] +1)
            tRNA_END       = str(row['end'])
            codon_START    = str(row['AntC_start'])
            codon_END      = str(row['AntC_end'])
            codon_POSITION = codon_START + ".." + codon_END
            if INTRON == True:
                INTRON_START   = str(int(row['intron_start'])+1)
                INTRON_END     = str(int(row['intron_end']))
            strand = "+"
        else: #strand[-]
            tRNA_START     = str(row['end'] +1)
            tRNA_END       = str(row['start'])
            codon_START    = str(row['AntC_end'])
            codon_END      = str(row['AntC_start'])
            out_STRAND     = "complement("
            out_STRAND_CLOSE = ")"
            codon_POSITION = out_STRAND + codon_START + ".." + codon_END +out_STRAND_CLOSE
            if INTRON == True:
                INTRON_START   = str(int(row['intron_end']+1))
                INTRON_END     = str(int(row['intron_start']))
            strand = "-"
###################gene        
        OUT_CHA += SOURCE + "\t" + \
                   "." + "\t" +\
                   "gene"  + "\t" +\
                   tRNA_START  + "\t" +\
                   tRNA_END  + "\t" +\
                   "."  + "\t" +\
                   strand  + "\t" +\
                   "."  + "\t" +\
                   "ID=" + "tRNA_" + str(ID_counter) + "_gene"\
                   ";Name=" + "tRNA-" + Type + "-" + Anticodon 
        if PSEUDO == True:
            OUT_CHA +=  ";pseudogene=unknown"
        OUT_CHA += "\n"
###################tRNA        
        OUT_CHA += SOURCE + "\t" +\
                   "." + "\t" +\
                   "tRNA"  + "\t" +\
                   tRNA_START  + "\t" +\
                   tRNA_END  + "\t" +\
                   "."  + "\t" +\
                   strand  + "\t" +\
                   "."  + "\t" +\
                   "Parent="+ "tRNA_" + str(ID_counter) + "_gene"\
                   ";ID=" + "tRNA_" + str(ID_counter) + "_tRNA"\
                   ";Name=" + "tRNA-" + Type + "-" + Anticodon
        if Type == "OTHER":
            OUT_CHA += ";note=" + "unknown AA"
        elif Type == "Sup":
            OUT_CHA += ";note=" + "Suppressor_tRNA"
        else:
            OUT_CHA += ";anticodon=" + "(pos:" +  codon_POSITION + ",aa:" + Type + ")"
        OUT_CHA += "\n"
###################exon        
        if INTRON == True:
            OUT_CHA += SOURCE + "\t" +\
                   "." + "\t" +\
                   "exon"  + "\t" +\
                   tRNA_START  + "\t" +\
                   INTRON_START  + "\t" +\
                   "."  + "\t" +\
                   strand  + "\t" +\
                   "."  + "\t" +\
                   "Parent=" + "tRNA_" + str(ID_counter) + "_tRNA;"\
                   "ID=" + "tRNA_" + str(ID_counter) + "_exon_1\n"
            OUT_CHA += SOURCE + "\t" +\
                   "." + "\t" +\
                   "exon"  + "\t" +\
                   INTRON_END  + "\t" +\
                   tRNA_END  + "\t" +\
                   "."  + "\t" +\
                   strand  + "\t" +\
                   "."  + "\t" +\
                   "Parent=" + "tRNA_" + str(ID_counter) + "_tRNA;"\
                   "ID=" + "tRNA_" + str(ID_counter) + "_exon_2\n"
        elif INTRON == False:
            OUT_CHA += SOURCE + "\t" +\
                   "." + "\t" +\
                   "exon"  + "\t" +\
                   tRNA_START  + "\t" +\
                   tRNA_END  + "\t" +\
                   "."  + "\t" +\
                   strand  + "\t" +\
                   "."  + "\t" +\
                   "Parent=" + "tRNA_" + str(ID_counter) + "_tRNA;"\
                   "ID=" + "tRNA_" + str(ID_counter) + "_exon\n"
    return OUT_CHA


if __name__ == '__main__':
    args = GET_ARGS()
    in_file = args.inf
    out_file = args.out
    anno_DF=pd.read_csv(in_file,sep='\t')
    out_cha=DF_Extract(anno_DF.round())
    with open(out_file, mode='w') as f:
        f.write(out_cha)
