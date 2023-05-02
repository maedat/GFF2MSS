
#!/usr/bin/env python3
# coding: UTF-8

# 
# ======================================================================
# Project Name    : GFF2MSS
# File Name       : GFF2MSS.py
# Version       : 4.0.3
# Encoding        : python
# Creation Date   : 2019/08/30
# Author : Taro Maeda 
# This software is released under the MIT License, see LICENSE.
# Copyright (c) 2019-2020 Taro Maeda
# ======================================================================
# 
# MSS (Mass Submission System) on DDBJ requires uniq annotation format file for data submission. 
# For my convenience, I made a python script convert the standard gff3 gene model file to the MSS annotation file. 
# This scrpt make a MSS file from gff3, annotation file (tsv file), and genomic fasta file. 


import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Seq
from Bio.Data import CodonTable
from BCBio import GFF
import gffpandas.gffpandas as gffpd
import re
from distutils.util import strtobool

def GET_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',  help="File path to a genome sequence file", required=True)
    parser.add_argument('-g','--gff', help="gff3 file for gene modeling", required=True)
    parser.add_argument('-a','--ann',help="tsv file for gene annotation The 'ID' and 'Description' columns are mandatory. 'Locus_tag' is optional.", required=True)
    parser.add_argument('-l','--loc',help="locus_tag prefix", type=str, required=True)
    parser.add_argument('-n','--nam',help="organism name", type=str, required=True)
    parser.add_argument('-s','--stn',help="strain", type=str, default='', required=False)
    parser.add_argument('-o','--out',help="output MSS file path (default = out.mss.txt)", required=True)
    parser.add_argument('-m','--mol',help="mol_type value (default = genomic DNA)", type=str, required=False)
    parser.add_argument('-p','--pid',help="file for protein ID (Only for the genome version-up)", type=str, required=False)
    parser.add_argument('-t','--gty',help="type of linkage_evidence (default = paired-ends)", type=str, required=False)
    parser.add_argument('-c','--gct',help="number of Genetic Code Tables (default = 1)", type=str, required=False)
    parser.add_argument('--ifc',help="default=%(default)s: inferring the completeness of gene models by the presence of start and stop codons and add '>' or '<' to the output.", default='no', type=strtobool, required=False)
    parser.add_argument('--stc',help="default=%(default)s: comma-separated list of start codons", default='ATG', type=str, required=False)
    parser.add_argument('--iso',help="default=%(default)s: The 'isolate' value. See https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html", default='', type=str, required=False)
    parser.add_argument('--sex',help="default=%(default)s: The 'sex' value. See https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html", default='', type=str, required=False)
    parser.add_argument('--cou',help="default=%(default)s: The 'country' value. See https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html", default='', type=str, required=False)
    parser.add_argument('--cod',help="default=%(default)s: The 'collection_date' value. See https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html", default='', type=str, required=False)
    parser.set_defaults(mol='genomic DNA', stn='', pid="NOFILE", gty='paired-ends', out='out.mss.txt', gct='1')

    return parser.parse_args()
    
def FASTA_CHA_SET(length_int, contig_name, organism_name_in, strain_in, mol_type_in, country_in, isolate_in,
                  collection_date_in, sex_in):
    OUT_CHA = NowContig + "\t" + "source" + "\t" + str(1) + ".." + str(length_int)
    if isolate_in != '':
        OUT_CHA += "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@" + "\n"
    else:
        OUT_CHA += "\t" + "ff_definition" + "\t" + "@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@" + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "mol_type" + "\t" + mol_type_in + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "organism" + "\t" + organism_name_in + "\n"
    if strain_in != '':
        OUT_CHA += "\t" + "\t" + "\t" + "strain" + "\t" + strain_in + "\n"
    if isolate_in != '':
        OUT_CHA += "\t" + "\t" + "\t" + "isolate" + "\t" + isolate_in + "\n"
    if sex_in != '':
        OUT_CHA += "\t" + "\t" + "\t" + "sex" + "\t" + sex_in + "\n"
    if country_in != '':
        OUT_CHA += "\t" + "\t" + "\t" + "country" + "\t" + country_in + "\n"
    if collection_date_in != '':
        OUT_CHA += "\t" + "\t" + "\t" + "collection_date" + "\t" + collection_date_in + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "submitter_seqid" + "\t" +  "@@[entry]@@" + "\n"
    return OUT_CHA
    
def CDS_CHA_SET(JOIN_CDS, locus_tag_prefix, locus_tag_counter, mRNA_ID, product_name, custom_locus_tag, ZFILL, CODON_START):
    if custom_locus_tag is None:
        OUT_CHA = "\tCDS\t"+ JOIN_CDS + "\t"+ "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(ZFILL) + "\n"
    else:
        OUT_CHA = "\tCDS\t"+ JOIN_CDS + "\t"+ "locus_tag\t" + custom_locus_tag + "\n"
    OUT_CHA += "\t\t\t" + "note\t" + "transcript_id:" + mRNA_ID+"\n"
    OUT_CHA += "\t\t\t" + "product\t" + product_name+"\n"
    OUT_CHA += "\t\t\t" + "transl_table\t" + transl_table+"\n"
    OUT_CHA += "\t\t\t" + "codon_start\t" + str(CODON_START)+"\n"    
    return OUT_CHA
    
def GAP_CHA_SET(gap_start, gap_end, link_evi):
    if(gap_start == gap_end): #1塩基のNの場合
        OUT_CHA = "\t" + "assembly_gap" + "\t" + str(gap_start) + "\t" + "estimated_length" + "\t" + "known" + "\n"
        OUT_CHA += "\t" + "\t" + "\t" + "gap_type" + "\t" + "within scaffold" + "\n"
        OUT_CHA += "\t" + "\t" + "\t" + "linkage_evidence" + "\t" + link_evi + "\n"
    else:
        OUT_CHA = "\t" + "assembly_gap" + "\t" + str(gap_start) + ".." + str(gap_end) + "\t" + "estimated_length" + "\t" + "known" + "\n"
        OUT_CHA += "\t" + "\t" + "\t" + "gap_type" + "\t" + "within scaffold" + "\n"
        OUT_CHA += "\t" + "\t" + "\t" + "linkage_evidence" + "\t" + link_evi + "\n"
    return OUT_CHA

def rDNA_cahnger(in_name, position):
    if in_name == "18S":
        OUT = "\t" + "rRNA" + "\t" + position + "\t" +"product" + "\t" + "18S rRNA" + "\n" 
    elif in_name == "5.8S":
        OUT = "\t" +"rRNA" + "\t" + position + "\t" +"product" + "\t" + "5.8S rRNA"+ "\n" 
    elif in_name == "28S":
        OUT = "\t" +"rRNA" + "\t" + position + "\t" +"product" + "\t" + "28S rRNA"+ "\n" 
    elif in_name == "ITS1":
        OUT = "\t" +"misc_RNA" + "\t" + position + "\t" +"note" + "\t" + "internal transcribed spacer 1"+ "\n" 
    elif in_name == "ITS2":
        OUT = "\t" +"misc_RNA" + "\t" + position + "\t" +"note" + "\t" + "internal transcribed spacer 2"+ "\n" 
    else:
        print("Undetactable rDNA type")
        OUT = "NA"
    return OUT

def tRNA_CHA_SET(POSITION, locus_tag_prefix, locus_tag_counter, product, anticodon, note, ZFILL):
    OUT_CHA = "\ttRNA\t"+ POSITION + "\t" + "product" + "\t" + product + "\n"
    OUT_CHA += "\t\t\t" + "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(ZFILL) + "\n"
    if anticodon != "":
        OUT_CHA += "\t\t\t" +  "anticodon" "\t" + anticodon + "\n"
    if note != "":
        OUT_CHA += "\t\t\t" +  "note" "\t" + note + "\n"
    return OUT_CHA

def DF_EXTRACT(anno_DF, query):
    tmp_DF = anno_DF[anno_DF['ID'] == query]
    product_name = tmp_DF.iat[0,1]
    if tmp_DF.shape[1]==2: # if the column "Locus_tag" is missing.
        custom_locus_tag = None
    else:
        custom_locus_tag = tmp_DF.iat[0,2]
    return product_name, custom_locus_tag
    
def GAP_DETECT_NP(NowSeq, OUT_CHA, link_evi):
    GAP_DF = pd.DataFrame({'start':[], 'end':[], 'seq_id':[]}) #Gapデータを初期化する
    NowSeq_n = pd.DataFrame(record,columns=["seq"]) #sequenceデータをpandas_dataに変える
    N_base_index = list((NowSeq_n.query('seq == "N"').index)) #N-baseである場所indexを取り出す
    N_base_index = list(map(lambda x: x+1, N_base_index))
    if N_base_index != []:
        N_base_index_pd = pd.DataFrame(N_base_index, columns=["n_base"])
    
        first_start = (N_base_index[0])
        last_end = (N_base_index[-1])
        #N-base連続部分はまとめる
    
        N_base_index_pd['n_base_2'] = N_base_index_pd['n_base'].shift(-1)
        N_base_index_pd = N_base_index_pd.fillna(0)
        N_base_index_pd["n_base_cal"] = N_base_index_pd["n_base"] - N_base_index_pd["n_base_2"]
        N_base_index_pd_end_start  = N_base_index_pd.query('n_base_cal < -1') 
        N_base_index_pd_end_start_S = N_base_index_pd_end_start["n_base_2"]
        N_base_index_pd_end_start_E = N_base_index_pd_end_start["n_base"]
        N_base_index_pd_end_start_S_list = N_base_index_pd_end_start_S.values.tolist()
        N_base_index_pd_end_start_E_list = N_base_index_pd_end_start_E.values.tolist()

        N_base_index_pd_end_start_S_list.insert(0, first_start)
        N_base_index_pd_end_start_E_list.append(last_end)
        
        GAP_DF = pd.DataFrame({ 'start' : N_base_index_pd_end_start_S_list,
                                'end' : N_base_index_pd_end_start_E_list}) #Save gap (N) region
        
        GAP_DF['seq_id'] = NowSeq.id #get contig name
                                   
        for start_list, end_list in zip(N_base_index_pd_end_start_S_list, N_base_index_pd_end_start_E_list):
            OUT_CHA += GAP_CHA_SET(int(start_list), int(end_list), link_evi)
    return OUT_CHA, GAP_DF

def EXON_GAP_DF_COMPA(EXON_DF_iterated, TEST_GAP_DF, strand):
    NEW_START_LIST = [] #output initialise 出力を初期化
    NEW_END_LIST = [] #output initialise 出力を初期化
    ART_LOC_FLAG = False #output initialise. If N-bases detected, True. 出力を初期化  N-base gapがある場合 Trueになる
    GAP_START_FLAG = False #output initialise. If Exon started with N-base gap, True. 出力を初期化  ExonがN-base gapから開始している場合（開始コドンがない可能性がある場合）Trueになる
    GAP_END_FLAG = False #output initialise. If Exon finished with N-base gap, True. 出力を初期化  ExonがN-base gapで終わっている場合（終始コドンがない可能性がある場合）Trueになる

    TEST_GAP_DF = TEST_GAP_DF.sort_values('start') #Sort the TEST_GAP_DF

    gap_count = 0 #新たなexonが始まるとリセット
    EXON_START = EXON_DF_iterated.start
    EXON_END = EXON_DF_iterated.end


    NEW_START_LIST = []
    NEW_END_LIST = [] 
    for GAP_index, GAP in TEST_GAP_DF.iterrows(): #Round for all gap region 
        GAP_START = GAP['start']
        GAP_END = GAP['end']
        if GAP_END <  EXON_START or EXON_END < GAP_START: #exon without gap (no operation)
            pass
        else: #exon with gap
            gap_count += 1 #Increment "count" if new gap was found in a exon.
            if EXON_START < GAP_START  : #no gap on upstream 上流にgapがない
                EXON_START = EXON_START #start point is Exon_start
                NEW_START_LIST.append(EXON_START)
                if EXON_START <  GAP_START and GAP_END < EXON_END: #GAP is present in the exon
                    if(gap_count == 1):
                        if (GAP_END - GAP_START+1)%3 == 1: #If codon are broken by the GAP, fix them.
                            if strand == '-': #Reverse strand
                                My_EXON_END = GAP_START-3 #GAP_START is the new temporary END
                                EXON_START = GAP_END+1 #GAP_END is now a new START (temporary)
                            else:
                                My_EXON_END = GAP_START-1 
                                EXON_START = GAP_END+3 
                        elif (GAP_END-GAP_START+1)%3 == 2:
                            if strand == '-':
                                My_EXON_END = GAP_START-2 
                                EXON_START = GAP_END+1 
                            else:
                                My_EXON_END = GAP_START-1 
                                EXON_START = GAP_END+2 
                        else:
                            My_EXON_END = GAP_START-1 
                            EXON_START = GAP_END+1 
                        NEW_END_LIST.append(My_EXON_END)
                        NEW_START_LIST.append(EXON_START)
                    elif(gap_count > 1):
                        NEW_END_LIST.pop()
                        NEW_START_LIST.pop()
                        if (GAP_END - GAP_START+1)%3 == 1:
                            if strand == '-':
                                My_EXON_END = GAP_START-3
                                EXON_START = GAP_END+1 
                            else:
                                My_EXON_END = GAP_START-1 
                                EXON_START = GAP_END+3 
                        elif (GAP_END-GAP_START+1)%3 == 2:
                            if strand == '-':
                                My_EXON_END = GAP_START-2
                                EXON_START = GAP_END+1 
                            else:
                                My_EXON_END = GAP_START-1 
                                EXON_START = GAP_END+2 
                        else:
                            My_EXON_END = GAP_START-1 
                            EXON_START = GAP_END+1 
                        NEW_END_LIST.append(My_EXON_END)
                        NEW_START_LIST.append(EXON_START)
                elif EXON_START <  GAP_START and EXON_END <= GAP_END : #GAP exists out of the downstream of exon
                    if(gap_count == 1):
                        EXON_END = GAP_START-1 #GAPstart is a new END.
                    elif(gap_count > 1):
                        NEW_END_LIST.pop()
                        NEW_START_LIST.pop()
                        EXON_END = GAP_START-1 #GAPstart is a new END.
            elif GAP_START <= EXON_START and GAP_END < EXON_END: #GAP exist upstream
                EXON_START = GAP_END+1 #GAP_END will be a new START point of exon
                NEW_START_LIST.append(EXON_START)#Added to the START list, which is a variable for output after processing
                EXON_END = EXON_END #Exon_END  become the new END
                GAP_START_FLAG = True #flag indicating a possible absence of a start codon
                print("Gap protrudes upstream of the exon")
            NEW_END_LIST.append(EXON_END) 
    if(gap_count == 0):
        NEW_START_LIST.append(EXON_START)
        NEW_END_LIST.append(EXON_END)
    elif(gap_count >=1):
        #Set up an artificial_location flag if you've ever been in a Gap process
        ART_LOC_FLAG = True
    
    
    OUT_EXON_DF = pd.DataFrame({'start' : list(map(int, NEW_START_LIST)), 'end' : list(map(int, NEW_END_LIST))}) #Save the EXON mock area
    return OUT_EXON_DF, ART_LOC_FLAG, GAP_START_FLAG, GAP_END_FLAG
    
def POSITION_SET(sub_features, COUNT, POSITION, GAP_DF, strand):
    JOINT, JOINT_CLOSE = "", ""
    ART_LOC_FLAG = False
    GAP_det, ART_LOC_FLAG_tmp, GAP_START_FLAG, GAP_END_FLAG  = EXON_GAP_DF_COMPA(sub_features, GAP_DF, strand) #It detects if the  area is covered by the CDS, cuts it off and sends it to GAP_det.
    if COUNT==1: #The first CDS/exon in the relevant RNA
        for row in GAP_det.itertuples():
            COUNT+=1
            CDS_START = row.start
            CDS_END = row.end
            if COUNT==2:
                if CDS_START==CDS_END:
                    POSITION += str(CDS_START)
                elif CDS_START>CDS_END:
                    pass
                else:
                    POSITION += str(CDS_START) + ".." + str(CDS_END)
            elif  COUNT>=3:
                if CDS_START==CDS_END:
                    POSITION += "," + str(CDS_START)
                elif CDS_START>CDS_END:
                    pass
                else:
                    POSITION += "," + str(CDS_START) + ".." + str(CDS_END)
                JOINT = "join("
                JOINT_CLOSE = ")"
    else: #The second and subsequent CDS/exon in the relevant RNA
        for row in GAP_det.itertuples():
            CDS_START = row.start
            CDS_END = row.end
            if CDS_START==CDS_END:
                POSITION += "," + str(CDS_START)
            elif CDS_START>CDS_END:
                pass
            else:
                POSITION += "," + str(CDS_START) + ".." + str(CDS_END)
            JOINT = "join("
            JOINT_CLOSE = ")"            
    ART_LOC_FLAG = ART_LOC_FLAG or ART_LOC_FLAG_tmp # If the flag is standing on either side, stand the flag (to tell the outside that the N-base is there).
    return POSITION, JOINT, JOINT_CLOSE, ART_LOC_FLAG
    
    
def INCOMP_DETECT(sub_features_STRAND, COUNT, GAP_DF, out_STRAND, CODON_START):
    INCOMP5 = False
    REVERSE_INCOMP5 = False 
    if COUNT==1: #The first CDS/exon in the relevant RNA.In the case of the reverse strand, the order is already aligned according to STRAND.
        CODON_START = sub_features_STRAND.phase +1
        if CODON_START != 1:
            if out_STRAND !="":
                REVERSE_INCOMP5 = True
            else:
                INCOMP5 =True
    else: 
        pass
    return INCOMP5, REVERSE_INCOMP5, CODON_START

def mRNA_MAKE_NP(gff_df_col, RNA_f, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF,
                 record, infer_boundary, start_codons, stop_codons):
    COUNT = 0 #Set count to 0 when entering a new mRNA.
    INCOMP5 = False #初期化
    REVERSE_INCOMP5 = False #初期化
    CODON_START = 1
    out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
    strand = RNA_f.strand
    if strand == '-':
        out_STRAND = "complement("
        out_STRAND_CLOSE = ")"
    ####GENE_INFORMATIONS
    mRNA_ID = RNA_f.ID
    product_name, custom_locus_tag = DF_EXTRACT(anno_DF, mRNA_ID) #DF_EXTRACT関数で該当するアノテーション情報を取り出す
    ####    
    gff_df_col_F_sub_sub = gff_df_col[gff_df_col['Parent'] == mRNA_ID] #mRNAに対応するsubfeatureを取り出す
    gff_df_col_F_sub_sub = gff_df_col_F_sub_sub[gff_df_col_F_sub_sub['type'] == "CDS"] #CDSだけ取り出す
    gff_df_col_F_sub_sub = gff_df_col_F_sub_sub.astype({'phase': int})
    
    if strand == '-':
        gff_df_col_F_sub_sub_STRAND = gff_df_col_F_sub_sub.sort_values(by ='start', ascending=False)
    else:
        gff_df_col_F_sub_sub_STRAND = gff_df_col_F_sub_sub.sort_values(by ='start')

    for rec_sub_sub, rec_sub_sub_strand in zip(gff_df_col_F_sub_sub.itertuples(), gff_df_col_F_sub_sub_STRAND.itertuples()):
        COUNT += 1
        POSITION, out_JOINT, out_JOINT_CLOSE, OUT_ART_LOC_FLAG = POSITION_SET(rec_sub_sub, COUNT, POSITION, GAP_DF, strand)        
        INCOMP5_tmp, REVERSE_INCOMP5_tmp, CODON_START = INCOMP_DETECT(rec_sub_sub_strand, COUNT, GAP_DF, out_STRAND, CODON_START)
        INCOMP5 = INCOMP5==True or INCOMP5_tmp==True
        REVERSE_INCOMP5 = REVERSE_INCOMP5==True or REVERSE_INCOMP5_tmp==True

    if infer_boundary:
        if strand == '+':
            start_codon_start = gff_df_col_F_sub_sub_STRAND.iloc[0].start
            start_codon_end = start_codon_start + 2
            stop_codon_end = gff_df_col_F_sub_sub_STRAND.iloc[-1].end
            stop_codon_start = stop_codon_end - 2
            start_codon_seq = record.seq[start_codon_start-1:start_codon_end]
            stop_codon_seq = record.seq[stop_codon_start-1:stop_codon_end]
            if start_codon_seq not in start_codons:
                print('Start codon not found for {}'.format(mRNA_ID), flush=True)
                INCOMP5 = True
            if stop_codon_seq not in stop_codons:
                print('Stop codon not found for {}'.format(mRNA_ID), flush=True)
                REVERSE_INCOMP5 = True
        elif strand == '-':
            start_codon_start = gff_df_col_F_sub_sub_STRAND.iloc[0].end
            start_codon_end = start_codon_start - 2
            stop_codon_end = gff_df_col_F_sub_sub_STRAND.iloc[-1].start
            stop_codon_start = stop_codon_end + 2
            start_codon_seq = record.seq[start_codon_end-1:start_codon_start].reverse_complement()
            stop_codon_seq = record.seq[stop_codon_end-1:stop_codon_start].reverse_complement()
            if start_codon_seq not in start_codons:
                print('Start codon not found for {}'.format(mRNA_ID), flush=True)
                REVERSE_INCOMP5 = True
            if stop_codon_seq not in stop_codons:
                print('Stop codon not found for {}'.format(mRNA_ID), flush=True)
                INCOMP5 = True

    if INCOMP5:
        POSITION = re.sub(r'^', '<', POSITION)
    if REVERSE_INCOMP5:
        POSITION = re.sub(r'\.([^.]*$)', r'.>\1', POSITION)
        
    JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
    
    
    OUT_CHA += CDS_CHA_SET(JOIN, locus_tag_prefix, locus_tag_counter, mRNA_ID, product_name, custom_locus_tag, 9, CODON_START)
    if OUT_ART_LOC_FLAG: #If FLAG is standing, add artificial_location
        OUT_CHA += "\t\t\t" + "artificial_location" + "\t"+ "low-quality sequence region" + "\n"
    if pid_DF != False:
        tmp_DF = ""
        tmp_DF = pid_DF[pid_DF['ID'] == mRNA_ID]
        if (len(tmp_DF)>0):
            pid_find= tmp_DF.iat[0,1]
            print("-> protein_id = " + pid_find)
            OUT_CHA += "\t" + "\t" + "\t" + "protein_id" + "\t" +  pid_find + "\n"
    print(mRNA_ID + " end")
    return(OUT_CHA) 
    
def rRNA_MAKE_NP(gff_df_col, RNA_f, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF):
    COUNT = 0 #When you enter the new RNA, set count to 0.
    INCOMP5 = False #initialization
    REVERSE_INCOMP5 = False #initialization
    CODON_START = "1"
    out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #Initialize each output item
    strand = RNA_f.strand
    if strand == '-':
        out_STRAND = "complement("
        out_STRAND_CLOSE = ")"
    ####rRNA_INFORMATIONS
    rRNA_ID = RNA_f.ID
    rRNA_name = RNA_f.Name
    rRNA_type = RNA_f.Type
    ####
    gff_df_col_F_sub_sub = gff_df_col[gff_df_col['Parent'] == rRNA_ID] #Extract the subfeature corresponding to the rRNA
    for rec_sub_sub in gff_df_col_F_sub_sub.itertuples():
        if rec_sub_sub.type == 'exon':
            COUNT += 1
            POSITION, out_JOINT, out_JOINT_CLOSE, OUT_ART_LOC_FLAG = POSITION_SET(rec_sub_sub, COUNT, POSITION, GAP_DF, strand)       
    JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
    OUT_CHA += rDNA_cahnger(rRNA_type,JOIN)
    OUT_CHA += "\t\t\t" + "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(9) + "\n"
    OUT_CHA += "\t\t\t" + "note" + "\t" + "transcript_id:" + rRNA_name + "\n"
    if OUT_ART_LOC_FLAG: #If FLAG is standing, add artificial_location
        OUT_CHA += "\t\t\t" + "artificial_location" + "\t"+ "low-quality sequence region" + "\n"
    print(rRNA_name + " end")
    return(OUT_CHA) 

def tRNA_MAKE_NP(gff_df_col, RNA_f, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF):                
    COUNT = 0 #When you enter the new RNA, set count to 0.
    INCOMP5 = False #initialization
    REVERSE_INCOMP5 = False #initialization
    CODON_START = "1"
    out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #Initialize each output item
    tRNA_ID, tRNA_name, tRNA_anticodon, tRNA_note = "", "", "", "" #initialization
    strand = RNA_f.strand
    if strand == '-':
        out_STRAND = "complement("
        out_STRAND_CLOSE = ")"
    ####rRNA_INFORMATIONS
    tRNA_ID = RNA_f.ID
    tRNA_name = RNA_f.Name
    tRNA_type = RNA_f.Type    
    
    if len(RNA_f.anticodon) != 0:
        tRNA_anticodon = RNA_f.anticodon
    else:
        tRNA_anticodon = "" 
    
    gff_df_col_F_sub_sub = gff_df_col[gff_df_col['Parent'] == tRNA_ID] #Extract the subfeatures corresponding to tRNAs
    ####
    for rec_sub_sub in gff_df_col_F_sub_sub.itertuples():
        if rec_sub_sub.type == 'exon':
            COUNT += 1
            POSITION, out_JOINT, out_JOINT_CLOSE, OUT_ART_LOC_FLAG = POSITION_SET(rec_sub_sub, COUNT, POSITION, GAP_DF, strand)          
    JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
    OUT_CHA += tRNA_CHA_SET(JOIN, locus_tag_prefix, locus_tag_counter, tRNA_name, tRNA_anticodon, tRNA_note, 9)
    if OUT_ART_LOC_FLAG: #If FLAG is standing, add artificial_location
        OUT_CHA += "\t\t\t" + "artificial_location" + "\t"+ "low-quality sequence region" + "\n"
    print(tRNA_name + " end") 
    return(OUT_CHA) 

def GFF_TO_CDS(gff_df_col, gh, NowContig, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF, record, infer_boundary, start_codons, stop_codons):
    if gff_df_col.empty:
        pass
    else:
        gff_df_col_F = gff_df_col[gff_df_col['seq_id'] == NowContig] #Select only the gff data corresponding to the contig
        for rec in gff_df_col_F.itertuples():
            if(rec.type == "gene"):
                locus_tag_counter += 100 #To make the splicing variants the same locus_tag, we updated locus_tag_counter on here. 
                Now_gene_ID = rec.ID
                gff_df_col_F_sub = gff_df_col_F[gff_df_col_F['Parent'] == Now_gene_ID] #Extract the subfeature corresponding to the gene
                for rec_sub in gff_df_col_F_sub.itertuples():
                    if rec_sub.type == 'mRNA':
                        OUT_CHA = mRNA_MAKE_NP(gff_df_col_F, rec_sub, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF, record, infer_boundary, start_codons, stop_codons)
                    elif rec_sub.type == 'rRNA':
                        OUT_CHA = rRNA_MAKE_NP(gff_df_col_F, rec_sub, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF)
                    elif rec_sub.type == 'tRNA':
                        OUT_CHA = tRNA_MAKE_NP(gff_df_col_F, rec_sub, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA, GAP_DF)
    return locus_tag_counter, OUT_CHA


def get_stop_codons(genetic_code):
    stop_codons = []
    for codon in CodonTable.unambiguous_dna_by_id[int(genetic_code)].stop_codons:
        stop_codons.append(str(codon))
    return stop_codons


if __name__ == '__main__':
    args = GET_ARGS()
    fasta_in = args.fasta
    in_file = args.gff
    anno_in = args.ann
    anno_DF = pd.read_csv(anno_in, sep='\t') #Finish loading the annotation file in one run
    if 'Locus_tag' in anno_DF.columns:
        print('The "Locus_tag" column exists in the annotation file. User-provided custom locus_tag values will be used.')
        anno_DF = anno_DF.loc[:,['ID','Description','Locus_tag']]
    else:
        txt = 'The "Locus_tag" column does not exist in the annotation file. locus_tag values will be generated with the provided prefix: {}.'
        print(txt.format(args.loc))
        anno_DF = anno_DF.loc[:,['ID','Description']]
    transl_table = args.gct
    
    out_path = args.out
    locus_tag_prefix = args.loc

    mol_type_in = args.mol
    protein_id = args.pid
    if protein_id != "NOFILE":
        pid_DF=pd.read_csv(protein_id, sep='\t')
    else:
        pid_DF=False
    organism_name_in = args.nam
    strain_in = args.stn
    country_in = args.cou
    isolate_in = args.iso
    collection_date_in = args.cod
    sex_in = args.sex
    
    link_evi = args.gty

    locus_tag_counter = 0  #initialization
    PreContig = "" #initialization
    Contig_Count = 0 #initialization
    OUT_CHA="" #initialization
    start_codons = args.stc.split(',')
    stop_codons = get_stop_codons(genetic_code=args.gct)
    
    gff_df = gffpd.read_gff3(in_file)
    gff_df_col = gff_df.attributes_to_columns()
    gff_df_col = gff_df_col.sort_values('start')
    with open(out_path, mode='w') as f:
        f.write(OUT_CHA)
    for record in SeqIO.parse(fasta_in, 'fasta'):
        record=record.upper() #All bases should be capitalized.    
        print("new_contig")
        features=[] #Initialize the features (sometimes it's not there)
        length = len(record) #Get array length
        NowContig = record.id #Get the array name.
        print("Processing " +  NowContig)
        OUT_CHA += FASTA_CHA_SET(length, NowContig, organism_name_in, strain_in, mol_type_in, country_in, isolate_in,
                                 collection_date_in, sex_in)
        ##Detect and describe the gap region from fasta.
        print("Gap finding")
        OUT_CHA, GAP_DF = GAP_DETECT_NP(record, OUT_CHA,link_evi)
        print("Gap find end")
        ##Read the features of the corresponding sequences from gff in order
        print("GFF Processing")            
        locus_tag_counter, OUT_CHA = GFF_TO_CDS(gff_df_col, in_file, NowContig, locus_tag_counter, anno_DF, pid_DF,
                                                OUT_CHA, GAP_DF, record, args.ifc, start_codons, stop_codons)
        print("GFF Processing end")
        ##Output string after MSS return
        with open(out_path, mode='a') as f:
            f.write(OUT_CHA)
        print("Processing " +NowContig+" end")
        OUT_CHA=""

