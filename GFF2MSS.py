#!/usr/bin/env python3
# coding: UTF-8

# 
# ======================================================================
# Project Name    : GFF2MSS
# File Name       : GFF2MSS.py
# Version       : 3.0.2
# Encoding        : python
# Creation Date   : 2019/08/30
# Author : Taro Maeda 
# license     MIT License (http://opensource.org/licenses/mit-license.php)
# Copyright (c) 2019 Taro Maeda
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
from BCBio import GFF
import gffpandas.gffpandas as gffpd


def DF_EXTRACT(anno_DF, query):
    tmp_DF = anno_DF[anno_DF['ID'] == query]
    out_DATA= tmp_DF.iat[0,1]
    return out_DATA


def GET_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',  help="File path to a genome sequence file", required=True)
    parser.add_argument('-g','--gff', help="gff3 file for gene modeling", required=True)
    parser.add_argument('-a','--ann',help="txt file for gene annotation (header = ID, Description)", required=True)
    parser.add_argument('-l','--loc',help="locus_tag prefix", type=str, required=True)
    parser.add_argument('-n','--nam',help="organism name", type=str, required=True)
    parser.add_argument('-s','--stn',help="strain", type=str, required=False)
    parser.add_argument('-o','--out',help="output MSS file path (default = out.mss.txt)", required=True)
    parser.add_argument('-m','--mol',help="mol_type value (default = genomic DNA)", type=str, required=False)
    parser.add_argument('-p','--pid',help="file for protein ID (Only for the genome version-up)", type=str, required=False)
    parser.add_argument('-t','--gty',help="type of linkage_evidence (default = paired-ends)", type=str, required=False)
    parser.set_defaults(mol='genomic DNA', stn='', pid="NOFILE", gty='paired-ends', out='out.mss.txt')
    return parser.parse_args()


def FASTA_CHA_SET(length_int, contig_name, organism_name_in, strain_in, mol_type_in):
    OUT_CHA = NowContig + "\t" + "source" + "\t" + str(1) + ".." + str(length_int) + "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@" + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "mol_type" + "\t" + mol_type_in + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "organism" + "\t" + organism_name_in + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "strain" + "\t" + strain_in + "\n"
    OUT_CHA += "\t" + "\t" + "\t" + "submitter_seqid" + "\t" +  "@@[entry]@@" + "\n"
    return OUT_CHA
    
    
def CDS_CHA_SET(JOIN_CDS, locus_tag_prefix, locus_tag_counter, mRNA_ID, product_name, transl_table, ZFILL):
    OUT_CHA = "\tCDS\t"+ JOIN_CDS + "\tcodon_start\t1"+"\n"
    OUT_CHA += "\t\t\t" + "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(ZFILL) + "\n"
    OUT_CHA += "\t\t\t" + "note\t" + "transcript_id:" + mRNA_ID+"\n"
    OUT_CHA += "\t\t\t" + "product\t" + product_name+"\n"
    OUT_CHA += "\t\t\t" + "transl_table\t" + transl_table+"\n"    
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

    
    
def GAP_DETECT_NP(NowSeq, OUT_CHA, link_evi):
    NowSeq_n = pd.DataFrame(record,columns=["seq"]) #pandas_dataに変える
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
    
        for start_list, end_list in zip(N_base_index_pd_end_start_S_list, N_base_index_pd_end_start_E_list):
            OUT_CHA += GAP_CHA_SET(int(start_list), int(end_list), link_evi)
    return OUT_CHA

    
def POSITION_SET(sub_features, COUNT, POSITION):
    JOINT, JOINT_CLOSE = "", ""
    if COUNT==1: #該当RNAにおける最初のCDS/exon
        CDS_START = sub_features.start +1
        CDS_END = sub_features.end
        POSITION += str(CDS_START) + ".." + str(CDS_END)
    else: #該当RNAにおける最初のCDS/exon
        CDS_START = sub_features.start +1
        CDS_END = sub_features.end
        POSITION += ","+ str(CDS_START) + ".." + str(CDS_END)
        JOINT = "join("
        JOINT_CLOSE = ")"
    return POSITION, JOINT, JOINT_CLOSE
    
    
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


def mRNA_MAKE_NP(gff_df_col, RNA_f, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA):
    COUNT = 0 #新しいmRNAに入ったらcountを０にする
    out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
    strand = RNA_f.strand
    if strand == '-':
        out_STRAND = "complement("
        out_STRAND_CLOSE = ")"
    ####GENE_INFORMATIONS
    mRNA_ID = RNA_f.ID
    product_name = DF_EXTRACT(anno_DF, mRNA_ID) #DF_EXTRACT関数で該当するアノテーション情報を取り出す
    ####    
    gff_df_col_F_sub_sub = gff_df_col[gff_df_col['Parent'] == mRNA_ID] #mRNAに対応するsubfeatureを取り出す
    for rec_sub_sub in gff_df_col_F_sub_sub.itertuples():
        if rec_sub_sub.type == 'CDS':
            COUNT += 1
            transl_table = "1"
            POSITION, out_JOINT, out_JOINT_CLOSE  = POSITION_SET(rec_sub_sub, COUNT, POSITION)
    JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
    OUT_CHA += CDS_CHA_SET(JOIN, locus_tag_prefix, locus_tag_counter, mRNA_ID, product_name, transl_table, 9)
    if pid_DF != False:
        tmp_DF = ""
        tmp_DF = pid_DF[pid_DF['ID'] == mRNA_ID]
        if (len(tmp_DF)>0):
            pid_find= tmp_DF.iat[0,1]
            print("-> protein_id = " + pid_find)
            OUT_CHA += "\t" + "\t" + "\t" + "protein_id" + "\t" +  pid_find + "\n"
    print(mRNA_ID + " end")
    return(OUT_CHA) 
    
    

def rRNA_MAKE_NP(gff_df_col, RNA_f, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA):
    COUNT = 0 #新しいrRNAに入ったらcountを０にする
    out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
    strand = RNA_f.strand
    if strand == '-':
        out_STRAND = "complement("
        out_STRAND_CLOSE = ")"
    ####rRNA_INFORMATIONS
    rRNA_ID = RNA_f.ID
    rRNA_name = RNA_f.Name
    rRNA_type = RNA_f.Type
    ####
    gff_df_col_F_sub_sub = gff_df_col[gff_df_col['Parent'] == rRNA_ID] #rRNAに対応するsubfeatureを取り出す
    for rec_sub_sub in gff_df_col_F_sub_sub.itertuples():
        if rec_sub_sub.type == 'exon':
            COUNT += 1
            POSITION, out_JOINT, out_JOINT_CLOSE  = POSITION_SET(rec_sub_sub, COUNT, POSITION)
    JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
    OUT_CHA += rDNA_cahnger(rRNA_type,JOIN)
    OUT_CHA += "\t\t\t" + "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(9) + "\n"
    OUT_CHA += "\t\t\t" + "note" + "\t" + "transcript_id:" + rRNA_name + "\n"
    print(rRNA_name + " end")
    return(OUT_CHA) 

    
def tRNA_MAKE_NP(gff_df_col, RNA_f, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA):                
    COUNT = 0 #新しいrRNAに入ったらcountを０にする
    out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
    tRNA_ID, tRNA_name, tRNA_anticodon, tRNA_note = "", "", "", "" #各出力項目を初期化
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
    
    gff_df_col_F_sub_sub = gff_df_col[gff_df_col['Parent'] == tRNA_ID] #tRNAに対応するsubfeatureを取り出す
    ####
    for rec_sub_sub in gff_df_col_F_sub_sub.itertuples():
        if rec_sub_sub.type == 'exon':
            COUNT += 1
            POSITION, out_JOINT, out_JOINT_CLOSE  = POSITION_SET(rec_sub_sub, COUNT, POSITION)
    JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
    OUT_CHA += tRNA_CHA_SET(JOIN, locus_tag_prefix, locus_tag_counter, tRNA_name, tRNA_anticodon, tRNA_note, 9)
    print(tRNA_name + " end") 
    return(OUT_CHA) 



def GFF_TO_CDS(gff_df_col,gh, NowContig, locus_tag_counter, anno_DF, pid_DF, OUT_CHA):
    gff_df_col_F = gff_df_col[gff_df_col['seq_id'] == NowContig]
    for rec in gff_df_col_F.itertuples():
    
    
        if(rec.type == "gene"):
            locus_tag_counter += 100 #ここでlocus_tagを更新することで、スプライシングバリアントを同じlocus_tagになるようにする
            Now_gene_ID = rec.ID
            gff_df_col_F_sub = gff_df_col_F[gff_df_col_F['Parent'] == Now_gene_ID] #geneに対応するsubfeatureを取り出す
            
            
            
            for rec_sub in gff_df_col_F_sub.itertuples():
                if rec_sub.type == 'mRNA':
                    OUT_CHA = mRNA_MAKE_NP(gff_df_col_F, rec_sub, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA)
                elif rec_sub.type == 'rRNA':
                    OUT_CHA = rRNA_MAKE_NP(gff_df_col_F, rec_sub, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA)
                elif rec_sub.type == 'tRNA':
                    OUT_CHA = tRNA_MAKE_NP(gff_df_col_F, rec_sub, locus_tag_prefix, locus_tag_counter, anno_DF, pid_DF, OUT_CHA)
    return locus_tag_counter, OUT_CHA




if __name__ == '__main__':
    args = GET_ARGS()
    fasta_in = args.fasta
    in_file = args.gff
    anno_in = args.ann
    anno_DF = pd.read_csv(anno_in, sep='\t') #annotationファイルの読み込みを１回で終わらせる
    
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
    
    link_evi = args.gty

    locus_tag_counter = 0  #初期化
    PreContig = "" #初期化
    Contig_Count = 0 #初期化
    OUT_CHA="" #初期化 
    
    gff_df = gffpd.read_gff3(in_file)
    gff_df_col = gff_df.attributes_to_columns()
    gff_df_col = gff_df_col.sort_values('start')
    with open(out_path, mode='w') as f:
        f.write(OUT_CHA)
    for record in SeqIO.parse(fasta_in, 'fasta'):
        record=record.upper() #塩基は全て大文字にしておく    
        print("new_contig")
        features=[] #featuresを初期化しておく(無い時があるので）
        length = len(record) #配列長を取得する
        NowContig = record.id #配列名を取得する
        print("Processing " +  NowContig)
        OUT_CHA += FASTA_CHA_SET(length, NowContig, organism_name_in, strain_in, mol_type_in)
        print("Preprocess end")
        ##fastaからgap領域を検出し記述する
        print("Gap finding")
        OUT_CHA = GAP_DETECT_NP(record, OUT_CHA,link_evi)
        print("Gap find end")
        ##gffから該当配列に関するfeatureを順番に読み込む
        print("GFF preocess")
        OUT_tmp = GFF_TO_CDS(gff_df_col,in_file, NowContig, locus_tag_counter, anno_DF, pid_DF, OUT_CHA)
        print("GFF preocess end")
        locus_tag_counter = OUT_tmp[0]
        OUT_CHA = OUT_tmp[1]
        ##MSS返還後文字列を出力
        with open(out_path, mode='a') as f:
            f.write(OUT_CHA)
            OUT_CHA=""

