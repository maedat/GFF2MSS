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
from Bio import SeqIO
from Bio import Seq
from BCBio import GFF


def DF_Extract(in_file, query):
    anno_DF=pd.read_csv(in_file, sep='\t')
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
    parser.add_argument('-o','--out',help="output MSS file path", required=True)
    parser.add_argument('-m','--mol',help="mol_type value (default = genomic DNA)", type=str, required=False)
    parser.add_argument('-p','--pid',help="file for protein ID (Only for the genome version-up)", type=str, required=False)
    parser.set_defaults(mol='genomic DNA', stn='', pid="NOFILE")
    return parser.parse_args()


def FASTA_CHA_SET(length_int, contig_name, organism_name_in, strain_in, mol_type_in):
    print("Processing " +  contig_name)
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

def POSITION_SET(sub_features, COUNT, POSITION):
    JOINT, JOINT_CLOSE = "", ""
    if COUNT==1: #該当RNAにおける最初のCDS/exon
        CDS_START = sub_features.location.start +1
        CDS_END = sub_features.location.end
        POSITION += str(CDS_START) + ".." + str(CDS_END)
    else: #該当RNAにおける最初のCDS/exon
        CDS_START = sub_features.location.start +1
        CDS_END = sub_features.location.end
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




if __name__ == '__main__':
    args = GET_ARGS()
    fasta_in = args.fasta
    in_file = args.gff
    anno_in = args.ann
    out_path = args.out
    locus_tag_prefix = args.loc

    mol_type_in = args.mol
    protein_id = args.pid
    
    organism_name_in = args.nam
    strain_in = args.stn


    locus_tag_counter = 0  #初期化
    PreContig = "" #初期化
    Contig_Count = 0 #初期化
    OUT_CHA="" #初期化
    with open(out_path, mode='w') as f:
        f.write(OUT_CHA)

    entries = [(fasta_in, in_file)]
    for index, (fasta_filename, gff_filename) in enumerate(entries):
    ##fastaから配列を順番に読み込む
        with open(fasta_filename, "rt") as fh:
            for record in SeqIO.parse(fh, "fasta"):
                features=[] #featuresを初期化しておく(無い時があるので）
                length = len(record) #配列長を取得する
                NowContig = record.id #配列名を取得する
                OUT_CHA += FASTA_CHA_SET(length, NowContig, organism_name_in, strain_in, mol_type_in)
    ##gffから該当配列に関するfeatureを順番に読み込む
                with open(gff_filename, "rt") as gh:
                    for rec in GFF.parse(gh):
                        if rec.id == NowContig:
                            for gene_f in rec.features:
                                locus_tag_counter += 100 #スプライシングバリアントを同じlocus_tagになるようにする
                                for RNA_f in gene_f.sub_features:
                                    if RNA_f.type == 'mRNA':
                                        COUNT = 0 #新しいmRNAに入ったらcountを０にする
                                        out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
                                        strand = RNA_f.strand
                                        if strand == -1:
                                            out_STRAND = "complement("
                                            out_STRAND_CLOSE = ")"
                                        ####GENE_INFORMATIONS
                                        mRNA_ID = RNA_f.qualifiers["ID"][0]
                                        product_name = DF_Extract(anno_in, mRNA_ID)
                                        ####
                                        for CDS_f in RNA_f.sub_features:
                                            if CDS_f.type == 'CDS':
                                                COUNT += 1
                                                transl_table = CDS_f.qualifiers.get("transl_table", ["1"])
                                                POSITION, out_JOINT, out_JOINT_CLOSE  = POSITION_SET(CDS_f, COUNT, POSITION)
                                        JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
                                        OUT_CHA += CDS_CHA_SET(JOIN, locus_tag_prefix, locus_tag_counter, mRNA_ID, product_name, transl_table[0], 9)
                                        if protein_id != "NOFILE":
                                            pid_DF=pd.read_csv(protein_id, sep='\t')
                                            tmp_DF = ""
                                            tmp_DF = pid_DF[pid_DF['ID'] == mRNA_ID]
                                            if (len(tmp_DF)>0):
                                                pid_find= tmp_DF.iat[0,1]
                                                print("-> protein_id = " + pid_find)
                                                OUT_CHA += "\t" + "\t" + "\t" + "protein_id" + "\t" +  pid_find + "\n"
                                        print(mRNA_ID + " end") 
                                    elif RNA_f.type == 'rRNA':
                                        COUNT = 0 #新しいrRNAに入ったらcountを０にする
                                        out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
                                        strand = RNA_f.strand
                                        if strand == -1:
                                            out_STRAND = "complement("
                                            out_STRAND_CLOSE = ")"
                                        ####rRNA_INFORMATIONS
                                        rRNA_ID = RNA_f.qualifiers["ID"][0]
                                        rRNA_name = RNA_f.qualifiers["Name"][0]
                                        rRNA_type = RNA_f.qualifiers["Type"][0]
                                        ####
                                        for Exon_f in RNA_f.sub_features:
                                            if Exon_f.type == 'exon':
                                                COUNT += 1
                                                POSITION, out_JOINT, out_JOINT_CLOSE  = POSITION_SET(Exon_f, COUNT, POSITION)
                                        JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
                                        OUT_CHA += rDNA_cahnger(rRNA_type,JOIN)
                                        OUT_CHA += "\t\t\t" + "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(9) + "\n"
                                        OUT_CHA += "\t\t\t" + "note" + "\t" + "transcript_id:" + rRNA_name + "\n"
                                        print(rRNA_name + " end") 
                                    elif RNA_f.type == 'tRNA':
                                        COUNT = 0 #新しいrRNAに入ったらcountを０にする
                                        out_STRAND, out_STRAND_CLOSE, POSITION, out_JOINT, out_JOINT_CLOSE="", "", "", "", "" #各出力項目を初期化
                                        tRNA_ID, tRNA_name, tRNA_anticodon, tRNA_note = "", "", "", "" #各出力項目を初期化
                                        strand = RNA_f.strand
                                        if strand == -1:
                                            out_STRAND = "complement("
                                            out_STRAND_CLOSE = ")"
                                        ####rRNA_INFORMATIONS
                                        tRNA_ID = RNA_f.qualifiers["ID"][0]
                                        tRNA_name = RNA_f.qualifiers["Name"][0]
                                        if "anticodon" in RNA_f.qualifiers:
                                            tRNA_anticodon = RNA_f.qualifiers["anticodon"][0]+"," +RNA_f.qualifiers["anticodon"][1]
                                        else:
                                            tRNA_anticodon = ""
                                        if "note" in RNA_f.qualifiers:
                                            tRNA_note = RNA_f.qualifiers["note"][0]
                                        else:
                                            tRNA_note = ""
                                        ####
                                        for Exon_f in RNA_f.sub_features:
                                            if Exon_f.type == 'exon':
                                                COUNT += 1
                                                POSITION, out_JOINT, out_JOINT_CLOSE  = POSITION_SET(Exon_f, COUNT, POSITION)
                                        JOIN = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
                                        OUT_CHA += tRNA_CHA_SET(JOIN, locus_tag_prefix, locus_tag_counter, tRNA_name, tRNA_anticodon, tRNA_note, 9)
                                        print(tRNA_name + " end") 
                with open(out_path, mode='a') as f:
                    f.write(OUT_CHA)
                    OUT_CHA=""
    




