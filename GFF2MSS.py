#!/usr/bin/python
# coding: UTF-8

# 
# ======================================================================
# Project Name    : GFF2MSS
# File Name       : GFF2MSS.py
# Version       : 1.0.8
# Encoding        : python
# Creation Date   : 2019/08/22
# Author : Taro Maeda 
# license     CC BY 4.0 
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
    parser.set_defaults(mol='genomic DNA', stn='')
    return parser.parse_args()



if __name__ == '__main__':
    args = GET_ARGS()
    fasta_in = args.fasta
    in_file = args.gff
    anno_in = args.ann
    out_path = args.out
    locus_tag_prefix = args.loc

    mol_type_in = args.mol
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
                print("Processing " +  NowContig)
                OUT_CHA = NowContig + "\t" + "source" + "\t" + str(1) + ".." + str(length) + "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, contig: " + NowContig + "\n"
                OUT_CHA += "\t" + "\t" + "\t" + "mol_type" + "\t" + mol_type_in + "\n"
                OUT_CHA += "\t" + "\t" + "\t" + "organism" + "\t" + organism_name_in + "\n"
                OUT_CHA += "\t" + "\t" + "\t" + "strain" + "\t" + strain_in + "\n"
                OUT_CHA += "\t" + "\t" + "\t" + "submitter_seqid" + "\t" +  NowContig + "\n"
    ##gffから該当配列に関するfeatureを順番に読み込む
                with open(gff_filename, "rt") as gh:
                    for rec in GFF.parse(gh):
                        if rec.id == NowContig:
                            for gene_f in rec.features:
                                for mRNA_f in gene_f.sub_features:
                                    if mRNA_f.type == 'mRNA':
                                        locus_tag_counter += 100
                                        COUNT = 0 #新しいmRNAに入ったらcountを０にする
                                        out_STRAND=""
                                        out_STRAND_CLOSE=""
                                        POSITION="" #各出力項目を初期化
                                        out_JOINT = ""
                                        out_JOINT_CLOSE=""
                                        strand = mRNA_f.strand
                                        if strand == -1:
                                            out_STRAND = "complement("
                                            out_STRAND_CLOSE = ")"
                                        ####GENE_INFORMATIONS
                                        mRNA_ID = mRNA_f.qualifiers["ID"][0]
                                        product_name = DF_Extract(anno_in, mRNA_ID)
                                        ####
                                        for CDS_f in mRNA_f.sub_features:
                                            if CDS_f.type == 'CDS':
                                                COUNT += 1
                                                transl_table = CDS_f.qualifiers.get("transl_table", ["1"])
                                                if COUNT==1: #該当mRNAにおける最初のCDS
                                                    CDS_START = CDS_f.location.start +1
                                                    CDS_END = CDS_f.location.end
                                                    POSITION = POSITION + str(CDS_START) + ".." + str(CDS_END)
                                                else: #該当mRNAにおけ二番目以降のCDS
                                                    CDS_START = CDS_f.location.start +1
                                                    CDS_END = CDS_f.location.end
                                                    POSITION = POSITION + ","+ str(CDS_START) + ".." + str(CDS_END)
                                                    out_JOINT = "join("
                                                    out_JOINT_CLOSE = ")"
                                        JOIN_CDS = out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE 
                                        OUT_CHA += "\tCDS\t"+ JOIN_CDS + "\tcodon_start\t1"+"\n"
                                        OUT_CHA += "\t\t\t" + "locus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(9) + "\n"
                                        OUT_CHA += "\t\t\t" + "note\t" + "transcript_id:" + mRNA_ID+"\n"
                                        OUT_CHA += "\t\t\t" + "product\t" + product_name+"\n"
                                        OUT_CHA += "\t\t\t" + "transl_table\t" + transl_table[0]+"\n"
                                        print(mRNA_ID + " end") 
                with open(out_path, mode='a') as f:
                    f.write(OUT_CHA)
                    OUT_CHA=""
    





