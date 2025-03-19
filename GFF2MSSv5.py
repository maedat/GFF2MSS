#!/usr/bin/env python3
# coding: UTF-8

# =============================================================================
# Project Name    : GFF2MSS
# File Name       : GFF2MSS.py
# Version         : 5.0.
# Encoding        : python
# Creation Date   : 2019/08/30
# Author          : Taro Maeda
# This software is released under the MIT License, see LICENSE.
# Copyright (c) 2019-2025 Taro Maeda
# =============================================================================
#
# MSS (Mass Submission System) on DDBJ requires a unique annotation format file for data submission.
# For convenience, this script converts a standard GFF3 gene model file into an MSS annotation file.
# This script generates an MSS file from a GFF3 file, an annotation TSV file, and a genomic FASTA file.
#
# -----------------------------------------------------------------------------
# 【変更点】
#  1) BCBio(GFF)のインポートを削除し、gffpandasのみを使用してGFF3を処理します。
#  2) 日本語でのコメントを追加し、コードの可読性を向上させています。
# -----------------------------------------------------------------------------

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import gffpandas.gffpandas as gffpd  # GFF3 を pandas DataFrame に変換するために使用
import re
from distutils.util import strtobool

def get_args():
    """
    コマンドライン引数をパースし、必要な情報を取得します。
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',  help="ゲノム配列ファイルへのパス", required=True)
    parser.add_argument('-g','--gff', help="GFF3ファイル (遺伝子モデル)", required=True)
    parser.add_argument('-a','--ann',help="アノテーションTSVファイル。'ID'と'Description'列が必須で、'Locus_tag'が任意。", required=True)
    parser.add_argument('-l','--loc',help="locus_tag のプレフィックス", type=str, required=True)
    parser.add_argument('-n','--nam',help="生物種名 (organism name)", type=str, required=True)
    parser.add_argument('-s','--stn',help="株名 (strain)", type=str, default='', required=False)
    parser.add_argument('-o','--out',help="MSS形式の出力先ファイルパス (default = out.mss.txt)", required=True)
    parser.add_argument('-m','--mol',help="mol_type (default = genomic DNA)", type=str, required=False)
    parser.add_argument('-p','--pid',help="protein ID ファイル (ゲノムのバージョンアップ時のみ)", type=str, required=False)
    parser.add_argument('-t','--gty',help="linkage_evidence の種類 (default = paired-ends)", type=str, required=False)
    parser.add_argument('-c','--gct',help="遺伝暗号表(Genetic Code Table)の番号 (default = 1)", type=str, required=False)
    parser.add_argument('--ifc',
                        help="default=%(default)s: start, stopコドンから遺伝子モデルの完全性を推定し、末端に '>' や '<' を付加するかどうか",
                        default='no', type=strtobool, required=False)
    parser.add_argument('--stc',
                        help="default=%(default)s: startコドンのリスト(カンマ区切り)",
                        default='ATG', type=str, required=False)
    parser.add_argument('--iso',
                        help="default=%(default)s: 'isolate' の値(DDJBフォーマット参照)",
                        default='', type=str, required=False)
    parser.add_argument('--sex',
                        help="default=%(default)s: 'sex' の値(DDJBフォーマット参照)",
                        default='', type=str, required=False)
    parser.add_argument('--cou',
                        help="default=%(default)s: 'country' の値(DDJBフォーマット参照)",
                        default='', type=str, required=False)
    parser.add_argument('--cod',
                        help="default=%(default)s: 'collection_date' の値(DDJBフォーマット参照)",
                        default='', type=str, required=False)
    parser.add_argument('--mag',
                        help="default=%(default)s: 'gap_assembly' として扱う最小のサイズ(Nの連続長)を指定",
                        default=0, type=int, required=False)
    parser.add_argument('--gel',
                        help="default=%(default)s: 'gap_assembly' の推定サイズが既知(known)か未知(unknown)か",
                        default='known', choices=['known','unknown'], metavar='known|unknown', type=str, required=False)
    parser.add_argument('--fwg',
                        help="default=%(default)s: 'assembly_gap' をまたぐアノテーションをどのように出力するか(asis|misc_feature)",
                        default='asis', choices=['asis','misc_feature'], metavar='asis|misc_feature', type=str, required=False)
    parser.add_argument('--mis',
                        help="default=%(default)s: 指定値より小さいイントロンを'artificial_location'とみなす閾値",
                        type=int, metavar='INT', default=0, required=False)

    parser.set_defaults(mol='genomic DNA', stn='', pid="NOFILE", gty='paired-ends', out='out.mss.txt', gct='1')
    return parser.parse_args()

def fasta_cha_set(length_int, contig_name, organism_name_in, strain_in, mol_type_in,
                  country_in, isolate_in, collection_date_in, sex_in):
    """
    FASTA配列(コンティグ)全体に付与すべきソース由来のアノテーション行を生成します。
    """
    out_cha = contig_name + "\t" + "source" + "\t" + str(1) + ".." + str(length_int)
    if isolate_in != '':
        out_cha += "\t" + "ff_definition" + "\t" + "@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@\n"
    else:
        out_cha += "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@\n"

    out_cha += "\t\t\tmol_type\t" + mol_type_in + "\n"
    out_cha += "\t\t\torganism\t" + organism_name_in + "\n"
    if strain_in != '':
        out_cha += "\t\t\tstrain\t" + strain_in + "\n"
    if isolate_in != '':
        out_cha += "\t\t\tisolate\t" + isolate_in + "\n"
    if sex_in != '':
        out_cha += "\t\t\tsex\t" + sex_in + "\n"
    if country_in != '':
        out_cha += "\t\t\tcountry\t" + country_in + "\n"
    if collection_date_in != '':
        out_cha += "\t\t\tcollection_date\t" + collection_date_in + "\n"
    out_cha += "\t\t\tsubmitter_seqid\t@@[entry]@@\n"
    return out_cha

def cds_cha_set(join_cds, locus_tag_prefix, locus_tag_counter, mrna_id,
                product_name, custom_locus_tag, zfill_val, codon_start, transl_table):
    """
    MSS用のCDS行を生成するヘルパー関数。
    """
    if custom_locus_tag is None:
        out_cha = "\tCDS\t" + join_cds + "\tlocus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(zfill_val) + "\n"
    else:
        out_cha = "\tCDS\t" + join_cds + "\tlocus_tag\t" + custom_locus_tag + "\n"

    out_cha += "\t\t\tnote\t" + "transcript_id:" + mrna_id + "\n"
    out_cha += "\t\t\tproduct\t" + product_name + "\n"
    out_cha += "\t\t\ttransl_table\t" + transl_table + "\n"
    out_cha += "\t\t\tcodon_start\t" + str(codon_start) + "\n"
    return out_cha

def gap_cha_set(gap_start, gap_end, link_evi, min_assembly_gap_size, gap_estimated_length):
    """
    FASTA中の連続したN塩基(Nのgap領域)に対してアノテーション行を生成する関数。
    gap_size が指定した閾値以上の場合のみアノテーションとして出力します。
    """
    gap_size = gap_end - gap_start + 1
    out_cha = ''
    if gap_size >= min_assembly_gap_size:
        if gap_start == gap_end:  # Nが1塩基のみ
            out_cha += "\tassembly_gap\t" + str(gap_start)
        else:
            out_cha += "\tassembly_gap\t" + str(gap_start) + ".." + str(gap_end)
        out_cha += "\testimated_length\t" + gap_estimated_length + "\n"
        out_cha += "\t\t\tgap_type\twithin scaffold\n"
        out_cha += "\t\t\tlinkage_evidence\t" + link_evi + "\n"
    return out_cha

def rdna_changer(in_name, position):
    """
    rRNAやITSなどの領域に応じてアノテーション行を出し分けるための補助関数。
    """
    if in_name == "18S":
        out = "\trRNA\t" + position + "\tproduct\t18S rRNA\n"
    elif in_name == "5.8S":
        out = "\trRNA\t" + position + "\tproduct\t5.8S rRNA\n"
    elif in_name == "28S":
        out = "\trRNA\t" + position + "\tproduct\t28S rRNA\n"
    elif in_name == "ITS1":
        out = "\tmisc_RNA\t" + position + "\tnote\tinternal transcribed spacer 1\n"
    elif in_name == "ITS2":
        out = "\tmisc_RNA\t" + position + "\tnote\tinternal transcribed spacer 2\n"
    else:
        print("Undetectable rDNA type")
        out = "NA"
    return out

def trna_cha_set(position, locus_tag_prefix, locus_tag_counter, product,
                 anticodon, note, zfill_val):
    """
    tRNAのアノテーション行を生成するための関数。
    """
    out_cha = "\ttRNA\t" + position + "\tproduct\t" + product + "\n"
    out_cha += "\t\t\tlocus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(zfill_val) + "\n"
    if anticodon != "":
        out_cha += "\t\t\tanticodon\t" + anticodon + "\n"
    if note != "":
        out_cha += "\t\t\tnote\t" + note + "\n"
    return out_cha

def df_extract(anno_df, query):
    """
    アノテーションTSVファイル(anno_df)から、指定ID(query)に合致する
    product_name と custom_locus_tag を取得する。
    """
    tmp_df = anno_df[anno_df['ID'] == query]
    product_name = tmp_df.iat[0, 1]
    if tmp_df.shape[1] == 2:  # 'Locus_tag' 列が無い場合
        custom_locus_tag = None
    else:
        custom_locus_tag = tmp_df.iat[0, 2]
    return product_name, custom_locus_tag

def gap_detect_np(record, out_cha, link_evi, min_assembly_gap_size, gap_estimated_length):
    """
    FASTA配列からN塩基の連続領域を探索し、対応するアノテーション行を生成する。
    gap情報は DataFrame形式でも返す。
    """
    gap_df = pd.DataFrame({'start':[], 'end':[], 'seq_id':[]})
    nowseq_n = pd.DataFrame(list(record.seq), columns=["seq"])  # 1文字ずつをDataFrameに

    # N塩基のインデックスを取得(0-based)
    n_base_index = list((nowseq_n.query('seq == "N"').index))
    # 1-basedに変換
    n_base_index = [x + 1 for x in n_base_index]

    if n_base_index:
        n_base_index_pd = pd.DataFrame(n_base_index, columns=["n_base"])
        first_start = n_base_index[0]
        last_end   = n_base_index[-1]

        # シフトして連続領域を検出
        n_base_index_pd['n_base_2'] = n_base_index_pd['n_base'].shift(-1).fillna(0)
        n_base_index_pd["n_base_cal"] = n_base_index_pd["n_base"] - n_base_index_pd["n_base_2"]

        # gap終了位置を割り出す
        n_base_index_pd_end_start  = n_base_index_pd.query('n_base_cal < -1')
        n_base_index_pd_end_start_s_list = n_base_index_pd_end_start["n_base_2"].values.tolist()
        n_base_index_pd_end_start_e_list = n_base_index_pd_end_start["n_base"].values.tolist()

        # 連続領域の開始・終了をリストに追加
        n_base_index_pd_end_start_s_list.insert(0, first_start)
        n_base_index_pd_end_start_e_list.append(last_end)

        gap_df = pd.DataFrame({
            'start': n_base_index_pd_end_start_s_list,
            'end':   n_base_index_pd_end_start_e_list
        })
        gap_df['seq_id'] = record.id

        # gap領域をアノテーションとして出力行に追加
        for start_list, end_list in zip(n_base_index_pd_end_start_s_list,
                                        n_base_index_pd_end_start_e_list):
            out_cha += gap_cha_set(int(start_list), int(end_list),
                                   link_evi, min_assembly_gap_size,
                                   gap_estimated_length)
    return out_cha, gap_df

def exon_gap_df_compa(exon_df_iterated, test_gap_df, strand):
    """
    Exon (CDSエリア)とgap領域(N塩基)が重複しているかを調べ、
    重複があった場合はアノテーション範囲を調整する（例: codon境界の修正）。
    GAPが含まれるかどうかのフラグも返す。
    """
    new_start_list = []
    new_end_list   = []
    gap_flag       = False
    gap_start_flag = False
    gap_end_flag   = False

    test_gap_df = test_gap_df.sort_values('start')
    gap_count   = 0
    exon_start  = exon_df_iterated.start
    exon_end    = exon_df_iterated.end

    for _, gap in test_gap_df.iterrows():
        gap_s = gap['start']
        gap_e = gap['end']
        # exon範囲外なら無視
        if gap_e < exon_start or exon_end < gap_s:
            continue
        else:
            gap_count += 1
            # exon開始部とgapが重ならない場合
            if exon_start < gap_s:
                # exon内部にgapがある場合
                if exon_start < gap_s and gap_e < exon_end:
                    if gap_count == 1:
                        # gapサイズがcodonを崩す場合、3の倍数に調整
                        gap_len = gap_e - gap_s + 1
                        if gap_len % 3 == 1:
                            if strand == '-':
                                my_exon_end = gap_s - 3
                                exon_start  = gap_e + 1
                            else:
                                my_exon_end = gap_s - 1
                                exon_start  = gap_e + 3
                        elif gap_len % 3 == 2:
                            if strand == '-':
                                my_exon_end = gap_s - 2
                                exon_start  = gap_e + 1
                            else:
                                my_exon_end = gap_s - 1
                                exon_start  = gap_e + 2
                        else:
                            my_exon_end = gap_s - 1
                            exon_start  = gap_e + 1
                        new_start_list.append(exon_start)
                        new_end_list.append(my_exon_end)
                    else:
                        # 既に分割済の場合の再調整
                        new_end_list.pop()
                        new_start_list.pop()
                        gap_len = gap_e - gap_s + 1
                        if gap_len % 3 == 1:
                            if strand == '-':
                                my_exon_end = gap_s - 3
                                exon_start  = gap_e + 1
                            else:
                                my_exon_end = gap_s - 1
                                exon_start  = gap_e + 3
                        elif gap_len % 3 == 2:
                            if strand == '-':
                                my_exon_end = gap_s - 2
                                exon_start  = gap_e + 1
                            else:
                                my_exon_end = gap_s - 1
                                exon_start  = gap_e + 2
                        else:
                            my_exon_end = gap_s - 1
                            exon_start  = gap_e + 1
                        new_end_list.append(my_exon_end)
                        new_start_list.append(exon_start)

                # gapがexon下流部を超える場合
                elif exon_start < gap_s and exon_end <= gap_e:
                    if gap_count == 1:
                        exon_end = gap_s - 1
                    else:
                        new_end_list.pop()
                        new_start_list.pop()
                        exon_end = gap_s - 1

            # exon開始付近がgapと重なる場合
            elif gap_s <= exon_start and gap_e < exon_end:
                exon_start = gap_e + 1
                new_start_list.append(exon_start)
                gap_start_flag = True
                print("Gap protrudes upstream of the exon")
            new_end_list.append(exon_end)

    if gap_count == 0:
        new_start_list.append(exon_start)
        new_end_list.append(exon_end)
    else:
        gap_flag = True

    out_exon_df = pd.DataFrame({
        'start': list(map(int, new_start_list)),
        'end':   list(map(int, new_end_list))
    })
    return out_exon_df, gap_flag, gap_start_flag, gap_end_flag

def position_set(sub_features, count_val, position_str, gap_df, strand):
    """
    各exon(CDS) の位置情報を文字列として結合し、MSS用の表現を作る。
    gap処理後の座標修正およびjoin句の作成を担当。
    """
    joint, joint_close = "", ""
    gap_det, gap_flag, _, _ = exon_gap_df_compa(sub_features, gap_df, strand)
    if count_val == 1:
        for row in gap_det.itertuples():
            count_val += 1
            cds_start = row.start
            cds_end   = row.end
            if count_val == 2:
                if cds_start == cds_end:
                    position_str += str(cds_start)
                elif cds_start > cds_end:
                    pass
                else:
                    position_str += str(cds_start) + ".." + str(cds_end)
            else:
                if cds_start == cds_end:
                    position_str += "," + str(cds_start)
                elif cds_start > cds_end:
                    pass
                else:
                    position_str += "," + str(cds_start) + ".." + str(cds_end)
                joint = "join("
                joint_close = ")"
    else:
        for row in gap_det.itertuples():
            cds_start = row.start
            cds_end   = row.end
            if cds_start == cds_end:
                position_str += "," + str(cds_start)
            elif cds_start > cds_end:
                pass
            else:
                position_str += "," + str(cds_start) + ".." + str(cds_end)
            joint = "join("
            joint_close = ")"

    return position_str, joint, joint_close, gap_flag

def incomp_detect(sub_features_strand, count_val, gap_df, out_strand, codon_start):
    """
    5' 側または3' 側が不完全なCDS(exon)を判定するための準備を行う。
    """
    incomp5 = False
    rev_incomp5 = False
    if count_val == 1:
        codon_start = sub_features_strand.phase + 1
        if codon_start != 1:
            if out_strand != "":
                rev_incomp5 = True
            else:
                incomp5 = True
    return incomp5, rev_incomp5, codon_start

def mrna_make_np(gff_df_col, rna_f, locus_tag_prefix, locus_tag_counter, anno_df,
                 pid_df, out_cha, gap_df, record, infer_boundary,
                 start_codons, stop_codons, feature_with_gap, minimum_intron_size_cutoff, transl_table):
    """
    mRNA (CDS) を MSS形式に変換するためのメイン関数。
    不完全な開始/終止コドンやgapを検出し、必要に応じて座標や注釈を補正する。
    """
    count_val = 0
    incomp5 = False
    rev_incomp5 = False
    incomp3 = False
    rev_incomp3 = False
    artificial_location_flag = False
    codon_start = 1

    out_strand, out_strand_close, position_str = "", "", ""
    out_joint, out_joint_close = "", ""
    strand = rna_f.strand
    if strand == '-':
        out_strand = "complement("
        out_strand_close = ")"

    mrna_id = rna_f.ID
    product_name, custom_locus_tag = df_extract(anno_df, mrna_id)

    # mRNA に対応する CDS feature を抽出
    gff_df_sub = gff_df_col[gff_df_col['Parent'] == mrna_id]
    gff_df_sub = gff_df_sub[gff_df_sub['type'] == "CDS"].astype({'phase': int})

    # ストランドに応じて sort 順を変える
    if strand == '-':
        gff_df_sub_strand = gff_df_sub.sort_values(by='start', ascending=False)
    else:
        gff_df_sub_strand = gff_df_sub.sort_values(by='start')

    # すべてのexon(CDS)を走査し、位置を設定
    for rec_sub, rec_sub_strand in zip(gff_df_sub.itertuples(), gff_df_sub_strand.itertuples()):
        count_val += 1
        position_str, out_joint, out_joint_close, out_gap_flag = position_set(rec_sub, count_val, position_str,
                                                                              gap_df, strand)
        incomp5_tmp, rev_incomp5_tmp, codon_start = incomp_detect(rec_sub_strand, count_val, gap_df,
                                                                  out_strand, codon_start)
        incomp5 = incomp5 or incomp5_tmp
        rev_incomp5 = rev_incomp5 or rev_incomp5_tmp

    # 開始コドンや終止コドンをチェックして、不完全なら < や > を付加
    if infer_boundary:
        subfeatures = gff_df_sub_strand[gff_df_sub_strand["type"] == "CDS"]
        spliced_cds = ""
        for exon in subfeatures.itertuples():
            exon_seq = record.seq[exon.start-1 : exon.end]  # 0-based indexing
            spliced_cds += str(exon_seq)
        if strand == '-':
            spliced_cds = str(Seq(spliced_cds).reverse_complement())

        start_codon_seq = spliced_cds[:3]
        stop_codon_seq  = spliced_cds[-3:]

        # startコドン不在
        if start_codon_seq not in start_codons:
            print('Start codon not found for {}'.format(mrna_id), flush=True)
            if strand == '+':
                incomp5 = True
            else:
                rev_incomp5 = True

        # stopコドン不在
        if stop_codon_seq not in stop_codons:
            print('Stop codon not found for {}'.format(mrna_id), flush=True)
            if strand == '+':
                incomp3 = True
            else:
                rev_incomp3 = True

    if incomp5 or rev_incomp3:
        # 5'または3'が不完全な場合、先頭に '<'
        position_str = re.sub(r'^', '<', position_str)
    if rev_incomp5 or incomp3:
        # 5'または3'が不完全な場合、末尾に '>'
        position_str = re.sub(r'([^.,]+$)', r'>\1', position_str)

    join_str = out_strand + out_joint + position_str + out_joint_close + out_strand_close

    # intronサイズをチェックし、非常に短い場合に'artificial_location'を付与する
    intron_sizes = []
    # 位置文字列のうち、..区切りの部分からintronサイズを推定(簡易版)
    # ※ 実際はCDS座標の連続間を計算する方が精密
    for end_start in position_str.split('..'):
        if ',' in end_start:
            end_val, start_val = end_start.split(',')
            end_val_clean   = end_val.strip('><')
            start_val_clean = start_val.strip('><')
            end_int   = int(end_val_clean)
            start_int = int(start_val_clean)
            intron_size = start_int - end_int - 1
            intron_sizes.append(intron_size)

    intron_sizes = np.array(intron_sizes)
    out_cha_tmp = cds_cha_set(join_str, locus_tag_prefix, locus_tag_counter, mrna_id,
                              product_name, custom_locus_tag, 9, codon_start, transl_table)

    if intron_sizes.shape[0] > 0:
        min_intron_size = intron_sizes.min()
        if min_intron_size < minimum_intron_size_cutoff:
            txt = 'Introns too small ({:,} bp). "artificial_location" will be added to: {}'
            print(txt.format(min_intron_size, mrna_id), flush=True)
            artificial_location_flag = True

    # gapがある場合の処理
    if out_gap_flag:
        if feature_with_gap == 'asis':
            print('Gap found. "artificial_location" will be added to: {}'.format(mrna_id), flush=True)
            artificial_location_flag = True
        elif feature_with_gap == 'misc_feature':
            # CDSではなくmisc_featureとして出力
            print('Gap found. {} will be described with "misc_feature".'.format(mrna_id), flush=True)
            out_cha_tmp = re.sub(r'^\tCDS\t', '\tmisc_feature\t', out_cha_tmp)
            # product, transl_table, codon_startはmisc_featureでは不要なので削除
            out_cha_tmp = re.sub(r'\n\t\t\tproduct\t.*\n', '\n', out_cha_tmp)
            out_cha_tmp = re.sub(r'\n\t\t\ttransl_table\t.*\n', '\n', out_cha_tmp)
            out_cha_tmp = re.sub(r'\n\t\t\tcodon_start\t.*\n', '\n', out_cha_tmp)
            artificial_location_flag = False

    if artificial_location_flag:
        out_cha_tmp += "\t\t\tartificial_location\tlow-quality sequence region\n"

    # protein_idファイルがある場合のみ処理
    if pid_df is not False:
        tmp_df = pid_df[pid_df['ID'] == mrna_id]
        if len(tmp_df) > 0:
            pid_find = tmp_df.iat[0, 1]
            print("-> protein_id = " + pid_find)
            out_cha_tmp += "\t\t\tprotein_id\t" + pid_find + "\n"

    out_cha += out_cha_tmp
    print(mrna_id + " end")
    return out_cha

def rrna_make_np(gff_df_col, rna_f, locus_tag_prefix, locus_tag_counter, anno_df,
                 pid_df, out_cha, gap_df):
    """
    rRNA featureを MSS 形式に変換する関数。
    """
    count_val = 0
    out_strand, out_strand_close, position_str = "", "", ""
    strand = rna_f.strand
    if strand == '-':
        out_strand = "complement("
        out_strand_close = ")"

    rrna_id   = rna_f.ID
    rrna_name = rna_f.Name
    rrna_type = rna_f.Type

    gff_df_sub = gff_df_col[gff_df_col['Parent'] == rrna_id]
    out_gap_flag = False

    for rec_sub in gff_df_sub.itertuples():
        if rec_sub.type == 'exon':
            count_val += 1
            position_str, _, _, out_gap_flag = position_set(rec_sub, count_val, position_str,
                                                            gap_df, strand)

    join_str = out_strand + position_str + out_strand_close
    out_cha += rdna_changer(rrna_type, join_str)
    out_cha += "\t\t\tlocus_tag\t" + locus_tag_prefix + str(locus_tag_counter).zfill(9) + "\n"
    out_cha += "\t\t\tnote\ttranscript_id:" + rrna_name + "\n"
    if out_gap_flag:
        out_cha += "\t\t\tartificial_location\tlow-quality sequence region\n"
    print(rrna_name + " end")
    return out_cha

def trna_make_np(gff_df_col, rna_f, locus_tag_prefix, locus_tag_counter, anno_df,
                 pid_df, out_cha, gap_df):
    """
    tRNA featureを MSS 形式に変換する関数。
    """
    count_val = 0
    out_strand, out_strand_close, position_str = "", "", ""
    strand = rna_f.strand
    if strand == '-':
        out_strand = "complement("
        out_strand_close = ")"

    trna_id   = rna_f.ID
    trna_name = rna_f.Name
    trna_type = rna_f.Type

    # GFF上のカスタム属性(anticodonなど)を取得
    if len(rna_f.anticodon) != 0:
        trna_anticodon = rna_f.anticodon
    else:
        trna_anticodon = ""

    gff_df_sub = gff_df_col[gff_df_col['Parent'] == trna_id]
    out_gap_flag = False

    for rec_sub in gff_df_sub.itertuples():
        if rec_sub.type == 'exon':
            count_val += 1
            position_str, _, _, out_gap_flag = position_set(rec_sub, count_val, position_str,
                                                            gap_df, strand)

    join_str = out_strand + position_str + out_strand_close
    out_cha += trna_cha_set(join_str, locus_tag_prefix, locus_tag_counter,
                            trna_name, trna_anticodon, "", 9)
    if out_gap_flag:
        out_cha += "\t\t\tartificial_location\tlow-quality sequence region\n"
    print(trna_name + " end")
    return out_cha

def gff_to_cds(gff_df_col, now_contig, locus_tag_counter, anno_df, pid_df, out_cha,
               gap_df, record, infer_boundary, start_codons, stop_codons,
               feature_with_gap, minimum_intron_size_cutoff, transl_table):
    """
    GFFのレコードを遺伝子単位で読込み、それぞれのサブfeature(mRNA, rRNA, tRNAなど)を
    MSS形式に変換して out_cha文字列に追記する。
    """
    if gff_df_col.empty:
        return locus_tag_counter, out_cha

    gff_df_sub = gff_df_col[gff_df_col['seq_id'] == now_contig]

    # gene単位でループ
    for rec in gff_df_sub.itertuples():
        if rec.type == "gene":
            # splice バリアント(アイソフォーム)で同じ locus_tag を使いたいので
            # geneごとに 100 ずつカウントを進める(例: gene1: 100, gene2: 200, ...)
            locus_tag_counter += 100
            now_gene_id = rec.ID

            # geneの子要素(mRNA, rRNA, tRNAなど)を抽出
            sub_df = gff_df_sub[gff_df_sub['Parent'] == now_gene_id]
            for rec_sub in sub_df.itertuples():
                if rec_sub.type == 'mRNA':
                    out_cha = mrna_make_np(gff_df_sub, rec_sub, locus_tag_prefix, locus_tag_counter,
                                           anno_df, pid_df, out_cha, gap_df, record, infer_boundary,
                                           start_codons, stop_codons, feature_with_gap,
                                           minimum_intron_size_cutoff, transl_table)
                elif rec_sub.type == 'rRNA':
                    out_cha = rrna_make_np(gff_df_sub, rec_sub, locus_tag_prefix, locus_tag_counter,
                                           anno_df, pid_df, out_cha, gap_df)
                elif rec_sub.type == 'tRNA':
                    out_cha = trna_make_np(gff_df_sub, rec_sub, locus_tag_prefix, locus_tag_counter,
                                           anno_df, pid_df, out_cha, gap_df)
    return locus_tag_counter, out_cha

def get_stop_codons(genetic_code):
    """
    指定した遺伝暗号表 ID に対応するストップコドン一覧を取得する。
    """
    table = CodonTable.unambiguous_dna_by_id[int(genetic_code)]
    return list(table.stop_codons)

def main():
    """
    メイン関数。コマンドライン引数を受け取り、GFF3・FASTA・アノテーションTSVを読み込み、
    MSS形式のファイルを出力する。
    """
    args = get_args()
    
    fasta_in   = args.fasta
    gff_in     = args.gff
    anno_in    = args.ann
    out_path   = args.out
    locus_tag_prefix = args.loc
    mol_type_in      = args.mol
    protein_id       = args.pid
    organism_name_in = args.nam
    strain_in        = args.stn
    country_in       = args.cou
    isolate_in       = args.iso
    collection_date_in = args.cod
    sex_in           = args.sex
    min_assembly_gap_size = args.mag
    gap_estimated_length  = args.gel
    link_evi         = args.gty
    transl_table     = args.gct

    # protein_id ファイルが指定されている場合のみ読み込み
    if protein_id != "NOFILE":
        pid_df = pd.read_csv(protein_id, sep='\t')
    else:
        pid_df = False

    # アノテーションTSVファイル読み込み
    anno_df = pd.read_csv(anno_in, sep='\t')
    if 'Locus_tag' in anno_df.columns:
        print('The "Locus_tag" column exists in the annotation file. Custom locus_tag values will be used.')
        anno_df = anno_df.loc[:, ['ID','Description','Locus_tag']]
    else:
        print('The "Locus_tag" column does not exist in the annotation file. locus_tag values will be generated.')
        anno_df = anno_df.loc[:, ['ID','Description']]

    # GFF3をDataFrame化
    gff_df = gffpd.read_gff3(gff_in)
    gff_df_col = gff_df.attributes_to_columns()
    gff_df_col = gff_df_col.sort_values('start')

    # stopコドン配列を取得
    stop_codons = get_stop_codons(genetic_code=transl_table)
    # startコドンをカンマ区切りで受け取る
    start_codons = args.stc.split(',')

    # 出力ファイルを初期化
    with open(out_path, mode='w') as f:
        f.write("")

    locus_tag_counter = 0
    out_cha = ""

    # FASTAファイルを1配列ずつ処理
    for record in SeqIO.parse(fasta_in, 'fasta'):
        record = record.upper()  # 塩基配列を大文字化
        length = len(record)
        now_contig = record.id
        print("Processing contig:", now_contig)

        # ソース情報 (source feature) を書き出し
        out_cha += fasta_cha_set(length, now_contig, organism_name_in, strain_in,
                                 mol_type_in, country_in, isolate_in, collection_date_in, sex_in)

        # ギャップ検出
        out_cha, gap_df = gap_detect_np(record, out_cha, link_evi,
                                        min_assembly_gap_size, gap_estimated_length)

        # GFFを用いてCDS等を処理
        locus_tag_counter, out_cha = gff_to_cds(gff_df_col, now_contig, locus_tag_counter,
                                                anno_df, pid_df, out_cha, gap_df, record,
                                                args.ifc, start_codons, stop_codons,
                                                feature_with_gap=args.fwg,
                                                minimum_intron_size_cutoff=args.mis,
                                                transl_table=transl_table)

        # コンティグ単位の結果をファイルに追記
        with open(out_path, mode='a') as f:
            f.write(out_cha)
        out_cha = ""

        print("Finished contig:", now_contig)

if __name__ == '__main__':
    main()