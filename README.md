# GFF2MSS: GFF3 converter for DDBJ submission via MSS

====

MSS (Mass Submission System) on DDBJ requires uniq annotation format file for data submission. I here made a python script converting the standard gff3 gene model file to the MSS annotation.  This script makes an MSS file from a gff3 file for gene modeling, a tsv file for annotation file, and a fasta file containing genomic sequence. We recommend pre-processing your gff3 data via [GFF3sort.pl](https://github.com/billzt/gff3sort), before your conversion. After the making of MSS file, you should fill "COMMON" entries (SUBMITTER, REFERENCE, etc.) before the submission for DDBJ. This software is a non-official converter for MSS. We do not guarantee that DDBJ accepts the generated files.

## v.4.2
our package is now available for installation via pip. 
You can easily install it using the following command:

```txt
pip install gff2mss
```

## v.4.1
Several new features have been added. Thanks to @kfuku52.
Modifications improve flexibility while ensuring compliance with DDBJ and annotation validation tools (e.g., transChecker).
    •    Added country and collection_date fields to the source feature.
            Example: country: Singapore:Tempines, collection_date: 2019-03-26
    •    Added isolate field to the source feature.
            Example: isolate: SING2019-196
    •    Modified strain field to sex.
            Example: strain: male → sex: male
    •    Updated ff_definition format to:
            @@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@
    •    Introduced new options in GFF2MSS:
            --iso <isolate> for setting isolate
            --sex <sex> for setting sex
            --cou <country> for setting country
            --cod <collection_date> for setting collection_date
            --mag <minimum size of gap_assembly> (sets the minimum size of assembly_gap)
            --gel <gap_assembly size known/unknown> (sets whether the gap size is known or unknown,)
    •    User-Defined locus_tag Support
            Allows users to specify their own locus_tag values instead of GFF2MSS’s numbering system.
            Uses the third column of the TSV file provided via --ann.
            If the TSV has only two columns (ID & Description), the behavior remains unchanged.
    •    Adds > or < symbols to indicate incomplete start or stop codons.


## v.4.0


- Support for ambiguous "N"-base in fasta. ("N/n"-base will be converted to "assembly_gap" feature in MSS). Please install "gffpandas" library. 
- According to the "N"-base supporting, "artificial_location" qualifier was adopted to mark the modified exon/CDS by GFF2MSS. To fix frameshift caused by "N"-base gapping, V.4.0 cut 1-2 bases of the exon as necessary. These artificially modified mRNA data will be marked with this qualifier on the MSS file. 
- Processing speed was improved by using pandas (e.g., v3 = 14.8s, v4 = 4.2s).



## v.3.0

- --pid option is available to marge the previous protein ID to the new submission.   
- "@@[entry]@@", and "@@[submitter_seqid]@@" were used in submitter_seqid. 


## v.2.0
- rRNA and tRNA gene models were supported.
- A utility for tRNAscan result processing was applied. 
- The license was changed from CC BY to the MIT License

## Requirement
Python 3.7. (Biopython, numpy, pandas, argparse, bcbio-gff, gffpandas)

## Usage
```sh
usage: GFF2MSS.py [-h] -f FASTA -g GFF -a ANN -l LOC -n NAM [-s STN] -o OUT [-m MOL] [-p PID] [-t GTY] [-c GCT]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        File path to a genome sequence file
  -g GFF, --gff GFF     gff3 file for gene modeling
  -a ANN, --ann ANN     txt file for gene annotation (header = ID, Description)
  -l LOC, --loc LOC     locus_tag prefix
  -n NAM, --nam NAM     organism name
  -s STN, --stn STN     strain
  -o OUT, --out OUT     output MSS file path (default = out.mss.txt)
  -m MOL, --mol MOL     mol_type value (default = genomic DNA)
  -p PID, --pid PID     file for protein ID (Only for the genome version-up)
  -t GTY, --gty GTY     type of linkage_evidence (default = paired-ends)
  -c GCT, --gct GCT     number of Genetic Code Tables (default = 1)

```

## Demo

Example 1: Plastid DNA

```sh
python3 GFF2MSS.py \
-f example/Lj3.0_Chloroplastl.fna \
-g example/Lj3.0_cp_gene_models.gff3  \
-a example/Lj3.0_anno.txt \
-l "PRE_TEST_" \
-n "Demo japonicus" \
-s "MG-20" \
-c "11" \
-o mss.ex1.out.txt 

```

Example 2: Eukaryotic nuclear DNA with rDNA and tRNA data

```sh:rDNA_tDNA.sh
python3 GFF2MSS.py \
-f example2/test.fa \
-g example2/test.gff  \
-a example2/annot.list \
-l "PRE_TEST_" \
-n "Demo japonicus" \
-s "DAOM100" \
-o mss.ex2.out.txt 
```

Example 3: Eukaryotic nuclear DNA with "N"-base gapping

```sh
python3 GFF2MSS.py \
-f example2/test.fa \
-g example2/test.gff  \
-a example2/annot.list \
-l "PRE_TEST_" \
-n "Demo japonicus" \
-s "BB2" \
-t "paired-ends" \
-o mss.ex2.out.txt 
```



## rDNA data

To distinct the type of rDNA sequence, please add an attribute, "Type=", for each rRNA sub-features. 

- 18S ribosomal RNA; 18S
- internal transcribed spacer 1; ITS1
- 5.8S ribosomal RNA; 5.8S
- internal transcribed spacer 2; ITS2
- 28S ribosomal RNA; 28S

```txt
e.g.,
chr1	.	gene	2472097	2473907	.	+	.	ID=rRNA_001;Name=chr1_1
chr1	.	rRNA	2472097	2473907	.	+	.	Parent=rRNA_001;ID=rRNA_0011;Name=chr1_1-18S;Type=18S
chr1	.	exon	2472097	2473907	.	+	.	Parent=rRNA_0011;ID=rRNA_00111;Name=rRNA_00111
chr1	.	gene	2473908	2474013	.	+	.	ID=rRNA_002;Name=chr1_1-ITS1;Type=ITS1
chr1	.	rRNA	2473908	2474013	.	+	.	Parent=rRNA_002;ID=rRNA_0021;Name=chr1_1-ITS1;Type=ITS1
chr1	.	exon	2473908	2474013	.	+	.	Parent=rRNA_0021;ID=rRNA_00211;Name=rRNA_00211
chr1	.	gene	2474014	2474165	.	+	.	ID=rRNA_003;Name=chr1-1-5.8S;Type=5.8S
chr1	.	rRNA	2474014	2474165	.	+	.	Parent=rRNA_003;ID=rRNA_0031;Name=chr1-1-5.8S;Type=5.8S
chr1	.	exon	2474014	2474165	.	+	.	Parent=rRNA_0031;ID=rRNA_00311;Name=rRNA_00311
```

## tDNA data
Please use the following GFF3 format for  tRNA data. 

```txt
e.g.,
chr1	.	gene	2	74	.	+	.	ID=t91_gene
chr1	.	tRNA	2	74	.	+	.	Parent=t91_gene;ID=t91_tRNA;Name=tRNA-Lys;anticodon=(pos:45..47,aa:Lys)
chr1	.	exon	2	51	.	+	.	Parent=t91_tRNA;ID=t91_exon_1
chr1	.	exon	60	74	.	+	.	Parent=t91_tRNA;ID=t91_exon_2
```

"Name" attribute is used as "product"  and "anticodon" attribute will be applied as "attribute" in MSS annotation file, as below. 

```txt
e.g.,
	tRNA	join(2..51,60..74)   	product	tRNA-Lys
			locus_tag	TES_000011100
			anticodon	(pos:45..47,aa:Lys)
```

./utl/tRNA2gff3.py converts the following type of tab-separated file, which is manually modified from tRNAscan structure prediction, to the GFF3

```txt
e.g.,
source	start	end	Type	Anticodon	AntC_start	AntC_end	intron_start	intron_end	Possible_pseudogene
chr1	491917	492013	Ile	TAT	491952	491954	491956	491977	
chr1	917911	917808	Ile	TAT	917876	917874	917872	917845	pseudogene
chr1	917525	917429	Asn	ATT	917490	917488	917486	917463	
chr1	915945	915827	Pro	CGG	915909	915907	915905	915869	pseudogene
chr1	899982	899904	Ile	GAT	899952	899950		
```



## Licence
[MIT License](http://opensource.org/licenses/mit-license.php)

## Author
[Taro Maeda](https://github.com/maedat)
