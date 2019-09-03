# GFF2MSS; GFF3 converter for DDBJ submission via MSS

====

MSS (Mass Submission System) on DDBJ requires Uniq annotation format file for data submission. I here made a python script converting the standard gff3 gene model file to the MSS annotation. 
This script makes an MSS file from a gff3 file for gene modeling, a tsv file for annotation file, and a fasta file containing genomic sequence. I recommend preprocessing [GFF3sort.pl](https://github.com/billzt/gff3sort). After the making of MSS file, you should fill "COMMON" entries (SUBMITTER, REFERENCE, etc.) before the submission for DDBJ. 


## v.2.0
- rRNA and tRNA gene models were supported.
- A utility for tRNAscan result processing was applied. 
- The license was changed from CC BY to the MIT License

## Requirement
Python 3.7. (Biopython, pandas, argparse, bcbio-gff)

## Usage
```sh
usage: GFF2MSS.py [-h] -f FASTA -g GFF -a ANN -l LOC -n NAM [-s STN] -o OUT
                  [-m MOL]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        File path to a genome sequence file
  -g GFF, --gff GFF     gff3 file for gene modeling
  -a ANN, --ann ANN     txt file for gene annotation (header = ID,
                        Description)
  -l LOC, --loc LOC     locus_tag prefix
  -n NAM, --nam NAM     organism name
  -s STN, --stn STN     strain
  -o OUT, --out OUT     output MSS file path
  -m MOL, --mol MOL     mol_type value (default = genomic DNA)
```
  
## Demo
```sh
python3 GFF2MSS.py \
-f example/Lj3.0_Chloroplastl.fna \
-g example/Lj3.0_cp_gene_models.gff3  \
-a example/Lj3.0_anno.txt \
-l "PRE_TEST_" \
-n "Lotus japonicus" \
-s "MG-20" \
-o mss.out.txt 

```

## rDNA data
To distinct the type of rDNA sequence, please add an attribute "Type=" for each rRNA sub-features. 

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
chr1	.	tRNA	2	74	.	+	.	Parent=t91_gene;ID=t91_tRNA;Name=tRNA-Lys-CTT;anticodon=(pos:45..47,aa:Lys)
chr1	.	exon	2	51	.	+	.	Parent=t91_tRNA;ID=t91_exon_1
chr1	.	exon	60	74	.	+	.	Parent=t91_tRNA;ID=t91_exon_2
```

"Name" attribute is used as "product"  and "anticodon" attribute will be applied as "attribute" in MSS annotation file, as below. 

```txt
e.g.,
	tRNA	join(2..51,60..74)   	product	tRNA-Lys-CTT;
			locus_tag	TES_000011100
			anticodon	(pos:45..47,aa:Lys)
```

./utl/tRNA2gff3.py converts  the following type of tab-separated file, which is manually modified from tRNAscan structure prediction, to the GFF3

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
