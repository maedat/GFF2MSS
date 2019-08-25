# GFF2MSS; GFF3 converter for DDBJ submission via MSS

====


MSS (Mass Submission System) on DDBJ requires Uniq annotation format file for data submission. I here made a python script converting the standard gff3 gene model file to the MSS annotation file. 
This script makes an MSS file from a gff3 file for gene modeling, a tsv file for annotation file, and a fasta file containing genomic sequence. After the making of MSS file, you should fill "COMMON" entries (SUBMITTER, REFERENCE, etc.) before the submission for DDBJ. 


## Requirement
Python 3.7. (Biopython, pandas, argparse)

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


## Licence
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.en)

## Author
[Taro Maeda](https://github.com/maedat)
