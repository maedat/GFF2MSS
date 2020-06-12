---
ftitle: 'GFF2MSS: GFF3 converter for DDBJ submission via MSS'
tags:
  - Python
  - biology
  - genome
authors:
  - name: Taro Maeda
    orcid: 0000-0003-4185-0135
    affiliation: "1"
affiliations:
 - name: Division of Symbiotic Systems, National Institute for Basic Biology
   index: 1
date: 6 June 2020
bibliography: paper.bib

---

# Summary 

For the open, reproducible research in biology, almost all journals require authors to deposit the decoded nucleotide sequences in their articles to a public database. The amount of sequence data in a paper has expanded, and dealing with multiple genomic data in a single report is becoming common. Hence, a user-friendly procedure of genome data registration is required in recently. @Benson2018-aj

The difference in the adopted data format between general research software and registration in the public database has been a barrier to data registration. For example, although Generic Feature Format Version 3 (GFF3) has become the de-facto reference format for annotations, the DNA Data Bank of Japan (DDBJ), which is a member of International Nucleotide Sequence Database Collaboration (INSDC), does not accept this data format for genome data registration but requires original flat format file (MSS format).

The difference between these formats are often overcome using manual operations and in the scripts developed by each user's own. These ad hoc solutions, however, are problematic in terms of consistency, i.e., it may cause discrepancies between the registered data and the data actually used for the analysis. To solve this problem, several of the public databases have developed tools to convert GFF3 to their annotation formats. Genbank, which uses the sqn format, has released a tool to convert GFF3 files to sqn format (table2asn). European Nucleotide Archive, which uses the EMBL flat-file format, has no official tools, but there are several third-party tools available (e.g., EMBLmyGFF3, ). However, DDBJ does not provide a tool to convert GFF3 to MSS format, either officially or as a third party.

Here we developed a tool filling this shortage call GFF2MSS, which is a GFF3 to MSS converter. To facilitate data management, the "product" information is separated from the GFF3 by a separate tab-delimited file. It supports the inheritance of "protein_id" from existing genome data. Scripts for the integration of structural RNA (ribosomal and transfer RNA genes) data are also included.

# References

