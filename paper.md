---
title: 'GFF2MSS: GFF3 converter for DDBJ submission via MSS'
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

For the open and reproducible research in biology, almost all journals require the deposittion of the decoded nucleotide sequences to a public database (e.g., International Nucleotide Sequence Database, INSD) [@Blaxter2016-aw]. Because it is becoming common to deal with multiple genomic data in a single report, there is a special necessity for developling of a user-friendly procedure of genome data registration. 

One of the large varrier to register is the difference in the adopted data format between general research software and the public database. For example, although Generic Feature Format Version 3 (GFF3) has become the de-facto reference format for annotations, the DNA Data Bank of Japan (DDBJ), which is a member of INSD network, does not accept GFF3 format for genome data registration but requires original flat format file (MSS format).

The difference between these formats are often overcome using manual operations and in the scripts developed by each user's own. These *ad hoc* solutions, however, are problematic in terms of consistency, i.e., it may cause discrepancies between the registered data and the data actually used for the analysis. To solve this problem, several of the public databases have developed tools to convert GFF3 to their annotation formats. Genbank, which uses the sqn format, has released a tool for the conversion, table2asn [@Benson2018-aj]. European Nucleotide Archive, which uses the EMBL flat-file format, has no official tools, but there are several third-party tools available (e.g., EMBLmyGFF3 and GFF3toEMBL) [@Norling2018-wj; @Page2016-uv]. However, DDBJ provide no conversion tools from GFF3 to MSS format, either officially or as a third party.

Here we developed a tool filling this shortage, GFF2MSS, which is a GFF3 to MSS converter. To facilitate data management, the "product" information is separated from the GFF3 by a separate tab-delimited file. It supports the inheritance of "protein_id" from existing genome data. Scripts for the integration of structural RNA (ribosomal and transfer RNA genes) data are also included.

# Acknowledgements

 This work was supported by the MEXT/JSPS KAKENHI (grant numbers 25128713 and 19K22269). 

# References

