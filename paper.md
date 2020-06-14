---
title: 'GFF2MSS: GFF3 converter for submission to DNA Data Bank of Japan'
tags:
  - Python
  - biology
  - genome
authors:
  - name: Taro Maeda
    orcid: 0000-0003-4185-0135
    affiliation: "1"
affiliations:
 - name: Faculty of Agriculture, Ryukoku University
   index: 1
date: 13 June 2020
bibliography: paper.bib

---

# Summary 

For the open and reproducible research in biology, almost all journals require the deposition of the decoded nucleotide sequences to a public database (e.g., International Nucleotide Sequence Database, INSD) [@Blaxter2016-aw]. Because it is becoming common to deal with multiple genomic data in a single article, there is a particular necessity for the development of a user-friendly tools for genomic data registration. 

One of the substantial barriers to register is the difference in the adopted data format between general research software and the public database. For example, although Generic Feature Format Version 3 (GFF3) has become the de-facto reference format for annotations, the DNA Data Bank of Japan (DDBJ), which is a member of INSD network, does not accept GFF3 format for genome data registration but requires original flat format file (MSS format).

The difference between these formats are often overcome using manual operations and in the scripts developed by each user's own. However, these *ad hoc* solutions are problematic in terms of consistency; i.e., it may cause discrepancies between the registered data and the analyzed data. Several public databases have developed tools to convert GFF3 to their annotation formats to solve this problem. Genbank, which uses the sqn format, has released a tool for the conversion, table2asn [@Benson2018-aj]. European Nucleotide Archive, which uses the EMBL flat-file format, has no official tools, but several third-party tools are available (e.g., EMBLmyGFF3 and GFF3toEMBL) [@Norling2018-wj; @Page2016-uv]. However, DDBJ provides no conversion tools from GFF3 to MSS format, either officially or as a third party.

Here we developed a tool filling this shortage, GFF2MSS, which is a GFF3 to MSS converter. The "product" information is separated from the GFF3 by a separate tab-delimited file to facilitate data management. It supports the inheritance of "protein_id" from existing genome data. Scripts for the integration of structural RNA (ribosomal and transfer RNA genes) data are also included.

# Acknowledgements

 This work was supported by the MEXT/JSPS KAKENHI (grant numbers 25128713 and 19K22269). 

# References

