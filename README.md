Before running the pipeline scripts, please create the required Conda environments to ensure all dependencies are installed correctly.

Run the following commands in your terminal:

```bash
conda env create -f py2.yaml    
conda env create -f py3.yml     



# Downloading circRNA detection tools:

## circsplice:
git clone https://github.com/GeneFeng/CircSplice.git

## find_circ
git clone https://github.com/marvin-jens/find_circ.git

## CIRI2
wget https://downloads.sourceforge.net/project/ciri/CIRI2/CIRI_v2.0.6.zip
unzip CIRI_v2.0.6.zip -d CIRI2

## circRNA_finder
conda install -c bioconda circrna_finder
which circrna_finder
cp $(which circrna_finder) /mnt/Data/research/tools/


This file lists example SRA studies and accession numbers retrieved from the NCBI GEO database. Each entry includes the exact advanced search filters used. 

# N.B. paste the filters as they are in the search bar


https://www.ncbi.nlm.nih.gov/geo/

Filters and accession numbers:

SRA study: SRP195418   SRR9016165  
"Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing" AND ("cellular senescence"[Title] OR "senescence signature"[Title]) AND ("Illumina HiSeq 2500" OR "Illumina HiSeq 4000") AND ("2019/05/01"[Publication Date] : "2019/06/30"[Publication Date])

SRA study: SRP043644   SRR1485146  
((homo sapiens[Organism]) AND oxidative-stress induced senescence[Description]) AND ("2014-01-01"[Publication Date] : "2014-12-31"[Publication Date])

SRA study: SRP447730   SRR25177406  
(((homo spaiens[Organism]) AND ("2023-01-01"[Publication Date] : "3000"[Publication Date])) AND senescence[Description]) AND GSE236738[GEO Accession]

SRA study: SRP389281   SRR20746222  
(((homo spaiens[Organism]) AND GSE210285[GEO Accession]) AND Integrated multi-omics approach revealed cellular senescence landscape[Title]) AND Cellular senescence[Description]

SRA study: SRP438520   SRR24652356  
((homo sapiens[Organism]) AND GSE232857[GEO Accession]) AND snoRNA[Description]

SRA study: SRP121031   SRR6205681  
((homo sapiens[Organism]) AND GSE105951[GEO Accession]) AND Histone acetyltransferase p300[Title]

SRA study: SRP165928   SRR8956965  
((homo sapiens[Organism]) AND senescence[Description]) AND GSE121276[GEO Accession]
