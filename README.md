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