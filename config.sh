# Reference Files
REFERENCE_GENOME="/mnt/Data/research/reference/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="/mnt/Data/research/reference/human/Homo_sapiens.GRCh38.111.gtf"

# Indexes
BWA_INDEX="/mnt/Data/research/reference/human/bwa/grch38_v111"
BOWTIE2_INDEX="/mnt/Data/research/reference/human/bowtie2/bowtie2"
HISAT2_INDEX="/mnt/Data/research/reference/human/hisat2/grch38/genome"
STAR_INDEX="/mnt/Data/research/reference/human/star/grch38_v111"

# Trimmomatic
TRIMMOMATIC_JAR="/home/anenelab/App/Trimmomatic-0.39/trimmomatic-0.39.jar"
TRIM_ADAPTER_PE="/home/anenelab/App/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
TRIM_ADAPTER_SE="/home/anenelab/App/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

# Alignment output
STAR_OUT="/mnt/Data/research/circRNA_project/hc_circRNA/alignment/circRNA_finder_SRR8060845/star"

# circRNA tools directories
CIRI2_PATH="/mnt/Data/research/circRNA_project/hc_circRNA"
FIND_CIRC_PATH="/mnt/Data/research/circRNA_project/hc_circRNA"
CIRCRNA_FINDER_PATH="/mnt/Data/research/circRNA_project/hc_circRNA"
CIRCTOOLS_PATH="/mnt/Data/research/circRNA_project/hc_circRNA/github/circtools/circtools"
CIRC_DETECT_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/alignment/circtools_SRR8060845"
CIRCTOOLS_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/github/circtools"
CIRCTOOLS_DEST="/mnt/Data/research/circRNA_project/hc_circRNA/circtools_input"
PREP_CIRCTOOLS="/mnt/Data/research/circRNA_project/hc_circRNA/github/circtools_detect/prep_circtools.sh"

PYTHON_BIN=$(which python)

# Output Directories
FASTQ_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/fastq"
TRIMMED_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/trimmed"
FASTQC_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/fastqc"
ALIGN_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/alignment"
CIRCRNA_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/circRNA"
LOG_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/logs"
CIRCTOOLS_ALIGN_DIR="/mnt/Data/research/circRNA_project/hc_circRNA/alignment/circtools_${SRA_ID}/star"      
CIRCTOOLS_OUTPUT_DIR="$CIRCRNA_DIR/circtools"    

# Thread Count
THREADS=8
