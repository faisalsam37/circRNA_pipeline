
#!/bin/bash
set -euo pipefail

# Input SRA list
SRA_LIST="sra_list.txt"

while read -r SRA_ID; do
  echo "Processing $SRA_ID..."

# Load config
source config.sh

echo "Activating conda environment 'circRNA_pipeline'..."

set +u

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    echo "Could not find Conda. Please ensure it is installed and in your PATH."
    exit 1
fi

conda activate circRNA_pipeline || {
    echo "Failed to activate Conda environment 'circRNA_pipeline'"
    exit 1
}

set -u

export PATH="$(conda info --base)/envs/circRNA_pipeline/bin:$PATH"


echo "Checking if Python can import dcc module..."
python -c "import dcc" || {
  echo "ERROR: Python cannot import dcc module. Make sure it is installed in 'circRNA_pipeline' environment."
  exit 1 
}

# Check key Python 3 packages
for pkg in pybedtools pysam; do
  echo "Checking Python package: $pkg ..."
  python -c "import $pkg" || {
    echo "$pkg is not installed in this environment.";
    exit 1;
  }
done

echo "Python 3 environment is correctly set up."

# Create output directories
mkdir -p "$FASTQ_DIR" "$TRIMMED_DIR" "$FASTQC_DIR" "$ALIGN_DIR" "$CIRCRNA_DIR" "$LOG_DIR"

# Tools list
CIRC_TOOLS=("circsplice" "circRNA_finder" "find_circ" "ciri2")

  echo "\n===== Processing $SRA_ID ====="

  # Setup logging
  TIMESTAMP=$(date +%Y%m%d_%H%M%S)
  LOG_FILE="$LOG_DIR/${SRA_ID}_$TIMESTAMP.log"
  exec > >(tee -a "$LOG_FILE") 2>&1
  echo "Logging to: $LOG_FILE"

 # Exit on error
set -e

# Check for existing gzipped FASTQ files (paired or single)
if [[ -f "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" && -f "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" ]] || [[ -f "$FASTQ_DIR/${SRA_ID}.fastq.gz" ]]; then
  echo "FASTQ files for $SRA_ID already exist. Skipping download."
else
  echo "Downloading $SRA_ID with fasterq-dump"
  if command -v fasterq-dump >/dev/null 2>&1; then
    fasterq-dump "$SRA_ID" -O "$FASTQ_DIR" --split-files --skip-technical
    # gzip only if uncompressed FASTQ files exist
    for f in "$FASTQ_DIR/${SRA_ID}"_*.fastq; do
      if [[ -f "$f" ]]; then
        gzip "$f"
      fi
    done
  else
    echo "fasterq-dump not found, falling back to fastq-dump"
    fastq-dump --gzip --split-files --skip-technical -O "$FASTQ_DIR" "$SRA_ID"
  fi
fi


 # Standardize filenames
FASTQ_FILES=($(ls "$FASTQ_DIR"/*"${SRA_ID}"*.fastq* 2>/dev/null || true))
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $NUM_FILES -eq 0 ]]; then
  echo "No FASTQ files found for $SRA_ID."
  continue

elif [[ $NUM_FILES -eq 1 ]]; then
  echo "Single-end detected"
  SRC="${FASTQ_FILES[0]}"
  DEST="$FASTQ_DIR/${SRA_ID}.fastq.gz"
  if [[ "$SRC" != "$DEST" ]]; then
    mv -f "$SRC" "$DEST"
  fi
  MODE="SE"

elif [[ $NUM_FILES -eq 2 ]]; then
  echo "Paired-end detected"
  SRC1="${FASTQ_FILES[0]}"
  DEST1="$FASTQ_DIR/${SRA_ID}_1.fastq.gz"
  if [[ "$SRC1" != "$DEST1" ]]; then
    mv -f "$SRC1" "$DEST1"
  fi

  SRC2="${FASTQ_FILES[1]}"
  DEST2="$FASTQ_DIR/${SRA_ID}_2.fastq.gz"
  if [[ "$SRC2" != "$DEST2" ]]; then
    mv -f "$SRC2" "$DEST2"
  fi
  MODE="PE"

else
  echo "Warning: more than 2 FASTQ files found, using first two."
  SRC1="${FASTQ_FILES[0]}"
  DEST1="$FASTQ_DIR/${SRA_ID}_1.fastq.gz"
  if [[ "$SRC1" != "$DEST1" ]]; then
    mv -f "$SRC1" "$DEST1"
  fi

  SRC2="${FASTQ_FILES[1]}"
  DEST2="$FASTQ_DIR/${SRA_ID}_2.fastq.gz"
  if [[ "$SRC2" != "$DEST2" ]]; then
    mv -f "$SRC2" "$DEST2"
  fi
  MODE="PE"
fi

# Run FastQC
if [[ "$MODE" == "PE" ]]; then
  SUMMARY1="$FASTQC_DIR/${SRA_ID}_1_fastqc/summary.txt"
  SUMMARY2="$FASTQC_DIR/${SRA_ID}_2_fastqc/summary.txt"

  if [[ -f "$SUMMARY1" && -f "$SUMMARY2" ]]; then
    echo "FastQC reports already exist for $SRA_ID. Skipping."
  else
    fastqc "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" -o "$FASTQC_DIR"
  fi
else
  SUMMARY1="$FASTQC_DIR/${SRA_ID}_fastqc/summary.txt"

  if [[ -f "$SUMMARY1" ]]; then
    echo "FastQC report already exists for $SRA_ID. Skipping."
  else
    fastqc "$FASTQ_DIR/${SRA_ID}.fastq.gz" -o "$FASTQC_DIR"
  fi
fi

  # Unzip FastQC
  for zipfile in "$FASTQC_DIR"/*.zip; do
    unzip -o "$zipfile" -d "$FASTQC_DIR" >/dev/null
  done

  # Determine trimming strategy
  if grep -q 'FAIL\|WARN' "$SUMMARY1" || { [[ "$MODE" == "PE" ]] && grep -q 'FAIL\|WARN' "$SUMMARY2"; }; then
    TRIM_PARAMS="SLIDINGWINDOW:4:10 LEADING:3 TRAILING:3 MINLEN:20"
  else
    TRIM_PARAMS="SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:36"
  fi

  # Trimming
  if [[ "$MODE" == "PE" ]]; then
    if [[ -f "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" && -f "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" ]]; then
      echo "Trimmed PE FASTQ exists, skipping."
    else
      java -jar "$TRIMMOMATIC_JAR" PE -threads 8 -phred33 \
        "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" \
        "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_1_unpaired.fastq.gz" \
        "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$TRIM_ADAPTER_PE":2:30:10 $TRIM_PARAMS
    fi
  else
    if [[ -f "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" ]]; then
      echo "Trimmed SE FASTQ exists, skipping."
    else
      java -jar "$TRIMMOMATIC_JAR" SE -threads 8 -phred33 \
        "$FASTQ_DIR/${SRA_ID}.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" \
        ILLUMINACLIP:"$TRIM_ADAPTER_SE":2:30:10 $TRIM_PARAMS
    fi
  fi

  echo "CIRC_TOOLS contains: ${CIRC_TOOLS[*]}"

for CIRC_TOOL in "${CIRC_TOOLS[@]}"; do
  echo -e "\n--- Running $CIRC_TOOL on $SRA_ID ---"
  
  TOOL_SUCCESS="false"

  case "$CIRC_TOOL" in
  circsplice)

ALIGN_OUT="$ALIGN_DIR/circsplice_${SRA_ID}/star"
CIRC_OUT="$CIRCRNA_DIR/circsplice_${SRA_ID}"

mkdir -p "$ALIGN_OUT" "$CIRC_OUT"

echo "Running circsplice for $SRA_ID..."

if [[ -f "$ALIGN_OUT/Chimeric.out.sam" ]]; then
  echo "🔁 Skipping STAR alignment for $SRA_ID — output already exists."
else
  echo "🔄 Performing STAR alignment for $SRA_ID..."
if [[ "$MODE" == "PE" ]]; then
  STAR --runThreadN 8 \
       --genomeDir "$STAR_INDEX" \
       --readFilesIn "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" \
       --readFilesCommand /usr/bin/gunzip -c \
       --chimSegmentMin 20 \
       --chimScoreMin 1 \
       --alignIntronMax 100000 \
       --outFilterMismatchNmax 4 \
       --alignTranscriptsPerReadNmax 100000 \
       --outFilterMultimapNmax 2 \
       --chimOutType Junctions SeparateSAMold \
       --outFileNamePrefix "$ALIGN_OUT/"
else
  STAR --runThreadN 8 \
       --genomeDir "$STAR_INDEX" \
       --readFilesIn "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" \
       --readFilesCommand /usr/bin/gunzip -c \
       --chimSegmentMin 20 \
       --chimScoreMin 1 \
       --alignIntronMax 100000 \
       --outFilterMismatchNmax 4 \
       --alignTranscriptsPerReadNmax 100000 \
       --outFilterMultimapNmax 2 \
       --chimOutType Junctions SeparateSAMold \
       --outFileNamePrefix "$ALIGN_OUT/"
fi

STAR_EXIT=$?
if [[ $STAR_EXIT -ne 0 ]]; then
  echo "❌ STAR alignment failed for $SRA_ID"
  continue
fi
fi

# After alignment, run circsplice detection if STAR output present
if [[ -f "$ALIGN_OUT/Chimeric.out.sam" ]]; then
  echo "Running circsplice on $SRA_ID..."
  perl "$CIRC_SPLICE_SCRIPT" "$ALIGN_OUT/Chimeric.out.sam" \
                         "$REFERENCE_GENOME" \
                         "$CIRC_SPLICE_REF_FLAT" \
                         > "$CIRC_OUT/${SRA_ID}_circsplice.txt"

else
  echo "❌ Missing STAR Chimeric.out.sam for circsplice. Skipping."
fi

;;


   circRNA_finder)

  ALIGN_OUT="$ALIGN_DIR/circRNA_finder_${SRA_ID}/star"
  CIRC_OUT="$CIRCRNA_DIR/circRNA_finder_${SRA_ID}"
  
  STAR_OUT="$ALIGN_DIR/circRNA_finder_${SRA_ID}/star"

  mkdir -p "$ALIGN_OUT" "$CIRC_OUT" "$STAR_OUT"

  echo "Running circRNA_finder for $SRA_ID..."
if [[ -f "$STAR_OUT/STAR_Chimeric.out.junction" ]]; then
  echo "🔁 Skipping STAR alignment for $SRA_ID — output already exists."
else
  echo "🔄 Performing STAR alignment for $SRA_ID..."

    if [[ "$MODE" == "PE" ]]; then
      STAR_CMD="$CIRCRNA_FINDER_PATH/runStar.pl \
        --inFile1 $TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz \
        --inFile2 $TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz \
        --genomeDir $STAR_INDEX \
        --outPrefix $STAR_OUT/STAR_"
    else
      STAR_CMD="$CIRCRNA_FINDER_PATH/runStar.pl \
        --inFile1 $TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz \
        --genomeDir $STAR_INDEX \
        --outPrefix $STAR_OUT/STAR_"
    fi

    echo "Executing: $STAR_CMD"
    eval "$STAR_CMD" | tee "$STAR_OUT/${SRA_ID}_runStar.log"
    STAR_EXIT=${PIPESTATUS[0]}

    if [[ $STAR_EXIT -ne 0 ]]; then
      echo "❌ STAR failed for $SRA_ID"
      continue
    fi

    if [[ ! -f "$STAR_OUT/STAR_Chimeric.out.junction" ]]; then
  echo "❌ Missing output: STAR_Chimeric.out.junction"
  continue
    fi

  fi
echo "✅ STAR alignment ready. Checking files in $STAR_OUT ..."

# List expected files
echo "Looking for files:"
echo "  $STAR_OUT/STAR_Chimeric.out.sam"
echo "  $STAR_OUT/STAR_Chimeric.out.junction"
echo "  $STAR_OUT/STAR_SJ.out.tab"

# Check files existence
missing=0
for file in STAR_Chimeric.out.sam STAR_Chimeric.out.junction STAR_SJ.out.tab; do

  if [ ! -f "$STAR_OUT/$file" ]; then
    echo "❌ Missing expected file: $STAR_OUT/$file"
    missing=1
  fi
done

if [ $missing -eq 1 ]; then
  echo "❌ Required STAR output files are missing. Skipping sample $SRA_ID."
  continue
fi

echo "STAR_OUT = '$STAR_OUT'"
echo "Files in STAR_OUT folder:"
ls -l "$STAR_OUT"


# Ensure STAR_OUT ends with a slash
case "$STAR_OUT" in
  */) ;; # already ends with slash, do nothing
  *) STAR_OUT="${STAR_OUT}/" ;;
esac

# Then call the postProcessStarAlignment.pl script
"$CIRCRNA_FINDER_PATH/postProcessStarAlignment.pl" \
  --starDir "$STAR_OUT" \
  --outDir "$CIRC_OUT" \
  --minLen 100 || {
    echo "❌ postProcessStarAlignment.pl failed"
    continue
}


TOOL_SUCCESS="true"

  ;;


ciri2)
  ALIGN_OUT="$ALIGN_DIR/ciri2_${SRA_ID}"
  CIRC_OUT="$CIRCRNA_DIR/ciri2_${SRA_ID}"
  echo "Running CIRI2 for $SRA_ID..."

  SAM_FILE="$ALIGN_OUT/${SRA_ID}_ciri2.sam"
  RESULT_FILE="$CIRC_OUT/ciri2_output.txt"
  mkdir -p "$CIRC_OUT" "$ALIGN_OUT"

  if [[ -f "$RESULT_FILE" ]]; then
    echo "CIRI2 output already exists for $SRA_ID at $RESULT_FILE. Skipping."
  else
    if [[ -f "$SAM_FILE" ]]; then
      echo "SAM file already exists for $SRA_ID at $SAM_FILE. Skipping alignment."
    else
      echo "Aligning reads with BWA..."
      if [[ "$MODE" == "PE" ]]; then
        bwa mem "$BWA_INDEX" <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz") <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz") > "$SAM_FILE"
      else
        bwa mem "$BWA_INDEX" <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz") > "$SAM_FILE"
      fi
    fi

    echo "Running CIRI2 detection..."
    perl "$CIRI2_PATH/CIRI2.pl" \
      -I "$SAM_FILE" \
      -O "$RESULT_FILE" \
      -F "$REFERENCE_GENOME" \
      -A "$GTF_FILE" \
      -M "$CIRI2_PATH/CIRI_Full_v2.1.1.jar" \
      -T 4 \
      -high
  fi
  ;;

find_circ)

  ALIGN_OUT="$ALIGN_DIR/find_circ_${SRA_ID}"
  CIRC_OUT="$CIRCRNA_DIR/find_circ_${SRA_ID}"
  
echo "Activating Python 2 environment 'circRNA_py2' for find_circ..."
ENV_NAME="circRNA_py2"
CONDA_SH=""

if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "Could not find conda.sh. Ensure Conda is installed."
  exit 1
fi

set +u
source "$CONDA_SH"
conda activate "$ENV_NAME" || {
  echo "Failed to activate conda environment '$ENV_NAME'"
  exit 1
}
set -u
if [[ "$(uname)" == "Darwin" ]]; then
  export DYLD_LIBRARY_PATH="$CONDA_PREFIX/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi


python -c "import pysam" 2>/dev/null || {
  echo "pysam not found in circRNA_py2 environment. Exiting."
  exit 1
}

  echo "Running find_circ for $SRA_ID..."


  SORTED_BAM="$ALIGN_OUT/${SRA_ID}_sorted.bam"
  ANCHORS_FASTQ="$ALIGN_OUT/${SRA_ID}_anchors.fastq.gz"
  RESULT_FILE="$CIRC_OUT/splice_sites.bed"
  mkdir -p "$ALIGN_OUT" "$CIRC_OUT"

  if [[ -f "$RESULT_FILE" ]]; then
    echo "find_circ result already exists for $SRA_ID at $RESULT_FILE. Skipping."
    continue
  fi

  if [[ ! -f "$SORTED_BAM" ]]; then
    echo "Aligning with Bowtie2 and sorting BAM..."
    if [[ "$MODE" == "PE" ]]; then
      bowtie2 -p 8 -x "$BOWTIE2_INDEX" \
        -1 <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz") \
        -2 <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz") \
        | samtools view -bS - | samtools sort -n -o "$SORTED_BAM"
    else
      bowtie2 -p 8 -x "$BOWTIE2_INDEX" \
        -U <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz") \
        | samtools view -bS - | samtools sort -n -o "$SORTED_BAM"
    fi
  else
    echo "Sorted BAM already exists for $SRA_ID. Skipping alignment."
  fi

  if [[ ! -f "$ANCHORS_FASTQ" ]]; then
    echo "Extracting unmapped reads and generating anchors..."
    samtools view -hf 4 "$SORTED_BAM" | python2 "$FIND_CIRC_PATH/unmapped2anchors.py" - | gzip > "$ANCHORS_FASTQ"
  else
    echo "Anchors FASTQ already exists. Skipping anchor generation."
  fi

  echo "Running find_circ.py..."
  bowtie2 -p 8 -x "$BOWTIE2_INDEX" -U "$ANCHORS_FASTQ" \
    | "$FIND_CIRC_PATH/find_circ.py" --genome="$REFERENCE_GENOME" --name "$SRA_ID" > "$RESULT_FILE"
  ;;

    esac 
    echo "$CIRC_TOOL complete for $SRA_ID"
  done  # ← closes the for CIRC_TOOL in ... loop

  echo "Completed all tools for $SRA_ID"
done < "$SRA_LIST"  # ← closes the first (and only) while loop



