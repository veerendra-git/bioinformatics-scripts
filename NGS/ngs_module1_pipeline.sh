#!/usr/bin/env bash
set -euo pipefail

# NGS Module 1 pipeline
# Processes each fastq_seq folder separately.
#
# Expected structure:
#   data/
#   ├── person1/
#   │   └── fastq_seq/
#   │       ├── sample_1.fastq
#   │       └── sample_2.fastq
#   └── person2/
#       └── fastq_seq/
#           ├── sample_1.fastq
#           └── sample_2.fastq

ROOT="${1:-$HOME/Desktop}"
THREADS="${THREADS:-4}"
OUTNAME="NGS_Module1_Output"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

msg() {
    echo
    echo "============================================================"
    echo "$*"
    echo "============================================================"
}

# Check required programs
command -v fastqc >/dev/null 2>&1 || \
    die "FastQC not found."

command -v java >/dev/null 2>&1 || \
    die "Java not found."

command -v cutadapt >/dev/null 2>&1 || \
    die "Cutadapt not found."

# Find Trimmomatic
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-}"

if [[ -z "$TRIMMOMATIC_JAR" ]]; then
    for file in \
        "$CONDA_PREFIX/share/trimmomatic-0.41-0/trimmomatic.jar" \
        "$ROOT/FASTQC/Seq2/trimmomatic.jar"
    do
        if [[ -f "$file" ]]; then
            TRIMMOMATIC_JAR="$file"
            break
        fi
    done
fi

[[ -f "$TRIMMOMATIC_JAR" ]] || \
    die "Trimmomatic JAR not found."

# Find Trimmomatic adapter file
ADAPTER_FILE="${TRIMMOMATIC_ADAPTER:-}"

if [[ -z "$ADAPTER_FILE" ]]; then
    for file in \
        "$CONDA_PREFIX/share/trimmomatic-0.41-0/adapters/TruSeq3-PE.fa" \
        "$ROOT/FASTQC/Seq2/TruSeq3-PE.fa"
    do
        if [[ -f "$file" ]]; then
            ADAPTER_FILE="$file"
            break
        fi
    done
fi

[[ -f "$ADAPTER_FILE" ]] || \
    die "TruSeq3-PE.fa not found."

# Adapter used in the Cutadapt practical
CUTADAPT_ADAPTER="${CUTADAPT_ADAPTER:-AAAAAAAAAAAA}"

# Find all fastq_seq directories
mapfile -t FASTQ_DIRS < <(
    find "$ROOT" -type d -name "fastq_seq" -print | sort
)

((${#FASTQ_DIRS[@]} > 0)) || \
    die "No fastq_seq directories found under $ROOT."

echo "Found ${#FASTQ_DIRS[@]} dataset(s)."

for FASTQ_DIR in "${FASTQ_DIRS[@]}"; do

    PERSON_DIR="$(dirname "$FASTQ_DIR")"
    PERSON="$(basename "$PERSON_DIR")"
    OUT="$PERSON_DIR/$OUTNAME"

    mkdir -p \
        "$OUT/Practical_1_SRA" \
        "$OUT/Practical_2_FastQC/Before_Trimming" \
        "$OUT/Practical_2_FastQC/After_Trimming" \
        "$OUT/Practical_3_Trimmomatic/Paired" \
        "$OUT/Practical_3_Trimmomatic/Unpaired" \
        "$OUT/Practical_3_Trimmomatic/FastQC" \
        "$OUT/Practical_4_Cutadapt/Trimmed" \
        "$OUT/Practical_4_Cutadapt/FastQC" \
        "$OUT/logs"

    LOG="$OUT/logs/${PERSON// /_}_pipeline.log"

    msg "Processing: $PERSON"

    # Find paired-end FASTQ files
    R1="$(
        find "$FASTQ_DIR" -maxdepth 1 -type f \( \
            -name "*_1.fastq" -o \
            -name "*_1.fq" -o \
            -name "*_1.fastq.gz" -o \
            -name "*_1.fq.gz" -o \
            -name "*_R1.fastq" -o \
            -name "*_R1.fq" -o \
            -name "*_R1.fastq.gz" -o \
            -name "*_R1.fq.gz" \
        \) | sort | head -n 1
    )"

    R2="$(
        find "$FASTQ_DIR" -maxdepth 1 -type f \( \
            -name "*_2.fastq" -o \
            -name "*_2.fq" -o \
            -name "*_2.fastq.gz" -o \
            -name "*_2.fq.gz" -o \
            -name "*_R2.fastq" -o \
            -name "*_R2.fq" -o \
            -name "*_R2.fastq.gz" -o \
            -name "*_R2.fq.gz" \
        \) | sort | head -n 1
    )"

    [[ -n "$R1" && -n "$R2" ]] || \
        die "$PERSON: paired FASTQ files not found."

    echo "R1: $R1" | tee -a "$LOG"
    echo "R2: $R2" | tee -a "$LOG"

    # Practical 1
    cat > "$OUT/Practical_1_SRA/README.txt" <<EOF
PRACTICAL 1 - SRA DATA RETRIEVAL STATUS

Dataset: $PERSON

Input FASTQ files:
R1: $R1
R2: $R2

The FASTQ files were already available locally, so this pipeline
does not download the SRA dataset again.

SRA download commands can be documented separately if required
for the practical record.
EOF

    # Practical 2 - FastQC before trimming
    msg "$PERSON: FastQC before trimming"

    fastqc \
        -t "$THREADS" \
        -o "$OUT/Practical_2_FastQC/Before_Trimming" \
        "$R1" "$R2" \
        2>&1 | tee -a "$LOG"

    # Get sample name
    SAMPLE="$(basename "$R1")"
    SAMPLE="${SAMPLE%%_R1*}"
    SAMPLE="${SAMPLE%%_1*}"
    SAMPLE="${SAMPLE%.fastq}"
    SAMPLE="${SAMPLE%.fq}"
    SAMPLE="${SAMPLE%.gz}"

    [[ -n "$SAMPLE" ]] || SAMPLE="${PERSON// /_}"

    TP1="$OUT/Practical_3_Trimmomatic/Paired/${SAMPLE}_R1_paired.fastq"
    TP2="$OUT/Practical_3_Trimmomatic/Paired/${SAMPLE}_R2_paired.fastq"
    TU1="$OUT/Practical_3_Trimmomatic/Unpaired/${SAMPLE}_R1_unpaired.fastq"
    TU2="$OUT/Practical_3_Trimmomatic/Unpaired/${SAMPLE}_R2_unpaired.fastq"

    # Practical 3 - Trimmomatic
    msg "$PERSON: Trimmomatic"

    java -jar "$TRIMMOMATIC_JAR" PE \
        -threads "$THREADS" \
        "$R1" "$R2" \
        "$TP1" "$TU1" "$TP2" "$TU2" \
        "ILLUMINACLIP:${ADAPTER_FILE}:2:30:10" \
        "SLIDINGWINDOW:4:20" \
        "MINLEN:36" \
        2>&1 | tee "$OUT/logs/trimmomatic.log"

    # FastQC after Trimmomatic
    msg "$PERSON: FastQC after Trimmomatic"

    fastqc \
        -t "$THREADS" \
        -o "$OUT/Practical_2_FastQC/After_Trimming" \
        "$TP1" "$TP2" \
        2>&1 | tee -a "$LOG"

    cp -f \
        "$OUT/Practical_2_FastQC/After_Trimming/"* \
        "$OUT/Practical_3_Trimmomatic/FastQC/" \
        2>/dev/null || true

    # Practical 4 - Cutadapt
    # This practical uses the original FASTQ files.
    msg "$PERSON: Cutadapt"

    CA1="$OUT/Practical_4_Cutadapt/Trimmed/${SAMPLE}_R1_cutadapt.fastq"
    CA2="$OUT/Practical_4_Cutadapt/Trimmed/${SAMPLE}_R2_cutadapt.fastq"

    cutadapt \
        -j "$THREADS" \
        -g "$CUTADAPT_ADAPTER" \
        -G "$CUTADAPT_ADAPTER" \
        -o "$CA1" \
        -p "$CA2" \
        "$R1" "$R2" \
        2>&1 | tee "$OUT/logs/cutadapt.log"

    # FastQC after Cutadapt
    msg "$PERSON: FastQC after Cutadapt"

    fastqc \
        -t "$THREADS" \
        -o "$OUT/Practical_4_Cutadapt/FastQC" \
        "$CA1" "$CA2" \
        2>&1 | tee -a "$LOG"

    # Output summary
    cat > "$OUT/output_manifest.txt" <<EOF
NGS MODULE 1 OUTPUT

Dataset: $PERSON

Practical 1:
  Practical_1_SRA/

Practical 2:
  Practical_2_FastQC/Before_Trimming/
  Practical_2_FastQC/After_Trimming/

Practical 3 - Trimmomatic:
  Paired:
    $TP1
    $TP2
  Unpaired:
    $TU1
    $TU2
  FastQC:
    Practical_3_Trimmomatic/FastQC/

Trimmomatic parameters:
  ILLUMINACLIP:${ADAPTER_FILE}:2:30:10
  SLIDINGWINDOW:4:20
  MINLEN:36

Practical 4 - Cutadapt:
  $CA1
  $CA2

Cutadapt adapter:
  $CUTADAPT_ADAPTER

Logs:
  logs/
EOF

    echo "Completed: $PERSON" | tee -a "$LOG"
done

msg "ALL DATASETS COMPLETED"
echo "Results are stored inside each person's NGS_Module1_Output folder."
