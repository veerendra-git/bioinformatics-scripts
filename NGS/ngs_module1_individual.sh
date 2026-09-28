#!/usr/bin/env bash
set -e

# NGS Module 1 - Single dataset pipeline
#
# Usage:
#   bash ngs_module1_individual.sh /path/to/fastq_seq
#
# If no directory is given, fastq_seq in the current directory is used.

FASTQ_DIR="${1:-./fastq_seq}"
OUTPUT_DIR="${2:-./NGS_Module1_Output}"
THREADS="${THREADS:-4}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

command -v fastqc >/dev/null 2>&1 || die "FastQC not found."
command -v java >/dev/null 2>&1 || die "Java not found."
command -v cutadapt >/dev/null 2>&1 || die "Cutadapt not found."

[[ -d "$FASTQ_DIR" ]] || die "FASTQ directory not found: $FASTQ_DIR"

# Find paired-end FASTQ files
R1="$(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*_1.fastq" | sort | head -n 1)"
R2="$(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*_2.fastq" | sort | head -n 1)"

[[ -n "$R1" && -n "$R2" ]] || \
    die "Paired-end FASTQ files not found in $FASTQ_DIR"

# Find Trimmomatic
TRIM_JAR="${TRIMMOMATIC_JAR:-}"

if [[ -z "$TRIM_JAR" ]]; then
    TRIM_JAR="$CONDA_PREFIX/share/trimmomatic-0.41-0/trimmomatic.jar"
fi

[[ -f "$TRIM_JAR" ]] || die "Trimmomatic JAR not found."

# Find adapter file
ADAPTER_FILE="${TRIMMOMATIC_ADAPTER:-$CONDA_PREFIX/share/trimmomatic-0.41-0/adapters/TruSeq3-PE.fa}"

[[ -f "$ADAPTER_FILE" ]] || \
    die "TruSeq3-PE.fa not found."

mkdir -p \
    "$OUTPUT_DIR/Practical_1_SRA" \
    "$OUTPUT_DIR/Practical_2_FastQC/Before_Trimming" \
    "$OUTPUT_DIR/Practical_3_Trimmomatic/Paired" \
    "$OUTPUT_DIR/Practical_3_Trimmomatic/Unpaired" \
    "$OUTPUT_DIR/Practical_3_Trimmomatic/FastQC" \
    "$OUTPUT_DIR/Practical_4_Cutadapt/Trimmed" \
    "$OUTPUT_DIR/Practical_4_Cutadapt/FastQC"

echo "R1: $R1"
echo "R2: $R2"

# Practical 1 - SRA input
cat > "$OUTPUT_DIR/Practical_1_SRA/README.txt" <<EOF
SRA Toolkit practical input:

R1 = $R1
R2 = $R2

The paired FASTQ files were already available locally,
so SRA retrieval/conversion was not repeated.
EOF

# Practical 2 - FastQC before trimming
echo "Running FastQC before trimming..."

fastqc \
    "$R1" "$R2" \
    -o "$OUTPUT_DIR/Practical_2_FastQC/Before_Trimming" \
    -t "$THREADS"

# Sample name
BASE="$(basename "$R1" _1.fastq)"

# Practical 3 - Trimmomatic
P1="$OUTPUT_DIR/Practical_3_Trimmomatic/Paired/${BASE}_1_paired.fastq"
U1="$OUTPUT_DIR/Practical_3_Trimmomatic/Unpaired/${BASE}_1_unpaired.fastq"
P2="$OUTPUT_DIR/Practical_3_Trimmomatic/Paired/${BASE}_2_paired.fastq"
U2="$OUTPUT_DIR/Practical_3_Trimmomatic/Unpaired/${BASE}_2_unpaired.fastq"

echo "Running Trimmomatic..."

java -jar "$TRIM_JAR" PE \
    -threads "$THREADS" \
    "$R1" "$R2" \
    "$P1" "$U1" "$P2" "$U2" \
    "ILLUMINACLIP:${ADAPTER_FILE}:2:30:10" \
    "SLIDINGWINDOW:4:20" \
    "MINLEN:36"

# FastQC after Trimmomatic
echo "Running FastQC after Trimmomatic..."

fastqc \
    "$P1" "$P2" \
    -o "$OUTPUT_DIR/Practical_3_Trimmomatic/FastQC" \
    -t "$THREADS"

# Practical 4 - Cutadapt
CA1="$OUTPUT_DIR/Practical_4_Cutadapt/Trimmed/${BASE}_1_cutadapt.fastq"
CA2="$OUTPUT_DIR/Practical_4_Cutadapt/Trimmed/${BASE}_2_cutadapt.fastq"

echo "Running Cutadapt..."

cutadapt \
    -j "$THREADS" \
    -g AAAAAAAAAAAA \
    -G AAAAAAAAAAAA \
    -o "$CA1" \
    -p "$CA2" \
    "$R1" "$R2"

# FastQC after Cutadapt
echo "Running FastQC after Cutadapt..."

fastqc \
    "$CA1" "$CA2" \
    -o "$OUTPUT_DIR/Practical_4_Cutadapt/FastQC" \
    -t "$THREADS"

echo
echo "NGS Module 1 analysis completed."
echo "Results: $OUTPUT_DIR"
