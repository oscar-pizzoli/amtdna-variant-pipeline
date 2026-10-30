#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="/home/pippo/amtdna-variant-pipeline/data/raw_fastq"
mkdir -p "${OUT_DIR}"

# 5 campioni moderni 1000 Genomes con aplotipi differenti (CEU, YRI, CHS, GBR, ASW)
declare -A SAMPLES=(
    ["NA12878"]="CEU:NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
    ["NA19240"]="YRI:NA19240.mapped.ILLUMINA.bwa.YRI.low_coverage.20121211.bam"
    ["HG00514"]="CHS:HG00514.mapped.ILLUMINA.bwa.CHS.low_coverage.20121211.bam"
    ["HG00096"]="GBR:HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20121211.bam"
    ["NA19700"]="ASW:NA19700.mapped.ILLUMINA.bwa.ASW.low_coverage.20121211.bam"
)

BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"

echo "================================================================================"
echo "📥 DOWNLOAD IN STREAMING DEI FASTQ MITOCONDRIALI (5 CAMPIONI MODERNI)"
echo "================================================================================"

for SAMPLE in "${!SAMPLES[@]}"; do
    INFO="${SAMPLES[$SAMPLE]}"
    POP="${INFO%%:*}"
    BAM_NAME="${INFO##*:}"
    REMOTE_BAM="${BASE_URL}/${SAMPLE}/alignment/${BAM_NAME}"
    
    TMP_BAM="${OUT_DIR}/${SAMPLE}_mt.bam"
    FASTQ_R1="${OUT_DIR}/${SAMPLE}_R1.fastq.gz"
    FASTQ_R2="${OUT_DIR}/${SAMPLE}_R2.fastq.gz"

    echo ""
    echo "[*] Processing ${SAMPLE} (${POP})..."
    echo "    -> Streaming reads da: ${REMOTE_BAM}"

    # 1. Estrazione remota selettiva del solo contig MT (o chrM)
    if ! samtools view -b -h "${REMOTE_BAM}" MT -o "${TMP_BAM}" 2>/dev/null; then
        echo "    [i] Retry con 'chrM'..."
        samtools view -b -h "${REMOTE_BAM}" chrM -o "${TMP_BAM}"
    fi

    # 2. Conversione da BAM a coppie di FASTQ ordinate per nome
    echo "    -> Conversione in FASTQ compressi..."
    samtools collate -u -O "${TMP_BAM}" | \
        samtools fastq \
            -1 "${FASTQ_R1}" \
            -2 "${FASTQ_R2}" \
            -0 /dev/null \
            -s /dev/null \
            -n -

    # 3. Pulizia file temporaneo
    rm -f "${TMP_BAM}"
    echo "    ✅ ${SAMPLE} completato: $(du -h "${FASTQ_R1}" | cut -f1) (R1) / $(du -h "${FASTQ_R2}" | cut -f1) (R2)"
done

echo ""
echo "================================================================================"
echo "🎉 TUTTI I 5 CAMPIONI MODERNI SONO STATI SCARICATI IN: ${OUT_DIR}"
echo "================================================================================"
ls -lh "${OUT_DIR}"/*.fastq.gz
