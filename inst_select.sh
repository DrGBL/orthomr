#!/bin/bash

# inst_select.sh - Select genetic instruments for non-linear MR analysis with dynamic column detection + optional variant restriction

# Function to display usage
usage() {
    cat << EOF
Usage: inst_select.sh [OPTIONS]

Required arguments:
  --path_summ_stats PATH        Path to directory containing summary statistics files
  --summ_stats_prefix PREFIX    Prefix for summary statistics files
  --path_tmp PATH               Path to temporary directory
  --path_ref PATH               Path to reference genotype data
  --ref_prefix PREFIX           Prefix for reference plink files
  --path_clumps PATH            Path to output clumps directory
  --max_degree INT              Maximum polynomial degree to process
  --col_chr COLNAME             Column name for chromosome in summary statistics
  --col_snp COLNAME             Column name for SNP (variant ID) in summary statistics
  --col_bp COLNAME              Column name for base pair position
  --col_a1 COLNAME              Column name for effect allele
  --col_p COLNAME               Column name for p-value

Optional arguments:
  --restrict_variants FILE      Path to a file containing variant IDs (1 column, no header) to restrict to
  --clump_p1 FLOAT              Clumping p-value threshold 1 (default: 1e-8)
  --clump_p2 FLOAT              Clumping p-value threshold 2 (default: 1e-8)
  --clump_r2 FLOAT              Clumping r² threshold (default: 0.0001)
  --clump_kb INT                Clumping kb window (default: 500)
  --sig_threshold FLOAT         Significance threshold for variant selection (default: 5e-8)
  --threads INT                 Number of threads for plink (default: 2)
  --memory INT                  Memory for plink in MB (default: 20000)
  --path_plink PATH             Path to directory containing plink binary (optional)
  -h, --help                    Show this help message

EOF
    exit 1
}

# Default values
clump_p1=1e-8
clump_p2=1e-8
clump_r2=0.0001
clump_kb=500
sig_threshold=5e-8
threads=2
memory=20000
path_plink=""
restrict_variants=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --path_summ_stats) path_summ_stats="$2"; shift 2 ;;
        --summ_stats_prefix) summ_stats_prefix="$2"; shift 2 ;;
        --path_tmp) path_tmp="$2"; shift 2 ;;
        --path_ref) path_ref="$2"; shift 2 ;;
        --ref_prefix) ref_prefix="$2"; shift 2 ;;
        --path_clumps) path_clumps="$2"; shift 2 ;;
        --max_degree) max_degree="$2"; shift 2 ;;
        --col_chr) col_chr="$2"; shift 2 ;;
        --col_snp) col_snp="$2"; shift 2 ;;
        --col_bp) col_bp="$2"; shift 2 ;;
        --col_a1) col_a1="$2"; shift 2 ;;
        --col_p) col_p="$2"; shift 2 ;;
        --restrict_variants) restrict_variants="$2"; shift 2 ;;
        --clump_p1) clump_p1="$2"; shift 2 ;;
        --clump_p2) clump_p2="$2"; shift 2 ;;
        --clump_r2) clump_r2="$2"; shift 2 ;;
        --clump_kb) clump_kb="$2"; shift 2 ;;
        --sig_threshold) sig_threshold="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --memory) memory="$2"; shift 2 ;;
        --path_plink) path_plink="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Check required arguments
if [[ -z "$path_summ_stats" || -z "$summ_stats_prefix" || -z "$path_tmp" || -z "$path_ref" || -z "$ref_prefix" || -z "$path_clumps" || -z "$max_degree" || -z "$col_chr" || -z "$col_snp" || -z "$col_bp" || -z "$col_a1" || -z "$col_p" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Create directories
mkdir -p "$path_tmp" "$path_clumps"

echo "Starting instrument selection..."
echo "Restriction file: ${restrict_variants:-None}"

# Create header file
printf "CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP\n" > ${path_tmp}header.txt

##############################################
# Detect column indices from first file
##############################################

first_file=$(ls ${path_summ_stats}/${summ_stats_prefix}_chr1_deg_1.tsv.gz 2>/dev/null | head -n 1)
if [[ -z "$first_file" ]]; then
  echo "Error: No input files found for dynamic column detection."
  exit 1
fi

zcat "$first_file" | head -n 1 > ${path_tmp}tmp_header.txt

col_chr_idx=$(awk -v col="$col_chr" '{for(i=1;i<=NF;i++) if($i==col) print i}' ${path_tmp}tmp_header.txt)
col_snp_idx=$(awk -v col="$col_snp" '{for(i=1;i<=NF;i++) if($i==col) print i}' ${path_tmp}tmp_header.txt)
col_bp_idx=$(awk -v col="$col_bp" '{for(i=1;i<=NF;i++) if($i==col) print i}' ${path_tmp}tmp_header.txt)
col_a1_idx=$(awk -v col="$col_a1" '{for(i=1;i<=NF;i++) if($i==col) print i}' ${path_tmp}tmp_header.txt)
col_p_idx=$(awk -v col="$col_p" '{for(i=1;i<=NF;i++) if($i==col) print i}' ${path_tmp}tmp_header.txt)

rm ${path_tmp}tmp_header.txt

if [[ -z "$col_chr_idx" || -z "$col_snp_idx" || -z "$col_bp_idx" || -z "$col_a1_idx" || -z "$col_p_idx" ]]; then
  echo "Error: One or more required columns not found."
  exit 1
fi

echo "Detected column indices:"
echo "CHR=$col_chr_idx SNP=$col_snp_idx BP=$col_bp_idx A1=$col_a1_idx P=$col_p_idx"

##############################################
# If restricting variants, preprocess the file
##############################################

if [[ -n "$restrict_variants" ]]; then
    echo "Preparing variant restriction set..."
    sort -u "$restrict_variants" > ${path_tmp}restrict_variants_sorted.txt
fi

##############################################
# Filter significant variants (+ optional restrict)
##############################################

echo "Filtering significant variants..."
for chr in {1..22}; do
  for deg in $(seq 1 $max_degree); do

    input_file="${path_summ_stats}/${summ_stats_prefix}_chr${chr}_deg_${deg}.tsv.gz"

    if [[ -f "$input_file" ]]; then

      zcat "$input_file" \
      | awk -v p_idx=$col_p_idx -v thresh=$sig_threshold 'NR==1{next} $p_idx < thresh' \
      | awk -v chr_idx=$col_chr_idx -v snp_idx=$col_snp_idx -v bp_idx=$col_bp_idx -v a1_idx=$col_a1_idx -v p_idx=$col_p_idx \
            'BEGIN{OFS="\t"} {print $chr_idx,$snp_idx,$bp_idx,$a1_idx,"TEST","NMISS","BETA","STAT",$p_idx}' \
      > ${path_tmp}tmp_preclump_unsorted.txt

      # Apply restriction if requested
      if [[ -n "$restrict_variants" ]]; then
          awk 'NR==FNR {keep[$1]; next} $2 in keep' \
              ${path_tmp}restrict_variants_sorted.txt \
              ${path_tmp}tmp_preclump_unsorted.txt \
              > ${path_tmp}tmp_preclump_filtered.txt
          mv ${path_tmp}tmp_preclump_filtered.txt ${path_tmp}tmp_preclump_unsorted.txt
      fi

      cat ${path_tmp}header.txt ${path_tmp}tmp_preclump_unsorted.txt \
        | gzip > ${path_tmp}${summ_stats_prefix}_chr${chr}_deg_${deg}_pre_clump.txt.gz

      rm ${path_tmp}tmp_preclump_unsorted.txt
    fi
  done
done

rm ${path_tmp}header.txt

##############################################
# Clumping
##############################################

echo "Performing clumping..."
for chr in {1..22}; do
  for deg in $(seq 1 $max_degree); do

    pre="${path_tmp}${summ_stats_prefix}_chr${chr}_deg_${deg}_pre_clump.txt.gz"

    if [[ -f "$pre" ]]; then
      zcat "$pre" | tail -n+2 | awk '{print $2}' > ${path_clumps}tmp_variants.txt

      n_line=$(wc -l < ${path_clumps}tmp_variants.txt)

      if (( n_line > 0 )); then
        ${path_plink}/plink \
          --bfile ${path_ref}/${ref_prefix}_chr${chr} \
          --extract ${path_clumps}tmp_variants.txt \
          --clump "$pre" \
          --clump-p1 ${clump_p1} \
          --clump-p2 ${clump_p2} \
          --clump-r2 ${clump_r2} \
          --clump-kb ${clump_kb} \
          --threads ${threads} \
          --memory ${memory} \
          --out ${path_clumps}${summ_stats_prefix}_chr${chr}_deg_${deg}_clumps
      fi

      rm ${path_clumps}tmp_variants.txt
    fi
  done
done

##############################################
# Collect lead variants
##############################################

echo "Collecting lead variants..."
cat ${path_clumps}${summ_stats_prefix}_chr*_deg_*_clumps.clumped 2>/dev/null \
  | grep -v '^$' | awk '!/SNP/' | awk '{print $3}' \
  | sort -t ':' -k1.4n,1.5 -k2,2n | uniq \
  > ${path_clumps}all_lead_variants.txt

##############################################
# Final LD pruning
##############################################

echo "Final LD pruning..."
> ${path_clumps}final_all_lead_variants.txt

for chr in {1..22}; do
  grep "chr${chr}:" ${path_clumps}all_lead_variants.txt > ${path_clumps}tmp_lead_variants.txt

  if [[ -s ${path_clumps}tmp_lead_variants.txt ]]; then
    ${path_plink}/plink \
      --bfile ${path_ref}/${ref_prefix}_chr${chr} \
      --extract ${path_clumps}tmp_lead_variants.txt \
      --indep-pairwise 500 1000 ${clump_r2} \
      --out ${path_clumps}tmp_lead_variants

    cat ${path_clumps}tmp_lead_variants.prune.in >> ${path_clumps}final_all_lead_variants.txt
    rm ${path_clumps}tmp_lead_variants*
  fi
done

##############################################
# DONE
##############################################

echo "Instrument selection complete!"
echo "Final variants: ${path_clumps}final_all_lead_variants.txt"
wc -l < ${path_clumps}final_all_lead_variants.txt
