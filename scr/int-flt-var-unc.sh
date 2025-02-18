#!/bin/bash

## header  --------------------------------------------------------------------

### the one that intersects the unmatched small variants 

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=24

v_dir="${base_dir}/var-calls/norm-vardict-unc"
g_dir="${base_dir}/var-calls/norm-gatk-unc"

out_dir="${base_dir}/var-calls/inter-filt-unc"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

### the reference
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")
fai_path="${ref_path}.fai"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that intersects the unmatched small variants..."

cd "${out_dir}"

pll_check=$((pll_runs + 1))
for ind_v in $(find "${v_dir}" -name "*vcf.gz"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  (
  ind_g=$(echo "${ind_v}" | sed 's|vardict|gatk|')
  # ind_s=$(echo "${ind_v}" | sed 's|vardict|strelka|')
  sample_id=$(basename "${ind_v}" | sed 's|-norm.vcf.gz||')
  ### extract and write records from vardict
  ### shared by both vardict and gatk using exact allele match
  bcftools isec "${ind_v}" "${ind_g}" \
  -n =2 -w 1 -O z -o "${sample_id}-temp.vcf.gz"
  mv "${sample_id}-temp.vcf.gz" "${sample_id}-isec.vcf.gz"
  
  # ### make the header
  # zgrep "^#" "${sample_id}-isec.vcf.gz" > "${sample_id}-isec-flt.vcf"
  # ### keep only high-impact variants
  # for ind_e in ${variants_kept}; do
  #   zgrep -v "^#" "${sample_id}-isec.vcf.gz" | \
  #   grep "STATUS=${ind_e}" \
  #   >> "${sample_id}-isec-flt.vcf"
  # done
  # rm -f "${sample_id}-isec.vcf.gz"
  
  ### if no filter is applied
  mv "${sample_id}-isec.vcf.gz" "${sample_id}-isec-flt.vcf.gz"
  
  ### sorting and indexing
  # bgzip "${sample_id}-isec-flt.vcf"
  bcftools sort "${sample_id}-isec-flt.vcf.gz" \
  -O z -o "${sample_id}-isec-flt-srt.vcf.gz"
  rm -f "${sample_id}-isec-flt.vcf.gz"
  tabix -f -p vcf "${sample_id}-isec-flt-srt.vcf.gz"
  ) &
done

wait
