#!/bin/bash

## header  --------------------------------------------------------------------

### the one that annotates the variants with snpeff

## settings  ------------------------------------------------------------------

snpeff_path="${1}"
snpsift_path="${2}"
full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=24
### the reference genome MUST match the reference used to map short-reads
ref_gen_ver="GRCh38.105"
### path to the dbs
dbnsfp_path="/home/shared/dbs/grch38/dbnsfp-4.1.0/dbnsfp.txt.gz"
cosmcod_path="/home/shared/dbs/grch38/cosmic-99.0.0/cosmic-cod-norm.vcf.gz"
cosmnonc_path="/home/shared/dbs/grch38/cosmic-99.0.0/cosmic-noncod-norm.vcf.gz"
dbsnp_path="/home/shared/dbs/grch38/dbsnp-155.0.0/dbsnp.vcf.gz"
### input folder
in_dir="${base_dir}/var-calls/inter-filt-unc"
### output folder
out_dir="${base_dir}/var-calls/anno-snpeff-unc"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that annotates the variants with snpeff..."

vcf_files=$(find "${in_dir}" -name "*-isec-flt-srt.vcf.gz")
pll_check=$((pll_runs + 1))
for ind_f in ${vcf_files}; do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  (
  ### get the ID of the control and the tumour samples
  file_name=$(basename "${ind_f}")
  control_id="unc"
  tumour_id=$(echo "${file_name}" | cut -d "-" -f 1)
  ### set output file and error file
  out_eff=$(echo "${file_name/%.vcf.gz/-anno.vcf}")
  err_eff=$(echo "${file_name/%.vcf.gz/-anno.log}")
  ### go to sample output folder since snpeff is dumb
  ### and does not have an output folder option
  ### (the report files are created in the current folder)
  ### jeez
  samp_out_dir="${out_dir}/${tumour_id}-vs-${control_id}"
  if [[ -d "${samp_out_dir}" ]]; then rm -rf "${samp_out_dir}"; fi
  mkdir -p "${samp_out_dir}"
  cd "${samp_out_dir}"
  ##### run the annotator
  java -Xmx8g -jar "${snpeff_path}" -v \
  "${ref_gen_ver}" \
  "${ind_f}" \
  > "${out_dir}/${out_eff}" \
  2> "${err_eff}"
  bgzip "${out_dir}/${out_eff}"
  tabix -f -p vcf "${out_dir}/${out_eff}.gz"
  
  ### set output file and error file for snpsift
  out_sift=$(echo "${file_name/%.vcf.gz/-anno-dbs.vcf}")
  err_sift=$(echo "${file_name/%.vcf.gz/-anno-dbs.log}")

  ### add snpsift annotations: dbnsfp
  java -jar "${snpsift_path}" dbnsfp -v -db "${dbnsfp_path}" \
  "${out_dir}/${out_eff}.gz" > "${out_dir}/${out_sift}.tmp1" \
  2> "${err_sift}"
  bgzip "${out_dir}/${out_sift}.tmp1"
  tabix -f -p vcf "${out_dir}/${out_sift}.tmp1.gz"

  ### a note on snpsift annotate: the manual entry for -info is not correct,
  ### the defualt is "none" while "all" is used 
  ### if the option is not invoked;
  ### the manual states:
  ### -info <list>: annotate using a list of info fields (list is a comma 
  ### separated list of fields). default: all.

  ### add snpsift annotations: cosmic coding
  java -jar "${snpsift_path}" annotate -a "${cosmcod_path}" \
  "${out_dir}/${out_sift}.tmp1.gz" > "${out_dir}/${out_sift}.tmp2" \
  2>> "${err_sift}"
  bgzip "${out_dir}/${out_sift}.tmp2"
  tabix -f -p vcf "${out_dir}/${out_sift}.tmp2.gz"

  ### add snpsift annotations: cosmic noncoding
  java -jar "${snpsift_path}" annotate -a "${cosmnonc_path}" \
  "${out_dir}/${out_sift}.tmp2.gz" > "${out_dir}/${out_sift}.tmp3" \
  2>> "${err_sift}"
  bgzip "${out_dir}/${out_sift}.tmp3"
  tabix -f -p vcf "${out_dir}/${out_sift}.tmp3.gz"

  ### add snpsift annotations: dbsnp
  java -jar "${snpsift_path}" annotate "${dbsnp_path}" \
  "${out_dir}/${out_sift}.tmp3.gz" > "${out_dir}/${out_sift}.tmp4" \
  2>> "${err_sift}"
  bgzip "${out_dir}/${out_sift}.tmp4"
  tabix -f -p vcf "${out_dir}/${out_sift}.tmp4.gz"

  mv "${out_dir}/${out_sift}.tmp4.gz" "${out_dir}/${out_sift}.gz"
  mv "${out_dir}/${out_sift}.tmp4.gz.tbi" "${out_dir}/${out_sift}.gz.tbi"

  ### cleaning temporary files
  rm -f "${out_dir}/${out_sift}.tmp1.gz" "${out_dir}/${out_sift}.tmp1.gz.tbi" 
  rm -f "${out_dir}/${out_sift}.tmp2.gz" "${out_dir}/${out_sift}.tmp2.gz.tbi" 
  rm -f "${out_dir}/${out_sift}.tmp3.gz" "${out_dir}/${out_sift}.tmp3.gz.tbi" 

  ) &
done

wait
