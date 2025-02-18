#!/bin/bash

## header  --------------------------------------------------------------------

### the one that makes the statistics of the coverage
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=8
intervals_file="$(find "${base_dir}/cpk" -name "*pad100.bed")"

### output folder
out_dir="${base_dir}/cov"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that makes the statistics of the coverage..."

cd "${base_dir}"

pll_check=$((pll_runs + 1))
for ind_m in $(ls "map-sr/"*"bam"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  (
  map_name=$(echo "${ind_m}" | cut -d "/" -f 2 | cut -d "-" -f 1,2)
  samtools depth "${ind_m}" > "${out_dir}/${map_name}-depth.txt"
  ref_len=$(samtools view -H "${ind_m}" | \
  grep "^@SQ" | cut -f 3 | cut -d ":" -f 2 | \
  awk '{sum += $1} END {print sum}')
  awk -v RL="${ref_len}" 'BEGIN {FS="\t"} {sum+=$3; sumsq+=$3*$3} \
  END {print "Covered sites = " NR; \
  print "Contig covered (%) = " 100*NR/RL; \
  print "Contig size = " RL; \
  print "Mean (depth) =", sum/NR; \
  print "SD (depth) =", sqrt(sumsq/NR - (sum/NR)**2)}' \ 
  "${out_dir}/${map_name}-depth.txt" > "${out_dir}/${map_name}-wg-mean.txt"
  
  samtools depth "${ind_m}" \
  -b "${intervals_file}" | \
  awk 'BEGIN {FS="\t"} {sum+=$3; sumsq+=$3*$3} \
  END {print "Covered sites: " NR; \
  print "Mean (depth): " sum/NR; \
  print "SD (depth): " sqrt(sumsq/NR - (sum/NR)**2)}' \
  > "${out_dir}/${map_name}-ex-mean.txt"
  ) &
done

wait

rm -f "cov/"*"depth.txt"
