#!/bin/bash

## header  --------------------------------------------------------------------

### the one that calls small variants with vardict unmatched

### This implementation of vardict performs 
### the following filters with var2vcf_paired.pl:
### - minimum allele frequency (min_af) 0.20
### - only one variant per locus is allowed (-A option is missing)
### - variants that didn't pass filters will not be present in the vcf file (-S)
### - strict somatic candidates (with -M)
### - minimal coverage (-d) 10 X
### - no strand bias for low-frequency variants (details are explained below)
###
### Another level of filtering may be added by checking the strand bias
### for ALL the variants
### using the SBF tag or the BIAS tag (better the latter).
### The p-value is already calculated and encoded in the SBF tag of the FORMAT field.
### We can use the BIAS tag of the FORMAT field which
### can take values [0-2];[0-2] (i.e. "0;2", "2;1" and separator can be another in paired and single VCF).
### The first value refers to reads that support the reference allele
### and the second to reads that support the variant allele. Tha values means:
### - 0: small total count of reads (less than 12 for the sum of forward and reverse reads) 
### - 1: strand bias 
### - 2: no strand bias
###
### According to the manual a variant is emitted if:
### - the frequency of the variant exceeds the threshold set by the -f option (default = 1%), 
### - the minimum number of high-quality reads supporting variant is larger than the threshold set by the -r option (default = 2)
### - the mean position of the variant in reads is larger than the value set by the -P option (default = 5)
### - the mean base quality (phred score) for the variant is larger than the threshold set by the -q option (default = 22.5)
### - variant frequency is more than 25% or reference allele does not have much better mapping quality than the variant
### - deletion variants are not located in the regions where the reference genome is missing
### - the ratio of high-quality reads to low-quality reads is larger than the threshold specified by -o option (default=1.5)
### - variant frequency exceeds 30%; if so, next steps won't be checked and variant considered as "good", 
### otherwise the following steps will be also checked
### - the mean mapping quality exceeds the threshold set by the -O option (default: no filtering)
### - in the case of an MSI (microsatellite) region, the variant size is less than 12 nucleotides for the non-monomer MSI
### or 15 for the monomer MSI; variant frequency is more than 10% for the non-monomer MSI (or set by --nmfreq option) 
### and 25% for the monomer MSI (or set by --mfreq option)
### - the variant has not "2;1" bias or variant frequency more than 20%; if both conditions aren't met, 
### then variant mustn't be SNV and any of variants REF allele 
### or ALT allele lengths must be more than 3 nucleotides

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=24
pll_runs=2
ref_name="${1}"
read -a popu_samp <<< "${2}"

### dev off
# full_dir="/home/ltattini/ems/plato/prog/zebruno/rerio-matched/scr"
# base_dir="/home/ltattini/ems/plato/prog/zebruno/rerio-matched"
# n_threads=22
# pll_runs=2
# ref_name="rerio"
# popu_samp=('t6' 't12' 't13')
# popu_samp=('t6')

### output folder
log_dir="${base_dir}/scr/logs"
out_dir="${base_dir}/var-calls/vardict-unc"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

### folder for regions bed file [only for WGS data]
# dir_bed="${base_dir}/var-calls/reg-bed"
# if [[ -d "${dir_bed}" ]]; then rm -rf "${dir_bed}"; fi
# mkdir -p "${dir_bed}"
### the reference
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")
fai_path="${ref_path}.fai"

### parameters
min_af="0.20"

### increase memory for java
export _JAVA_OPTIONS="-Xms16g -Xmx48g"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that calls small variants with vardict unmatched..."

## make the bed (only 3 columns) for callable regions -------------------------

### if the bed file has more than 3 column vardict exits with an error

### remarkably, if "${dir_bed}/call-regions.bed" ends with an empty line (ergo also "\n")
### we get an error:
### Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length 1
###         at com.astrazeneca.vardict.RegionBuilder.buildRegions(RegionBuilder.java:65)
###         at com.astrazeneca.vardict.VarDictLauncher.initResources(VarDictLauncher.java:102)
###         at com.astrazeneca.vardict.VarDictLauncher.start(VarDictLauncher.java:49)
###         at com.astrazeneca.vardict.Main.main(Main.java:15)

### raw bed file
raw_bed=$(find "${base_dir}/cpk" -name "*intervals.bed")
tmp_bed="${raw_bed}.tmp"
if [[ -f "${tmp_bed}" ]]; then rm -f "${tmp_bed}"; fi

### search for bed contigs in the reference fasta
ref_chr=$(grep "^>" "${ref_path}" | cut -d " " -f 1 | sed 's|^.||')
for ind_c in ${ref_chr}; do
  grep -w "^${ind_c}" "${raw_bed}" >> "${tmp_bed}"
done

### make bed file with padding (100 bp)
intervals_file="${raw_bed%.bed}-pad100.bed"
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2-100,$3+100}' "${tmp_bed}" \
> "${intervals_file}"
rm -f "${tmp_bed}"

## unmatched tumor calls ------------------------------------------------------

cd "${base_dir}/map-sr"

seq_dim=$(echo "${#popu_samp[@]}")
pll_check=$((pll_runs + 1))
for (( ind_i=0; ind_i<seq_dim; ind_i++ )); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  (
  ### call an unmatched sample
  popu_bam="${popu_samp[ind_i]}-${ref_name}-srt-mdp.bam"
  tmp_vcf="${out_dir}/${popu_samp[ind_i]}-tmp.vcf"
  out_vcf_gz="${out_dir}/${popu_samp[ind_i]}.vcf.gz"
  vardict -G "${ref_path}" -f "${min_af}" \
  -th "${n_threads}" \
  -b "${popu_bam}" \
  -N "${popu_samp[ind_i]}" --nosv \
  --fisher \
  -C -c 1 -S 2 -E 3 "${intervals_file}" | \
  var2vcf_valid.pl -N "${popu_samp[ind_i]}" -E -f "${min_af}" \
  > "${tmp_vcf}"
  ### add the contig info in the header, sort and compress
  bcftools reheader --fai "${ref_path}.fai" "${tmp_vcf}" \
  > "${tmp_vcf}.antani"
  bcftools sort "${tmp_vcf}.antani" \
  -O z -o "${out_vcf_gz}"
  ### index for the normalisation step
  tabix -f -p vcf "${out_vcf_gz}"
    
  ### clean temporary files
  rm -f "${tmp_vcf}"
  rm -f "${tmp_vcf}.antani"
  ) &
done

wait
