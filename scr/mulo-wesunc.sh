#!/bin/bash

## header ---------------------------------------------------------------------

### this script is the mulo-wesunc runner

## user's settings ------------------------------------------------------------

ref_name="grch38-p14"

### run-1
popu_samp="SBS01 \
SBS02 \
SBS03 \
SBS04 \
SBS05 \
SBS06 \
SBS07 \
SBS08 \
SBS09 \
SBS10 \
SBS11 \
SBS12 \
SBS13 \
SBS14 \
SBS15 \
SBS16 \
SBS17 \
SBS18 \
SBS19 \
SBS20 \
SBS21 \
SBS22 \
SBS23 \
SBS24 \
SBS25 \
SBS26 \
SBS27 \
SBS28 \
SBS29 \
SBS30"

## system's settings ----------------------------------------------------------

### check logs folder
if [[ ! -d "logs" ]]; then mkdir "logs"; fi

### set path to snpeff.jar
snpeff_path="/home/tools/lib/snpeff/snpeff-5.2.0/snpEff.jar"
snpsift_path="/home/tools/lib/snpeff/snpeff-5.2.0/SnpSift.jar"

### path to gatk panel of normal samples
pon_path="/home/shared/dbs/grch38/gatk-pon/1000g-pon-hg38.vcf.gz"

## clmnt ----------------------------------------------------------------------

# ### quality check
# /usr/bin/time -v bash fq-check.sh \
# > "logs/fq-check.out" 2> "logs/fq-check.err" &

# ### reference indexing
# /usr/bin/time -v bash index-ref.sh "${ref_name}" \
# > "logs/index-ref.out" 2> "logs/index-ref.err"

# ### mapping
# /usr/bin/time -v bash map-sr.sh "${ref_name}" "${popu_samp}" \
# > "logs/map-sr.out" 2> "logs/map-sr.err"

# ### coverage statistics
# /usr/bin/time -v bash depth-stats.sh \
# > "logs/depth-stats.out" 2> "logs/depth-stats.err" &

# ### vardict sets the filters and reads the padded bed done before
# /usr/bin/time -v bash call-vardict-unc.sh "${ref_name}" "${popu_samp}" \
# > "logs/call-vardict-unc.out" 2> "logs/call-vardict-unc.err"

# ### calling with gatk
# /usr/bin/time -v bash call-gatk-unc.sh "${ref_name}" \
# "${popu_samp}" "${pon_path}" \
# > "logs/call-gatk-unc.out" 2> "logs/call-gatk-unc.err"

# ### normalise and de-duplicate the variants
# /usr/bin/time -v bash norm-var-unc.sh \
# > "logs/norm-var-unc.out" 2> "logs/norm-var-unc.err"

### intersect and filter
/usr/bin/time -v bash int-flt-var-unc.sh \
> "logs/int-flt-var-unc.out" 2> "logs/int-flt-var-unc.err"

### annotation with snpeff
/usr/bin/time -v bash anno-snpeff-unc.sh "${snpeff_path}" "${snpsift_path}" \
> "logs/anno-snpeff-unc.out" 2> "logs/anno-snpeff-unc.err"

### call copy-number variants

echo "Prematura la supercazola o scherziamo?"
echo "[Conte Lello Mascetti]"
