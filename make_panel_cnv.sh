#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=common.cfg.sh
. "${DIR_SCRIPT}"/common.cfg.sh


function usage() {
  echo "usage: miracum_pipe.sh -d dir [-h]"
  echo "  -d  dir             specify relative folder of patient"
  echo "  -p                  computing as parallel process"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:ph option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  p) readonly PARALLEL_PROCESSES=2 ;;
  h) usage ;;
  \?)
    echo "Unknown option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Missing option argument for -$OPTARG" >&2
    exit 1
    ;;
  *)
    echo "Unimplemented option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# if no patient is defined
if [[ -z "${PARAM_DIR_PATIENT}" ]]; then
  echo "no patient defined."
  echo "--"
  usage
fi


# load patient yaml
readonly CFG_SEX=$(get_config_value sex "${PARAM_DIR_PATIENT}")
if [[ "$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")" = "panel" ]]; then
  readonly CFG_CASE=panelTumor
fi

echo ${CFG_CASE}

# check inputs
readonly VALID_SEXES=("XX XY")

if [[ ! " ${VALID_SEXES[@]} " =~ " ${CFG_SEX} " ]]; then
  echo "unknown sex: ${CFG_SEX}"
  echo "use one of the following values: $(join_by ' ' ${VALID_SEXES})"
  exit 1
fi

##################################################################################################################

## load programs
# shellcheck source=programs.cfg.sh
. "${DIR_SCRIPT}/programs.cfg.sh"

##################################################################################################################

# CNV Calling with ControlFREEC

# names
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td
readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}

readonly DIR_CNV_OUTPUT="${DIR_WES}/CNV"
readonly DIR_SEQUENZA="${DIR_WES}/${NameD}_sequenza"

[[ -d "${DIR_CNV_OUTPUT}" ]] || mkdir -p "${DIR_CNV_OUTPUT}"

# create sequenza directory if necessary
[[ -d "${DIR_SEQUENZA}" ]] || mkdir -p "${DIR_SEQUENZA}"

# keep
readonly bam=${DIR_WES}/${NameTD}_output.sort.rmdup.realigned.fixed.recal.bam
readonly HRD_OUTPUT=${DIR_WES}/${NameD}.seqz.gz
readonly HRD_OUTPUT_SMALL=${DIR_WES}/${NameD}.small.seqz.gz
readonly HRD_MODEL=${DIR_WES}/${NameD}_sequenza/${NameD}_alternative_solutions.txt
touch ${HRD_MODEL}

# HRD
if [ ! -f "${HRD_REF_WIG}" ]; then
    echo "${HRD_REF_WIG} does not exist. Generating ..."
    ${SEQUENZA_UTILS} gc_wiggle --fasta "${FILE_GENOME}" -w "${SEQUENZA_WINDOW}" -o "${HRD_REF_WIG}"
fi

${SEQUENZA_UTILS} bam2seqz -S "${BIN_SAMTOOLS}" -gc "${HRD_REF_WIG}" --fasta "${FILE_GENOME}" -n "${bam}" --tumor "${bam}" --normal2 "${SEQUENZA_NON_MATCHING_NORMAL}" \
  -C ${SEQUENZA_CHROMOSOMES} -o "${HRD_OUTPUT}"
${SEQUENZA_UTILS} seqz_binning -s "${HRD_OUTPUT}" -w "${SEQUENZA_WINDOW}" -o "${HRD_OUTPUT_SMALL}"
${BIN_RSCRIPT} "${DIR_RSCRIPT}/HRD.R" "${NameD}" "${DIR_WES}" "${HRD_MODEL}"

# Extract ploidy and purity from sequenza
PLOIDY=$(${BIN_RSCRIPT} "--vanilla" "-e" "cat(round(read.delim('${HRD_MODEL}', header = T)[1, 2]))")
CONTAMINATION=$(${BIN_RSCRIPT} "--vanilla" "-e" "cat(1-read.delim('${HRD_MODEL}', header = T)[1, 1])")

# build config file
cat >"${DIR_WES}"/CNV_config.txt <<EOI
[general]

chrFiles = ${DIR_CHROMOSOMES}
chrLenFile = ${CFG_REFERENCE_LENGTH}
breakPointType = 4
breakPointThreshold = 1.2
forceGCcontentNormalization = 1
gemMappabilityFile = ${FILE_REFERENCE_MAPPABILITY}
intercept = 1
minCNAlength = 3
maxThreads = ${CFG_COMMON_CPUCORES}
noisyData = TRUE
outputDir = ${DIR_CNV_OUTPUT}
ploidy = 1,2,3,4
contamination = 0
contaminationAdjustment = TRUE
printNA = FALSE
readCountThreshold = 50
samtools = ${BIN_SAMTOOLS}
bedtools = ${DIR_TOOLS}/bedtools2/bin/bedtools
sex = ${CFG_SEX}
step = 0
window = 0
uniqueMatch = TRUE
BedGraphOutput = TRUE

[sample]

mateFile = ${bam}
inputFormat = BAM
mateOrientation = FR

[BAF]

makePileup = ${CFG_REFERENCE_DBSNP}
fastaFile = ${FILE_GENOME}
minimalCoveragePerPosition = ${CFG_PANEL_VARSCAN_MPILEUP2SNP_MINCOVERAGE}
minimalQualityPerPosition = ${CFG_GENERAL_MINBASEQUAL}
SNPfile = ${CFG_REFERENCE_DBSNP}

[target]

captureRegions = ${CFG_REFERENCE_CAPTUREREGIONS}
EOI

# run CNV calling
export PATH="${BIN_SAMTOOLS}:${PATH}"
${BIN_FREEC} -conf "${DIR_WES}"/CNV_config.txt

# nun Uebertragen von purity und ploidy in alternative_file:
if [ $(wc -l < ${HRD_MODEL}) -eq 1 ]; then
	readonly info=${DIR_CNV_OUTPUT}/${NameTD}_output.sort.rmdup.realigned.fixed.recal.bam_info.txt
	readonly output_ploidy=$(grep "Output_Ploidy" ${info} | cut -f 2)
	readonly sample_purity=$(grep "Sample_Purity" ${info} | cut -f 2)
	printf "${sample_purity}\t${output_ploidy}\tNA" >> ${HRD_MODEL}
fi
