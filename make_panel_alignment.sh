#!/usr/bin/env bash

# required params: task

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=common.cfg.sh
. "${DIR_SCRIPT}"/common.cfg.sh

# for now only "tumor only" panels
readonly VALID_TASKS=("td")
readonly VALID_SEXES=("XX XY")

function usage() {
  echo "usage: miracum_pipe.sh -d dir [-h]"
  echo "  -d  dir             specify relative folder of patient"
  echo "  -t  task            specify task: $(join_by ' ' ${VALID_TASKS})"
  echo "  -p                  computing as parallel process"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:t:ph option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  t) readonly PARAM_TASK=$OPTARG ;;
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


# check inputs
if [[ ! " ${VALID_TASKS[@]} " =~ " ${PARAM_TASK} " ]]; then
  echo "unknown task: ${PARAM_TASK}"
  echo "use one of the following values: $(join_by ' ' ${VALID_TASKS})"
  exit 1
fi

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

[[ -d "${DIR_WES}" ]] || mkdir -p "${DIR_WES}"

# SAMPLE
readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_${PARAM_TASK}

# temp files
readonly fastq_o1_p_t=${DIR_TMP}/${NameD}_output1_paired_trimmed.fastq.gz
readonly fastq_o1_u_t=${DIR_TMP}/${NameD}_output1_unpaired_trimmed.fastq.gz
readonly fastq_o2_p_t=${DIR_TMP}/${NameD}_output2_paired_trimmed.fastq.gz
readonly fastq_o2_u_t=${DIR_TMP}/${NameD}_output2_unpaired_trimmed.fastq.gz

readonly bam=${DIR_TMP}/${NameD}_output.bam
readonly prefixsort=${DIR_TMP}/${NameD}_output.sort
readonly sortbam=${DIR_TMP}/${NameD}_output.sort.bam
readonly rmdupbam=${DIR_TMP}/${NameD}_output.sort.rmdup.bam
readonly bai=${DIR_TMP}/${NameD}_output.sort.rmdup.bai
readonly bamlist=${DIR_TMP}/${NameD}_output.sort.rmdup.bam.list
readonly realignedbam=${DIR_TMP}/${NameD}_output.sort.rmdup.realigned.bam
readonly realignedbai=${DIR_TMP}/${NameD}_output.sort.rmdup.realigned.bai
readonly fixedbam=${DIR_TMP}/${NameD}_output.sort.rmdup.realigned.fixed.bam
readonly fixedbai=${DIR_TMP}/${NameD}_output.sort.rmdup.realigned.fixed.bai
readonly csv=${DIR_TMP}/${NameD}_output.sort.rmdup.realigned.fixed.recal_data.csv

# keep
readonly recalbam=${DIR_WES}/${NameD}_output.sort.rmdup.realigned.fixed.recal.bam
readonly statstxt=${DIR_WES}/${NameD}_stats.txt
readonly coveragetxt=${DIR_WES}/${NameD}_coverage.all.txt
readonly coverageexons=${DIR_WES}/${NameD}_coverage.exons.txt

# Panel Files are spread over multiple files, e.g. 4
# filenames without extension in yaml file!
if [[ ${CFG_PANEL_FILE_NUMBER} -gt 1 ]]
then
  file_numbers=$(seq -s "" 1 ${CFG_PANEL_FILE_NUMBER})
  # tmp merged fastq
  readonly merged_fq1=${DIR_TMP}/${NameD}_merged_R1.fastq.gz
  readonly merged_fq2=${DIR_TMP}/${NameD}_merged_R2.fastq.gz

  if [[ ! -f ${merged_fq1} ]] 
  then
    cat "${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_PANEL_FILE_TUMOR}"_R1.fastq.gz > "${merged_fq1}" || exit 1;
  fi

  if [[ ! -f ${merged_fq2} ]]
  then
    cat "${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_PANEL_FILE_TUMOR}"_R2.fastq.gz > "${merged_fq2}" || exit 1;
  fi
else
  readonly merged_fq1=${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_PANEL_FILE_TUMOR}_R1.fastq.gz
  readonly merged_fq2=${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_PANEL_FILE_TUMOR}_R2.fastq.gz
fi

# fastqc zip to WES
${BIN_FASTQC} "${merged_fq1}" -o "${DIR_WES}"
${BIN_FASTQC} "${merged_fq2}" -o "${DIR_WES}"

# trim fastq
${BIN_TRIM} "${merged_fq1}" "${merged_fq2}" "${fastq_o1_p_t}" "${fastq_o1_u_t}" "${fastq_o2_p_t}" "${fastq_o2_u_t}" \
ILLUMINACLIP:"${DIR_TRIMMOMATIC_ADAPTER}"/TruSeq3-PE-2.fa:2:30:10 HEADCROP:3 TRAILING:10 MINLEN:70

# fastqc
${BIN_FASTQC} "${fastq_o1_p_t}" -o "${DIR_WES}"
${BIN_FASTQC} "${fastq_o2_p_t}" -o "${DIR_WES}"

# make bam
${BIN_BWAMEM} -R "@RG\tID:${NameD}\tSM:${NameD}\tPL:illumina\tLB:lib1\tPU:unit1" -t "${CFG_COMMON_CPUCORES}" "${FILE_GENOME}" \
"${fastq_o1_p_t}" "${fastq_o2_p_t}" | ${BIN_SAMVIEW} -bS - >"${bam}"

# stats
${BIN_STATS} "${bam}" >"${statstxt}"

# sort bam
${BIN_SAMSORT} "${bam}" -T "${prefixsort}" -o "${sortbam}"

# rmdup bam
${BIN_SAMVIEW} -b -f 0x2 -q "${CFG_SAMTOOLS_MPILEUP_MINMQ}" "${sortbam}" | ${BIN_SAMRMDUP} - "${rmdupbam}"

# make bai
${BIN_SAMINDEX} "${rmdupbam}" "${bai}"

# make bam list
${BIN_REALIGNER_TARGER_CREATOR} -o "${bamlist}" -I "${rmdupbam}"

# realign bam
${BIN_INDEL_REALIGNER} -I "${rmdupbam}" -targetIntervals "${bamlist}" -o "${realignedbam}" # not working with panels due to the many reads per position

# fix bam
${BIN_FIX_MATE} -INPUT "${realignedbam}" -OUTPUT "${fixedbam}" -SO coordinate -VALIDATION_STRINGENCY LENIENT -CREATE_INDEX true

# make csv
${BIN_BASE_RECALIBRATOR} -I "${fixedbam}" \
-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o "${csv}"

# recal bam
${BIN_PRINT_READS} -I "${fixedbam}" -BQSR "${csv}" -o "${recalbam}"

# coverage
${BIN_COVERAGE} -b "${recalbam}" -a "${CFG_REFERENCE_CAPTUREREGIONS}" | grep '^all' >"${coveragetxt}"

# advanced qc / coverage of exonic regions
${BIN_COVERAGE} -b "${recalbam}" -a "${CFG_REFERENCE_COVERED_EXONS}" | grep '^all' > "${coverageexons}"

# fastqc
${BIN_FASTQC} "${recalbam}" -o "${DIR_WES}"
