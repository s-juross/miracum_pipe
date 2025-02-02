#!/usr/bin/env bash

# This function gets configuration parameters from yaml files
# first trying a look-up inside the patient specific config and
# falls back to the global config if requested config key not available
#
#   usage: get_config_value some.key relative_patient_path
#
function get_config_value()
{
  SCRIPT_PATH=$(
    cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
    pwd -P
  )
  #echo ${DIR_SCRIPT}
  #echo ${2}
  #echo ${1}
  # first look inside the patient specific config
  patient_conf="${SCRIPT_PATH}"/assets/input/"${2}"/patient.yaml
  if [[ -f "${patient_conf}" ]]; then
    value=$(shyaml get-value "${1}" < "${patient_conf}" 2> /dev/null)
  else
    echo "no patient config in ${2}"
    exit 1
  fi

  # if value not available in patient specific config, look into custom config
  if [[ -z "${value}" ]]; then
    value=$(shyaml get-value "${1}" < "${SCRIPT_PATH}"/conf/custom.yaml 2> /dev/null)

# if value not available in custom config, take global config's value
    if [[ -z "${value}" ]]; then
      value=$(shyaml get-value "${1}" < "${SCRIPT_PATH}"/conf/default.yaml 2> /dev/null)
    fi
  fi

  echo "${value}"
}
export -f get_config_value

function join_by { local IFS="$1"; shift; echo "$*"; }

export -f join_by

## Directories
## General
readonly DIR_ASSETS="${DIR_SCRIPT}/assets"
readonly DIR_TOOLS="${DIR_SCRIPT}/tools"
readonly DIR_INPUT="${DIR_ASSETS}/input"            # folder containing the raw data (.fastq files)
readonly DIR_OUTPUT="${DIR_ASSETS}/output"
readonly DIR_REF="${DIR_ASSETS}/references"         # reference genome
readonly DIR_CHROMOSOMES="${DIR_REF}/chromosomes"
readonly DIR_SEQUENCING="${DIR_REF}/sequencing"
readonly DIR_RSCRIPT="${DIR_SCRIPT}/RScripts"
readonly DIR_DATABASE="${DIR_SCRIPT}/databases"
readonly DIR_DBSNP="${DIR_DATABASE}/dbSNP"


export DIR_ASSETS
export DIR_TOOLS
export DIR_INPUT
export DIR_OUTPUT
export DIR_REF
export DIR_CHROMOSOMES
export DIR_DBSNP
export DIR_DATABASE
export DIR_RSCRIPT
export DIR_SEQUENCING
