#!/usr/bin/env bash
#SCRIPT_PATH$(
#    cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
#    pwd -P
#)
## load settings
# shellcheck source=common.cfg.sh
#. "${DIR_SCRIPT}"/common.cfg.sh
##################################################################################################################
## load programs
# shellcheck source=programs.cfg.sh
#. "${DIR_SCRIPT}/programs.cfg.sh"
##################################################################################################################
# fill covered exon file
tools/bedtools2/bin/bedtools intersect -u -a assets/references/sequencing/V5UTR.bed -b assets/references/sequencing/V5UTR_Targets.txt > assets/references/sequencing/exons_5UTR.bed