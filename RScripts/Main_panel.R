#####################
# Variant Analysis #
####################

# Packages
library(OmicCircos)
library(foreach)
library(doMC)
library(rtracklayer)
library(openxlsx)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(RMySQL)
library(biomaRt)
library(RColorBrewer)
library(stringr)
library(Rsamtools)
library(gtrellis)
library(circlize)
library(ComplexHeatmap)
library(YAPSA)
library(SomaticSignatures)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringi)
library(tidyr)
library(dplyr)
library(magrittr)
library(pracma)
library(scarHRD)
library(gtools)

args <- commandArgs()

#############
# Parameter #

sample <- paste(args[6:7], collapse = "_")
protocol <- args[6]
id <- args[7]
germline <- args[8]
tumor <- args[9]
base_path <- args[10]
path_output <- paste(base_path, "Analyses/", sep = "/")
path_input <- paste(base_path, "WES/", sep = "/")
path_script <- args[11]
path_data <- args[12]
targets_txt <- args[13]
#covered_region <- args[14]
author <- args[15]
center <- args[16]
bed_file <- args[17]
sureselect_type <- args[18]
ref_genome <- args[19]
target_capture_cor_factors <- args[20]
vaf <- as.numeric(args[21])*100
min_var_count <- as.numeric(args[22])
maf_cutoff <- as.numeric(args[23])
actionable_genes <- ifelse(args[24]=="/opt/MIRACUM-Pipe/databases/", NA, args[24])
covered_exons <- args[25]
entity <- as.character(args[26])
gender <- as.character(args[27])
fusion_genes <- args[28]
ampl_genes_txt <- args[29]
ucsc_server <- args[30]
cnv_region_annotation <- args[31]

print(ref_genome)

#############
# Functions #
# Mutation Analyses #
print("Load functions.")
source(paste(path_script, "filtering.R", sep = "/"))
source(paste(path_script, "filtering_tools.R", sep = "/"))
source(paste(path_script, "mutationAnalysis.R", sep = "/"))
source(paste(path_script, "MutAna_tools.R", sep = "/"))
source(paste(path_script, "stats.R", sep = "/"))
source(paste(path_script, "stats_tools.R", sep = "/"))
# Copy Number Analyses #
source(paste(path_script, "cnv_analysis.R", sep = "/"))
source(paste(path_script, "cnv_ana_tools.R", sep = "/"))
# Mutation Signature Analysis #
source(paste(path_script, "mutationSignatureAnalysisFunction.R", sep = "/"))
source(paste(path_script, "Mut_sig_tools.R", sep = "/"))
source(paste(path_script, "mut_sig_ana.R", sep = "/"))

##################
# FASTQC Reports #

tumor <- paste0(
    path_input, strsplit(x = tumor, split = ".", fixed = T)[[1]][1],
    "_fastqc/Images/per_base_quality.png"
  )
  tumor_bsqr <- paste0(
    path_input, sample,
    "_td_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png"
  )

# DEFINE FILES
print("Preparations.")
# TUMOR
snp_file_td <- paste0(
	path_input, sample,
	"_vc.output.snp.Sample1.avinput.hg19_multianno.csv"
	)
indel_file_td <- paste0(
	path_input, sample,
	"_vc.output.indel.Sample1.avinput.hg19_multianno.csv"
	)
snpefffile_snp_td <- paste0(
	path_input, sample,
	"_vc.output.snp.SnpEff.vcf"
	)
snpefffile_indel_td <- paste0(
	path_input, sample,
	"_vc.output.indel.SnpEff.vcf"
	)
filter_out_td <- paste0(
	path_output, sample,
	"_VC_TUMOR.xlsx"
	)
mutect2_vcf <- paste0(
	path_input, sample,
	"_td_gatk4_mutect2_filtered.hg19_multianno.txt"
	)
mutect2_snpEff_vcf <- paste0(
	path_input, sample,
	"_td_gatk4_mutect2_filtered_SnpEff.vcf"
	)
msi_file <- paste0(
	path_input, sample,
	"_td_MSI"
	)
# Results
outfile_circos <- paste0(path_output, sample, "_TD_circos.pdf")
outfile_circos_mutect2 <- paste0(path_output, sample, "_TD_circos_mutect2.pdf")
coverage_out <- paste0(path_output, sample, "_coverage.pdf")
# MAFs
maf_td <- paste0(path_output, sample, "_TUMOR.maf")
maf_complete <- paste0(path_output, sample, ".maf")
maf_td_mutect2 <- paste0(path_output, sample, "_Mutect2_TUMOR.maf")
# Biomarker
biomarker_out <- paste0(path_output, sample, "_Biomarker.txt")

##############
# STATISTICS #
## Input Files
print("Statistics.")
stats_td <- paste0(path_input, sample, "_td_stats.txt")

stats <- stats_func(path = path_input, outfile_pdf = coverage_out,
	stats_td = stats_td, protocol = protocol, sureselect_type = sureselect_type)

#####################
# Mutation Analyses
#####################

print("Filtering for Tumor.")
filt_result_td <- filtering(
    snpfile = snp_file_td,
    indelfile = indel_file_td,
    snpefffile_snp = snpefffile_snp_td,
    snpefffile_indel = snpefffile_indel_td,
    outfile = filter_out_td,
    path_data = path_data,
    path_script = path_script,
    mode = "T",
    center = center,
    id = id,
    protocol = protocol,
    sureselect = bed_file,
    vaf = vaf,
    min_var_count = min_var_count,
    maf = maf_cutoff,
    covered_exons = covered_exons,
    cov_t = stats$cover_exons$perc[[1]][1],
    sureselect_type = sureselect_type
  )
# GATK4 Mutect2
filt_result_td_mutect2 <- filtering_mutect2(
    snpfile = mutect2_vcf,
    snpefffile = mutect2_snpEff_vcf,
    id = id,
    path_data = path_data,
    path_script = path_script,
    mode = "T",
    center = center,
    protocol = protocol,
    sureselect = bed_file,
    vaf = vaf,
    min_var_count = min_var_count,
    maf = maf_cutoff,
    covered_exons = covered_exons,
    cov_t = stats$cover_exons$perc[[1]][1],
    sureselect_type = sureselect_type
  )

filt_result_loh <- list(table = NULL, tmb = NULL)

# Analyses
print("Variant Analyses.")
mutation_analysis_result <- mutation_analysis(
    loh = filt_result_loh$table,
    somatic = filt_result_td$table,
    tumbu = filt_result_td$tmb,
    outfile_circos = outfile_circos,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    protocol = protocol,
    sureselect = bed_file,
    sureselect_type = sureselect_type,
    msi_file = msi_file
  )
# GATK4 Mutect2
mutation_analysis_result_mutect2 <- mutation_analysis(
    loh = filt_result_loh$table,
    somatic = filt_result_td_mutect2$table,
    tumbu = filt_result_td_mutect2$tmb,
    outfile_circos = outfile_circos_mutect2,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    protocol = protocol,
    sureselect = bed_file,
    sureselect_type = sureselect_type,
    msi_file = msi_file
  )

# Combine MAF files to obtain one complete maf per patient
maf_comb <- filt_result_td_mutect2$maf
write.table(
x = maf_comb,
file = maf_complete,
append = F,
quote = F,
sep = "\t",
col.names = T,
row.names = F
)


########################
# Copy Number Analysis #
print("CNV Analyses.") #
########################
## Input/Output Files
ratio_file <- paste0(
  path_input, "CNV/", sample,
  "_td_output.sort.rmdup.realigned.fixed.recal.bam_ratio.txt"
)
cnvs_file <- paste0(
  path_input, "CNV/", sample,
  "_td_output.sort.rmdup.realigned.fixed.recal.bam_CNVs"
)
cpn_file <- paste0(
  path_input, "CNV/", sample,
  "_td_output.sort.rmdup.realigned.fixed.recal.bam_sample.cpn"
)
# HRD/Purity
purity_file <- paste0(
  path_input,
  sample,
  "_sequenza/",
  sample,
  "_alternative_solutions.txt"
)
hrd_file <- paste0(
  path_input,
  sample,
  "_HRD.txt"
)

# Results
cnv_pvalue_txt <- paste0(cnvs_file, ".p.value.txt")
cnv_ideogram_plot <- paste0(path_output, sample, "_CNV_Plot_Ideogram.pdf")
outfile_cnvs_cbioportal <- paste0(path_output, sample, "_CNV_cbioportal.txt")
outfile_cnvs_seg <- paste0(path_output, sample, "_CNV.seg")

cnv_analysis_results <- cnv_analysis(
  ratio_file = ratio_file,
  cnvs_file = cnvs_file,
  cpn_file = cpn_file,
  cnv_pvalue_txt = cnv_pvalue_txt,
  outfile_ideogram = cnv_ideogram_plot,
  path_data = path_data,
  path_script = path_script,
  targets_txt = targets_txt,
  ampl_genes_txt = ampl_genes_txt,
  outfile_cbioportal = outfile_cnvs_cbioportal,
  outfile_seg = outfile_cnvs_seg,
  id = id,
  protocol = protocol,
  sureselect_type = sureselect_type,
  gender = gender,
  purity_file = purity_file,
  hrd_file = hrd_file,
  ucsc_server = ucsc_server,
  cnv_region_annotation = cnv_region_annotation
)

###############################
# Mutation Signature Analysis #
print("Mutation Signature Analysis.")

vcf <- paste0(path_input, sample, "_td_gatk4_mutect2_filtered.vcf")

outfile_mutsig_cbioportal <- paste0(path_output, sample, "_mutsig_cbioportal")
mut_sig_analysis <- mutation_signature_analysis(
vcf_file = vcf,
cutoff = 0.01,
sample_name = NULL,
only_coding = FALSE,
path_data = path_data,
path_output = path_output,
outfile_cbioportal = outfile_mutsig_cbioportal
)
mut_sig <- data.frame(Signature = mut_sig_analysis$CosmicValid_cutoffGen_LCDlist$out_sig_ind_df$sig,
                Process = mut_sig_analysis$CosmicValid_cutoffGen_LCDlist$out_sig_ind_df$process,
                Frequency = mut_sig_analysis$CosmicValid_cutoffGen_LCDlist$norm_exposures[,1])
rownames(mut_sig) <- mut_sig$Signature

# Write Excel File
print("Write Excel Table.")

print("3 - panelTumor")
output <- list(
Mutations = filt_result_td_mutect2$table,
CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated,
Mutation_Signatures = mut_sig
)
write.xlsx(
x = output,
file = paste0(path_output, sample, "_results.xlsx"),
overwrite = TRUE,
na.string = "."
)


# Export TMB, MSI, HRD and BRCAness as text file for cBioPortal
# das kann ja auch eigentlich erstmal raus
print("cBioPortal Export.")
if(protocol == "panelTumor" & sureselect_type == "TSO500") {
  if (mutation_analysis_result_mutect2$msi < 20) {
    msi_helper <- "Non-MSI-H"
  } else {
    msi_helper <- "Instable"
  }
  brca_helper <- which(mut_sig ==  "AC3")
  if (length(brca_helper) == 1 & mut_sig["AC3", 3] * 100 > 1) {
    brca_helper <- paste0(round(mut_sig["AC3", 3] * 100, digits = 1))
  } else {
    brca_helper <- "<1%"
  }
  biomarker <- data.frame(
    Tumor_Sample_Barcode = paste(as.character(id),"TD",sep = "_"),
    MSI_SCORE = mutation_analysis_result_mutect2$msi,
    MSI_TYPE = msi_helper,
    CVR_TMB_SCORE = filt_result_td_mutect2$tmb,
    BRCAness = brca_helper,
    HRD = cnv_analysis_results$hrd$sum,
    Purity = cnv_analysis_results$purity$purity,
    Ploidity = cnv_analysis_results$purity$ploidy
  )
  write.table(x = biomarker, file = biomarker_out , append = F, quote = F,
              sep = '\t', col.names = T, row.names = F)
}

#######################
####### FUSIONS #######
#######################
print("Fusion Analysis.")
fus_file <- paste(
  path_input,
  "../RNA/fusioncatcher/final-list_candidate-fusion-genes.hg19.txt",
  sep = "/"
)
outfile_fusions_cbioportal <- paste0(
  path_output, sample, "_fusions_cbioportal.txt"
)
if (file.exists(fus_file)) {
  source(paste(path_script, "fusions.R", sep = "/"))
  fusions <- fusions_ana(
    fus_file = fus_file,
    path_data = path_data,
    fusion_genes = fusion_genes
  )
  cat(dim(fusions$Table))
  # print Circosplot with Fusions
  if (!is.null(fusions$Table)) {
    if (dim(fusions$Table)[1] != 0) {
      if (length(which(duplicated(fusions$Table[, c(1,3)]))) != 0) {
        fus_tab <- fusions$Table[-which(
          duplicated(fusions$Table[, c(1,3)])
        ), c(1:4)]
      } else {
        fus_tab <- fusions$Table[, c(1:5) ]
      }
      id_1 <- grep(pattern = "not-converted", x = fus_tab[, 2])
      id_2 <- grep(pattern = "not-converted", x = fus_tab[, 4])
      if(length(union(id_1, id_2)) > 0){
        fus_tab <- fus_tab[-union(id_1, id_2), ]
      }
      # Prepare Fusionstable for CircosPlot
      sep <- strsplit(x = as.character(fus_tab$Bruch1), split = ":")
      fus_tab$Chr1 <- unlist(lapply(sep, function(x){return( x[[1]])}))
      fus_tab$Pos1 <- unlist(lapply(sep, function(x){return( x[[2]])}))
      sep <- strsplit(x = as.character(fus_tab$Bruch2), split = ":")
      fus_tab$Chr2 <- unlist(lapply(sep, function(x){return( x[[1]])}))
      fus_tab$Pos2 <- unlist(lapply(sep, function(x){return( x[[2]])}))
      id_genes <- union(
        which(
          nchar(fus_tab$Chr1) > 2
        ),
        which(
          nchar(fus_tab$Chr2) > 2
        ))
      if (length(id_genes) > 0){
        fus_tab <- fus_tab[-id_genes, c(
          "Gen1",
          "Chr1",
          "Pos1",
          "Gen2",
          "Chr2",
          "Pos2"
        )]
      } else {
        fus_tab <- fus_tab[, c("Gen1", "Chr1", "Pos1", "Gen2", "Chr2", "Pos2")]
      }
      if (sureselect_type == "TSO500") {
        sub <- div(
          filt_result_td_mutect2$table,
          NULL,
          TRUE,
          protocol = protocol,
          sureselect_type = "V5UTR"
        )
      } else {
        sub <- div(
          filt_result_td$table,
          NULL,
          TRUE,
          protocol = protocol,
          sureselect_type = "V5UTR"
        )
      }

      cc <- circos_colors(
        x_s_snp = sub$x_s_snp,
        x_s_indel = sub$x_s_indel,
        x_l_snp = sub$x_l_snp,
        x_l_indel = sub$x_l_indel, no_loh = sub$no_loh,
        no_indel_somatic = sub$no_indel_somatic,
        no_snp = sub$no_snp,
        no_indel_loh = sub$no_indel_loh,
        no_snp_loh = sub$no_snp_loh
      )
      omicCircosFus2(
        listOfMap = as.matrix(cc$map_mat),
        fusions = fus_tab,
        label = NULL, minR = 125,
        outfile = outfile_circos_mutect2,
        circosColors = as.vector(cc$circoscolors),
        mode = sureselect_type,
        protocol = "tumorOnly",
        path_data = path_data,
        trgt = bed_file
      )
      if (any(is.na(fusions$Plots$file))) {
        fus_38_files <- paste(
          fusions$Plots$paths[!is.na(fusions$Plots$file)],
          fusions$Plots$file[!is.na(fusions$Plots$file)], sep = "/"
        )
        for (i in 1:length(fus_38_files)) {
          system(paste0("cp ", fus_38_files[i], " ", path_output))
        }
      } else {
        fus_38_files <- paste(
          fusions$Plots$paths,
          fusions$Plots$file,
          sep = "/"
        )
        for (i in 1:length(fus_38_files)) {
          system(paste0("cp ", fus_38_files[i], " ", path_output))
        }
      }
      fusions2cbioportal(fusions, sample, outfile_fusions_cbioportal)
    } else {
      fusions <- NULL
    }
  } else {
    fusions <- NULL
  } 
} else {
  fusions <- NULL
}
save.image(file = "MTB.RData")

# Prepare Report
print("Report Preparation.")
source(paste(path_script, "Report_tools.R", sep = "/"))

tmb_med <- med_tmb(as.character(entity))
key_results <- keys(
	mut_sig = mut_sig,
	mutation_analysis_result = mutation_analysis_result_mutect2,
	mutation_analysis_result_gd = NULL,
	filt_result_td = filt_result_td_mutect2,
	cnv_analysis_results = cnv_analysis_results,
	filt_result_gd = NULL,
	med_tmb = tmb_med,
	protocol = protocol,
	fusions = fusions
	)

highlight_table <- highlight(
	muts_tab = mutation_analysis_result_mutect2$som_mut_tab, protocol = protocol
	)

sq_tab <- summary_quality(stats = stats, protocol = protocol)

som_muts <- mutation_analysis_result_mutect2$mut_tab
sum_mut_cg <- highlight_detail(
	muts_tab = mutation_analysis_result_mutect2$ts_og,
	Mode = "Tumor",
	protocol = protocol
	)
sum_loh_cg <- NULL

cnvs <- cnv_panel(cnv_results = cnv_analysis_results$out)

som_mut_pthw <- pthws_mut(
	df = mutation_analysis_result_mutect2$important_pathways,
	protocol = protocol
	)
som_mut_topart <- topart_mut(
	df = mutation_analysis_result_mutect2$important_pathways,
	protocol = protocol
	)
cnvs_pthws <- pathws_cnv(df = cnv_analysis_results$impa)

germ_mut_cg <- NULL
germ_mut_pthw <- NULL

som_all <- highlight_detail(
	muts_tab = mutation_analysis_result_mutect2$som_mut_tab,
	Mode = "Tumor",
	protocol = protocol
	)
germ_all <- NULL
loh_all <- NULL

save.image("Report.RData")
