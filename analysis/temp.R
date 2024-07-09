library(data.table)
ad_studies <- c("AD_Bellenguez_2022","AD_Bellenguez_EADB_2022",
                "AD_Bellenguez_EADI_2022","AD_Bellenguez_GRACE_2022",
                "AD_Jansen_2021","AD_Kunkle_Stage1_2019", 
                "AD_Wightman_Excluding23andMe_2021",
                "AD_Wightman_ExcludingUKBand23andME_2021", 
                "AD_Wightman_Full_2021")
gwas <- fread(file.path("/project2/mstephens/fungen_xqtl",
                        "ftp_fgc_xqtl/analysis_result",
                        "Fungen_xQTL_allQTL.overlapped.gwas.export.csv.gz"),
              sep = ",",verbose = TRUE,stringsAsFactors = FALSE)
class(gwas) <- "data.frame"
gwas <- subset(gwas,
               cs_coverage_0.95 > 0 &
               is.element(study,ad_studies))
gwas <- gwas[c("variant_id","pip","study","method","region","betahat",
               "sebetahat","maf")]

ids <- subset(fsusie_mqtl$cs,pip > 0.25)$id
ids <- unique(intersect(ids,gwas$variant_id))
gwas_mqtl <- subset(gwas,is.element(variant_id,ids))

ids <- subset(fsusie_haqtl$cs,pip > 0.25)$id
ids <- unique(intersect(ids,gwas$variant_id))
gwas_haqtl <- subset(gwas,is.element(variant_id,ids))

temp <- subset(fsusie_mqtl$cs,
               get_chr_from_id(id) == "chr20" &
               abs(get_pos_from_id(id) - 56412112) < 2e4)
head(sort(abs(get_pos_from_id(temp$id) -
              subset(genes,gene_name == "CASS4")$start)))

min(abs(get_pos_from_id(subset(gwas,chr == "chr20")$variant_id) - 56409008))
