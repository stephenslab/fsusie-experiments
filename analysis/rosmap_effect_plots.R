library(fsusieR)
highconf_overlap <- read.csv("../outputs/rosmap_highconf_overlap.csv",
                             stringsAsFactors = FALSE)
highconf_overlap <- transform(highconf_overlap,chr = factor(chr))

# Candidate regions identified:
# 
#                 id   chr       pos     region_eqtl region_haqtl region_mqtl
# chr1:170663243:A:G  chr1 170663243 ENSG00000116132      TADB_75     TADB_75
# chr17:63440644:T:C chr17  63440644 ENSG00000008283    TADB_1205   TADB_1205
# chr19:46435226:C:T chr19  46435226 ENSG00000169515    TADB_1261   TADB_1261
#  chr3:44721338:A:G  chr3  44721338 ENSG00000163808     TADB_258    TADB_258
#
# *Start with chr3:44721338:A:G.*
#
