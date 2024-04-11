library(tidyverse)
library(ggplot2)    
group1 <- read.table('../application/enrichment_result/Enrichment_snp_postive_set_lifted_methy_gregor_output_enrichment_results.txt', header=T)
group2 <- read.table('../application/enrichment_result/Enrichment_snp_postive_set_lifted_ha_gregor_output_enrichment_results.txt', header=T)
group1$group <- 'methy'
group2$group <- 'ha'
res_all <- rbind(group1, group2)
res_all$feature <- gsub(".bed","",res_all$Bed_File)
res_all <- res_all %>%filter(!(str_detect(feature,"hg38")))%>%filter(!(p_fisher == 1))%>% na.omit  # Because the database is hg19
p <- res_all%>%
    arrange(odds)%>%ggplot()+geom_point(aes(x = odds, y = reorder(feature,-odds), color = group))+
        geom_vline(aes( xintercept = 0 ))+theme_bw()+theme(text = element_text(size = 20))+xlab("Odds Ratio")+
            geom_linerange(aes(xmin = low, xmax= high , y = feature ), color= 'grey40')+ylab("Functional Annotations")
ggsave(plot = p, filename = '../plot/methy_vs_ha_enrichment.pdf', height = 24, width = 20)


