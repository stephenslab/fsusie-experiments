source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/scripts_plot/cases_study/ggplot_HSP90AA1_case_study.R")
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/scripts_plot/cases_study/ggplot_AHCYL1_case_study.R")
  source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/scripts_plot/cases_study/ABHD17A_case_study.R")
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/scripts_plot/cases_study/SCD_case_study.R" )



library(gridExtra)
fout= grid.arrange(sCD,ABH,AHCL,HSP, ncol=2)

ggsave(fout, height = 16, width = 20, filename="combined_plot.pdf")


library(cowplot)

# Combine plots in 2x2 layout with labels A–D
fout <- plot_grid(
  sCD, ABH, AHCL, HSP,
  labels = c("A", "B", "C", "D"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 2
)

# Save to PDF
ggsave("combined_plot_GTeX.pdf", fout, height = 16, width = 20)
