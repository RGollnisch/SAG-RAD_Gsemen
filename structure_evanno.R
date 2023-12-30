## Analysis to determine the best value for K in Structure runs using the Evanno method

# remotes::install_github('royfrancis/pophelper')

library(pophelper)
library(gridExtra)
library(readr)
library(dplyr)
library(ggplot2)


# EUR populations Structure analysis path
structure_path_eur = "./analysis/9.structure/EUR_50k_p13_r70_coraf/allK_EUR_50k_p13_r70_coraf/"
structure_files_eur <- list.files(path = structure_path_eur)

# USA populations Structure analysis path
structure_path_usa = "./analysis/9.structure/USA_noSMLLCC_50k_p7_r70_coraf/allK_USA_noSMLLCC_50k_p7_r70_coraf/"
structure_files_usa <- list.files(path = structure_path_usa)


## EUR evanno analysis
setwd(structure_path_eur)

# Convert STRUCTURE run files to qlist
structure_list_eur <- readQ(files = structure_files_eur, filetype = "structure", indlabfromfile = TRUE, readci = FALSE) 

structure_table_eur <- tabulateQ(structure_list_eur)
structure_summary_eur <- summariseQ(structure_table_eur)

evanno_eur <- evannoMethodStructure(structure_summary_eur, returnplot = TRUE, basesize = 10, linesize = 0.7)
evanno_eur$data
evanno_eur.plot <- grid.arrange(evanno_eur$plot)
ggsave(filename = "./analysis/manuscript_figures/EUR_structure_evanno.png", plot = evanno_eur.plot, width = 210, height = 162, unit = "mm", dpi = 300)


## USA evanno analysis
setwd(structure_path_usa)

# Convert STRUCTURE run files to qlist
structure_list_usa <- readQ(files = structure_files_usa, filetype = "structure", indlabfromfile = TRUE, readci = FALSE) 

structure_table_usa <- tabulateQ(structure_list_usa)
structure_summary_usa <- summariseQ(structure_table_usa)

evanno_usa <- evannoMethodStructure(structure_summary_usa, returnplot = TRUE, basesize = 10, linesize = 0.7)
evanno_usa$data
evanno_usa.plot <- grid.arrange(evanno_usa$plot)
ggsave(filename = "./analysis/USA_structure_evanno.png", plot = evanno_usa.plot, width = 210, height = 162, unit = "mm", dpi = 300)