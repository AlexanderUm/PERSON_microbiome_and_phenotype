#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls$General$Seed)

#-------------------------------------------------------------------------------
# Summary of the reads 
#-------------------------------------------------------------------------------
source("R/long_otu_table.R")
source("R/create_subdir.R")

create_subdir(main_dir = prm.ls$General$root_dir, 
              sub_dir_ls = list("miscelleneous", c("reads_summ")))

OtuTabLong <- long_otu_table(data.ls$with_contr$PS$ASV$Raw) %>% 
                summarise(Sum_Abundace = sum(Abundance), .by = SeqID) %>% 
                left_join(data.ls$with_contr$meta, by = "SeqID") %>% 
                mutate(Sum_Abundace = Sum_Abundace/1000)

ReadsDistrPlot <- ggplot(OtuTabLong, aes(x = Sum_Abundace)) + 
                      geom_histogram(fill = "steelblue", color = "white") + 
                      facet_wrap(~Phenotype) + 
                      theme_bw() + 
                      xlab(label = "Reads (x1000)") + 
                      ylab(label = "Samples Count") + 
                      scale_x_log10()

ggsave(filename = paste0(prm.ls$General$root_dir, 
                       "/miscelleneous/reads_summ/distr.png"), 
       plot = ReadsDistrPlot, width = 5.5, height = 2.5, units = "in")


# Summary of reads 
OtuTabLong %>% 
  summarise(Min = min(Sum_Abundace*1000), 
            Max = max(Sum_Abundace*1000), 
            Median = median(Sum_Abundace*1000), 
            Mean = round(mean(Sum_Abundace*1000), 0), .by = "Phenotype") %>% 
  write.csv(paste0(prm.ls$General$root_dir, 
                   "/miscelleneous/reads_summ/summary.csv"), 
            row.names = FALSE)
