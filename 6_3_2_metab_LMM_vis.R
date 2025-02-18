#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")
load("out/supp/6_3_1_metab_LMM.Rdata")

set.seed(prm.ls$General$Seed)

source("R/create_subdir.R")


################################################################################
# Deferential abundance with LMM
################################################################################
#-------------------------------------------------------------------------------
# Variables and data to use 
#-------------------------------------------------------------------------------
# Parameter variables 
MainEff <- prm.ls$metab$da_main_effect

FixEff <- c(prm.ls$metab$da_main_effect, prm.ls$metab$da_covar) 

RandEff <- prm.ls$metab$da_rand

MetabTabs <- prm.ls$metab$da_tabs

TransCountForm <- prm.ls$metab$da_trans_count 

ParticipIDcol <- prm.ls$General$Part_id_col

MaxQval <- prm.ls$metab$da_max_qval

DirOut <- prm.ls$metab$dir_out

NcolPlots <- prm.ls$metab$da_ncol_plots


#-------------------------------------------------------------------------------
# Visualize significant differences if any 
#-------------------------------------------------------------------------------
# Keep only main therm and add adjusted pval 
StatResForm <- StatResDf %>% 
                filter(grepl(MainEff, term)) %>% 
                group_by(Set, Method) %>% 
                arrange(p.value, .by_group = TRUE) %>% 
                mutate(qval = p.adjust(p.value, 
                                       method = prm.ls$metab$da_p_adj_method)) %>% 
                ungroup() %>% 
                mutate(Significant = ifelse(qval <= prm.ls$metab$da_max_qval, 
                                            "Yes", "No"))

# Write filtered results per set 
for(i in unique(StatResForm$Set)) {
  
  create_subdir(DirOut, 
                list("Set" = i, 
                     "Type" = c("plots", "tabs"), 
                     "SubType" = "LMM"))
  
  StatResForm %>% 
    filter(Set == i) %>% 
    write.csv(paste0(DirOut, "/", i, "/tabs/LMM/", 
                     MainEff, "_metabol_lmm.csv"))
    
  
}

# Keep only significant 
StatResFormSig <- StatResForm  %>% 
                    filter(qval <= MaxQval, 
                           Method == "LMM")

MetPlots <- list()

# Visualize per 
if(nrow(StatResFormSig) > 1) {
  
  for(i in unique(StatResFormSig[["Set"]])) { 
    
    iSigMet <- StatResFormSig %>% 
                    filter(Set == i)
    
    MetabDfSigP <- MetabSetsForm[[i]] %>% 
                        pivot_longer(cols = starts_with("met_"), 
                                     names_to = "Metabolite") %>% 
                        filter(Metabolite %in% iSigMet$Metabol) %>% 
                        mutate(Metabolite = gsub("met_", "", Metabolite)) %>% 
                        mutate(Metabolite = gsub("_", " ", Metabolite))
    
    iViolin <- ggplot(MetabDfSigP, 
                      aes(y = value, 
                          x = .data[[MainEff]],
                          color = .data[[MainEff]])) +
                  geom_jitter(width = 0.25, 
                              height = 0, 
                              alpha = 0.1) + 
                  geom_violin(fill = NA) +
                  theme_bw() +
                  scale_color_manual(values = aes.ls$col[["phen"]]) +
                  facet_wrap("Metabolite", scales = "free", 
                             ncol = 4) +
                  theme(strip.text = element_text(face = "italic", size = 7), 
                  panel.grid = element_blank()) 
    
    # Assumptions 
    FitResFilt <- FitResDf %>% 
                    filter(Set == i, 
                           Index %in% iSigMet$Metabol, 
                           Method == "LMM")
    
    iResidPlot <- ggplot(FitResFilt, 
                         aes(y = Residuals, x = Fitted)) + 
                      geom_point(color = "steelblue", alpha = 0.5) + 
                      geom_abline(intercept = 0, slope = 0) + 
                      facet_wrap(c("Method", "Index"), 
                                 scales = "free", 
                                 ncol = 4) + 
                      theme_bw() 
    
    iQqPlot <- ggplot(FitResFilt, aes(sample = Residuals)) + 
                  stat_qq(color = "steelblue", alpha = 0.5) + 
                  stat_qq_line() + 
                  facet_wrap(c("Method", "Index"), 
                             scales = "free", 
                             ncol = 4) + 
                  theme_bw()
    
    MetPlots[[i]][["violin"]] <- iViolin
    
    MetPlots[[i]][["resid"]] <- iResidPlot 
    
    MetPlots[[i]][["qq"]] <- iQqPlot
    
    
    #---------------------------------------------------------------------------
    # Save plots 
    #---------------------------------------------------------------------------
    
    PlotWidth <- NcolPlots
    
    if(NcolPlots > nrow(iSigMet)) {
      
      PlotWidth <- nrow(iSigMet)
      
    }
    
    for(j in c("iViolin", "iResidPlot", "iQqPlot")) {
      
      ggsave(filename = paste0(DirOut, "/", i, "/plots/LMM/", 
                                j, MaxQval, ".png"), 
             plot = get(j), 
             width = (PlotWidth*2 + 1), 
             height = ceiling(PlotWidth/nrow(iSigMet))*2.5)
    }
    
  }
  
  StatResForm %>% 
    filter(Set == i) %>% 
    write.csv(paste0(DirOut, "/", i, "/tabs/LMM/", 
                     i, "_stat_res.csv"))
  
}


#-------------------------------------------------------------------------------
# General overview - Volcano plot 
#-------------------------------------------------------------------------------
for(i in MetabTabs) {
  
  create_subdir(DirOut, 
                list("Set" = i, "Type" = c("plots", "tabs"), "SubType" = "LMM"))
  
  iVilcanoDf <- StatResForm %>% 
                    filter(Method == "LMM", Set == i) %>% 
                    rename("pval" = "p.value") %>% 
                    pivot_longer(cols = c("pval", "qval")) %>% 
                    mutate(!!"p/q-value" := ifelse(value < 0.001, 
                                                   "<0.001", 
                                                   ifelse(value <= 0.01, 
                                                          "0.001-0.01", 
                                                          ifelse(value > 0.05, 
                                                                 ">0.05", 
                                                                 ">0.01-0.05")))) %>% 
                    mutate(!!"p/q-value" := factor(.data[["p/q-value"]], 
                                                   levels = c("<0.001", 
                                                              "0.001-0.01", 
                                                              ">0.01-0.05", 
                                                              ">0.05")))
  
  iVulcanoPlot <- ggplot(iVilcanoDf, 
                         aes(x = estimate, 
                             y = -log(value), 
                             colour = .data[["p/q-value"]])) + 
                      geom_point(size = 2, alpha = 0.5) + 
                      facet_wrap("name", scales = "free") + 
                      theme_bw() +
                      scale_color_manual(values = c("<0.001" = "#E41A1C", 
                                                    "0.001-0.01" = "#377EB8", 
                                                    ">0.01-0.05" = "#4DAF4A", 
                                                    ">0.05" = "grey40")) + 
                      geom_vline(xintercept=0) + 
                      theme(panel.grid = element_blank())
  
  
  ggsave(paste0(DirOut, "/", i, "/plots/LMM/Volcano_overview.png"), 
         plot = iVulcanoPlot, width = 8, height = 4)
  
}

# Clean environment
rm(list = ls())
gc()