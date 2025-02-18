#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls[["General"]][["Seed"]])

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_taxa_filter.R")
source("R/fix_taxa_names_for_plot.R")
source("R/long_otu_table.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
OutDir <- prm.ls$DA$out_dir

DataSets <- prm.ls$DA$data_set_ps

TaxaLvl <- prm.ls$DA$tax_lvl

TaxaNorm <- prm.ls$DA$ps_norm

TaxaNormPlot <- prm.ls$DA$plot_norm

TaxaNormSumm <- prm.ls$DA$summary_norm

TaxaMinPrev <- prm.ls$DA$feature_min_prev

MaasPadjMethod <- prm.ls$DA$p_adj_method

MaasRandEff <- prm.ls$DA$maas_rand

MaasMainEff <- prm.ls$DA$maas_main_effect

MaasFixCovarEff <- prm.ls$DA$maas_covar

MaasPrmSets <- prm.ls$DA$maas_method_prm

MaxQvalue <- prm.ls$DA$max_qval


#-------------------------------------------------------------------------------
# Make directories
#-------------------------------------------------------------------------------
DirsOut <- expand.grid("Set" = DataSets, 
                       "Dir" = c("plots", "tabs"), 
                       stringsAsFactors = FALSE)

for(i in 1:nrow(DirsOut)) {
  
 dir.create(paste0(OutDir, "/", 
                   paste(DirsOut[i, ], collapse = "/")),
            recursive = TRUE, showWarnings = FALSE) 
  
}

#-------------------------------------------------------------------------------
# Data and effects grid 
#-------------------------------------------------------------------------------
PrmSets <- expand.grid("Set" = DataSets, 
                       "Lvl" = TaxaLvl,
                       "Norm" = TaxaNorm, 
                       "MaasPrmSet" = names(MaasPrmSets),
                       stringsAsFactors = FALSE)

# Fixed effects variables to use for maaslin 
MaasFixed <- c(MaasMainEff, prm.ls$DA$maas_covar)

# Run MaasLin2
MaasResLs <- list()

MaasResDf <- NULL

for(i in 1:nrow(PrmSets)) { 
  
  # Variables instances 
  iSet <- PrmSets[["Set"]][[i]]
  
  iLvl <- PrmSets[["Lvl"]][[i]]
  
  iNorm <- PrmSets[["Norm"]][[i]]
  
  iMaasPrm <- MaasPrmSets[[PrmSets[["MaasPrmSet"]][[i]]]]
  
  # MaasLin parameters name
  iMaasPrmName <- paste(iMaasPrm, collapse = "-")
  
  # Make dir for maaslin
  MaasDirInst <- paste0(OutDir, "/maas_res/", 
                        paste(PrmSets[i, ], collapse = "-"))
  
  dir.create(MaasDirInst, recursive = TRUE, showWarnings = FALSE)
  
  # Extract data 
  MetaInst <- data.ls[[iSet]]$meta
  
  OtuInst <- data.ls[[iSet]]$PS[[iLvl]][[iNorm]] %>% 
                  otu_table() %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame()
  
  # Run MaasLin2 
  InstMaasRes <- Maaslin2(
                  input_data =  OtuInst,
                  input_metadata = MetaInst,
                  output = MaasDirInst,
                  fixed_effects = MaasFixed, 
                  random_effects = MaasRandEff,
                  correction = MaasPadjMethod,
                  cores = prm.ls$General$ncores,
                  min_abundance = 0,
                  min_prevalence = TaxaMinPrev,
                  min_variance = 0,
                  normalization = iMaasPrm[["norm"]],
                  transform = iMaasPrm[["trans"]],
                  analysis_method = iMaasPrm[["method"]],
                  max_significance = 1,
                  plot_scatter = FALSE,
                  plot_heatmap = FALSE,
                  save_scatter = FALSE)
  
  MaasResLs[[iSet]][[iNorm]][[iMaasPrmName]][[iLvl]] <- InstMaasRes
  
  MaasResDf <- InstMaasRes$results %>% 
                          bind_cols(PrmSets[i, ], 
                                    data.frame(MaasParm = iMaasPrmName)) %>% 
                          bind_rows(MaasResDf, .)
}


#-------------------------------------------------------------------------------
# Plot results 
#-------------------------------------------------------------------------------
# Variables 
FacetNcol <- 4

# Empty objects 
CombPlotsLs <- list()

PlotsLs <- list()

SigMainTaxLs <- list()

for(i in 1:nrow(PrmSets)) {
 
  # Variables instances 
  iSet <- PrmSets[["Set"]][[i]]
  
  iLvl <- PrmSets[["Lvl"]][[i]]
  
  iNorm <- PrmSets[["Norm"]][[i]]
  
  iMaasPrm <- MaasPrmSets[[PrmSets[["MaasPrmSet"]][[i]]]]
  
  # MaasLin parameters 
  iMaasPrmName <- paste(iMaasPrm, collapse = "-")
  
  
  # Add information about abundance and write out results 
  #-----------------------------------------------------------------------------
  OtuSummInst <- long_otu_table(data.ls[[iSet]]$PS[[iLvl]][[TaxaNormSumm]], 
                                names_to = "SeqID", 
                                values_to = "Abundance") %>% 
                      left_join(., MetaInst, by = "SeqID") %>% 
                      reframe(Mean = mean(Abundance),
                              Median = median(Abundance),
                              SD = sd(Abundance),
                              IQR = IQR(Abundance),
                              .by = all_of(c("feature", MaasMainEff))) %>% 
                      pivot_wider(names_from = MaasMainEff, 
                                  values_from = all_of(c("Mean", "Median", 
                                                         "SD", "IQR")))
  
  # MaasLin results 
  MaasResSummInst <- MaasResLs[[iSet]][[iNorm]][[iMaasPrmName]][[iLvl]]$results %>% 
                      left_join(OtuSummInst, by = "feature") %>% 
                      mutate(Summary_norm = TaxaNormSumm, 
                             Maas_Prm = iMaasPrmName, 
                             Taxa_Lvl = iLvl,
                             Samples_Set = iSet, 
                             feature = fix_taxa_names_for_plot(feature))
 
  write.csv(MaasResSummInst, 
            paste0(OutDir, "/", iSet, "/tabs/MaasAll_", 
                   iLvl, "-", iMaasPrmName, ".csv"))
  
  # Plot data 
  #-----------------------------------------------------------------------------
  # Extract and transform OTU table 
  OtuPlotInst <- long_otu_table(data.ls[[iSet]]$PS[[iLvl]][[TaxaNormPlot]], 
                                names_to = "SeqID", 
                                values_to = "Abundance") %>% 
                      left_join(., MetaInst, by = "SeqID")
  
  # Significant taxa 
  SigTaxAllInst <- MaasResLs[[iSet]][[iNorm]][[iMaasPrmName]][[iLvl]]$results %>% 
                          filter(qval <= MaxQvalue)
  
  # Significant taxa corresponding to the main effect 
  SigTaxMainInst <- SigTaxAllInst %>% 
                        filter(metadata == MaasMainEff)
  
  # Plot taxa if any are significant 
  if(nrow(SigTaxAllInst) > 0) {
    
      # Plot main effect taxa as a boxplot (if any)
    if(nrow(SigTaxMainInst) > 0) {
      
      CoefDf <- SigTaxMainInst %>% 
                    mutate(Coefficient = ifelse(coef > 0, "Positive", "Negative")) %>% 
                    mutate(Coefficient = factor(Coefficient, 
                                                levels = c("Positive", "Negative")))
      
      # Write out significant taxa as separate csv
      MaasResSummInst %>% 
        filter(metadata == MaasMainEff, 
               feature %in% SigTaxMainInst$feature) %>% 
        write.csv(paste0(OutDir, "/", iSet, "/tabs/MaasSigPlot_", 
                         iLvl, "-", iMaasPrmName, ".csv"))
      
      # Extract and format abundance data 
      OtuPlotSigInst <- OtuPlotInst %>% 
                          filter(feature %in% SigTaxMainInst$feature) %>% 
                          left_join(., CoefDf, by = "feature") %>% 
                          mutate(feature = fix_taxa_names_for_plot(feature))
      
      # Plot results
      SigPlot <- ggplot(OtuPlotSigInst, 
                        aes(x = Abundance, 
                            y = .data[[MaasMainEff]],
                            color = .data[[MaasMainEff]], 
                            linetype = Coefficient)) +
                      geom_jitter(width = 0, 
                                  height = 0.25, 
                                  alpha = 0.1) + 
                      geom_violin(fill = NA) +
                      theme_bw() +
                      scale_color_manual(values = aes.ls$col[[iSet]]) +
                      facet_wrap("feature", scales = "free", 
                                 ncol = FacetNcol) +
                      theme(axis.title.y = element_blank(), 
                            axis.ticks.y = element_blank(),
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(),
                            strip.text = element_text(face = "italic", 
                                                      size = 7), 
                            panel.grid = element_blank()) + 
                      scale_x_log10()
      
      PlotsLs[[iSet]][[iNorm]][[iMaasPrmName]][[iLvl]][["v1"]] <- 
        SigPlot                 
      
      PlotsLs[[iSet]][[iNorm]][[iMaasPrmName]][[iLvl]][["ntax"]] <- 
        nrow(SigTaxMainInst)
      
    }
    
    # Heatmap overview of the significant taxa 
    SigTaxHeatDf <- MaasResLs[[iSet]][[iNorm]][[iMaasPrmName]][[iLvl]]$results %>% 
                        filter(feature %in% SigTaxAllInst$feature) %>% 
                        select(feature, metadata, qval) %>% 
                        mutate(feature = fix_taxa_names_for_plot(feature)) %>% 
                        pivot_wider(names_from = metadata, values_from = qval) %>% 
                        column_to_rownames(var = "feature") %>% 
                        select(all_of(MaasFixed))
    
  
    SigTaxHeat <- Heatmap(SigTaxHeatDf, 
                          name = "Q-value",
                          row_names_side = "left", 
                          row_names_gp = gpar(fontface = "italic", 
                                              fontfamily = "serif"), 
                          column_names_rot = 45, 
                          height = unit(nrow(SigTaxHeatDf)*0.25, "in"), 
                          width = unit(ncol(SigTaxHeatDf)*0.4, "in"),
                          cluster_rows = TRUE, 
                          show_row_dend = FALSE,
                          cluster_columns = FALSE, 
                          rect_gp = gpar(col = "gray25", lwd = 0.5),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if(SigTaxHeatDf[i, j] <= MaxQvalue) {
                              grid.text(sprintf("%.3f", SigTaxHeatDf[i, j]), 
                                        x, y, gp = gpar(fontsize = 8, 
                                                        col = "white", 
                                                        fontface = "bold"))
                            } else {
                              grid.text(sprintf("%.3f", SigTaxHeatDf[i, j]), 
                                        x, y, gp = gpar(fontsize = 6, 
                                                        col = "gray40"))}})
    
    
    png(filename = paste0(OutDir, "/", iSet, "/plots/", 
                           paste(PrmSets[i, ], collapse = "_"), "_",
                          iMaasPrmName, ".png"), 
        width = ncol(SigTaxHeatDf)*0.4 + 5, 
        height = nrow(SigTaxHeatDf)*0.25 + 2, 
        units = "in", 
        res = 600)
    
    draw(SigTaxHeat)
    
    dev.off()
    
    CombPlotsLs[[iSet]][[iNorm]][[iMaasPrmName]][["Heat"]][[iLvl]] <- SigTaxHeat
    
  }
  
}


#-------------------------------------------------------------------------------
# Combine plots from different taxonomic levels 
#-------------------------------------------------------------------------------
PlotsGap <- 0.25

PlotsLegendAdd <- 0.5

PrmSetNoLvl <- PrmSets %>% 
                  select(-Lvl) %>% 
                  distinct()

for(i in 1:nrow(PrmSetNoLvl)) {
  
  iSet <- PrmSetNoLvl[["Set"]][[i]]
  
  iNorm <- PrmSets[["Norm"]][[i]]
  
  iMaasPrm <- MaasPrmSets[[PrmSetNoLvl[["MaasPrmSet"]][[i]]]]
  
  iMaasPrmName <- paste(iMaasPrm, collapse = "-")
  
  PlotsInst <- PlotsLs[[iSet]][[iNorm]][[iMaasPrmName]]
  
  # Objects for combining plots 
  PlotsProp <- c()
  
  PlotsLabs <- c("A")
  
  PlotsList <- list()
  
  CountVar <- 1
  
  for(j in names(PlotsInst)) {
    
    CountVar <- CountVar + 1
    
    PlotsList[paste0("null_", j)] <- list(NULL)
    
    if(j != last(names(PlotsInst))) {
      
      jPlot <- PlotsInst[[j]]$v1 + 
                    theme(legend.position = "none")
      
      jProp <- ceiling(PlotsInst[[j]]$ntax/FacetNcol)
      
      PlotsLabs <- c(PlotsLabs, "", LETTERS[CountVar])
      
    } else {
      
      jPlot <- PlotsInst[[j]]$v1 + 
                    theme(legend.position = "bottom")
      
      jProp <- (ceiling(PlotsInst[[j]]$ntax/FacetNcol) + PlotsLegendAdd)
      
      PlotsLabs <- c(PlotsLabs, "")
      
    }
    
    PlotsList[[j]] <- jPlot
    
    PlotsProp <- c(PlotsProp, PlotsGap, jProp)
    
  }
  
  iCombPlots <- plot_grid(plotlist = PlotsList, 
                          ncol = 1, 
                          rel_heights = PlotsProp, 
                          labels = PlotsLabs)
  
  ggsave(filename = paste0(OutDir, "/", iSet, "/plots/Main-", 
                           paste(PrmSetNoLvl[i, ], collapse = "_"), "_",
                           iMaasPrmName, ".png"), 
         plot = iCombPlots, 
         width = FacetNcol*2, 
         height = sum(PlotsProp)*1 + 1, 
         dpi = 600)
  
  CombPlotsLs[[iSet]][[iNorm]][[iMaasPrmName]][["Main"]] <- iCombPlots

}


#-------------------------------------------------------------------------------
# Write and clean up 
#-------------------------------------------------------------------------------
save(file = "out/supp/DA_res.Rdata", list = c("CombPlotsLs", 
                                              "MaasResDf", 
                                              "SigMainTaxLs", 
                                              "MaasResLs"))
# Clean environment
rm(list = ls())
gc()
