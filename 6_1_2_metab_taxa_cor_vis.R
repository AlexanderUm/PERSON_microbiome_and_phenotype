#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")
load("out/supp/6_1_1_met_corr.Rdata")

set.seed(prm.ls$General$Seed)


#-------------------------------------------------------------------------------
# Custom functions
#-------------------------------------------------------------------------------
source("R/create_subdir.R")
source("R/fix_taxa_names_for_plot.R")

# Calculate prevalence per taxa from a phyloseq 
phy_tax_pervalence <- function(phy_obj, per_group = NULL) {
  
  TransPhy <- transform_sample_counts(phy_obj, 
                             function(x){ifelse(x==0, 0, 1)}) %>%
                             suppressMessages()
  
  TaxSums <- TransPhy %>% 
                    taxa_sums() 
                  
  TaxPrev <- TaxSums/nsamples(phy_obj) 
  
  TaxPrevDf <- data.frame(Overall = TaxPrev)
    
  
  if(!is.null(per_group)) {
    
    StrataVec <- sample_data(TransPhy)[[per_group]]
    
    for(i in unique(StrataVec)) {
      
      TransPhyFilt <- prune_samples(StrataVec == i, TransPhy)  
      
      TaxPrevStr <- taxa_sums(TransPhyFilt)/nsamples(TransPhyFilt) 
      
      TaxPrevDf <- data.frame(TaxPrevStr) %>% 
                        setNames(i) %>% 
                        bind_cols(TaxPrevDf, .)
      
    }
    
  }
  
  return(TaxPrevDf)
  
}


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Plot as a heat map 
SigVal <- prm.ls$metab$cor_heat_val

SigValCut <- prm.ls$metab$cor_heat_val_cut

QvalMethod <- prm.ls$metab$cor_p_adj_method

CorStrata <- prm.ls$metab$cor_strata

CorFiltPrm <- prm.ls$metab$cor_filt_plot_prm

MinTaxPrev <- prm.ls$metab$cor_tax_prev

MetabTabs <- prm.ls$metab$cor_metab_tab

CorMethod <- prm.ls$metab$cor_method

NetWorkMinEst <- prm.ls$metab$cor_network_min_est

NetWorkMaxSig <- prm.ls$metab$cor_network_max_sig

NetHeatValCols <- prm.ls$metab$cor_netheat_val

FullHeatMinEst <- prm.ls$metab$cor_full_heat_min_est

NetJitter <- prm.ls$metab$cor_network_jitter

MinEstPvalAdj <- prm.ls$metab$cor_padj_filt_est


#-------------------------------------------------------------------------------
# Directories 
#-------------------------------------------------------------------------------
create_subdir(prm.ls$metab$dir_out, 
              list("Set" = prm.ls$metab$cca_tabs, 
                   "Type" = c("plots", "tabs")))

#-------------------------------------------------------------------------------
# Plot significant correlations 
#-------------------------------------------------------------------------------
# Correlation between metabolites and individual taxa 
CorGrid <- prm.ls$metab[c("cor_tax_lvl", 
                          "cor_tax_norm", 
                          "cor_metab_tab")] %>% 
                  expand.grid(stringsAsFactors = FALSE)

SigPlotsLs <- list()

for(i in 1:nrow(CorGrid)) {
  
  # Variables 
  iTaxLvl <- CorGrid[i, "cor_tax_lvl"]
  
  iTaxNorm <- CorGrid[i, "cor_tax_norm"]
  
  iMetabTabId <- CorGrid[i, "cor_metab_tab"]
  
  # Retrieve filtered tabs
  iOtuTabF <- FiltTabs[["otu"]][[paste0(iTaxLvl, "_", iTaxNorm)]] 
  
  iMetabTabF <- FiltTabs[["met"]][[iMetabTabId]]
  
  # Calculate prevalence of taxa
  iTaxPrev <- data.ls$phen$PS[[iTaxLvl]][[iTaxNorm]] %>% 
                    prune_samples(paste0(rownames(iOtuTabF), "CIW1"), .) %>% 
                    phy_tax_pervalence(per_group = CorStrata)
  
  # Taxa only with prevalence higher than cutoff. 
  iTaxPrevFilt <- iTaxPrev %>% 
                    select(-Overall) %>% 
                    filter(if_any(where(is.numeric), ~ . >= MinTaxPrev))
  
  # Retrieve correlations results table and add q-values and additional info
  iCorRes <- CorResLs[[paste0(iTaxLvl, "_", iTaxNorm)]][[iMetabTabId]] %>% 
                filter(x %in% rownames(iTaxPrevFilt)) %>% 
                mutate(Cor_Id = paste0(x, "--", y), 
                       y_Strata = paste0(Strata, "--", y), 
                       Cor_Str = ifelse(abs(estimate) >= MinEstPvalAdj,
                                       "Strong", "Weak")) %>% 
                arrange(Strata, p.value) %>% 
                mutate(qval_overall = p.adjust(p.value, method = QvalMethod), 
                       .by = Strata) %>% 
                arrange(Strata, x, p.value) %>% 
                mutate(qval_per_tax = p.adjust(p.value, method = QvalMethod),
                       .by = c(Strata, x)) %>% 
                arrange(Strata, Cor_Str, p.value) %>% 
                mutate(qval_min_est = p.adjust(p.value, method = QvalMethod),
                         .by = c(Strata, Cor_Str))
  
  # Write out correlation table 
  iCorRes %>% 
    filter(abs(estimate) >= MinEstPvalAdj) %>%
    mutate(x = fix_taxa_names_for_plot(x), 
           y = gsub("_", " ", y)) %>% 
    write.csv(paste0(prm.ls$metab$dir_out, "/", iMetabTabId, "/tabs/", 
                     iMetabTabId, "--", CorMethod, "_more", 
                     gsub("\\.", "", MinEstPvalAdj), "_cor.csv"))
  
  
  #-----------------------------------------------------------------------------
  # Full correlation heat map 
  #-----------------------------------------------------------------------------
  iCompleteHeatDf <- iCorRes %>% 
                        select(estimate, x, y_Strata) %>% 
                        pivot_wider(names_from = y_Strata, 
                                    values_from = estimate) %>% 
                        column_to_rownames(var = "x") %>% 
                        mutate(across(everything(), 
                                      function(x){replace(x, is.na(x), 0)}))
  
  iPhenotypeVec <- gsub("--.*", "", names(iCompleteHeatDf)) %>% 
                      gsub(".*_", "", .)
  
  iTopAnot <- HeatmapAnnotation(Phenothype = iPhenotypeVec, 
                                col = list(Phenothype = aes.ls$col$phen))
  
  iFullHeat <- Heatmap(iCompleteHeatDf, 
                        name = "Estimate",
                        show_column_names = FALSE, 
                        show_row_names = FALSE, 
                        top_annotation = iTopAnot)
  
  png(filename = paste0(prm.ls$metab$dir_out, "/", iMetabTabId, "/plots/", 
                        iMetabTabId, "--", CorMethod, "_all_cor.png"), 
      width = 16, 
      height = 8, 
      units = "in", 
      res = 600)
  
  draw(iFullHeat)
  
  dev.off()
    
  
  #-----------------------------------------------------------------------------
  # Selected correlations heatmap
  #-----------------------------------------------------------------------------
  iCorResFilt <- iCorRes
  
  for(j in CorFiltPrm) {
    
    iCorResFilt <- iCorResFilt %>% 
                       filter(eval(parse(text=j)))
    
  }
  
  
  if(nrow(iCorResFilt) > 0) {
    
    iCorResFilt2 <- iCorRes %>% 
                      filter(x %in% iCorResFilt$x) %>% 
                      filter(y %in% iCorResFilt$y)
    
    # Heatmap data frames 
    iHeatDf <- list()
    
    for(j in c("estimate", SigVal)) {
      
      iHeatDf[[j]] <- iCorResFilt2 %>% 
                        select(all_of(c(j, "x", "y_Strata"))) %>% 
                        arrange(y_Strata) %>% 
                        mutate(x = fix_taxa_names_for_plot(x)) %>% 
                        pivot_wider(values_from = all_of(j), 
                                    names_from = y_Strata) %>% 
                        column_to_rownames(var = "x")
    }
    
    iHeatSplit <- colnames(iHeatDf$estimate) %>% 
                        sub(paste0(CorStrata, "_"), "", .) %>% 
                        sub("--.*", "", .)
    
    iHeatColsId <- colnames(iHeatDf$estimate) %>%
                        sub(".*--", "", .) %>% 
                        gsub("_", " ", .)
    
    iSigTaxHeat <- Heatmap(iHeatDf$estimate, 
                           name = "Estimate",
                           row_names_side = "left", 
                           row_names_gp = gpar(fontface = "italic", 
                                               fontfamily = "serif"), 
                           column_names_rot = 45, 
                           height = unit(nrow(iHeatDf$estimate)*0.3, "in"), 
                           width = unit(ncol(iHeatDf$estimate)*0.4, "in"),
                           column_split = iHeatSplit, 
                           column_labels = iHeatColsId,
                           cluster_rows = FALSE, 
                           show_row_dend = FALSE,
                           cluster_columns = FALSE, 
                           rect_gp = gpar(col = "gray25", lwd = 0.5),
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             if(iHeatDf[[SigVal]][i, j] <= SigValCut) {
                               grid.text(sprintf("%.3f", iHeatDf[[SigVal]][i, j]), 
                                         x, y, gp = gpar(fontsize = 8, 
                                                         col = "white", 
                                                         fontface = "bold"))
                             } else {
                               grid.text(sprintf("%.3f", iHeatDf[[SigVal]][i, j]), 
                                         x, y, gp = gpar(fontsize = 6, 
                                                         col = "gray40"))}})
    
    SigPlotsLs[[paste0(iTaxLvl, "_", iTaxNorm)]][[iMetabTabId]] <- iSigTaxHeat
    
    
    png(filename = paste0(prm.ls$metab$dir_out, "/", iMetabTabId, "/plots/", 
                          iMetabTabId, "--", CorMethod, "_cor.png"), 
        width = ncol(iHeatDf$estimate)*0.4 + 5, 
        height = nrow(iHeatDf$estimate)*0.25 + 2, 
        units = "in", 
        res = 600)
    
    draw(iSigTaxHeat)
    
    dev.off()
    
  }
  
  
  #-----------------------------------------------------------------------------
  # Network of strong correlations 
  #-----------------------------------------------------------------------------
  NetFiltGrid <- expand.grid("Est" = NetWorkMinEst, 
                             "MaxSig" = NetWorkMaxSig, 
                             stringsAsFactors = FALSE)
  
  for(FiltPrm in 1:nrow(NetFiltGrid)) { 
    
    MinEst <- NetFiltGrid[FiltPrm, "Est"]
    
    SigFilt <- NetFiltGrid[FiltPrm, "MaxSig"]
    
    # Prepare data 
    iCorResNet <- iCorRes %>% 
                    filter(abs(estimate) >= MinEst) %>%
                    filter(eval(parse(text = SigFilt))) %>% 
                    rename("T" = x, "M" = y) %>% 
                    mutate(Correlation = ifelse(estimate > 0, 
                                                "Positive", 
                                                "Negative"), 
                           EstAbs = abs(estimate)) 
    
    if(nrow(iCorResNet) > 0) {
      
      iVertCount <- iCorResNet %>% 
                        select(all_of(c("T", "M", "Strata"))) %>% 
                        pivot_longer(cols = c("T", "M"), 
                                     values_to = "vertex.names", 
                                     names_to = "Type") %>% 
                        mutate(Vert_Count = n(), 
                               .by = c(Strata, vertex.names)) %>% 
                        distinct()  
      
      iVertNames <- iVertCount %>% 
                        arrange(Type, desc(Vert_Count)) %>% 
                        select(Type, vertex.names) %>% 
                        distinct() %>% 
                        mutate(Vert_ID = paste0(Type, 1:n()), 
                               .by = Type) %>% 
                        select(-Type)
      
      iVertNameCount <- left_join(iVertCount, 
                                  iVertNames, 
                                  by = "vertex.names")
      
      iNetPlotLs <- list()
      
      for(j in unique(iCorResNet$Strata)) {
        
        #-----------------------------------------------------------------------
        # Network itself 
        #-----------------------------------------------------------------------
        jCorResNetFilt <- iCorResNet %>% 
                              filter(Strata == j)
        
        jVertNameCountFilt <- iVertNameCount %>% 
                                  filter(Strata == j)  
        
        # Convert into Network object (package "network")              
        jNetData <- jCorResNetFilt[, c("T", "M")] %>% 
                          network(directed = FALSE, multiple = TRUE)
        
        # Add parameters to edges
        set.edge.attribute(jNetData, 
                           c("Estimate", "Correlation"), 
                           jCorResNetFilt[, c("EstAbs", "Correlation")])
        
        # Convert into data frame for plotting 
        jNetDataDf <- ggnetwork(jNetData, 
                                layout = "fruchtermanreingold", 
                                niter = 5000, 
                                cell.jitter = NetJitter, 
                                cell.pointpointrad = 3) %>% 
                          left_join(jVertNameCountFilt, by = "vertex.names") %>% 
                          mutate(vertex.names = Vert_ID)
        
        # Plot 
        jNetPlot <- ggplot(jNetDataDf, 
                                  aes(x = x, y = y, 
                                      xend = xend, yend = yend)) +
                      geom_edges(aes(linewidth = Estimate, 
                                     color = Correlation)) +
                      geom_nodes(aes(size = Vert_Count, fill = Type), 
                                 shape = 21) +
                      geom_nodetext(aes(label = vertex.names),
                                    fontface = "bold", size = 2.75) +
                      theme_blank() + 
                      scale_color_manual(values = c("Negative" = "#E41A1C", 
                                                    "Positive" = "#4DAF4A")) + 
                      scale_linewidth(range = c(0.5, 3.5)) + 
                      scale_size(range = c(10, 20)) +
                      scale_fill_manual(values = c("M" = "#377EB8", 
                                                   "T" = "#984EA3")) + 
                      theme(legend.position = "none")
        
        # Collect/record output 
        iNetPlotLs[[j]] <- jNetPlot
        
        jDirPath <- paste0(prm.ls$metab$dir_out, "/", 
                           iMetabTabId, "/plots/Network_est", 
                           MinEst, "_", SigFilt) %>% 
                      gsub("<|\\=", "", .)
        
        dir.create(jDirPath, recursive = TRUE, showWarnings = FALSE)
        
        # Plot dimensions 
        jPlotSize <- nrow(jVertNameCountFilt)*0.1
        
        if(jPlotSize < 10) {jPlotSize <- 10}
        
        
        ggsave(filename = paste0(jDirPath, "/NetWork_", j, "--",
                                 CorMethod, ".png"), 
               plot = jNetPlot, 
               width = jPlotSize, 
               height = jPlotSize)
        
        jVertNameCountFilt %>% 
          arrange(Type, desc(Vert_Count)) %>% 
          mutate(vertex.names = fix_taxa_names_for_plot(vertex.names)) %>% 
          mutate(vertex.names = trimws(gsub("_", " ", vertex.names))) %>% 
          write.csv(file = paste0(jDirPath, "/NetWork_", j, "--",
                                  iMetabTabId, "--", CorMethod, ".csv"))
        
        write.csv(jCorResNetFilt, 
                  file = paste0(jDirPath, "/CorrNetWork_", j, "--",
                                iMetabTabId, "--", CorMethod, ".csv"))
        
      }
      
    }
  
  }

}


#-------------------------------------------------------------------------------
# Clean environment
rm(list = ls())
gc()


