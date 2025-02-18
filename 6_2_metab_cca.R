################################################################################
# CCA 
################################################################################
#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls$General$Seed)


#-------------------------------------------------------------------------------
# Custom functions
#-------------------------------------------------------------------------------
source("R/create_subdir.R")
source("R/cca_extract_for_plot.R")

# Function to format text 
sig_text <- function(col, n_dig, pref) {
  
  ColRound <- round(col, n_dig)
  
  TextIfLess <- paste0(pref, "<", 
                       paste(replicate((n_dig-1), "0"), 
                             collapse = ""), "1")
  
  TextOut <- ifelse(ColRound == 0, 
                    TextIfLess,
                    paste0(pref, "==", sprintf(paste0("%.", n_dig,"f"), 
                                               ColRound)))
  
  return(TextOut)
  
}


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
MetabSets <- prm.ls$metab$cca_tabs

EffMain <- prm.ls$metab$cca_main

EffCovar <- prm.ls$metab$cca_covar

TransFun <- prm.ls$metab$cca_trans_count


#-------------------------------------------------------------------------------
# Run analysis
#-------------------------------------------------------------------------------
# Samples overlap 
SamplesIDPs <- data.ls$phen$PS[[1]][[1]] %>% 
                  sample_names() %>% 
                  gsub("CIW1", "", .)

# Make sub directories for results 
create_subdir(prm.ls$metab$dir_out, 
              list("Set" = MetabSets, 
                   "Type" = c("plots", "tabs")))

# Formula for CCA
CcaForm <- paste(c(EffMain, EffCovar), collapse = "+") %>% 
                paste0("iMetabol~", .)

CcaResDf <- NULL

CcaResLs <- list()

for(i in MetabSets) {
  
  # Extract data metabolites 
  iMetabol <- data.ls$Metabol[[i]] %>% 
                    setNames(trimws(names(.))) %>% 
                    filter(ID %in% SamplesIDPs) %>% 
                    column_to_rownames(var = "ID")  %>% 
                    mutate(across(where(is.numeric), 
                                  function(x){replace(x, is.na(x), 
                                                      median(x, na.rm = TRUE))}))
  
  if(!is.null(TransFun)) {
    
    iMetabol <- iMetabol %>% 
                  mutate(across(everything(), TransFun))
    
  }
                
  # Extract metadata 
  iMeta <- data.ls$phen$meta
  
  rownames(iMeta) <- gsub("CIW1", "", rownames(iMeta))
  
  iMeta <- iMeta[rownames(iMetabol), ] %>% 
                  mutate(across(where(is.numeric), 
                                function(x){replace(x, is.na(x), 
                                                    median(x, na.rm = TRUE))}))
                
              
  # CCA analysis 
  iCcaRes <- cca(as.formula(CcaForm), iMeta)
  
  iCcaAnova <- anova.cca(iCcaRes, 
                         permutations = prm.ls$metab$cca_nperm, 
                         by = "margin")
  
  iCcaAnovaDf <- iCcaAnova %>% 
                    broom.mixed::tidy() %>% 
                    mutate(Set = i, 
                           Formula = CcaForm) %>% 
                    suppressWarnings()
  
  write.csv(iCcaAnovaDf, 
            file = paste0(prm.ls$metab$dir_out, 
                          "/", i, "/tabs/", i, "_stat_cca.csv"), 
            row.names = FALSE)
  
  # Collect data 
  CcaResDf <- bind_rows(CcaResDf, iCcaAnovaDf)
  
  CcaResLs[[i]][["cca"]] <- iCcaRes
  
  CcaResLs[[i]][["anova"]] <- iCcaAnova
  
  
  #-----------------------------------------------------------------------------
  # Plot Ordination
  #-----------------------------------------------------------------------------
  iPlotData <- cca_extract_for_plot(cca_obj = iCcaRes, 
                                    metadata = iMeta, 
                                    const_var = c(prm.ls$metab$cca_main, 
                                                  prm.ls$metab$cca_covar))
  
  iGroupStatText <- iCcaAnovaDf %>% 
                        filter(term == prm.ls$metab$cca_main) %>% 
                        mutate(Text = paste0(sig_text(p.value, 3, "P"), 
                                             "~~(", 
                                             sig_text(ChiSquare, 3, "X^2"), 
                                             ")"), 
                               x = min(iPlotData$main[[1]]), 
                               y = max(iPlotData$main[[2]]))  
  
  iPlotBase <- ggplot(iPlotData$main,
                      aes(x = .data[[colnames(iPlotData$main)[1]]], 
                          y = .data[[colnames(iPlotData$main)[2]]], 
                          colour = Phenotype)) +
                    geom_point(alpha = 0.5) +
                    coord_fixed() +
                    theme_bw() + 
                    theme(axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + 
                    xlab(iPlotData$var_expl[1]) +
                    ylab(iPlotData$var_expl[2]) +
                    theme(plot.title = element_text(size=11, 
                                                    face="italic")) + 
                    scale_color_manual(values = aes.ls$col$phen) + 
                    geom_text(data = iGroupStatText, 
                              aes(x = x, y = y, label = Text), 
                              inherit.aes = FALSE,
                              hjust = 0, parse = T)
  
  # Significant terms 
  iCcaAnovaDfSig <- iCcaAnovaDf %>% 
                        filter(p.value <= prm.ls$metab$cca_max_pval)
  
  # Extract significant vectors and centroids if any
  if(nrow(iCcaAnovaDfSig) > 0) {
    
    iVectSub <- iPlotData$vect %>% 
                    filter(label %in% iCcaAnovaDfSig$term) %>% 
                    rename(Vector = label)
    
    iCenterSub <- iPlotData$centr %>% 
                        filter(grepl(paste(iCcaAnovaDfSig$term, collapse = "|"), 
                                     label)) %>% 
                        rename(Centroid = label)
    
    # Add vectors to plot 
    if(nrow(iVectSub) > 0) {
      
      iPlotBase <- iPlotBase + 
                      geom_segment(data = iVectSub, 
                                   aes(x = x, y = y, 
                                       xend = 0, yend = 0), 
                                   inherit.aes = FALSE, 
                                   arrow = arrow(length=unit(0.10,"cm"), 
                                                 ends="first", 
                                                 type = "closed")) +
                      geom_text(data = iVectSub, 
                                aes(x = x, y = y, label = Vector), 
                                inherit.aes = FALSE, hjust = -0.1)
      
    }
    
    # Add centroids to plot 
    if(nrow(iCenterSub) > 0) { 
      
      iPlotBase <- iPlotBase + 
                      geom_point(data = iCenterSub, 
                                 aes(x = x, y = y, 
                                     shape = Centroid), 
                                 inherit.aes = FALSE, 
                                 size = 2)
    }
    
  }
  
  CcaResLs[[i]][["plot"]] <- iPlotBase
  
  ggsave(filename = paste0(prm.ls$metab$dir_out, 
                           "/", i, "/plots/", 
                           i, "_cca.png"), 
         plot = iPlotBase, 
         width = prm.ls$metab$cca_plot_dim$w, 
         height = prm.ls$metab$cca_plot_dim$h)
  
}


#-------------------------------------------------------------------------------
save(list = c("CcaResLs", "CcaResDf"), 
     file = "out/supp/6_2_met_cca.Rdata")

# Clean environment
rm(list = ls())
gc()
