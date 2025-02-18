################################################################################
# Capture RF permute output as a dataframe 
#-------------------------------------------------------------------------------
RF_catch_out <- function(rfPermut_results) {
  
  require(tidyverse)
  require(randomForest)
  
  # Capture output
  cons.out <- capture.output(rfPermut_results) 
  
  output.end <- length(rfPermut_results$rf$classes) + 12
  
  catch.out <- c(paste0("Group", cons.out[11]), cons.out[12:output.end]) %>% 
                  str_squish() %>% 
                  str_split(pattern = " ", simplify = TRUE) %>% 
                  as.data.frame() %>% 
                  select(-c("V2", "V3")) %>% 
                  setNames(.[1, ]) %>% 
                  .[-1, ] %>% 
                  mutate(across(2:4, as.numeric)) %>% 
                  setNames(c("Group", "Correct(%)", "CI(low)", "CI(high)")) 
  
  return(catch.out)
}


################################################################################
# Crud tuning of mtry
#-------------------------------------------------------------------------------
RF_get_mtry <- function(data_for_RF, 
                        group_col, 
                        ntreeTry, 
                        max_mtry_prop = 0.35) {
  
    require(tidyverse)
    require(randomForest)
    
    TuneDf <- tuneRF(x = select(data_for_RF,  -all_of(group_col)),
                      y = data_for_RF[[group_col]],
                      ntreeTry = ntreeTry,
                      improve = 1e-5,
                      trace = FALSE,
                      plot = FALSE) %>% 
                  as.data.frame()
  
  MtryOut <- TuneDf$mtry[which.min(TuneDf$OOBError)]
  
  if(any(c(MtryOut == 1, MtryOut >= ceiling((ncol(data_for_RF)-1)*max_mtry_prop)))) {
    
    MtryOut <- ceiling(sqrt((ncol(data_for_RF)-1)))
    
  }
  
  return(MtryOut)
  
}


################################################################################
# Extract features importance 
#-------------------------------------------------------------------------------
RF_extract_importance <- function(rfPermut_results, 
                                  col_arrange_by = "MeanDecreaseAccuracy_importance") {
  
  require(tidyverse)
  
  imp.df <- rfPermut_results$rf$importance %>% 
                as.data.frame() %>% 
                setNames(paste0(names(.), "_importance")) 
  
  imp.sd.df <- rfPermut_results$rf$importanceSD %>% 
                as.data.frame() %>% 
                setNames(paste0(names(.), "_SD")) 
              
  p.val.df <- rfPermut_results$pval %>% 
                as.data.frame() %>% 
                setNames(paste0(names(.), "_pval")) 
  
  comb.df <- bind_cols(imp.df, imp.sd.df, p.val.df) %>% 
                arrange(desc(.data[[col_arrange_by]])) %>% 
                rownames_to_column(var = "feature")
  
  return(comb.df) 
  
}


################################################################################
# Run RFpermute with parameters and create AUC plot 
#-------------------------------------------------------------------------------
RF_run_with_par <- function(data_in, 
                            RF_formula,
                            group_col, 
                            n_trees, 
                            n_perm, 
                            gr_size_prop = NA,
                            Mtry_auto_tune = FALSE,
                            n_cores = 4) {
  
  require(pROC)
  require(rfPermute)
  
  #-----------------------------------------------------------------------------
  # Initial Tuning 
  #-----------------------------------------------------------------------------
  # Run mtry parameter tuning. Custom function
  
  if(Mtry_auto_tune) {
    
    MtryInst <- RF_get_mtry(data_in, 
                            group_col = group_col, 
                            ntreeTry = n_trees)
    
  } else {
    
    MtryInst <- ceiling(sqrt((ncol(data_in)-1)))
    
  }

  
  # Run RF
  if(all(is.numeric(gr_size_prop), !is.na(gr_size_prop))) {
    
    # Balance sample size
    SampSizeRf <- balancedSampsize(data_in[[group_col]],
                                   pct = gr_size_prop)
    
    #-----------------------------------------------------------------------------
    # Run RF per for original and group column permuted data sets
    #-----------------------------------------------------------------------------
    RfRes <- rfPermute(RF_formula,
                       data = data_in,
                       ntree = n_trees,
                       mtry = MtryInst,
                       num.cores = n_cores,
                       num.rep = n_perm, 
                       strata = data_in[[group_col]],
                       sampsize = SampSizeRf,
                       proximity = TRUE)
    
  } else {
    
    message("No group balancing performed.")
    
    RfRes <- rfPermute(RF_formula,
                       data = data_in,
                       ntree = n_trees,
                       mtry = MtryInst,
                       num.cores = n_cores,
                       num.rep = n_perm, 
                       proximity = TRUE)
    
  }
  
  
  #---------------------------------------------------------------------------
  # Plot ROC and calculate AUC
  #---------------------------------------------------------------------------
  RfRoc <- roc(data_in[[group_col]], RfRes$rf$votes[,2])
  
  Auc <- round(auc(RfRoc), 2)
  
  RocPlot <- ggroc(RfRoc, 
                   colour = 'steelblue', 
                   size = 1) +
    ggtitle(paste0('ROC Curve ', '(AUC = ', Auc, ')')) + 
    geom_abline(intercept=1, 
                slope=1, 
                linetype="dashed") + 
    theme_classic()
  
  
  #-----------------------------------------------------------------------------
  # Capture output. Custom function
  #-----------------------------------------------------------------------------
  RfResDf <- RF_catch_out(RfRes) %>% 
                mutate(mtry = MtryInst, 
                       n_features = (ncol(data_in)-1), 
                       AUC = Auc) 
  
  
  #-----------------------------------------------------------------------------
  # Results out
  #-----------------------------------------------------------------------------
  RfResOut <- list("Full" = RfRes, 
                   "Summary_df" = RfResDf, 
                   "ROC_objects" = list("plot" = RocPlot, 
                                        "AUC" = Auc, 
                                        "roc" = RfRoc))
  
  return(RfResOut)
  
}


################################################################################
# Plot accuracy of RF
#-------------------------------------------------------------------------------
RF_plot_accuracy <- function(RF_res_main, 
                             x_col_name,
                             color_col, 
                             legend_name = NULL, 
                             accuracy_col = "Correct(%)", 
                             low_CI = "CI(low)", 
                             high_CI = "CI(high)") {
  
  SummaryPlot <- ggplot(RF_res_main, 
                        aes(color = .data[[color_col]])) +
                geom_hline(yintercept = 50,
                           size = 0.5, 
                           color = "black", 
                           linetype = "dashed") + 
                geom_point(aes(x = .data[[x_col_name]], 
                               y = .data[[accuracy_col]]), 
                           # stat="identity",
                           position=position_dodge(.5), 
                           size = 2.5) +
                geom_errorbar(aes(x = .data[[x_col_name]], 
                                  ymin = .data[[low_CI]],  
                                  ymax = .data[[high_CI]]),
                              width=.2,
                              position=position_dodge(.5), 
                              size = 0.75) + 
                ylim(c(0, 100)) + 
                theme_bw() + 
                xlab("") + 
                theme(axis.text.x = element_text(angle = 45, 
                                                 vjust = 1, 
                                                 hjust=1)) 
  
  return(SummaryPlot)
}


################################################################################
# Plot Features importance 
#-------------------------------------------------------------------------------
RF_plot_importance <- function(importance_tab, 
                               cols_to_plot = c("LIR_importance",
                                                "MIR_importance",
                                                "MeanDecreaseAccuracy_importance",
                                                "MeanDecreaseGini_importance"), 
                               col_arrange_by = "MeanDecreaseAccuracy_importance", 
                               pvals_to_use = "unscaled",
                               feature_col = "feature", 
                               remove_from_names = "_importance", 
                               extend_y_axis = 1.75, 
                               y_text_face_italic = TRUE) {
  
  require(tidyverse)
  
  pval.cols <- gsub(remove_from_names, "", cols_to_plot) %>% 
                paste0(., ".*", pvals_to_use) %>% 
                paste(., collapse = "|") %>% 
                grep(., colnames(importance_tab), value = TRUE)
  
  DfToPlot0 <- importance_tab %>% 
                  select(all_of(c(feature_col, cols_to_plot, pval.cols))) %>% 
                  setNames(gsub("\\.", "_", names(.))) %>% 
                  setNames(gsub(paste0("_", pvals_to_use), "", names(.))) %>% 
                  arrange(.data[[col_arrange_by]])
  
  DfToPlot <- DfToPlot0 %>% 
                  pivot_longer(-all_of(feature_col), 
                               names_to = c("Index", ".value"), 
                               names_sep = "_") %>% 
                  arrange(importance) %>% 
                  mutate(PvalShort = sprintf("%.3f", round(pval, 3)), 
                         Index = factor(Index, 
                                        levels = gsub(remove_from_names, "", 
                                                      cols_to_plot)), 
                         feature = factor(feature, levels = DfToPlot0$feature)) %>% 
                  mutate(PvalText = ifelse(PvalShort == 0, 
                                           "P<0.001", 
                                           paste0("P=", PvalShort)), 
                         !!"P-value" := ifelse(pval <= 0.05, pval, NA), 
                         PvalLabPos = ifelse(importance < 0, 0, importance))
  
  PlotOut <- ggplot(DfToPlot, aes(y = .data[[feature_col]], 
                                  x = importance, 
                                  fill = .data[["P-value"]])) + 
                geom_bar(stat = "identity") + 
                geom_point(aes(x = importance*extend_y_axis), shape = NA) +
                geom_text(aes(x = PvalLabPos, label = PvalText), 
                          hjust = -0.1, size = 3) +
                facet_grid(~Index, scales = "free") + 
                theme_bw() + 
                theme(axis.title = element_blank())
  
  
  if(y_text_face_italic) {
    
    PlotOut <-  PlotOut + 
                   theme(axis.text.y = element_text(face = "italic"))
    
  }
  
  return(PlotOut)
  
}



