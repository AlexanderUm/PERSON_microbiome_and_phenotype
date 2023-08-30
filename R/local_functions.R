###############################################################################
# General functions 
################################################################################
  #-----------------------------------------------------------------------------
  # Function: print a list of kables (tables)
  #-----------------------------------------------------------------------------
  print_kables <- function(tabs_named_ls) {
    
    for (i.tab in names(tabs_named_ls)) {
      
      tab.print <- tabs_named_ls[[i.tab]]
      
      print(knitr::kable(tab.print, caption = i.tab))
      
      cat('\n\n<!-- -->\n\n')
      
    }
    
  }
  
  
  #-----------------------------------------------------------------------------
  # Function: dynamically change size of figures
  #-----------------------------------------------------------------------------
  subchunkify <- function(g, fig_height=7, fig_width=5) {
    g_deparsed <- paste0(deparse(
      function() {g}
    ), collapse = '')
    
    sub_chunk <- paste0("`","``{r sub_chunk_", 
                        floor(runif(1) * 10000), 
                        ", fig.height=", fig_height, 
                        ", fig.width=", fig_width, ", echo=FALSE}",
                        "\n(", g_deparsed, ")()","\n`","``")
    
    cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
    
  }
  
  
  #-----------------------------------------------------------------------------
  # Function: recursive lapply for data frames in nested lists
  #------------------------------------------------------------------------------
  recurse_lapply <- function(df_nested_list, fun) {
    
    if(inherits(df_nested_list, "data.frame")){
      
      fun(df_nested_list)
      
    } else if(!is.list(df_nested_list) && length(df_nested_list) == 1) {
      
      fun(df_nested_list)
      
    }  else {lapply(df_nested_list, recurse_lapply, fun) }
    
  }
  
  
  
  #-----------------------------------------------------------------------------
  # Function: Filter otu table (taxa are columns) based on prevalence 
  #-----------------------------------------------------------------------------
  prev_filt_otu_tab <- function(otu_tab, meta_data, group_col = FALSE, prev_cutoff) {
    
    if(is.logical(group_col)) {
      
      group.col <- rep("a", nrow(otu_tab))
      
    } else {
      
      group.col <- pull(meta_data, group_col) }
    
    filt.vec <- otu_tab %>%
      mutate(across(where(is.numeric), function(x){ifelse(x == 0, 0, 1)})) %>%
      mutate(group_col = group.col)  %>%
      group_by(group_col) %>%
      summarise(across(where(is.numeric), 
                       function(x){sum(x)/length(x)}), .groups = "drop") %>%
      dplyr::select(-group_col) %>% 
      apply(., 2, function(x){any(x >= prev_cutoff)}, simplify = TRUE)
    
    tab.out <- otu_tab[, filt.vec]
    
    return(tab.out)
    
  }

  
################################################################################
# Analysis specific 
################################################################################
  #-------------------------------------------------------------------------------
  # Wilcox test and summary from the data 
  #-------------------------------------------------------------------------------
  wilcox_plus <- function(numer_df, meta_data, group_col) {
    
    if(!identical(rownames(numer_df), rownames(meta_data))) {
      stop("Row names in data frames are not matching")
    }
    
    comb.df <- bind_cols(numer_df, select(meta_data, group_col))
    
    all.res <- NULL
    
    for(i in colnames(numer_df)) {
      
      wil.form <- paste0("`", i, "`", " ~ ", "`", group_col, "`")
      
      wil.res <- wilcox.test(as.formula(wil.form), data = comb.df)
      
      all.res <- comb.df %>% 
        select(all_of(c(i, group_col))) %>% 
        group_by(get(group_col)) %>% 
        summarise(Mean = mean(get(i)), 
                  SD = sd(get(i))) %>% 
        pivot_wider(names_from = "get(group_col)", 
                    values_from = c("Mean", "SD")) %>% 
        mutate(p_value = wil.res$p.value, 
               stat_W = wil.res$statistic, 
               ID = i, 
               Formula = wil.form) %>% 
        bind_rows(all.res, .)
      
    }
    
    return(all.res)
  }
  
  
  #-------------------------------------------------------------------------------
  # Constrain ordination with vegan 
  #-------------------------------------------------------------------------------
  cca_test_plot <- function(data, 
                            metadata, 
                            formula, 
                            color_col,
                            bi_var_to_plot = color_col,
                            veg_method = "cca",
                            permut = how(nperm = 999), 
                            by = "terms") {
    
    if(!identical(rownames(data), rownames(metadata))) {
      stop("Row names are not identical.")
    }
    
    method.form <- paste0("data", "~", formula)
    
    ord.obj <- get(veg_method)(as.formula(method.form), data = metadata)
    
    anova.obj <- anova.cca(ord.obj, permutations = permut, by = by)
    
    summary.obj <- summary(ord.obj)
    
    # Extract data for plot-------------------------------------------------------
    var.to.plot <- paste(bi_var_to_plot, collapse = "|")
    
    sites.df <- summary.obj$sites[, 1:2] %>% 
      bind_cols(., select(metadata, all_of(color_col)))
    
    prop.expl <- round((summary.obj$cont$importance[2, ]*100), 1) 
    
    
    bi.p <- ggplot() + 
      geom_point(data = sites.df, aes_string(x = colnames(sites.df)[1], 
                                             y = colnames(sites.df)[2], 
                                             color = color_col)) + 
      theme_bw() + 
      xlab(label = paste0(colnames(sites.df)[1], " [", prop.expl[1], "%]")) + 
      ylab(label = paste0(colnames(sites.df)[2], " [", prop.expl[2], "%]"))
    
    
    # Centroids and arrows 
    if(!is.na(summary.obj$centroids)) {
      
      cent.df <- as.data.frame(summary.obj$centroids)[, 1:2] %>% 
        mutate(Centroids = rownames(.)) %>% 
        filter(grepl(var.to.plot, .$Centroids)) %>% 
        mutate(Centroids = gsub(var.to.plot, "", Centroids))
      
      bi.p <- bi.p + 
        geom_point(data = cent.df, 
                   aes_string(x = colnames(cent.df)[1], 
                              y = colnames(cent.df)[2], 
                              shape = "Centroids"),
                   size = 2.5) # + 
        # geom_text(data = cent.df, 
        #           aes_string(x = colnames(cent.df)[1], 
        #                      y = colnames(cent.df)[2], 
        #                      label = "Centroids"), 
        #           alpha = 0.5, 
        #           hjust = -0.25)
      
    }
    
    bip.d <- scores(ord.obj, choices = 1:2, display = "bp")
    
    bip.df <- as.data.frame(bip.d)
    
    if (any(!rownames(bip.df) %in% rownames(summary.obj$centroids))) {
      
      arrow.scale <- ordiArrowMul(bip.d, fill = 0.6)*10
      
      arrows.df <- bip.df[!rownames(bip.df) %in% 
                            rownames(summary.obj$centroids), ]*arrow.scale 
      
      arrows.df <- arrows.df %>% 
        mutate(Arrow_name = rownames(.)) %>% 
        filter(grepl(var.to.plot, .$Arrow_name)) 
      
      
      bi.p <- bi.p +
        geom_segment(data = arrows.df,
                     aes_string(xend = colnames(arrows.df)[1],
                                yend = colnames(arrows.df)[2], 
                                x = 0, y = 0), 
                     arrow = arrow(length = unit(0.2, "cm"))) +
        geom_text(data = arrows.df,
                  aes_string(x = colnames(arrows.df)[1],
                             y = colnames(arrows.df)[2],
                             label = "Arrow_name"),
                  alpha = 0.5,
                  hjust = -0.2)
      
    }
    
    
    # Triplot----------------------------------------------------------------------
    species.df <- summary.obj$species[, 1:2] %>% 
      as.data.frame(.) %>% 
      mutate(Spec = rownames(.))
    
    tri.p <- bi.p + 
      geom_text(data = species.df, 
                aes_string(x = colnames(species.df)[1], 
                           y = colnames(species.df)[2], 
                           label = "Spec"), 
                fontface = "italic", 
                alpha = 0.35, 
                color = "blue")
    
    obj.out <- list(anoca_cca = anova.obj, 
                    bi_plot = bi.p, 
                    tri_plot= tri.p)
    
    return(obj.out)
    
  }
  
  
  #-------------------------------------------------------------------------------
  # Correlation between a list of data frames (otu tables) and a data frame 
  #-------------------------------------------------------------------------------
  bulk_corr_matrix <- function(otu_tabs_list, 
                               meta_data, 
                               second_df, 
                               strata_cols_meta = NULL,
                               prev_cut_offs = c(0.25), 
                               corr_method = "spearman", 
                               p_adj_method = "BH") {
    
    
    # Create data frame for stratification ---------------------------------------
    if(is.null(strata_cols_meta)) {
      
      strata.data <- data.frame(Strata = rep("No Stratification", 
                                             nrow(meta_data)))
      
    } else {  
      
      strata.data <- meta_data %>% 
        dplyr::select(all_of(strata_cols_meta)) %>% 
        unite(Strata, all_of(strata_cols_meta))
    }
    #-----------------------------------------------------------------------------
    
    
    # Create a list of lists that will be used later 
    res.out <- list()
    
    #--- Loop through: otu tables ---
    for (i.otu in names(otu_tabs_list)) {
      
      # --- Loop through: strata ---
      for (i.strata in unique(strata.data$Strata)) {
        
        # Subset data ------------------------------------------------------------
        otu.tab.i.strat <- otu_tabs_list[[i.otu]] %>% 
          .[strata.data$Strata == i.strata, ]
        
        second.df.i.strat <- second_df[strata.data$Strata == i.strata, ]
        
        meta.data.strata <- meta_data[strata.data$Strata == i.strata, ]
        
        
        # Row names check --------------------------------------------------------
        if(!identical(rownames(otu.tab.i.strat), 
                      rownames(meta.data.strata))) {
          stop("Row names are not identical")}
        
        if(!identical(rownames(otu.tab.i.strat), 
                      rownames(second.df.i.strat))) {
          stop("Row names are not identical")}
        #-------------------------------------------------------------------------
        
        # Filter based on prevalence 
        otu.tab.i.strat.f <- prev_filt_otu_tab(otu_tab = otu.tab.i.strat, 
                                               meta_data = meta.data.strata, 
                                               prev_cutoff = min(prev_cut_offs))
        
        # Stop if no taxa left after filtering----------------------------------
        if(ncol(otu.tab.i.strat.f) < 1) {
          paste0("No taxa left after prevalece filtering. \n", 
                 "Cut Off :", i.cut.prev, "\n",
                 "Otu table: ", i.otu, "\n",
                 "Starta :", i.strata) %>% 
            stop(.)}
        #-----------------------------------------------------------------------
        
        #-----------------------------------------------------------------------
        # Correlation loop
        cor.res.comb <- NULL
        
        for (i.second in colnames(second.df.i.strat)) {
          
          cor.res <- NULL 
          
          for (i.tax in colnames(otu.tab.i.strat.f)) {
            
            cor.res <- cor.test(second.df.i.strat[, i.second],
                                otu.tab.i.strat.f[, i.tax],
                                method = corr_method) %>%
              tidy() %>%
              mutate(Taxa = i.tax) %>%
              bind_rows(cor.res, .)
            
          }
          
          cor.res.comb <- cor.res %>% 
            mutate(Corre_vector = i.second) %>% 
            bind_rows(cor.res.comb, .)
        }
        
        
        # Add strata and otu table names
        cor.res.comb <- cor.res.comb %>% 
          mutate(Strata = i.strata, 
                 OTU_table_name = i.otu)
        
        
        # Loop though prevalence cutoff and adjust q value per metabolite separately 
        for (i.prev in prev_cut_offs) {
          
          # Filter based on prevalence 
          sel.taxa <- prev_filt_otu_tab(otu_tab = otu.tab.i.strat.f, 
                                        meta_data = strata.data, 
                                        prev_cutoff = i.prev) %>% 
            colnames(.) 
          
          cor.res.comb.f <-  cor.res.comb %>% 
            group_by(Corre_vector) %>% 
            filter(Taxa %in% sel.taxa) %>% 
            mutate(Prevalence_cutoff = i.prev, 
                   qval = p.adjust(p.value, 
                                   method = p_adj_method)) %>% 
            ungroup()
          
          firs.ls.name <- paste0(i.otu, "--", i.prev)
          
          res.out[[firs.ls.name]][[i.strata]] <- cor.res.comb.f
          
          
        }
        
      }
      
    }
    
    return(res.out)
  }
  
  
  #-------------------------------------------------------------------------------
  # Plot scatter plot for correlations 
  #-------------------------------------------------------------------------------
  plot_bulk_corr <- function(otu_tab, 
                             second_df, 
                             meta_data, 
                             bulk_cor_res_long,
                             strata_cols_meta = NULL, 
                             qval_cutoff = 1, 
                             n_top_corr = nrow(bulk_cor_res_long),  
                             n_plots_per_line = 8, 
                             strip_text_size = 7) { 
    
    # Create data frame for stratification ---------------------------------------
    if(is.null(strata_cols_meta)) {
      
      strata.data <- data.frame(Strata = rep("No Stratification", 
                                             nrow(meta_data)))
      
    } else {  
      
      strata.data <- meta_data %>% 
        dplyr::select(all_of(strata_cols_meta)) %>% 
        unite(Strata, all_of(strata_cols_meta))
    }
    #-----------------------------------------------------------------------------
    
    
    # Filter data 
    #-----------------------------------------------------------------------------
    cor.res.f <- bulk_cor_res_long %>% 
      dplyr::filter(qval <= qval_cutoff) %>% 
      dplyr::arrange(desc(estimate)) %>% 
      slice_head(n = n_top_corr) %>% 
      mutate(combID = paste0(Taxa, "--", Corre_vector))
    
    if(nrow(cor.res.f) < 1) {
      
      return(NULL)
      
    } else {
      
      otu.tab.f <- otu_tab %>% 
        dplyr::select(all_of(unique(cor.res.f$Taxa)))
      
      second.df.f <- second_df %>% 
        dplyr::select(all_of(unique(cor.res.f$Corre_vector)))
      
      
      # Row names check-------------------------------------------------------------
      if(!identical(rownames(otu.tab.f), rownames(second.df.f)) & 
         identical(rownames(otu.tab.f), rownames(meta_data))) {
        stop("Row names in the input data are not matching.") }
      # ----------------------------------------------------------------------------
      
      plots.ls <- list()
      
      for(i.comb in unique(cor.res.f$combID)) {
        
        text.df <- bulk_cor_res_long %>% 
          mutate(combID = paste0(Taxa, "--", Corre_vector)) %>% 
          filter(combID == i.comb) %>% 
          mutate(lable = paste0("qval=", round(qval, 2), " / ", 
                                "est=", round(estimate, 2))) %>% 
          select(Strata, lable)
        
        i.tax <- gsub("--.*", "", i.comb)
        
        i.second <- gsub(".*--", "", i.comb)
        
        
        data.p <- bind_cols(dplyr::select(otu.tab.f, all_of(i.tax)), 
                            dplyr::select(second.df.f, all_of(i.second)), 
                            strata.data) %>% 
          mutate(Taxa = i.tax, 
                 y.text = max(get(i.tax)*1.1, na.rm = TRUE), 
                 x.text = (min(get(i.second), na.rm = TRUE) + 
                             max(get(i.second), na.rm = TRUE))/2 , 
                 y.border = max(get(i.tax)*1.2), na.rm = TRUE) %>% 
          left_join(., text.df, by = "Strata") %>% 
          mutate(lable = ifelse(is.na(lable), "Not tested.", lable))
        
        # Protect against weird names 
        x.lab <- paste0("`", colnames(data.p)[2], "`")
        
        y.lab <- paste0("`", colnames(data.p)[1], "`")
        
        p <- ggplot(data.p, aes_string(y = y.lab, 
                                       x = x.lab)) + 
          facet_grid(c("Taxa", "Strata")) +
          geom_point(alpha = 0.5) + 
          geom_smooth(method = "lm", se = FALSE) + 
          geom_text(aes(x = x.text, 
                        y = y.text, 
                        label = lable), size = 2.1, hjust = 0.5) + 
          geom_point(aes(x = x.text, y = y.border), shape = NA) + 
          theme_bw() + 
          theme(legend.position = "none", 
                axis.title.y = element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.text.x=element_blank(), 
                axis.ticks.y=element_blank(), 
                # axis.text.y=element_blank(), 
                panel.spacing = unit(0, "cm"), 
                strip.text.y = element_text(size = strip_text_size, 
                                            face = "italic")) + 
          xlab(label = colnames(data.p)[2]) 
        
        plots.ls[[paste0(i.tax, "_", i.second)]] <- p
        
      }
      
      
      # Calculate number of columns to plot-----------------------------------------
      n.strata <- length(unique(strata.data$Strata))
      
      if (n.strata < n_plots_per_line) {
        
        n.col <- floor(n_plots_per_line/n.strata)
        
      } else {n.col <- 1}
      
      if (n_plots_per_line > n.strata*nrow(cor.res.f)) {
        
        n.col <- nrow(cor.res.f)
        
      }
      
      plot.out <- arrangeGrob(grobs = plots.ls, ncol = n.col)
      
      return(plot.out)
      
    }
  }
  
# Funciton to plot #############################################################
cor_heat_plot <- function(cor_func_res, 
                            sup_df, 
                            he_hight = NULL, 
                            out_folder, 
                            qval_f_cut, 
                            estimate_cut) {
    
    dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
    
    
    
    sig.tax <- cor_func_res %>% 
      recurse_lapply(function(x){x %>% 
          filter(qval_f <= qval_f_cut, 
                 abs(estimate) >= estimate_cut) %>% 
          pull(Taxa)}) %>% 
      unlist() %>% unique()
    
    if(is.null(he_hight)) {
      he_hight <- length(sig.tax)*1
    }
    
    sig.vec <- cor_func_res %>% 
      recurse_lapply(function(x){x %>% 
          filter(qval_f <= qval_f_cut, 
                 abs(estimate) >= estimate_cut) %>% 
          pull(Corre_vector)}) %>% 
      unlist() %>% unique()
    
    cor_func_res.f <- cor_func_res %>% 
      recurse_lapply(function(x){ x %>% 
          filter(Taxa %in% sig.tax, 
                 Corre_vector %in% sig.vec)})
    
    for(i in names(cor_func_res.f)) {
      cor_func_res.f[[i]] %>% 
        bind_rows() %>% 
        write.csv(file = paste0(out_folder, 
                                "/heat_stat_", 
                                gsub("\\.", "_", trimws(i)), 
                                ".csv")) }
    
    est.mat.ls <- cor_func_res.f %>% 
      recurse_lapply(function(x){x %>% 
          select(Taxa, Corre_vector, 
                 estimate, Strata) %>% 
          pivot_wider(names_from = Corre_vector, 
                      values_from = estimate)})
    
    qp.mat.ls <- cor_func_res.f %>% 
      recurse_lapply(function(x){x %>% 
          select(Taxa, Corre_vector, 
                 qval_p, Strata) %>% 
          mutate(qval_p = ifelse(is.na(qval_p),
                                 "", qval_p)) %>% 
          pivot_wider(names_from = Corre_vector, 
                      values_from = qval_p)})
    
    
    for (i.mat in names(cor_func_res)) { 
      
      estim.inst <- est.mat.ls[[i.mat]] %>% 
        bind_rows() 
      
      if(all(dim(estim.inst) != 0)) {
        
        qval.inst <- qp.mat.ls[[i.mat]] %>% 
          bind_rows() 
 
        he.width <- (ncol(qval.inst) - 2)
          
          if(he.width < 2) {
            
            he.width <- 2}
          
        
        n.inst <- gsub("\\..*", "", i.mat)
        
        ht_opt$TITLE_PADDING = unit(c(5.5, 5.5), "points")
        
        heats.inst <- Heatmap(select(estim.inst, -Taxa, -Strata), 
                              column_title = n.inst, 
                              column_title_gp = gpar(fill = sup_df[n.inst, "color"], 
                                                     col = "black", 
                                                     border = "black", 
                                                     fontsize = 12, 
                                                     fontface = "bold"), 
                              show_row_dend = FALSE, 
                              name = "Estimate", 
                              width = unit(he.width, "cm"), 
                              height = unit(he_hight, "cm"), 
                              row_split = estim.inst$Strata, 
                              rect_gp = gpar(col = "gray45", lwd = 1),
                              column_names_rot = 45, 
                              cluster_rows = FALSE, 
                              show_column_names = TRUE, 
                              show_column_dend = FALSE, 
                              show_row_names = TRUE, 
                              row_labels = estim.inst$Taxa, 
                              cell_fun = function(j, i, x, y, w, h, fill) {
                                grid.text(select(qval.inst, -Taxa, -Strata)[i, j], x, y)})
        
        
        png(paste0(out_folder, "/heat_", n.inst, ".png"), 
            width = (he.width+15), 
            height = (he_hight+15), res = 600, units = "cm",)
        
        draw(heats.inst)
        
        dev.off()
        
      }
      
    }
    
  }
  ################################################################################
  
