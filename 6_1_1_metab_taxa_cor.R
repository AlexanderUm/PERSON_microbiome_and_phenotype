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

# Function: 
# Correlate each column of one data frame with each column of the second.
dfs_correlation <- function(df_1, df_2, 
                            method = "spearman", 
                            exact = FALSE) {
  
  if(!identical(rownames(df_1), rownames(df_2))) {
    
    stop("Row names are not identical between dataframes.")
    
  }
  
  
  ResDf <- NULL
  
  for(i in colnames(df_1)) {
    
    for(j in colnames(df_2)) {
      
      ResDf <-  cor.test(df_1[[i]], df_2[[j]], 
                         method = method, 
                         exact = exact) %>% 
                    tidy() %>% 
                    mutate(x = i, y = j) %>% 
                    bind_rows(ResDf, .)
      
    }
    
  }
  
  return(ResDf)
  
}


################################################################################
# Correlation analysis
################################################################################

# Variables 
CorStrata <- prm.ls$metab$cor_strata

CorMethod <- prm.ls$metab$cor_method

BalanceGroupSize <- prm.ls$metab$cor_balance_size

# Correlation between metabolites and individual taxa 
CorGrid <- prm.ls$metab[c("cor_tax_lvl", 
                          "cor_tax_norm", 
                          "cor_metab_tab")] %>% 
                  expand.grid(stringsAsFactors = FALSE)

CorResDf <- NULL

CorResLs <- list()

FiltTabs <- list()

for(i in 1:nrow(CorGrid)) {
  
  # Variables 
  iTaxLvl <- CorGrid[i, "cor_tax_lvl"]
  
  iTaxNorm <- CorGrid[i, "cor_tax_norm"]
  
  iMetabTabId <- CorGrid[i, "cor_metab_tab"]
  
  # Tables for correlation
  iPs <- data.ls$phen$PS[[iTaxLvl]][[iTaxNorm]]
  
  if(taxa_are_rows(iPs)) {
    
    iPs <- phyloseq::t(iPs)
    
  }
  
  iOtuTab <- iPs %>% 
              otu_table() %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              suppressMessages() # "phylo" class message 
  
  rownames(iOtuTab) <- gsub("CIW1", "", rownames(iOtuTab))
  
  # Metabolites table
  iMetabTab <- data.ls$Metabol[[iMetabTabId]] %>% 
                    as.data.frame() %>% 
                    column_to_rownames(var = "ID")
  
  # Select only overlapping rows 
  iOverSamp <- intersect(rownames(iOtuTab), rownames(iMetabTab))
  
  iMetabTabF <- iMetabTab[iOverSamp, ]
  
  iOtuTabF <- iOtuTab[iOverSamp, ]
  
  iMeta <- data.ls$phen$meta %>% 
              'rownames<-'(gsub("CIW1", "", rownames(.))) %>% 
              .[iOverSamp, ]
  
  
  # Collect filtered tabs
  FiltTabs[["otu"]][[paste0(iTaxLvl, "_", iTaxNorm)]] <- iOtuTabF
  
  FiltTabs[["met"]][[iMetabTabId]] <- iMetabTabF
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Shuffle correlation analysis 
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # iMeta[[CorStrata]] <- sample(iMeta[[CorStrata]])
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  #-----------------------------------------------------------------------------
  # Correlation analysis 
  #-----------------------------------------------------------------------------
  iCorRes <- NULL 
  
  # Extract levels of the grouping columns if define otherwise return NULL
  if(!is.null(CorStrata)) {
    
    jVec <- levels(iMeta[[CorStrata]]) 
    
  } else {
    
    jVec <- NULL
    
  }
  
  # Make correlation analysis per metadata group 
  for(j in jVec) { 
    
    # Extract samples per group if define otherwise take all
    if(is.null(jVec)) {
      
      jSamp <- rep(TRUE, nrow(iMeta))
      
    } else {
      
      jSamp <- (iMeta[[CorStrata]] == j)
      
    }
    
    # Subset data
    jMetDf <- iMetabTabF[jSamp, ]
    
    jOtuDf <- iOtuTabF[jSamp, ]
    
    # Test - balance group size between group
    #TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    if(!is.null(CorStrata) & BalanceGroupSize) {
      
      jSampSet <- sort(sample(1:sum(jSamp), min(table(iMeta[[CorStrata]]))))
      
      jMetDf <-  jMetDf[jSampSet, ]
      
      jOtuDf <-  jOtuDf[jSampSet, ]
      
    }
    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    
    # Correlation test 
    iCorRes <- dfs_correlation(df_1 = jOtuDf, 
                               df_2 = jMetDf, 
                               method = CorMethod) %>% 
                      mutate(Strata = paste0(CorStrata, "_", j)) %>% 
                      bind_rows(iCorRes, .)
    
  }
  
  CorResDf <- iCorRes %>% 
                bind_cols(CorGrid[i, ]) %>% 
                bind_rows(CorResDf, .)
  
  CorResLs[[paste0(iTaxLvl, "_", iTaxNorm)]][[iMetabTabId]] <- iCorRes
  
  dir.create(paste0(prm.ls$metab$dir_out, "/", iMetabTabId, "/tabs/"), 
             showWarnings = FALSE, recursive = TRUE)
  
  write.csv(iCorRes, 
            paste0(prm.ls$metab$dir_out, "/", iMetabTabId, "/tabs/", 
                   iMetabTabId, "--", CorMethod, "_cor.csv"))
  
}


#-------------------------------------------------------------------------------
save(list = c("CorResLs", "FiltTabs"), 
     file = "out/supp/6_1_1_met_corr.Rdata")

# Clean environment
rm(list = ls())
gc()


