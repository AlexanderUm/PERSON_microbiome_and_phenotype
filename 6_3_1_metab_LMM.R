#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls$General$Seed)


################################################################################
# Deferential abundance with LMM
################################################################################
#-------------------------------------------------------------------------------
# Variables and data to use 
#-------------------------------------------------------------------------------
# Parameter variables 
ModelType <- prm.ls$metab$model_type

MainEff <- prm.ls$metab$da_main_effect

FixEff <- c(prm.ls$metab$da_main_effect, prm.ls$metab$da_covar) 

RandEff <- prm.ls$metab$da_rand

RandEffSpecif <- prm.ls$metab$da_rand_specific

MetabTabs <- prm.ls$metab$da_tabs

TransCountForm <- prm.ls$metab$da_trans_count 

ParticipIDcol <- prm.ls$General$Part_id_col

DirOut <- prm.ls$metab$dir_out

# Data sets 
MetaDf <- data.ls$phen$meta

MetabSets <- data.ls$Metabol

# Empty objects for results 
StatResDf <- NULL

FitResDf <- NULL

StatResLs <- list()

MetabSetsForm <- list()


for(i in MetabTabs) { 
  
  # Data frame for results per data set 
  iStatResDf <- NULL
  
  # Extract metabolites set 
  iData0 <- MetabSets[[i]] 
  
  # Vector with metabolites names old and new names 
  iColNames <- colnames(iData0)[colnames(iData0) != ParticipIDcol]
  
  names(iColNames) <- paste0("met_", iColNames)
  
  # Format data for modeling 
  iData <- iData0 %>% 
              dplyr::rename(iColNames) %>% 
              left_join(., MetaDf, by = "ID") %>% 
              filter(!is.na(.data[[MainEff]])) %>% 
              droplevels() 
  
  iDataSum <- iData %>% 
                select(all_of(c(MainEff, FixEff, RandEff, RandEffSpecif, 
                                names(iColNames)))) %>% 
                pivot_longer(cols = names(iColNames))
  
  # Transform metabolite count if specified 
  if(!is.null(TransCountForm)) {
    
    iData <- iData %>% 
                mutate(across(all_of(names(iColNames)), TransCountForm))
    
  }
  
  MetabSetsForm[[i]] <- iData
  
  # Format random effect if present in parameters files 
  if(!is.null(RandEff)) {
    
    iRandEff <- RandEff
    
    if(i %in% names(RandEffSpecif)) {
      
      iRandEff <- c(RandEff, RandEffSpecif[i])
      
    }
      
    iRandEff <- iData[, iRandEff] %>% 
                    apply(2, function(x){length(unique(x))}) %>% 
                    .[. > 1] %>% 
                    names() %>% 
                    paste0("(1|", ., ")") %>% 
                    paste(collapse = " + ")
      
  } else {
    
    iRandEff <- NULL
    
  }
  
  # Right side of the formula
  iFormComb <- paste(c(FixEff, iRandEff), collapse = "+")
  
  for(j in names(iColNames)) { 
    
    jDataSum <- iDataSum %>% 
                  filter(name == j) %>% 
                  mutate(value_trans = TransCountForm(value)) %>% 
                  summarise(Median = median(value), 
                            IQR = IQR(value),
                            Median_Trans = median(value_trans), 
                            IQR_Trans = IQR(value_trans), 
                            .by = all_of(MainEff)) %>% 
                  pivot_wider(names_from = all_of(MainEff), 
                              values_from = c(Median, IQR, 
                                              Median_Trans, IQR_Trans,)) %>% 
                  mutate(Metabol = j)
    
    # Formula adjustment
    jForm <- paste0(j, "~", iFormComb)
    
    for(k in ModelType) {
      
      if(k == "LM") {
        
        jRes <- gsub("\\(1\\||\\)", "", jForm) %>% 
                  as.formula() %>% 
                  lm(iData) 
        
      } else {
        
        if(!is.null(iRandEff)) {
          
          jRes <- lmerTest::lmer(as.formula(jForm), iData) 
          
        }
        
      }
      
      # Model results in a single data frame
      iStatResDf <- broom.mixed::tidy(jRes) %>% 
                                mutate(Set = i, 
                                       Method = k,
                                       Metabol = j, 
                                       Formula = jForm) %>% 
                                left_join(., jDataSum, by = "Metabol") %>% 
                                bind_rows(iStatResDf, .)
      
      # Model results in a list 
      StatResLs[[i]][[j]][[k]] <- jRes
      
      # Fitted and residuals 
      FitResDf <- data.frame(Fitted = fitted(jRes), 
                             Residuals = resid(jRes), 
                             Method = k,
                             Set = i, 
                             Index = j) %>% 
                   bind_rows(FitResDf, .)
    }
    
  }
  
  StatResDf <- bind_rows(StatResDf, iStatResDf)
  
  dir.create(paste0(DirOut, "/", i, "/tabs/LMM/"), 
             showWarnings = FALSE,
             recursive = TRUE)
  
  write.csv(iStatResDf, file = paste0(DirOut, "/", i, "/tabs/LMM/metabol_lmm.csv"))
  
}

save(list = c("FitResDf", "StatResLs", "MetabSetsForm", "StatResDf"), 
     file = "out/supp/6_3_1_metab_LMM.Rdata")

# Clean environment
rm(list = ls())
gc()



