#-------------------------------------------------------------------------------
# Load data
#------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls[["General"]][["Seed"]])

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_taxa_filter.R")
source("R/phy_alpha.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Data sets to use
SampleSet <- prm.ls$Alpha$data_set_ps

TaxLvl <- prm.ls$Alpha$Tax_lvl

norms <- prm.ls$Alpha$Norm

AlphaCountTrans <- prm.ls$Alpha$measures_trans


# Alpha parameters
AlphaIndx <- prm.ls$Alpha$measures

p.text.size <- 3.1


################################################################################
# Run alpha diversity analysis
################################################################################
alpha.res.ls <- list()

StatResComb <- NULL

FitRes <- NULL

for(i in SampleSet) { 
  
  # Extract data
  PsInst <- data.ls[[i]][["PS"]][[TaxLvl]][[norms]]
  
  MetaInst <- data.ls[[i]][["meta"]]
  
  # Calculate alpha diversity
  AlphaDfInst <- phy_alpha(PsInst, measures = AlphaIndx) %>% 
                    bind_cols(., MetaInst)
  
  # Transform index count if specied 
  if(!is.null(AlphaCountTrans)) {
    
    for(j in 1:length(AlphaCountTrans)) {
      
    AlphaDfInst <- AlphaDfInst %>% 
                    mutate(across(AlphaCountTrans[[j]]$measures, 
                                  AlphaCountTrans[[j]]$trans))
    
    }
  
  }
  
  # Elongate the diversity table 
  AlphaDfLongInst <- AlphaDfInst %>% 
                        pivot_longer(cols = all_of(AlphaIndx), 
                                     names_to = "Index")
  
  # Create output directories
  dir.create(paste0(prm.ls$Alpha$out_dir, "/", i, "/tabs/"), 
             recursive = TRUE, showWarnings = FALSE)
  
  dir.create(paste0(prm.ls$Alpha$out_dir, "/", i, "/plots/"), 
             recursive = TRUE, showWarnings = FALSE)
  
  
  #-----------------------------------------------------------------------------
  # Statistics 
  #-----------------------------------------------------------------------------
  FormLmm <- paste(c(prm.ls$Alpha$main_var, 
                     prm.ls$Alpha$fix_covar[[i]], 
                     paste0("(1|", prm.ls$Alpha$rand_vars[[i]], ")")), 
                   collapse = " + ")
  
  StatResInst <- NULL
 
  for(j in AlphaIndx) { 
      
      FormInst <- paste0(j, " ~ ", FormLmm)
      
      ResLmm <- lmerTest::lmer(as.formula(FormInst), AlphaDfInst)
      
      # Collect results
      StatResInst <- broom.mixed::tidy(ResLmm) %>% 
                          mutate(across(where(is.numeric), 
                                        function(x){round(x, digits = 4)})) %>% 
                          mutate(Set = i, 
                                 Index = j, 
                                 Formula = FormInst) %>% 
                          bind_rows(StatResInst, .)
      
      # Assumptions data 
      FitRes <- data.frame(Fitted = fitted(ResLmm), 
                           Residuals = resid(ResLmm), 
                           Set = i, 
                           Index = j) %>% 
                      bind_rows(FitRes, .)
  }
  
  StatResComb <- bind_rows(StatResComb, StatResInst)
    
  write.csv(StatResInst, 
            paste0(prm.ls$Alpha$out_dir, "/", i, "/tabs/", "test_res.csv"), 
            row.names = FALSE)
  
  #-----------------------------------------------------------------------------
  # Plot results
  #-----------------------------------------------------------------------------
  # Base plot
  p.base <- ggplot(AlphaDfLongInst, 
                   aes(y = value, x = .data[[prm.ls$Alpha$main_var]])) + 
              geom_jitter(aes(colour = .data[[prm.ls$Alpha$main_var]]), 
                          height = 0, 
                          width = 0.15, 
                          alpha = 0.5) +
              geom_violin(fill = NA) + 
              facet_grid(c("Index"), scales = "free") + 
              theme_bw() + 
              theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.title.y = element_blank())
  
  # Plot annotation dataframe
  min.max.df <- AlphaDfLongInst %>% 
                    reframe(y.min = min(value), 
                            y.max = max(value), 
                            .by = Index)
  
  SigDf <- StatResInst %>% 
              filter(grepl(prm.ls$Alpha$main_var, term)) %>% 
              mutate(p.short = round(p.value, 3), 
                     End = gsub(prm.ls$Alpha$main_var, "", term), 
                     Start = levels(MetaInst[[prm.ls$Alpha$main_var]])[1]) %>% 
              mutate(p.text = ifelse(p.short == 0, 
                                     "P<0.001", 
                                     paste0("P=", sprintf("%.3f", 
                                                          p.short)))) %>% 
              left_join(., min.max.df, by = "Index") %>% 
              mutate(y.sig = y.max + ((y.max - y.min)*(1:n()*0.35)),
                     y.inv.point = y.max + ((y.max - y.min)*((n()+0.5)*0.35)), 
                     .by = Index)
  
  FinalPlot <-  p.base +  
                  geom_signif(data = SigDf,
                              aes(xmin = Start,
                                  xmax = End,
                                  annotations = p.text,
                                  y_position = y.sig),
                              textsize = p.text.size, 
                              vjust = -0.2,
                              manual = TRUE, 
                              margin_top = 1) + 
                  geom_point(data = SigDf,
                             aes(x = End, 
                                 y = y.inv.point), 
                             x=NA) + 
                  scale_color_manual(values = aes.ls$col[[i]])
                
  # Save plot 
  ggsave(filename = paste0(prm.ls$Alpha$out_dir, "/", i, "/plots/alpha.png"), 
         plot = FinalPlot, 
         width = length(levels(AlphaDfInst[[prm.ls$Alpha$main_var]]))*1 + 1.5, 
         height = length(AlphaIndx)*2 + 0.25, 
         dpi = 600)
  
  # Write out results 
  alpha.res.ls[[i]][["plot"]] <- FinalPlot
  
  alpha.res.ls[[i]][["res"]] <- StatResInst
    
}
 

#-------------------------------------------------------------------------------
# Visualize models assumptions
#-------------------------------------------------------------------------------
FitResPlot <- ggplot(FitRes, aes(y = Residuals, x = Fitted)) + 
                geom_point(color = "steelblue", alpha = 0.5) + 
                geom_abline(intercept = 0, slope = 0) + 
                facet_wrap(c("Set", "Index"), 
                           scales = "free", 
                           ncol = 4) + 
                theme_bw() 

ggsave(plot = FitResPlot, 
       filename = paste0(prm.ls$Alpha$out_dir, 
                         "/Full_Fitted_vs_Residual.png"), 
       width = 12, 
       height = 9) 


QqResPlot <- ggplot(FitRes, aes(sample = Residuals)) + 
                stat_qq(color = "steelblue", alpha = 0.5) + 
                stat_qq_line() + 
                facet_wrap(c("Set", "Index"), 
                           scales = "free", 
                           ncol = 4) + 
                theme_bw()

ggsave(plot = QqResPlot, 
       filename = paste0(prm.ls$Alpha$out_dir, 
                         "/Full_QQ_Residual.png"), 
       width = 12, 
       height = 9)

#-------------------------------------------------------------------------------
save(list = c("alpha.res.ls"), 
     file = "out/supp/2_alpha.Rdata")

# Clean environment
rm(list = ls())
gc()
