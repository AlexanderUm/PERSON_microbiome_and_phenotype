#-----------------------------------------------------------------------------
# Function to make short names. 
#-----------------------------------------------------------------------------

phy_shorten_tax_names <- function(phy_obj, 
                                  clarify_taxa = c("A2", "RF39", "V9D2013", 
                                                   "S5_A14a", "P5D1_392", "TM7x", 
                                                   "WPS_2", "CHKCI002", "CHKCI001")) {
  
  char.rep <- function(in_str) {
    
    in_str <- gsub(" ", ".", in_str) 
    
    in_str <- gsub("_", "__", in_str)
    
    in_str <- gsub("-", "_", in_str)
    
    in_str <- gsub("\\[|\\]", "", in_str)
    
    return(in_str)
    
    }
  
  tax.t <- as(tax_table(phy_obj), "matrix")
  
  tax.t[tax.t == ""] <- NA
  
  tax.out <- c()
  
  excle.names <- paste(c("__NA", "uncultured", "gut metagenome",
                         "metagenome", "human__gut", "uncultured__archaeon", 
                         "human gut",  "swine fecal",  "swine__fecal", 
                         "<empty>", "Incertae__Sedis", "unidentified"), 
                       collapse = "|")
  
  tax.lvl.pref <- substr(colnames(tax.t), 1, 1)
  
  
  #------------------------------------------------------------------------------
  # Uncultured genus fix
  uncult.pref <- c("__GCA", "__UCG", "__UBA", 
                   "__A2", "__RF39", "__NK4A214", 
                   "__DTU", "__CAG")
  
  ucg.pref <- NULL
  
  for(i.pref in uncult.pref) {
    
    ucg.inst <- paste(paste0(tax.lvl.pref, i.pref), collapse = "|")  
    
    ucg.pref <- c(ucg.pref, ucg.inst)
    
  }
  
  ucg.pref <- paste(ucg.pref, collapse = "|")
  
  
  #-----------------------------------------------------------------------------
  # Fix names 
  #-----------------------------------------------------------------------------
  for (r in 1:nrow(tax.t)) {
    
    ind.tax <- paste0(tax.lvl.pref, "_", tax.t[r, ])
    
    ind.tax <- char.rep(ind.tax)
    
    l.tax <-  tail(ind.tax[!grepl(excle.names, 
                                  ind.tax)], 1)
    
    if(length(l.tax) == 0) {
      
      l.tax <- NA
      
    }
    
    # Unidentified genus fix 
    if(grepl(ucg.pref, l.tax)) {
      
      no.ucg.tax <- ind.tax[!grepl(ucg.pref, ind.tax)] 
      
      l.tax.p <- tail(no.ucg.tax[!grepl(excle.names, 
                                        no.ucg.tax)], 1)
      
      l.tax <- paste0(l.tax.p, "_", gsub(".*_", "", l.tax))
      
    }
    
    #---------------------------------------------------------------------------
    # Clarify user specified taxa 
    #---------------------------------------------------------------------------
    if(l.tax %in% clarify_taxa) {

      ind.tax.c <- ind.tax[!grepl(excle.names, ind.tax)]
      
      l.tax <- paste0(tail(ind.tax.c[!ind.tax.c %in% l.tax], 1), 
                      "__", 
                      gsub(paste(tax.lvl.pref, collapse = "__|"), "", l.tax))

    }
    
    tax.out <- c(tax.out, l.tax)
    
  }
  
  return(tax.out)
  
}

