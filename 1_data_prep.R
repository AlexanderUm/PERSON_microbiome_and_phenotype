#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")

# Set Seed 
set.seed(prm.ls[["General"]][["Seed"]])

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_shorten_tax_names.R")
source("R/phy_css_norm.R")


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Paths to data
seq.path <- prm.ls$Paths$path_seq

picrust.path <- prm.ls$Paths$path_picrust2

meta.paths <- prm.ls$Paths$meta

metabol.paths <- prm.ls$Paths$metabol

# Filtering 
min.r <- prm.ls$Data$min_read_tax

metabol.prev <- prm.ls$Data$metabolite_min_prev

# Summary
tax.lvls <- prm.ls$Data$tax_lvls

# Columns 
seq.id.col <- prm.ls$General$SeqID_col

gr.col <- prm.ls$General$Group_col

sex.col <- prm.ls$General$Sex_col

center.col <- prm.ls$General$Center_col

antib.col <- prm.ls$General$Antib_col

id.col <- prm.ls$General$Part_id_col

# Object to fill up 
data.ls <- list()

###############################################################################
# Read data/ Create phyloseq
#-------------------------------------------------------------------------------
ps <- qza_to_phyloseq(features = paste0(seq.path, "asv_table.qza"), 
                      tree = paste0(seq.path, "tree/rooted-tree.qza"), 
                      taxonomy = paste0(seq.path, "taxonomy_07.qza"))

#-------------------------------------------------------------------------------
# Filter out taxa
#-------------------------------------------------------------------------------
# filter taxa with with less than X reads in total   
ps <- prune_taxa(taxa_sums(ps) >= min.r, ps)

# Remove ASVs: 
ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Mitochondria", ps)

ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Chloroplast", ps)

ps <- prune_taxa(!is.na(tax_table(ps)[, "Phylum"])[, "Phylum"], ps)

ps <- prune_taxa(tax_table(ps)[, "Kingdom"] %in% c("d__Bacteria", 
                                                   "d__Archaea"), ps)


################################################################################
# Glom to higher taxonomic level ->
# Normalize count ->
# Make strata
#-------------------------------------------------------------------------------
ps.ls <- list()

for(i.lvl in tax.lvls)  {
  
  if(i.lvl == "ASV") {
    
    ps.inst <- ps
  
    # Adjust taxa names
    taxa_names(ps.inst) <- phy_shorten_tax_names(ps.inst) %>% 
                              as.data.frame() %>% 
                              setNames("feature") %>% 
                              mutate(feature2 = if(n() > 1) {
                                paste0(feature, "__asv", row_number())} else {
                                  paste0(feature, "__asv")}, .by = "feature") %>% 
                              pull(feature2)
  
  } else { 
    
    ps.inst <- tax_glom(ps, i.lvl, NArm = FALSE) 
    
    # Adjust taxa names
    taxa_names(ps.inst) <- phy_shorten_tax_names(ps.inst) %>% 
                                as.data.frame() %>% 
                                setNames("feature") %>% 
                                mutate(feature2 = if(n() > 1) {
                                  paste0(feature, "__", 
                                         tolower(str_sub(i.lvl, 1, 1)), 
                                         row_number())} else {
                                           feature}, .by = "feature") %>% 
                                pull(feature2)
    }

  ps.ls[[i.lvl]] <- list("Raw" = ps.inst, 
                         "Rare" = rarefy_even_depth(ps.inst, rngseed = 347), 
                         "CSS" = phy_css_norm(ps.inst))
}


################################################################################
# Construct phyloseq with PICRUST data 
################################################################################
picr.tab <- read_tsv(picrust.path, show_col_types = FALSE) %>% 
                as.data.frame() %>% 
                column_to_rownames(var = "pathway")

# Construct phyloseq with PICRUST2 data 
picr.otu <-  otu_table(picr.tab, taxa_are_rows = TRUE)

picr.tax <-  as.matrix(rownames(picr.tab)) %>% 
              'rownames<-'(rownames(picr.tab)) %>% 
              'colnames<-'("Pathway") %>% 
              tax_table()

ps.picr <- phyloseq(picr.otu, picr.tax)

# Normolize count and Place into the same list as other phyloseqs
ps.ls[["Picrust"]] <- list("Raw" = ps.picr, 
                           "Rare" = rarefy_even_depth(ps.picr, rngseed = 347), 
                           "CSS" = phy_css_norm(ps.picr))


################################################################################
# Add full metadata 
################################################################################
# Metadata 
phys.grid <- expand.grid(names(ps.ls), 
                         names(ps.ls[[1]]), 
                         stringsAsFactors = FALSE)

for(i.meta.path in names(meta.paths)) {
  
  i.name <- gsub("meta_", "", i.meta.path)
  
  meta <- read.csv(meta.paths[[i.meta.path]]) 

  rownames(meta) <- meta[[seq.id.col]]
  
  # Overlapping samples - base phyloseq object is used
  over.samp <- intersect(sample_names(ps), rownames(meta))
  
  # Adjust metadata
  meta <- meta[over.samp, ] %>% 
            mutate(!!gr.col := as.factor(.[[gr.col]]), 
                   !!center.col := as.factor(.[[center.col]]), 
                   !!sex.col := as.factor(.[[sex.col]]), 
                   !!antib.col := as.factor(.[[antib.col]]))
  
  
  data.ls[[i.name]][["meta"]] <- meta
  
  # Add metadata for each phyloseq
  for(i.ps in 1:nrow(phys.grid)) {
    
    i.lvl <- phys.grid[i.ps, 1]
    
    i.norm <- phys.grid[i.ps, 2]
    
    ps.inst <- ps.ls[[i.lvl]][[i.norm]]
    
    if(!is.null(ps.inst)) {
      
    ps.inst <- prune_samples(over.samp, ps.inst)
  
    ps.inst <- prune_taxa(taxa_sums(ps.inst) > 0, ps.inst)
    
    sample_data(ps.inst) <- meta
    
    data.ls[[i.name]][["PS"]][[i.lvl]][[i.norm]] <- ps.inst
      
    }
    
  }
  
}

################################################################################
# Metabolites table
################################################################################

for(i.metabol.path in names(metabol.paths)) {
  
  i.path <- metabol.paths[[i.metabol.path]]
  
  metabol.inst <- read_csv(i.path, show_col_types = FALSE)
  
  # Filter out low abundant metabolites
  met.to.keep <- (colSums(!is.na(metabol.inst))/nrow(metabol.inst) >= metabol.prev)
  
  metabol.inst <- metabol.inst[, met.to.keep]
  
  metabol.inst[is.na(metabol.inst)] <- 0 
  
  data.ls[["Metabol"]][[i.metabol.path]] <- metabol.inst
  
}


#-------------------------------------------------------------------------------
# Aesthetics 
#-------------------------------------------------------------------------------
aes.ls <- list()

col.vec <- c("#7FC97F", "#BEAED4", "#386CB0", "#BF5B17","#FDC086",
             "#1B9E77", "#7570B3", "#66A61E", "#E6AB02",
             "#A6761D", "#A6CEE3", "#33A02C", "#FB9A99", "#FDBF6F",
             "#FF7F00", "#6A3D9A", "#FFFF99", "#B15928", "#CCEBC5",
             "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
             "#B3E2CD", "#FDCDAC", "#CBD5E8", "#E6F5C9", "#FFF2AE",
             "#F1E2CC", "#CCCCCC", "#E41A1C", "#984EA3", "#FF7F00",
             "#FFFF33", "#A65628", "#66C2A5", "#FC8D62", "#8DA0CB", 
             "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#8DD3C7", 
             "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
             "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#999999")


aes.ls[["col"]][["with_contr"]] <- c("No IR"= col.vec[1], 
                                     "LIR" = col.vec[18], 
                                     "MIR" = col.vec[16])

aes.ls[["col"]][["phen"]] <- c("LIR" = col.vec[18], 
                               "MIR" = col.vec[16])


#-------------------------------------------------------------------------------
save(list = c("data.ls", "aes.ls"), 
     file = "out/supp/data.Rdata")

# Clean environment
rm(list = ls())
gc()
