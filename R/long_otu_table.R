# Function to extract and elongate OTU table 
long_otu_table <- function(phy, 
                           names_to = "SeqID", 
                           values_to = "Abundance") {
  
  LongTable <- phy %>% 
    otu_table() %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    mutate(feature = rownames(.)) %>% 
    pivot_longer(cols = - feature, 
                 names_to = names_to, 
                 values_to = values_to)
  return(LongTable)
  
}