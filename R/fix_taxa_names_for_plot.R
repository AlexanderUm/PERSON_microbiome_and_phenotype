#-------------------------------------------------------------------------------
# Fix names of taxa for more elegant ploting 
#-------------------------------------------------------------------------------
fix_taxa_names_for_plot <- function(names_vector, 
                                    separator_to_replace = "__", 
                                    new_separator = " ",
                                    pref_to_replace = c("S__", "G__", "F__", 
                                                        "C__", "O__", "P__", 
                                                        "K__"), 
                                    new_pref = c("S. ", "G. ", "F. ", 
                                                 "C. ", "O. ", "P. ", 
                                                 "K. ")) {
 ReturnVec <- names_vector
   
 for(i in 1:length(pref_to_replace)) {
   
   ReturnVec <- sub(pref_to_replace[i], new_pref[i], ReturnVec)
   
 }
  
 ReturnVec <- gsub(separator_to_replace, new_separator, ReturnVec)
 
 return(ReturnVec)
  
}
