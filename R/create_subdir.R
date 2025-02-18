#-------------------------------------------------------------------------------
# Function to create subdirectories
#-------------------------------------------------------------------------------
create_subdir <- function(main_dir, sub_dir_ls) {
  
  DirsOut <- expand.grid(sub_dir_ls, stringsAsFactors = FALSE)
  
  for(i in 1:nrow(DirsOut)) {
    
    dir.create(paste0(main_dir, "/", 
                      paste(DirsOut[i, ], collapse = "/")),
               recursive = TRUE, showWarnings = FALSE) 
    
  }
  
}