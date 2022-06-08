

largeDfLoad <- function(filename){
FDRandSensitivity <- data.frame()


for(i in seq(runs)){

  f <- function(data, pos) {
    data <- data %>% filter(run == i)
    if(nrow(data) > 0){
      conf <- calculateConfusion(data)
      conf$method <- unlist(strsplit(conf$method, "_"))[!unlist(strsplit(conf$method, "_")) == "pvalsAdj"] # Remove .pvals suffix
    } else{
      conf <- NULL
    }
    # FDRandSensitivity <- calculateFDRandSensitivity(conf)
    return(conf)
  }



  data <- read_delim_chunked(file = filename,
                             callback = DataFrameCallback$new(f),
                             chunk_size = chunk_size, delim = c(","), col_names = T, na = c("NA", NA)
  ) # na = c("NA", "NA")

  stats <- calculateFDRandSensitivity(data)
  FDRandSensitivity <- rbind(FDRandSensitivity, stats)

}
return(FDRandSensitivity)
}
