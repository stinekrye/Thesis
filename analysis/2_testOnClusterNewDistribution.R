####################################################################################################
########################################## Load libraries ##########################################
####################################################################################################

library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(ANCOMBC)
library(tidyverse)




####################################################################################################
########################################### Set variables ##########################################
####################################################################################################

# Variables for simulation of data
change <- c(10)#c(1,3,5,7,10)#c(1,2,3,4,5,6,7,8,9,10)
sample_size <- c(9) #<- c(2,3,4,5,6,7,8,9)#c(2,3,4,5,6,7,8,9)
runs <- 10
n_taxa <- 10000
n_random <- c(150) #c(0,5,50,75,150)#c(0, 5, 50, 75, 150) #c(0, 3, 15, 30, 75, 150)
v <- 2
s <- 120

# Choose tests
tests <- c("ANCOM.BC", "DESEq2", "ALDEx2", "Student.t.test", "Wilcoxon")#c("ANCOM.BC", "DESEq2", "ALDEx2", "Student.t.test", "Wilcoxon")

# Declare distribution used for simulating data (Used in name generation)
distribution <- "LogNormalAndNegBinom3"


####################################################################################################
######################################## Helper functions ##########################################
####################################################################################################

createDfLogNormalAndNegBinom3 <- function(variables, sample_size, n_taxa, run, change, n_random, seed) { # DO NOT CHANGE! MAKE NEW ONE IF NEEDED
  set.seed(seed)
  data <- data.frame(matrix(nrow = n_taxa , ncol = variables*sample_size))
  colnames(data) <- c(sort(paste(as.vector(outer(seq(variables), seq(sample_size), paste, sep = ".")), sep = "")))
  zeroesProb <- rnorm(1, mean = 0.97, sd = 0.01)
  trueZero <- rbinom(n = n_taxa, size = 1, prob = zeroesProb)
  # marginalProb <- rlnorm(n_taxa, meanlog = log(1), sdlog = log(2.718282*7))
  marginalProb <- ifelse(trueZero==1, 0.01, rlnorm(sum(trueZero==0), meanlog = log(3), sdlog = log(2.718282*6)))


  data$taxon <- seq(n_taxa)
  data$sample_size <- sample_size
  data$n_taxa <- n_taxa
  data$n_random <- n_random
  data$run <- run
  data$change <- change
  data$diff <- FALSE

  data <- data %>% select(taxon, run, sample_size, n_taxa, n_random, change, diff, everything())

  # Counters
  uc <- which(startsWith(colnames(data), "1.1"))
  ec <- which(startsWith(colnames(data), "2.1"))

  # In order to make sure that the counts are consistant even when the sample size changes, it is important to first simulate 1.2 then 2.1 before 2.1 and 2.2 and so on. It is easier to just assign the data to the correct column.
  # The point is to get all the uneven numbers first and then all the even numbers.
  for (j in which(startsWith(colnames(data), "1.1")):ncol(data)) {
    # Generate count data for each taxon (SAMPLES)
    if(j %% 2 == 0){ # Even numbers
      data[ec] <- as.integer(rnbinom(n_taxa, mu = marginalProb, size = 100)) # Size is the shape parameter af the gamma distrib. The lower size the more variance is introduced
      ec <- ec + 1
    } else { # Uneven numbers
      data[uc] <- as.integer(rnbinom(n_taxa, mu = marginalProb, size = 100)) # Size is the shape parameter af the gamma distrib. The lower size the more variance is introduced
      uc <- uc + 1
    }}

  return(data)
}

# Makes n_taxa number of taxa have different abundance
diffAbundance <- function(df, change, random_taxon, variables) {
  if(length(random_taxon) == 0) {
    return(df)
  } #else {
  # for(i in 1:length(random_taxon)) {
  #     random_taxon[i] <- which.min(abs(as.integer(row.names(df))-random_taxon[i]))
  # }
  df[random_taxon, startsWith(colnames(df), "2")] <- df[random_taxon, startsWith(colnames(df), "2")]*change
  df[random_taxon,]$diff <- TRUE
  return(df)
}



relativeAbundanceParalel <- function(data) {
  for(i in colnames(data)[startsWith(colnames(data), "2") | startsWith(colnames(data), "1")]){
    data[[i]] <- data[[i]]/sum(data[[i]])
  }
  return(data)
}


cleanData <- function(data){
  d <- data[,startsWith(colnames(data), "2") | startsWith(colnames(data), "1")]
  boolean <- d > 0

  rsums <- rowSums(boolean, na.rm = T)
  rkeep <- rsums/ncol(d) > 0.05

  data <- data[rkeep,]


  csums <- colSums(boolean, na.rm = T)
  ckeep <- csums > 0

  data <- data[,ckeep]
  return(data)
}




########################################################
################# FUNCTION TO RUN TESTS ################
########################################################

# Recreate run parameter somehow
runTests <- function(tests, countData, meta, variables, seed) {
  set.seed(seed)

  # Define variables and phyloseq object
  n_taxa <- dim(countData)[1] # number of taxa
  n_samples_total <- dim(countData[startsWith(colnames(countData), "2") | startsWith(colnames(countData), "1")])[2]
  n_tests <- length(tests)
  phyloData <- phyloseq(otu_table(countData[startsWith(colnames(countData), "2") | startsWith(colnames(countData), "1")], taxa_are_rows = T), sample_data(meta)) # Bruges i ANCOM-BE

  # Make dataframes for p-values and effect sizes
  pvals <- cbind(data.frame("taxon" = countData$taxon,
                            "run" = countData$run,
                            "sample_size" = countData$sample_size,
                            "n_taxa" = countData$n_taxa,
                            "n_random" = countData$n_random,
                            "change" = countData$change,
                            "diff" = countData$diff),
                 data.frame(matrix(nrow = n_taxa, ncol = n_tests, dimnames = list(seq(n_taxa), tests))))
  pvalsAdj <- cbind(data.frame("taxon" = countData$taxon,
                               "run" = countData$run,
                               "sample_size" = countData$sample_size,
                               "n_taxa" = countData$n_taxa,
                               "n_random" = countData$n_random,
                               "change" = countData$change,
                               "diff" = countData$diff),
                    data.frame(matrix(nrow = n_taxa, ncol = n_tests, dimnames = list(seq(n_taxa), tests))))
  effect <- cbind(data.frame("taxon" = countData$taxon,
                             "run" = countData$run,
                             "sample_size" = countData$sample_size,
                             "n_taxa" = countData$n_taxa,
                             "n_random" = countData$n_random,
                             "change" = countData$change,
                             "diff" = countData$diff),
                  data.frame(matrix(nrow = n_taxa, ncol = n_tests, dimnames = list(seq(n_taxa), paste(tests, "_effect", sep = "")))))

  # Remove idenfication columns - keep only counts
  countData <- countData[startsWith(colnames(countData), "2") | startsWith(colnames(countData), "1")]
  relData <- relativeAbundanceParalel(countData)

  if ("DESEq2" %in% tests) {
    message("##### DESEq2 #####")
    DESEQ2 <- DESeqDataSetFromMatrix(countData = countData,
                                     colData = meta,
                                     design= ~ treatment)

    DESEQ2 <- DESeq(DESEQ2)
    DESEQ2_res <- results(DESEQ2)
    pvals$DESEq2 <- DESEQ2_res$pvalue
    pvalsAdj$DESEq2 <- DESEQ2_res$padj
    effect$DESEq2_effect <- DESEQ2_res$log2FoldChange

  }
  if ("ALDEx2" %in% tests) {
    # https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html#5_ALDEx2_outputs
    message("##### ALDEx2 #####")
    m <- as.vector(unlist(meta[1])) # Convert metaData data frame to vector

    resALDEX2 <- aldex(countData, m, mc.samples=64, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE) # Object contain reaslts from many different tests. See link for
    res <- merge(x = countData, y = resALDEX2, by = 0, all = T)
    res <- merge(x = pvals, y = res, by.x = "taxon", by.y =  "Row.names")
    pvals$ALDEx2 <- res$wi.ep
    pvalsAdj$ALDEx2 <- res$wi.eBH
    effect$ALDEx2_effect <- res$effect


  }

  xAll <- relData  %>% select(starts_with("1"))
  yAll <- relData  %>% select(starts_with("2"))

  if ("Student.t.test" %in% tests) {
    message("##### Student.t.test #####")
    r <- vector(length = n_taxa)
    eff <- vector(length = n_taxa)
    for (taxa in 1:n_taxa){
      x <- xAll[taxa,]
      y <- yAll[taxa,]
      t <- t.test(x = x, y = y)
      r[taxa] <- t$p.value
      eff[taxa] <- t$statistic


    }
    pvals$Student.t.test <- r
    pvalsAdj$Student.t.test <- p.adjust(r, method = "BH")
    effect$Student.t.test_effect <- eff
  }

  if ("Wilcoxon" %in% tests) {
    options(digits=8)
    message("##### Wilcoxon #####")
    r <- vector(length = n_taxa)
    eff <- vector(length = n_taxa)
    for (taxa in 1:n_taxa){
      x <- xAll[taxa,]
      y <- yAll[taxa,]
      x <- as.numeric(x)
      y <- as.numeric(y)
      t <- wilcox.test(x = x, y = y)
      r[taxa] <- t$p.value
      eff[taxa] <- t$statistic

    }
    pvals$Wilcoxon <- r
    pvalsAdj$Wilcoxon <- p.adjust(r, method = "BH")
    effect$Wilcoxon_effect <- eff
  }

  if ("ANCOM.BC" %in% tests) {
    # https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html#5_ALDEx2_outputs
    message("##### ANCOM.BC #####")
    resancombc <- ancombc(phyloseq = phyloData, formula = "treatment",
                          p_adj_method = "none", zero_cut = 0.90, tol = 1e-5,
                          max_iter = 100, conserve = TRUE, alpha = 0.05) # lib_cut = 1000, struc_zero = TRUE, neg_lb = TRUE, global = TRUE,group = "treatment"

    # res_global <- out$res_global # Only for more than one group
    resancombc <- resancombc$res

    pv <- merge(countData[1], resancombc$p_val, by = 0, all = T)
    pv <- merge(x = pvals, y = pv, by.x = "taxon", by.y = "Row.names", all = T)$treatment1
    ev <- merge(countData[1], resancombc$beta, by = 0, all = T)
    ev <- merge(x = effect, y = ev, by.x = "taxon", by.y = "Row.names", all = T)$treatment1
    pvals$ANCOM.BC <- pv
    pvalsAdj$ANCOM.BC <- p.adjust(pv, method = "BH")
    effect$ANCOM.BC_effect <- ev

  }

  # names(lst) <- c("pvalues", "effect-size")
  if("DESEq2" %in% tests){
    pvals <- pvals %>% mutate(DESEq2 = replace_na(DESEq2, 1)) # Some are removed due to the filtering step
  }
  result <- merge(pvals, pvalsAdj, by = c("taxon", "run", "sample_size", "n_taxa", "n_random", "change", "diff"), suffixes = c("_pvals","_pvalsAdj"))
  result <- merge(result, effect, by = c("taxon", "run", "sample_size", "n_taxa", "n_random", "change", "diff"))
  return(result)
}



####################################################################################
######################### FUNCTIONS TO RUN THE ACTUAL TEST #########################
####################################################################################

generateData <- function(variables, sample_size, n_taxa, run, change, n_random, seed){
  diff_data_list <- data.frame()
  for(m in n_taxa) {
    for(n in n_random) {
      for(j in sample_size){
        for(i in change){
          d <- createDfLogNormalAndNegBinom3(variables = v, sample_size = j, n_taxa = m, run, change = i, n_random = n, seed = seed)
          d <- cleanData(d)
          if(nrow(d) == 0){

            return(diff_data_list)
          }
          random_taxon <- sample(nrow(d), size = n, replace=F) # Randomly pick one taxon for differential abundance
          diff <- diffAbundance(df = d, change = i, random_taxon, variables = v)
          diff_data_list <- append(diff_data_list, list(diff))
        }
      }
    }
  }
  return(diff_data_list)
}

filterData <- function(data) {
  for(list in data){

    mintest <- min(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")]), na.rm = T) # Check min sample size
    maxtest <- max(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")]), na.rm = T) # Check max sample size
    meantest <- mean(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")]), na.rm = T)
    # if(mintest < 15000 | maxtest > 50000) {
    #   return(FALSE)
    # }
  }
  if(meantest < 20000 | meantest > 40000){
    return(FALSE)
  }
  print(mean(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")])), na.rm = T)
  return(TRUE)
}

getData <- function(variables, sample_size, n_taxa, run, change, n_random, seed){
  boolean <- FALSE
  s <- seed
  while(boolean == FALSE){

    # Generate all data
    data <- generateData(variables = variables, sample_size = sample_size, n_taxa = n_taxa, run = run, change = change, n_random = n_random, seed = s)
    if(length(data) > 0){
      # Test all data (each seperate data frame)
      boolean <- filterData(data)}
    else {boolean <- F}

    # Change seed to get another data frame
    s <- s + 100
  }
  return(data)
}

runAll <- function(variables = variables, sample_size = sample_size, n_taxa = n_taxa, runs = runs, change = change, n_random = n_random, seed = seed){
  diff_data_list <- data.frame()
  results <- data.frame()

  for(run in runs){
    s <- seed+run
    # set.seed(s)
    d <- getData(variables = variables, sample_size = sample_size, n_taxa = n_taxa, run = run, change = change, n_random = n_random, seed = s)
    for(df in d){
      n_samples <- sum(startsWith(colnames(df), "1"))
      metaData <- data.frame("treatment" = c(rep("1", times = n_samples),rep("0", times = n_samples)), row.names = colnames(df[,startsWith(colnames(df), "1") | startsWith(colnames(df), "2")]))
      res <- runTests(tests = tests, countData = df, meta = metaData, variables = v, seed = s)
      diff_data_list <- bind_rows(diff_data_list, df)
      results <- rbind(results, res)
    }

  }

  all_data <- full_join(diff_data_list, results)
  return(all_data)
}





####################################################################################################
########################### Main call: Simulate data and run test suite ############################
####################################################################################################

finalResults <- runAll(variables = variables, sample_size = sample_size, n_taxa = n_taxa, runs = runs, change = change, n_random = n_random, seed = s)


####################################################################################################
####################################### Write and save file ########################################
####################################################################################################


file_name <- paste(distribution,
                   "#n_taxa(", paste(n_taxa, collapse = ","), ")",
                   "#n_random(", paste(n_random, collapse = ","), ")",
                   "#n_runs(", paste(runs, collapse = ","), ")",
                   "#sample_size(", paste(sample_size, collapse = ","), ")",
                   "#change(", paste(change, collapse = ","), ")",
                   "#tests(", paste(tests, collapse = ","), ")",
                   "#seed(", s, ")", sep = "")

write_csv(finalResults, paste(file_name, ".csv", sep = ""))
