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
change <- c(2,3) #c(1,2,5,8,15)
sample_size <-3# c(3,5,10,15)
runs <- 2 # 5
n_taxa <- 10#100
n_random <- 2# c(2,5)#c(2,5)
v <- 2

# Choose tests
tests <- c("DESEq2", "ALDEx2", "Student.t.test", "Wilcoxon", "ANCOM.BC")

# Declare distribution used for simulating data (Used in name generation)
distribution <- "LogNormalAndNegBinom1"

# Set seed
seed <- 128
set.seed(seed)

####################################################################################################
######################################## Helper functions ##########################################
####################################################################################################

createDfLogNormalAndNegBinom1 <- function(variables, sample_size, n_taxa, run, change, n_random, random_taxon, seed) { # DO NOT CHANGE! MAKE NEW ONE IF NEEDED
  set.seed(seed)
  data <- data.frame(matrix(nrow = n_taxa , ncol = variables*sample_size))
  colnames(data) <- c(sort(paste(as.vector(outer(seq(variables), seq(sample_size), paste, sep = ".")), sep = "")))
  marginalProb <- rlnorm(n_taxa, meanlog = 4, sdlog = 2.5)

  data$taxon <- seq(n_taxa)
  data$sample_size <- sample_size
  data$n_taxa <- n_taxa
  data$n_random <- n_random
  data$run <- run
  data$change <- change
  data$diff <- data$taxon %in% random_taxon

  data <- data %>% select(taxon, run, sample_size, n_taxa, n_random, change, diff, everything())

  for (j in 8:ncol(data)) {
    # Generate count data for each taxon (SAMPLES)
    data[j] <- as.integer(rnbinom(n_taxa, mu = marginalProb, size = 100)) # Size is the shape parameter af the gamma distrib. The lower size the more variance is introduced

  }

  return(data)
}

# Makes n_taxa number of taxa have different abundance
diffAbundance <- function(n_taxa, df, change, random_taxon, n_dataframes, variables) {
  df[random_taxon, startsWith(colnames(df), "2")] <- df[random_taxon, startsWith(colnames(df), "2")]*change
  return(df)

}

# Transform counts to relative abundances
relativeAbundanceParalel <- function(data) {
  for(i in colnames(data)[startsWith(colnames(data), "2") | startsWith(colnames(data), "1")]){
    data[[i]] <- data[[i]]/sum(data[[i]])
  }
  return(data)
}




####################################################################################################
###################################### Test suite function #########################################
####################################################################################################

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

    resALDEX2 <- aldex(countData, m, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE) # Object contain reaslts from many different tests. See link for
    if (!dim(resALDEX2)[1] == n_taxa){
      n <- row.names(resALDEX2)
      m <- pvals$taxon
      diff <- setdiff(m,n)
      for (j in diff){
        resALDEX2 <- rbind(resALDEX2[1:j-1,], rep(NA, times = dim(resALDEX2)[2]), resALDEX2[-(1:j-1),])
      }


    }


    pvals$ALDEx2 <- resALDEX2$wi.ep
    pvalsAdj$ALDEx2 <- resALDEX2$wi.eBH
    effect$ALDEx2_effect <- resALDEX2$effect

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
      if (sd(x) != 0 & sd(x) != 0) {
        t <- t.test(x = x, y = y)
        r[taxa] <- t$p.value
        eff[taxa] <- t$statistic

      } else {
        r[taxa] <- NA
        eff[taxa] <- NA
      }

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

    if (!dim(resancombc$p_val)[1] == n_taxa){
      n <- as.numeric(substr(row.names(resancombc$p_val), 3,8))
      m <- pvals$taxon
      diff <- setdiff(m,n)
      pv <- resancombc$p_val$treatment1
      for (i in diff){
        pv <- c(pv[1:i-1], NA, pv[-(1:i-1)])
      }
      pvals$ANCOM.BC <- pv
      pvalsAdj$ANCOM.BC <- p.adjust(pv, method = "BH")
    }
    if (!dim(resancombc$beta)[1] == n_taxa){
      n <- as.numeric(substr(row.names(resancombc$beta), 3,8))
      m <- pvals$taxon
      diff <- setdiff(m,n)
      ev <- resancombc$beta$treatment1
      for (i in diff){
        ev <- c(ev[1:i-1], NA, ev[-(1:i-1)])

      }
      effect$ANCOM.BC_effect <- ev
    } else {
      pvals$ANCOM.BC <- resancombc$p_val$treatment1 # Will break if more treatments are added
      pvalsAdj$ANCOM.BC <- p.adjust(resancombc$p_val$treatment1, method = "BH")
      effect$ANCOM.BC_effect <- resancombc$beta$treatment1 # Will break

    }

  }

  # names(lst) <- c("pvalues", "effect-size")
  if("DESEq2" %in% tests){
    pvals <- pvals %>% mutate(DESEq2 = replace_na(DESEq2, 1)) # Some are removed due to the filterin step
  }
  result <- merge(pvals, pvalsAdj, by = c("taxon", "run", "sample_size", "n_taxa", "n_random", "change", "diff"), suffixes = c("_pvals","_pvalsAdj"))
  result <- merge(result, effect, by = c("taxon", "run", "sample_size", "n_taxa", "n_random", "change", "diff"))
  return(result)
}



####################################################################################################
######################### Main function: Simulate data and run test suite ##########################
####################################################################################################


diff_data_list <- data.frame()
results <- data.frame()

system.time({
  for(k in seq(runs)){
    s <- k # Is set here, so the abundance is the same for taxon # in entire run k (So we can compare the effect of change + the other variables)
    for(m in n_taxa) {
      for(n in n_random) {
        random_taxon <- sample(m, size = n, replace=F) # Randomly pick one taxon for differential abundance
        for(j in sample_size){
          for(i in change){
            d <- createDfLogNormalAndNegBinom1(variables = v, sample_size = j, n_taxa = m, run = k, change = i, n_random = n,random_taxon = random_taxon, seed = s)
            diff <- diffAbundance(n_taxa = m, df = d, change = i, random_taxon, variables = v)
            metaData <- data.frame("treatment" = c(rep("1", times = j),rep("0", times = j)), row.names = colnames(d[8:ncol(d)]))
            res <- runTests(tests = tests, countData = diff, meta = metaData, variables = v, seed = s)
            diff_data_list <- bind_rows(diff_data_list, diff)
            results <- rbind(results, res)
          }
        }
      }
    }
  }



  all_data <- merge(diff_data_list, results)
})


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
                   "#seed(", seed, ")", sep = "")

write_csv(all_data, paste("./data/compositionConstraintTest_data/", file_name, ".csv", sep = ""))
