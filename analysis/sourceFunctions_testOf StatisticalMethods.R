############################################################
############### FUNCTIONS TO SIMULATE DATA #################
############################################################

# Note: The relative abundance data function makes the central log ratio trans!!

createDfLogNormalAndNegBinom1 <- function(variables, sample_size, n_taxa, run, change, n_random, random_taxon, seed) { # DO NOT CHANGE! MAKE NEW ONE IF NEEDED
  # set.seed(seed)
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

createDfLogNormalAndNegBinom2 <- function(variables, sample_size, n_taxa, run, change, n_random, random_taxon, seed) { # DO NOT CHANGE! MAKE NEW ONE IF NEEDED
  # set.seed(seed)
  data <- data.frame(matrix(nrow = n_taxa , ncol = variables*sample_size))
  colnames(data) <- c(sort(paste(as.vector(outer(seq(variables), seq(sample_size), paste, sep = ".")), sep = "")))
  zeroesProb <- rnorm(1, mean = 0.97, sd = 0.01)
  trueZero <- rbinom(n = n_taxa, size = 1, prob = zeroesProb)
  # marginalProb <- rlnorm(n_taxa, meanlog = log(1), sdlog = log(2.718282*7))
  marginalProb <- ifelse(trueZero==1, 0, rlnorm(sum(trueZero==0), meanlog = log(2.5), sdlog = log(2.718282*5)))


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


createDfLogNormalAndNegBinom3 <- function(variables, sample_size, n_taxa, run, change, n_random, seed) { # DO NOT CHANGE! MAKE NEW ONE IF NEEDED
  print("createDf")
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
  uv <- which(startsWith(colnames(data), "1."))
  ev <- which(startsWith(colnames(data), "2."))
  uc <- 1
  ec <- 1

  # In order to make sure that the counts are consistant even when the sample size changes, it is important to first simulate 1.2 then 2.1 before 2.1 and 2.2 and so on. It is easier to just assign the data to the correct column.
  # The point is to get all the uneven numbers first and then all the even numbers.
  for (j in which(startsWith(colnames(data), "1.")):ncol(data)) {
    # Generate count data for each taxon (SAMPLES)
    print(j)
    if(j %% 2 == 0){ # Even numbers
      data[ev[ec]] <- as.integer(rnbinom(n_taxa, mu = marginalProb, size = 100)) # Size is the shape parameter af the gamma distrib. The lower size the more variance is introduced
      ec <- ec + 1
    } else { # Uneven numbers
      data[uv[uc]] <- as.integer(rnbinom(n_taxa, mu = marginalProb, size = 100)) # Size is the shape parameter af the gamma distrib. The lower size the more variance is introduced
      uc <- uc + 1
    }}

  return(data)
}

# Makes n_taxa number of taxa have different abundance
diffAbundance <- function(df, change, random_taxon, variables) {
  print("diffAbundance")
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
  rows <- rownames(data)
  cols <- colnames(data)
  data <- as.matrix(data)
  for(i in 1:ncol(data)){
    data[,i] <- data[,i]/sum(data[,i])
    data[,i] <- clr(data[,i])
  }

  data <- as.data.frame(data)
  colnames(data) <- cols
  rownames(data) <- rows

  return(data)
}


cleanData <- function(data, t){
  print("cleanData")
  d <- data[,startsWith(colnames(data), "2") | startsWith(colnames(data), "1")]
  boolean <- d > 0

  rsums <- rowSums(boolean, na.rm = T)
  rkeep <- rsums/ncol(d) >= t

  data <- data[rkeep,]


  # csums <- colSums(boolean, na.rm = T)
  # ckeep <- csums > 0
  # # print(ckeep)
  #
  # data <- data[,ckeep]
  print("return cleanData")
  return(data)
}




########################################################
################# FUNCTION TO RUN TESTS ################
########################################################

# Recreate run parameter somehow
runTests <- function(tests, countData, meta, variables, seed) {
  print("runTests")
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
        pvals <- pvals %>% mutate(DESEq2 = replace_na(DESEq2, 1)) # Some are removed due to the filterin step
      }
    result <- merge(pvals, pvalsAdj, by = c("taxon", "run", "sample_size", "n_taxa", "n_random", "change", "diff"), suffixes = c("_pvals","_pvalsAdj"))
    result <- merge(result, effect, by = c("taxon", "run", "sample_size", "n_taxa", "n_random", "change", "diff"))
  return(result)
}



####################################################################################
######################### FUNCTIONS TO RUN THE ACTUAL TEST #########################
####################################################################################

generateData <- function(variables, sample_size, n_taxa, run, change, n_random, seed){
  set.seed(seed)
  diff_data_list <- data.frame()
  for(m in n_taxa) {
    for(n in n_random) {
      random_taxon <- FALSE
      for(j in sample_size){
        for(i in change){
          d <- createDfLogNormalAndNegBinom3(variables = v, sample_size = j, n_taxa = m, run, change = i, n_random = n, seed = seed)
          d <- cleanData(d, 0.05)
          if(nrow(d) == 0){

            return(diff_data_list)
          }
          if(length(random_taxon) == 0){
            random_taxon <- sample(nrow(d), size = n, replace=F)
          }else if(random_taxon == FALSE){
            random_taxon <- sample(nrow(d), size = n, replace=F) # Randomly pick one taxon for differential abundance
          } else {
            random_taxon <- which(d$taxon %in% random_taxon_names)
          }
          diff <- diffAbundance(df = d, change = i, random_taxon, variables = v)
          random_taxon_names <- diff$taxon[diff$diff]
          diff_data_list <- append(diff_data_list, list(diff))
        }
      }
    }
  }
  return(diff_data_list)
}

filterData <- function(data) {
  print("filter data")
  for(list in data){

  mintest <- min(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")]), na.rm = T) # Check min sample size
  maxtest <- max(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")]), na.rm = T) # Check max sample size
  meantest <- mean(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")]), na.rm = T)
  print(mintest)
  print(mean(colSums(list[,startsWith(colnames(list), "1") | startsWith(colnames(list), "2")])), na.rm = T)
  print(maxtest)
  # if(mintest < 15000 | maxtest > 50000) {
  #   return(FALSE)
  # }
  }
  # if(meantest < 20000 | meantest > 40000){
  #   return(FALSE)
  # }
  print("return filterData")
  return(TRUE)
  }

getData <- function(variables, sample_size, n_taxa, run, change, n_random, seed){
  print("getData")
  boolean <- FALSE
  s <- seed
  while(boolean == FALSE){

    # Generate all data
    data <- generateData(variables = variables, sample_size = sample_size, n_taxa = n_taxa, run = run, change = change, n_random = n_random, seed = s)
  print("data < 0")
    if(length(data) > 0){
    # Test all data (each seperate data frame)
    boolean <- filterData(data)}
    else {boolean <- F}

    # Change seed to get another data frame
    s <- s + 100
  }
  print("return getData")
  return(data)
}

runAll <- function(variables = variables, sample_size = sample_size, n_taxa = n_taxa, runs = runs, change = change, n_random = n_random, seed = seed){
  print("runAll")
  diff_data_list <- data.frame()
  results <- data.frame()

  for(run in seq(runs)){
    s <- seed+run
    # set.seed(s)
    d <- getData(variables = variables, sample_size = sample_size, n_taxa = n_taxa, run = run, change = change, n_random = n_random, seed = s)
    for(df in d){
      n_samples <- sum(startsWith(colnames(df), "1"))
      print(n_samples)
      metaData <- data.frame("treatment" = c(rep("1", times = n_samples),rep("0", times = n_samples)), row.names = colnames(df[,startsWith(colnames(df), "1") | startsWith(colnames(df), "2")]))
      res <- runTests(tests = tests, countData = df, meta = metaData, variables = v, seed = s)
      diff_data_list <- bind_rows(diff_data_list, df)
      results <- rbind(results, res)
    }

  }
  # Reorder columns (They are in correct order for each single data frame when testing. meta data and the columns are matching!)
  diff_data_list <- diff_data_list %>% select(taxon, run, sample_size, n_taxa, n_random, change, diff, starts_with(c("1.")), starts_with(c("2.")))
  all_data <- full_join(diff_data_list, results)
  return(all_data)
}

########################################################################################
####################### FUNCTIONS TO RUN ALL WITH EXTRA FILTER STEP ####################
########################################################################################
extraClean <- function(data, t) {
  for(i in 1:length(data)){
    data[[i]] <- cleanData(data[[i]], t)
  }
  return(data)
}



runAllFilterChange <- function(variables = variables, sample_size = sample_size, n_taxa = n_taxa, runs = runs, change = change, n_random = n_random, seed = seed, preThresh = preThresh){
  print("runAllFilterChange")
  diff_data_list <- data.frame()
  results <- data.frame()

  for(run in seq(runs)){
    s <- seed+run
    # set.seed(s)
    d <- getData(variables = variables, sample_size = sample_size, n_taxa = n_taxa, run = run, change = change, n_random = n_random, seed = s)
    d <- extraClean(d, preThresh)
    for(df in d){
      n_samples <- sum(startsWith(colnames(df), "1"))
      print(n_samples)
      metaData <- data.frame("treatment" = c(rep("1", times = n_samples),rep("0", times = n_samples)), row.names = colnames(df[,startsWith(colnames(df), "1") | startsWith(colnames(df), "2")]))
      res <- runTests(tests = tests, countData = df, meta = metaData, variables = v, seed = s)
      diff_data_list <- bind_rows(diff_data_list, df)
      results <- rbind(results, res)
    }

  }
  # Reorder columns (They are in correct order for each single data frame when testing. meta data and the columns are matching!)
  diff_data_list <- diff_data_list %>% select(taxon, run, sample_size, n_taxa, n_random, change, diff, starts_with(c("1.")), starts_with(c("2.")))
  all_data <- full_join(diff_data_list, results)
  return(all_data)
}


####################################################################################
################ FUNCTIONS TO MODIFY TEST RESULTS AND MAKE FIGURES #################
####################################################################################


calculateConfusion <- function(data){
  p <- data %>% select(-ends_with(c("effect", "pvals"))) %>% pivot_longer(cols = !c("taxon", "run", "change", "sample_size", "n_taxa", "n_random", "diff", starts_with(c("1", "2"))), names_to = "method", values_to = "pvalueAdj")
  p <- p %>% mutate(pvalueAdj = replace_na(pvalueAdj, 1))

  p <- p %>% mutate(class = ifelse(pvalueAdj < 0.05 & diff == TRUE & !change == 1, "TP",
                            ifelse(pvalueAdj < 0.05 & diff == TRUE & change == 1, "FP",
                            ifelse(pvalueAdj < 0.05 & !diff == TRUE, "FP",
                            ifelse(pvalueAdj > 0.05 & diff == TRUE & !change == 1,"FN",
                            ifelse(pvalueAdj > 0.05 & diff == TRUE & change == 1,"TN",
                            ifelse(pvalueAdj > 0.05 & !diff == TRUE, "TN","N"
      )))))))
  p <- p %>% mutate(class = replace_na(class, "NA"))
  return(as.data.frame(p))

}

calculateConfusionAndEffect <- function(data){
  p <- data %>% select(-ends_with(c("pvals"))) %>% pivot_longer(cols = !c("taxon", "run", "change", "sample_size", "n_taxa", "n_random", "diff", starts_with(c("1", "2"))), names_to = "method", values_to = "pvalueAdj")
  pEffect <- p %>% filter(grepl("effect",method)) %>% rename(pvalueAdj = "effect")
  pEffect$method <- unlist(strsplit(pEffect$method, "_"))[!unlist(strsplit(pEffect$method, "_")) == "effect"] # Remove _effect suffix

  pPvals <- p %>% filter(grepl("pvalsAdj",method))
  pPvals$method <- unlist(strsplit(pPvals$method, "_"))[!unlist(strsplit(pPvals$method, "_")) == "pvalsAdj"] # Remove _pvalsAdj suffix
  p <- inner_join(pEffect, pPvals)


  p <- p %>% mutate(pvalueAdj = replace_na(pvalueAdj, 1))

  p <- p %>% mutate(class = ifelse(pvalueAdj < 0.05 & diff == TRUE & !change == 1, "TP",
                                   ifelse(pvalueAdj < 0.05 & diff == TRUE & change == 1, "FP",
                                          ifelse(pvalueAdj < 0.05 & !diff == TRUE, "FP",
                                                 ifelse(pvalueAdj > 0.05 & diff == TRUE & !change == 1,"FN",
                                                        ifelse(pvalueAdj > 0.05 & diff == TRUE & change == 1,"TN",
                                                               ifelse(pvalueAdj > 0.05 & !diff == TRUE, "TN","N"
                                                               )))))))
  p <- p %>% mutate(class = replace_na(class, "NA"))
  return(as.data.frame(p))

}

calculateFDRandSensitivity <- function(data) {
  stats <- data %>% group_by(method, run, change, sample_size,n_taxa, n_random) %>% summarise(FDR = sum(class == "FP", na.rm=TRUE)/(sum(class == "FP", na.rm=TRUE)+sum(class == "TP", na.rm=TRUE))*100,
                                                Sensitivity = sum(class == "TP", na.rm=TRUE)/(sum(class == "TP", na.rm=TRUE)+sum(class == "FN", na.rm=TRUE))*100,
                                                Specificity = sum(class == "TN", na.rm=TRUE)/(sum(class == "TN", na.rm=TRUE)+sum(class == "FP", na.rm=TRUE))*100)

  stats <- stats %>% mutate(FDR = replace(FDR, FDR == "NaN", 0)) # If no FP or TP are found then no false positives are found

  stats <- melt(stats, id.vars = c("method", "run", "sample_size", "n_taxa", "n_random", "change"), variable.name = "class")
  return(stats)
}

figureChange <- function(data, s, nt, nr) {
  data$run <- as.factor(data$run)
  # data$change <- as.factor(data$change)

  data <- data %>% filter(sample_size == s, n_taxa == nt, n_random == nr) %>% group_by(change, method, class) %>% mutate(mean_val = mean(value), sd_val = sd(value))



  p <- data %>% ggplot(aes(x = change, y = value)) +
    geom_point(aes(color = run), position = position_dodge(width = 0.3))+
    geom_line(aes(color = run), position = position_dodge(width = 0.3))+
    # geom_point(aes(x = change, y = mean_val ), position = position_dodge(width = 0.3))+
    # geom_line(aes(x = change, y = mean_val, color = "Mean"), position = position_dodge(width = 0.3))+
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "point", aes(color = "Mean"))+
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "line", aes(color = "Mean"))+
    geom_errorbar(aes(y = mean_val, x = change,
                      ymin=mean_val-sd_val,
                      ymax=mean_val+sd_val),
                      color = "black", width = 1)+
    facet_grid(cols = vars(method), rows = vars(class), scales = "free")+
    ylim(-20,120)+
    scale_y_continuous(breaks = c(0,20,40,60,80,100))+
    scale_x_continuous(breaks = change, labels = change)+
    labs(color='Run')+
    ylab("Percent %") +
    xlab("Fold change") +
    scale_colour_manual(values= append(brewer.pal(length(unique(data$run)), "Set1"), "#000000"))+
    theme_bw()

  return(p)

}

figureSample <- function(data, c, nt, nr) {
  data$run <- as.factor(data$run)
  # data$change <- as.factor(data$change)

  data <- data %>% filter(change == c, n_taxa == nt, n_random == nr) %>% group_by(sample_size, method, class) %>% mutate(mean_val = mean(value), sd_val = sd(value))



  p <- data %>% ggplot(aes(x = sample_size, y = value)) +
    geom_point(aes(color = run), position = position_dodge(width = 0.3))+
    geom_line(aes(color = run), position = position_dodge(width = 0.3))+
    # geom_point(aes(x = sample_size, y = mean_val ), position = position_dodge(width = 0.3))+
    # geom_line(aes(x = sample_size, y = mean_val, color = "Mean"), position = position_dodge(width = 0.3))+
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "point", aes(color = "Mean"))+
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "line", aes(color = "Mean"))+
    geom_errorbar(aes(y = mean_val, x = sample_size,
                      ymin=mean_val-sd_val,
                      ymax=mean_val+sd_val),
                  color = "black", width = 1)+
    facet_grid(cols = vars(method), rows = vars(class), scales = "free")+
    ylim(-20,120)+
    scale_y_continuous(breaks = c(0,20,40,60,80,100))+
    scale_x_continuous(breaks = sample_size, labels = sample_size)+
    labs(color='Run')+
    ylab("Percent %") +
    xlab("Sample size") +
    scale_colour_manual(values= append(brewer.pal(length(unique(data$run)), "Set1"), "#000000"))+
    theme_bw()

  return(p)

}

figureDiff <- function(data, c, s, nt) {
  data$run <- as.factor(data$run)
  # data$change <- as.factor(data$change)

  data <- data %>% filter(change == c, sample_size == s, n_taxa == nt) %>% group_by(n_random, method, class) %>% mutate(mean_val = mean(value), sd_val = sd(value))



  p <- data %>% ggplot(aes(x = n_random, y = value)) +
    geom_point(aes(color = run), position = position_dodge(width = 0.3))+
    geom_line(aes(color = run), position = position_dodge(width = 0.3))+
    # geom_point(aes(x = n_random, y = mean_val ), position = position_dodge(width = 0.3))+
    # geom_line(aes(x = n_random, y = mean_val, color = "Mean"), position = position_dodge(width = 0.3))+
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "point", aes(color = "Mean"))+
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "line", aes(color = "Mean"))+
    geom_errorbar(aes(y = mean_val, x = n_random,
                      ymin=mean_val-sd_val,
                      ymax=mean_val+sd_val),
                  color = "black", width = 1)+
    facet_grid(cols = vars(method), rows = vars(class), scales = "free")+
    ylim(-20,120)+
    scale_y_continuous(breaks = c(0,20,40,60,80,100))+
    scale_x_continuous(breaks = n_random, labels = n_random)+
    labs(color='Run')+
    ylab("Percent %") +
    xlab("Number of random taxon") +
    scale_colour_manual(values= append(brewer.pal(length(unique(data$run)), "Set1"), "#000000"))+
    theme_bw()

  return(p)

}

makeVenn2 <- function(conf_data,tests, t, r, c, j, m, n) {

  all_counts <- data.frame() # Mean over all runs

  for(n_run in 1:r){

  data <- conf_data %>% mutate(dummy = row_number()) %>% filter(class == t, run == n_run, change == c, sample_size == j, n_taxa == m, n_random == n) %>% pivot_wider(names_from = method, values_from = taxon) %>% select(-c(run, change, sample_size, n_taxa, n_random, diff, dummy, starts_with(c("1", "2"))))
  l <- list()

  if(length(tests[!tests %in% colnames(data)]) > 0){
    for(i in tests[!tests %in% colnames(data)]){
      df <- data.frame("X" = rep(NA, nrow(data)))
      colnames(df) <- i
      data <- cbind(data, df)

    }

  }

  for(i in 3:ncol(data)){
    l <- append(l,na.omit(data[i]))

  }
  if(length(l) == 1) {
    return("Only one test has this type of data. Venn diagram is not possible to make")
  }

  names(l) <- colnames(data)[3:ncol(data)]
  l <- l[order(names(l))]

  sum <- 0
  for(list in l){
  sum <- sum + length(list)
  }

  if(sum == 0){
    cs <- rep(0,32) #32 = number of elements in Venn diagram with 5 groups
  } else{

    p <- venn(l)
    cs <- p$counts
  }
    all_counts <- rbind(all_counts, cs)
  }

  mean_counts <- colSums(all_counts)/length(r) # round
  v_plot <-   venn(length(tests), snames = names(l)[order(names(l))], counts = mean_counts, zcolor = scales::viridis_pal(option = "D")(length(names(l))), opacity = 0.4, ellipse = T,
                   borders = TRUE, box = TRUE, par = TRUE, ggplot = TRUE) #c(ilcs = 10, sncs = 10)
# zcolor = brewer.pal(length(names(l)), "Dark2")
#zcolor = scales::viridis_pal(option = "D")(length(names(l))

  scales::viridis_pal()(4)
return(v_plot)
}

makeVenn3 <- function(conf_data, m, c, r) {

  all_counts <- data.frame() # Mean over all runs

  for(n_run in 1:r){

    data <- conf_data %>% mutate(dummy = row_number()) %>% filter(method == m, class == c, run == n_run) %>%  pivot_wider(names_from = filter, values_from = taxon) %>% select(-c(run, change, sample_size, n_taxa, n_random, diff, dummy, starts_with(c("1", "2"))))
    l <- list()


    if(ncol(data) <= 6){
      for(i in unique(conf_data$filter)[!unique(conf_data$filter) %in% colnames(data[startsWith(colnames(data),"0")])]){
        df <- data.frame("X" = rep(NA, nrow(data)))
        colnames(df) <- i
        data <- cbind(data, df)

      }}

    # Make df with results
    for(i in 4:ncol(data)){
      l <- append(l,na.omit(data[i]))

    }

    # make results numeric
    for(i in 1:length(l)){
      l[[i]] <- as.numeric(l[[i]])
    }

    # Rename results
    names(l) <- paste("Filter: ", colnames(data)[4:ncol(data)])
    l <- l[order(names(l))]

    sum <- 0
    for(list in l){
      sum <- sum + length(list)
    }

    if(sum == 0){
      cs <- rep(0,7) #7 = number of elements in Venn diagram with 3 groups
    } else{

      p <- venn(l)
      cs <- p$counts
    }
    all_counts <- rbind(all_counts, cs)
    }


  mean_counts <- colSums(all_counts)/r # round
  v_plot <-   venn(length(l), snames = names(l), counts = mean_counts, zcolor = scales::viridis_pal(option = "D")(length(names(l))), opacity = 0.4,
                   borders = TRUE, box = TRUE, par = TRUE, ggplot = TRUE) #c(ilcs = 10, sncs = 10)
  # zcolor = brewer.pal(length(names(l)), "Dark2")
  #zcolor = scales::viridis_pal(option = "D")(length(names(l))

  scales::viridis_pal()(4)
  return(v_plot)
}


getNA <- function(conf_data,tests, r, c, j, m, n){
  all_NA <- table(conf_data %>% filter(class == "NA", change == c, sample_size == j, n_taxa == m, n_random == n) %>% select(method))
  NAtable <- all_NA/length(r)
  return(NAtable)
}

performancePlot <- function(data, nr, c){
  data <- data %>% group_by(method, sample_size, n_taxa, n_random, change, class) %>% mutate(mean = mean(value),
                                                                                             std.dev = sd(value))
  data <- data %>% rename("class" = "Measure", "change" = "Fold change")
  data <- data %>% filter(!Measure == "Specificity")

  plot <- data %>% filter(n_random == nr, `Fold change` %in% c) %>%
    ggplot(aes(x = sample_size, y = std)) +
    geom_line(aes(color = method)) +
    geom_point(aes(shape = method, color = method), size = 3) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              colour = "gray", fill = NA) +
    facet_grid(cols = vars(Measure), rows = vars(`Fold change`),labeller = label_both)+
    labs(color='Method', shape = "Method")+
    xlab("Sample size")+
    ylab("Percentage %")+
    geom_hline(aes(yintercept = line), linetype = 2, color = "grey20")+
    # scale_color_viridis_d(option = "D")+
    scale_color_manual(values = c(viridis(4, alpha = 1, begin = 0, end = 0.80, direction = 1, option = "D"), "#f1c232")) +
    # scale_color_manual(values = c("#E95694", "#7E4EAC", "#3C93C2", "#40AD5A", "#FD8D3C"))+
    theme_minimal()+
    theme(strip.background=element_rect(colour="gray",
                                        fill="gray90"))
  return(plot)

}

performancePlotDESEq <- function(data, nr, c, fc){
  data <- data %>% group_by(method, sample_size, n_taxa, n_random, change, class, filter) %>% mutate(mean = mean(value),
                                                                                             std.dev = sd(value))

  data <- data %>% filter(!class == "Specificity", method == "DESEq2", change == fc)
  data <- data %>% rename("class" = "Measure", "change" = "Fold change")

  plot <- data %>% filter(n_random == nr, `Fold change` %in% c) %>%
    ggplot(aes(x = sample_size, y = mean)) +
    geom_line(aes(color = filter)) +
    geom_point(aes(shape = filter, color = filter), size = 3) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              colour = "gray", fill = NA) +
    facet_grid(cols = vars(Measure), rows = vars(`Fold change`),labeller = label_both)+
    labs(color='Filter', shape = "Filter")+
    xlab("Sample size")+
    ylab("Percentage %")+
    geom_hline(aes(yintercept = line), linetype = 2, color = "grey20")+
    # scale_color_viridis_d(option = "D")+
    scale_color_manual(values = c(viridis(6, alpha = 1, begin = 0, end = 0.80, direction = 1, option = "B"))) +
    # scale_color_manual(values = c("#E95694", "#7E4EAC", "#3C93C2", "#40AD5A", "#FD8D3C"))+
    theme_minimal()+
    theme(strip.background=element_rect(colour="gray",
                                        fill="gray90"))
  return(plot)

}
####################################################################################
######### FUNCTION TO LOAD LARGE DATA FRAME AND CALCULATE FDR/SENSITIVITY ##########
####################################################################################

largeDfLoad <- function(filename, runs, chunck_size){
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

largeDfLoadConf <- function(filename, runs, chunck_size, size){
  FDRandSensitivity <- data.frame()



    f <- function(data, pos) {
      data <- data %>% filter(change == 10, sample_size == size, n_random == 150)
      if(nrow(data) > 0){
        conf <- calculateConfusion(data)
        conf$method <- unlist(strsplit(conf$method, "_"))[!unlist(strsplit(conf$method, "_")) == "pvalsAdj"] # Remove .pvals suffix
      } else{
        conf <- NULL
      }
      return(conf)
    }



    data <- read_delim_chunked(file = filename,
                               callback = DataFrameCallback$new(f),
                               chunk_size = chunk_size, delim = c(","), col_names = T, na = c("NA", NA)
    ) # na = c("NA", "NA")

  return(data)
}

largeDfLoadConfAndEffect <- function(filename, runs, chunck_size, size){
  FDRandSensitivity <- data.frame()



  f <- function(data, pos) {
    data <- data %>% filter(change == 10, sample_size == size, n_random == 150)
    if(nrow(data) > 0){
      conf <- calculateConfusionAndEffect(data)
    } else{
      conf <- NULL
    }
    return(conf)
  }



  data <- read_delim_chunked(file = filename,
                             callback = DataFrameCallback$new(f),
                             chunk_size = chunk_size, delim = c(","), col_names = T, na = c("NA", NA)
  ) # na = c("NA", "NA")

  return(data)
}



###########################################################################
###########################################################################

###########################################################################
###########################################################################
################################ OBSOLETE #################################
###########################################################################
###########################################################################

###########################################################################
###########################################################################
#
#
#
# # Creates one data frame
# createDfPoisAndNegBinom <- function(variables, sample_size, sizes, n_taxa, seed) {
#   set.seed(seed)
#   data <- data.frame(matrix(nrow = n_taxa , ncol = variables*sample_size), row.names = paste("Taxon", seq(n_taxa), sep = ""))
#   colnames(data) <- c(sort(paste("Sample", as.vector(outer(seq(variables), seq(sample_size), paste, sep = ".")), sep = "")))
#   k = 1 # Counter to help fill table
#   # Get marginal probabilities (TRUE BACKGROUND)
#   marginalProb <- as.integer(rpois(n_taxa, lambda = sizes), seed = seed)
#
#   for (j in 1:(sample_size*variables)) {
#     # Generate count data for each taxon (SAMPLES) (seed in the loop)
#     data[k] <- as.integer(rnbinom(n_taxa, size = marginalProb, prob = 0.5), seed = seed+k)
#     k = k+1
#   }
#   return(data)
# }
#
# createDfLogNormalAndMultinorm <- function(variables, sample_size, sizes, n_taxa, seed) {
#   set.seed(seed)
#   data <- data.frame(matrix(nrow = n_taxa , ncol = variables*sample_size), row.names = paste("Taxon", seq(n_taxa), sep = ""))
#   colnames(data) <- c(sort(paste("Sample", as.vector(outer(seq(variables), seq(sample_size), paste, sep = ".")), sep = "")))
#   k = 1 # Counter to help fill table
#   marginalProb <- rlnorm(n_taxa, meanlog = 4, sdlog = 2.5)
#
#   for (j in 1:(sample_size*variables)) {
#     # Generate count data for each taxon (SAMPLES) (seed in the loop)
#     data[k] <- as.integer(rnbinom(n_taxa, mu = marginalProb, size = 100), seed = seed+k) # Size is the shape parameter af the gamma distrib. The lower size the more variance is introduced
#     k = k+1
#   }
#   return(data)
# }
#
#
# combineData <- function(data, confusion, n_taxa) {
#
#   t <- confusion %>% pivot_wider(names_from =  "method", values_from = c("pvalue", "class"))
#
#   for (i in 1:length(data)) {
#     data[[i]]$taxon = seq(n_taxa)
#   }
#   f <- bind_rows(data)
#
#   g <- cbind(t,f)
#   g <- g[1:(ncol(g)-2)]
#
#   return(g)
# }
#
# # Uses createDf to create a list of data frames
# createDfList <- function(n_dataframes,variables, sample_size, sizes, n_taxa, change, seed, n_random, run, mean_distribution){
#   df_list <- list()
#   if(mean_distribution == "poisson") {
#     for (i in 1:(n_dataframes*length(change))) {
#       data <- createDfPoisAndNegBinom(variables, sample_size, sizes, n_taxa, seed+i)
#       df_list[[i]] <- data
#     }
#   } else if (mean_distribution == "log_normal") {
#     for (i in 1:(n_dataframes*length(change))) {
#       data <- createDfLogNormalAndNegBinom(variables, sample_size, sizes, n_taxa,n_random, run, seed+i)
#       df_list[[i]] <- data
#     }
#
#   }
#   return(df_list)
# }
#
# constraint_summary <- function(relBackgroundData, relDiffAbundanceData){
#   backAb <- data.frame()
#   diffAb <- data.frame()
#   for (i in 1:length(relBackgroundData)){
#     n <- relBackgroundData[[i]] %>% rownames_to_column(var = "Taxon")
#     m <- relDiffAbundanceData[[i]] %>%  rownames_to_column(var = "Taxon")
#     backAb <- rbind(backAb, n)
#     diffAb <- rbind(diffAb, m)
#   }
#   diffAb <- diffAb %>% filter(!Taxon %in% paste("Taxon", random_taxon, sep = "")) # remove random taxon which are diff abundant
#   backAb <- backAb %>% filter(!Taxon %in% paste("Taxon", random_taxon, sep = "")) # remove random taxon which are diff abundant
#   change <- diffAb %>% select(change)
#   diffAb <- diffAb %>% select(!c(Taxon, change))
#   backAb <- backAb %>% select(!c(Taxon))
#
#
#   taxa_means <- rowMeans(backAb)
#
#   difference <- diffAb-backAb
#   difference <- difference/taxa_means*100
#
#   difference <- difference %>% select(!c(Sample1.1, Sample1.2, Sample1.3))
#   difference <- difference %>% mutate(mean = rowMeans(.))
#   difference <- cbind(difference, change)
#   difference <- na.omit(difference)
#   hist(difference$mean, breaks = 30)
#   table <- difference %>% group_by(change) %>% summarise(medians = median(mean))
#   table %>%
#     kable(caption="Median change in relative abundance") %>%
#     kable_styling()
#   print(table)
#   return(table)
# }
#
# makeVenn <- function(conf_data, t, r, c) {
#
#   data <- conf_data %>%  filter(class == t, run == r, change == c) %>% select(-c(run, change)) %>% pivot_wider(names_from = method, values_from = taxon)
#   m <- list()
#
#   for(i in 3:ncol(data)){
#     m <- append(m,na.omit(data[i]))
#
#   }
#   if(length(m) == 1) {
#     return("Only one test has this type of data. Venn diagram is not possible to make")
#   }
#
#   names(m) <- colnames(data)[3:ncol(data)]
#
#   ggvenn(m,
#          fill_color = brewer.pal(length(names(m)), "PiYG"),
#          stroke_size = 0.5, set_name_size = 3,
#          text_size = 4)
#
# }
#
#
#


largeDftest <- function(filename, runs, chunck_size){


  for(i in seq(runs)){

    f <- function(data, pos) {
      data <- data %>% filter(run == i)
      if(nrow(data) > 0){
        size <- length(unique(data$taxon))
      } else{
        size <- NULL
      }
      # FDRandSensitivity <- calculateFDRandSensitivity(conf)
      return(size)
    }



    data <- read_delim_chunked(file = filename,
                               callback = DataFrameCallback$new(f),
                               chunk_size = chunk_size, delim = c(","), col_names = T, na = c("NA", NA)
    ) # na = c("NA", "NA")
  }
  return(data)
}
