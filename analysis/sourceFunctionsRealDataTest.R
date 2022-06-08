###########################################
######## FUNCTION TO PREPARE DATA #########
###########################################

prepareData <- function(data){
  final <- data.frame()
  meta <- as.data.frame(as.matrix(sample_data(data)))
  counts <- as.data.frame(otu_table(data))


  meta$rownames <- rownames(meta)

  # Create column with rownames used to indicate treatment(days) and replicate
  meta <- meta %>% mutate(newNames = ifelse(days == "116", paste(2, replicate, sep = "."), paste(1, replicate, sep = ".")) )


  # Split data into individual data frames used in the comparison
  test <- meta %>% group_by(location, oil) %>% group_split()# Make list of df


  finalData <- list()
  for(i in 1:length(test)){

  # For each df get the respective counts (mapped by sample name)
  c <- counts %>% select(all_of(test[[i]]$rownames))

  # Transpose count df and sort
  c <- as.data.frame(t(c))
  c$rownames <- rownames(c)
  c <- c %>% arrange(rownames)

  # Sort meta df
  test[[i]] <- test[[i]] %>% arrange(rownames)

  # cbind count and meta (Can't be merged due to loss of precision)
  test3 <- cbind(test[[i]], c)

  # Manipulate
  test3 <- test3 %>% select(-rownames)

  # Reshape data frame
  test3 <- test3 %>% pivot_longer(cols = starts_with("ASV")) # Is correct
  test3 <- test3 %>% pivot_wider(id_cols = c(name, oil,location), names_from = newNames, values_from = value)

  # Rename column
  colnames(test3)[1] <- "taxon"

  # Sort columns
  test3 <- test3 %>% select(taxon, oil, location, starts_with("1."), starts_with("2."))

  # Add to final list
  finalData <- append(finalData, list(test3))
  }


  return(finalData)

}

########################################################
################# FUNCTION TO RUN TESTS ################
########################################################

# Recreate run parameter somehow
runTestsReal <- function(tests, countData, meta, variables, seed) {
  print("runTests")
  set.seed(seed)

  # Define variables and phyloseq object
  n_taxa <- dim(countData)[1] # number of taxa
  n_samples_total <- dim(countData[startsWith(colnames(countData), "2") | startsWith(colnames(countData), "1")])[2]
  n_tests <- length(tests)
  phyloData <- phyloseq(otu_table(countData[startsWith(colnames(countData), "2") | startsWith(colnames(countData), "1")], taxa_are_rows = T), sample_data(meta)) # Bruges i ANCOM-BE

  # Make dataframes for p-values and effect sizes
  pvals <- cbind(data.frame("taxon" = countData$taxon,
                            "oil" = countData$oil,
                            "location" = countData$location))
  pvalsAdj <- cbind(data.frame("taxon" = countData$taxon,
                               "oil" = countData$oil,
                               "location" = countData$location))
  effect <- cbind(data.frame("taxon" = countData$taxon,
                             "oil" = countData$oil,
                             "location" = countData$location))

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
    # res <- merge(x = countData, y = resALDEX2, by = 0, all = T)
    # res <- merge(x = pvals, y = res, by.x = "taxon", by.y =  "Row.names")
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

    # pv <- merge(countData[1], resancombc$p_val, by = 0, all = T)
    # pv <- merge(x = pvals, y = pv, by.x = "taxon", by.y = "Row.names", all = T)$treatment1
    # ev <- merge(countData[1], resancombc$beta, by = 0, all = T)
    # ev <- merge(x = effect, y = ev, by.x = "taxon", by.y = "Row.names", all = T)$treatment1
    pvals$ANCOM.BC <- remove_rownames(resancombc$p_val)$treatment1
    pvalsAdj$ANCOM.BC <- p.adjust(unlist(remove_rownames(resancombc$p_val)$treatment1), method = "BH")
    effect$ANCOM.BC_effect <- remove_rownames(resancombc$beta)$treatment1

  }

  # names(lst) <- c("pvalues", "effect-size")
  if("DESEq2" %in% tests){
    pvals <- pvals %>% mutate(DESEq2 = replace_na(DESEq2, 1)) # Some are removed due to the filterin step
  }
  result <- inner_join(pvals, pvalsAdj, by = c("taxon", "oil", "location"), suffix = c("_pvals","_pvalsAdj"))
  result <- inner_join(result, effect, by = c("taxon", "oil", "location"))
  return(result)
}





runAllReal <- function(data, tests, prevalence, s){
  print("runAll")
  diff_data_list <- data.frame()
  results <- data.frame()


  for(df in data){
    # Remove zeros and apply prevalence filter
    df <- cleanData(df, prevalence)

    n_samples <- sum(startsWith(colnames(df), "1"))
    metaData <- data.frame("treatment" = c(rep("1", times = n_samples),rep("0", times = n_samples)), row.names = colnames(df[,startsWith(colnames(df), "1") | startsWith(colnames(df), "2")]))
    res <- runTestsReal(tests = tests, countData = df, meta = metaData, variables = v, seed = s)
    diff_data_list <- bind_rows(diff_data_list, df)
    results <- rbind(results, res)

}
  # Reorder columns (They are in correct order for each single data frame when testing. meta data and the columns are matching!)
  diff_data_list <- diff_data_list %>% select(taxon, oil, location, starts_with(c("1.")), starts_with(c("2.")))
  all_data <- full_join(diff_data_list, results)
  return(all_data)
}

runAllRealDup <- function(data, tests, prevalence, s, pair){
  print("runAll")
  diff_data_list <- data.frame()
  results <- data.frame()


  for(df in data){
    # Subsample duplicates
    df <- df %>% select(taxon, oil, location, ends_with(pair))

    # Remove zeros and apply prevalence filter
    df <- cleanData(df, prevalence)

    n_samples <- sum(startsWith(colnames(df), "1"))
    metaData <- data.frame("treatment" = c(rep("1", times = n_samples),rep("0", times = n_samples)), row.names = colnames(df[,startsWith(colnames(df), "1") | startsWith(colnames(df), "2")]))
    res <- runTestsReal(tests = tests, countData = df, meta = metaData, variables = v, seed = s)
    diff_data_list <- bind_rows(diff_data_list, df)
    results <- rbind(results, res)

  }
  # Reorder columns (They are in correct order for each single data frame when testing. meta data and the columns are matching!)
  diff_data_list <- diff_data_list %>% select(taxon, oil, location, starts_with(c("1.")), starts_with(c("2.")))
  all_data <- full_join(diff_data_list, results)
  return(all_data)
}

###################################################
####### FUNCTIONS TO WORK WITH REAL DATA ##########
###################################################

calculateConfusionAndEffectReal2 <- function(data){
  p <- data %>% select(-ends_with(c("pvals"))) %>% pivot_longer(cols = !c("taxon", "oil", "location", starts_with(c("1", "2"))), names_to = "method", values_to = "pvalueAdj")
  pEffect <- p %>% filter(grepl("effect",method)) %>% rename(pvalueAdj = "effect")
  pEffect$method <- unlist(strsplit(pEffect$method, "_"))[!unlist(strsplit(pEffect$method, "_")) == "effect"] # Remove _effect suffix

  pPvals <- p %>% filter(grepl("pvalsAdj",method))
  pPvals$method <- unlist(strsplit(pPvals$method, "_"))[!unlist(strsplit(pPvals$method, "_")) == "pvalsAdj"] # Remove _pvalsAdj suffix
  p <- inner_join(pEffect, pPvals)

  p <- p %>% mutate(class = ifelse(pvalueAdj < 0.05, "P", "N"))
  p <- p %>% mutate(class = replace_na(class, "NA"))
  return(as.data.frame(p))

}

calculateConfusionAndEffectReal2pval <- function(data){
  p <- data %>% select(-ends_with(c("pvalsAdj"))) %>% pivot_longer(cols = !c("taxon", "oil", "location", starts_with(c("1", "2"))), names_to = "method", values_to = "pvalue")
  pEffect <- p %>% filter(grepl("effect",method)) %>% rename(pvalue = "effect")
  pEffect$method <- unlist(strsplit(pEffect$method, "_"))[!unlist(strsplit(pEffect$method, "_")) == "effect"] # Remove _effect suffix

  pPvals <- p %>% filter(grepl("pvals",method))
  pPvals$method <- unlist(strsplit(pPvals$method, "_"))[!unlist(strsplit(pPvals$method, "_")) == "pvals"] # Remove _pvalsAdj suffix
  p <- inner_join(pEffect, pPvals, by = c("taxon", "oil", "location", "1.1", "1.2", "1.3", "2.1", "2.2", "2.3", "method"))

  p <- p %>% mutate(class = ifelse(pvalue < 0.05, "P", "N"))
  p <- p %>% mutate(class = replace_na(class, "NA"))
  return(as.data.frame(p))

}


calculateConfusionAndEffectReal <- function(data){
  p <- data %>% select(-ends_with(c("pvals"))) %>% pivot_longer(cols = !c("taxon", "oil", "location", starts_with(c("1", "2"))), names_to = "method", values_to = "pvalueAdj")
  pEffect <- p %>% filter(grepl("effect",method)) %>% rename(pvalueAdj = "effect")
  pEffect$method <- unlist(strsplit(pEffect$method, "_"))[!unlist(strsplit(pEffect$method, "_")) == "effect"] # Remove _effect suffix

  pPvals <- p %>% filter(grepl("pvalsAdj",method))
  pPvals$method <- unlist(strsplit(pPvals$method, "_"))[!unlist(strsplit(pPvals$method, "_")) == "pvalsAdj"] # Remove _pvalsAdj suffix
  p <- inner_join(pEffect, pPvals)

  p <- p %>%  pivot_longer(cols = !c("taxon", "oil", "location","method", "effect", "pvalueAdj"), names_to = "sample", values_to = "counts")
  p <- p %>% mutate(pvalueAdj = replace_na(pvalueAdj, 1))

  p <- p %>% mutate(class = ifelse(pvalueAdj < 0.05, "P", "N"))
  p <- p %>% mutate(class = replace_na(class, "NA"))
  return(as.data.frame(p))

}

makeVennReal <- function(conf_data, oilt, loc, m) {

  all_counts <- data.frame() # Mean over all runs


for(k in loc){
  for(o in oilt){
    if(k == "C" & o == "blank"){
      next
    }
    data <- conf_data %>% mutate(row = row_number()) %>% filter(oil == o, location == k, method == m) %>%  pivot_wider(names_from = df, values_from = taxon) %>% select(-row)
    l <- list()


    if(ncol(data) <= 9){
      for(i in c("1+2+3","1+2","1+3","2+3")[!c("1+2+3","1+2","1+3","2+3") %in% colnames(data)]){
        df <- data.frame("X" = rep(NA, nrow(data)))
        colnames(df) <- i
        data <- cbind(data, df)

      }}

    # Make df with results
    for(i in 7:ncol(data)){
      l <- append(l,na.omit(data[i]))

    }


    # Rename results
    names(l) <- paste("Replicates: ", colnames(data)[7:ncol(data)])
    # l <- l[order(names(l))]

    sum <- 0
    for(list in l){
      sum <- sum + length(list)
    }

    if(sum == 0){
      cs <- rep(0,16) #12 = number of elements in Venn diagram with 4 groups
    } else{

      p <- venn(l)
      cs <- p$counts
    }
    all_counts <- rbind(all_counts, cs)
}}

  mean_counts <- round(colSums(all_counts)/(length(oilt)+length(loc))) # round
    # mean_counts <- colSums(all_counts)
  v_plot <-   venn(length(l), snames = names(l), counts = mean_counts, zcolor = scales::viridis_pal(option = "D")(length(names(l))), opacity = 0.4, box = F,
                   par = F, ggplot = F, ellipse = T) #c(ilcs = 10, sncs = 10) #ilcs = 1.5, sncs = 1.5,
  # zcolor = brewer.pal(length(names(l)), "Dark2")
  #zcolor = scales::viridis_pal(option = "D")(length(names(l))

  scales::viridis_pal()(4)
  return(v_plot)
}


makeVennReal2 <- function(conf_data, oilt, loc, m) {

  all_counts <- data.frame() # Mean over all runs


  for(k in loc){
    for(o in oilt){
      if(k == "C" & o == "blank"){
        next
      }
      data <- conf_data %>% mutate(row = row_number()) %>% filter(oil == o, location == k, method == m) %>%  pivot_wider(names_from = filter, values_from = taxon) %>% select(-row)
      l <- list()


      if(ncol(data) <= 9){
        for(i in c("0.05","0.25","0.75")[!c("0.05","0.25","0.75") %in% colnames(data)]){
          df <- data.frame("X" = rep(NA, nrow(data)))
          colnames(df) <- i
          data <- cbind(data, df)

        }}

      # Make df with results
      for(i in 7:ncol(data)){
        l <- append(l,na.omit(data[i]))

      }


      # Rename results
      names(l) <- paste("Filter: ", colnames(data)[7:ncol(data)])
      l <- l[order(names(l))]

      sum <- 0
      for(list in l){
        sum <- sum + length(list)
      }

      if(sum == 0){
        cs <- rep(0,8) #7 = number of elements in Venn diagram with 4 groups +1
      } else{

        p <- venn(l)
        cs <- p$counts
      }
      all_counts <- rbind(all_counts, cs)
    }}

  mean_counts <- round(colSums(all_counts)/(length(oilt)+length(loc))) # round
  # mean_counts <- colSums(all_counts)
  v_plot <-   venn(length(l), snames = names(l), counts = mean_counts, zcolor = scales::viridis_pal(option = "D")(length(names(l))), opacity = 0.4, box = F,
                   borders = T, par = T, ggplot = F, ellipse = F) #c(ilcs = 10, sncs = 10) #ilcs = 1.5, sncs = 1.5,

  # zcolor = brewer.pal(length(names(l)), "Dark2")
  #zcolor = scales::viridis_pal(option = "D")(length(names(l))

  # scales::viridis_pal()(4)
  return(v_plot)
}

makeVennReal3 <- function(conf_data, oilt, loc, m) {

  all_counts <- data.frame() # Mean over all runs


  for(k in loc){
    for(o in oilt){
      if(k == "C" & o == "blank"){
        next
      }
      data <- conf_data %>% mutate(row = row_number()) %>% filter(oil == o, location == k) %>%  pivot_wider(names_from = method, values_from = taxon) %>% select(-row)
      l <- list()


      if(ncol(data) <= 9){
        for(i in c("ANCOM.BC","DESEq2","ALDEx2", "Student.t.test", "Wilcoxon")[!c("ANCOM.BC","DESEq2","ALDEx2", "Student.t.test", "Wilcoxon") %in% colnames(data)]){
          df <- data.frame("X" = rep(NA, nrow(data)))
          colnames(df) <- i
          data <- cbind(data, df)

        }}

      # Make df with results
      for(i in 6:ncol(data)){
        l <- append(l,na.omit(data[i]))

      }


      # Rename results
      names(l) <- colnames(data)[6:ncol(data)]
      l <- l[order(names(l))]

      sum <- 0
      for(list in l){
        sum <- sum + length(list)
      }

      if(sum == 0){
        cs <- rep(0,32) #32 = number of elements in Venn diagram with 5 groups +1
      } else{

        p <- venn(l)
        cs <- p$counts
      }
      all_counts <- rbind(all_counts, cs)
    }}

  mean_counts <- round(colSums(all_counts)/(length(oilt)+length(loc))) # round
  # mean_counts <- colSums(all_counts)
  v_plot <-   venn(length(l), snames = names(l), counts = mean_counts, zcolor = scales::viridis_pal(option = "D")(length(names(l))), opacity = 0.4, box = F,
                   borders = T, par = T, ggplot = F, ellipse = T) #c(ilcs = 10, sncs = 10) #ilcs = 1.5, sncs = 1.5,

  # zcolor = brewer.pal(length(names(l)), "Dark2")
  #zcolor = scales::viridis_pal(option = "D")(length(names(l))

  # scales::viridis_pal()(4)
  return(v_plot)
}

# https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}
