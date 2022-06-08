
cleanData <- function(phylos){
phylos <- prune_samples(!sample_sums(phylos) < 19000, phylos)
phylos <- prune_samples(!sample_sums(phylos) > 40000, phylos)
phylos <- prune_taxa(!taxa_sums(phylos) == 0, phylos)
# phylos <- filter_taxa(phylos, function(x){sum(x > 0) > nsamples(phylos)*0.05}, prune = TRUE)
# phylos <- prune_taxa(names(sort(taxa_sums(phylos),TRUE)[1:1000]), phylos)
# phylos <- prune_samples(!sample_sums(phylos) > 50000, phylos)
return(phylos)
}
