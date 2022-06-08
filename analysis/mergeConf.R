library(tidyverse)
library(dplyr)
#######
# FDR #
#######
old <- read_csv("./data/dataFinal/noClr/ConfLogNormalAndNegBinom3#preThresh(0.75).csv")
new <- read_csv("./data/dataFinal/Clr/ConfLogNormalAndNegBinom3Clr#preThresh(0.75).csv")


# Seperate old into t-test and not
oldTtest <- old %>% filter(method %in% c("Wilcoxon", "Student.t.test"))
oldOther <- old %>% filter(!method %in% c("Wilcoxon", "Student.t.test"))

# Merge t-test with new
test <- inner_join(oldTtest, new, by = c("method", "run", "sample_size", "n_taxa", "n_random", "change", "class"), suffix = c("_old", "_new"))

# Drop old p-vals column
test <- test %>%  select(-value_old)
test <- as.data.frame(test)

# rename colmn
test <- rename(test, "value_new" = "value")

# rbind df

newFull <- rbind(oldOther, test)

# write file

write_csv(newFull,"./data/dataFinal/FDRLogNormalAndNegBinom3All#preThresh(0.75).csv")

#######
# CONF #
#######

old <- read_csv("./data/dataFinal/noClr/ConfLogNormalAndNegBinom3#preThresh(0.05)Samplesize2.csv")
new <- read_csv("./data/dataFinal/Clr/ConfLogNormalAndNegBinom3Clr#preThresh(0.05)Samplesize2.csv")


# Seperate old into t-test and not
oldTtest <- old %>% filter(method %in% c("Wilcoxon", "Student.t.test"))
oldOther <- old %>% filter(!method %in% c("Wilcoxon", "Student.t.test"))

# Merge t-test with new
test <- inner_join(oldTtest, new, by = c("taxon", "method", "run", "sample_size", "n_taxa", "n_random", "change", "diff", colnames(new)[startsWith(colnames(new), "2")], colnames(new)[startsWith(colnames(new), "1")]), suffix = c("_old", "_new"))
# When using inner_join I get only a data frame of size 2519, so I guess that the data is not perfect.



# Drop old p-vals column
test <- test %>%  select(-pvalueAdj_old, -class_old)

# Rename column
test <- test %>% rename(pvalueAdj_new = "pvalueAdj",
                        class_new = "class")
test <- as.data.frame(test)

# rbind df

newFull <- rbind(oldOther, test)

# write file

write_csv(newFull,"./data/dataFinal/simulatedMergedClrnoClr/ConfLogNormalAndNegBinom3Samplesize2#preThresh(0.05).csv")

###################
# CONF AND EFFECT #
###################

old <- read_csv("./data/dataFinal/noClr/ConfAndEffectLogNormalAndNegBinom3#preThresh(0.05)Samplesize3.csv")
new <- read_csv("./data/dataFinal/Clr/ConfAndEffectLogNormalAndNegBinom3Clr#preThresh(0.05)Samplesize3.csv")


# Seperate old into t-test and not
oldTtest <- old %>% filter(method %in% c("Wilcoxon", "Student.t.test"))
oldOther <- old %>% filter(!method %in% c("Wilcoxon", "Student.t.test"))

# Merge t-test with new
test <- inner_join(oldTtest, new, by = c("taxon", "method", "run", "sample_size", "n_taxa", "n_random", "change", "diff", colnames(new)[startsWith(colnames(new), "2")], colnames(new)[startsWith(colnames(new), "1")]), suffix = c("_old", "_new"))
# When using inner_join I get only a data frame of size 2519, so I guess that the data is not perfect.



# Drop old p-vals column
test <- test %>%  select(-ends_with("old"))#-pvalueAdj_old, -class_old)

# Rename column
test <- test %>% rename(pvalueAdj_new = "pvalueAdj",
                        class_new = "class",
                        effect_new = "effect")
test <- as.data.frame(test)

# rbind df

newFull <- rbind(oldOther, test)

# write file

write_csv(newFull,"./data/dataFinal/simulatedMergedClrnoClr/ConfAndEffectLogNormalAndNegBinom3Samplesize3#preThresh(0.05).csv")


