##############################################################
################### FIGURES FOR THESIS #######################
##############################################################
library(ggplot2)
set.seed(15456)
n <- 10000


TaxaSize <- data.frame(rnbinom(n, mu = 1000, size = 5))
colnames(TaxaSize)[1] <- "x"

ggplot(TaxaSize,aes(x = x))+
  geom_histogram(fill = "grey70", color = "grey10", bins = 40)+
  ylab("Count")+
  xlab("Taxa abundance")+
  geom_vline(xintercept = 1000, color = "firebrick", size = 1)+
  theme_classic()


test <- rlnorm(n, meanlog = log(1), sdlog = log(2.718282*7))
p <- hist(log(test))






##############################################################
################### EXPLORATORY ##############################
##############################################################

set.seed(15466)
n <- 8000
zeroesProb <- 0.96 # normal distributed!!
trueZero <- rbinom(n = n, size = 1, prob = zeroesProb)
test <- ifelse(trueZero==1, 0, rlnorm(sum(trueZero==0), meanlog = log(4), sdlog = log(2.718282*5)))
data <- rnbinom(n, mu = test, size = 100)
data
p <- density(log(data), bw = 0.2484)
plot(p, xlim = range(-1,9))

sum(test == 0)/n
sum(data == 0)/n

mean(data)
sd(data)
sum(data)
max(data)



######### TRY TO INTRODUCE ZEROES AFTER SIMULATION #############
set.seed(15456)
n <- 8000
zeroesProb <- 0.95 # normal distributed!!

test <- rlnorm(n, meanlog = log(1), sdlog = log(2.718282*7))
data <- rnbinom(n, mu = test, size = 100)

trueZero <- rbinom(n = n, size = 1, prob = zeroesProb)
data[trueZero == 1] <- 0


data
p <- density(log(data), bw = 0.2484)
plot(p, xlim = range(-1,9))

sum(test == 0)/n
sum(data == 0)/n

mean(data)
sd(data)
sum(data)
max(data)














##############################################################
################### SIMULATION ##############################
##############################################################

# Alle parametrene ligger i den lave ende. Hvordan ser det ud efter diff abundance?






df[random_taxon, startsWith(colnames(df), "2")] <- df[random_taxon, startsWith(colnames(df), "2")]*change







# Simular mange flere samples og se om de minder om figurerene for rigtigt data





















