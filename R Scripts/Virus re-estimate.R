library(ggplot2)
library(cowplot)
library(codependent)

### PART A. TOTAL VIRAL RICHNESS ###

# Read in data 
virus <- read.csv("../Github/brevity/Olival Nature 2017 Raw Data/associations.csv")
virus <- virus[,c(2:1)]
names(virus) <- c('Host','Parasite')
# Exclude humans
virus <- virus[!(virus$Host %in% c('Homo_sapiens')),]

# Number of hosts
length(unique(virus$Host))
# Number of viruses
length(unique(virus$Parasite))
# These get used to scale out estimates later

# Average host range of a virus
mean(data.frame(table(virus$Host))$Freq)

### Predict total richness ### 
raw.ests <- copredict(virus,5291,1000)
raw.ests

### Derive the rates of sampling ###

kevin <- read.csv("../Github/brevity/Olival Nature 2017 Raw Data/associations.csv")
mac.rare <- read.csv('../GitHub/GVP_Science/data/output/macaca_estimates_by_viral_family.csv')
mac.total = 184
mac.raw <- kevin[kevin$hHostNameFinal=='Macaca_mulatta',]$vVirusNameCorrected

mac.est <- sum(mac.rare$Estimator) + (mac.total-sum(mac.rare$Observed))
mac.lci <- sum(mac.rare$LCL) + (mac.total-sum(mac.rare$Observed))
mac.uci <- sum(mac.rare$UCL) + (mac.total-sum(mac.rare$Observed))

mac.rate <- length(mac.raw)/c(mac.est,mac.lci,mac.uci)


pt.rare <- read.csv('../GitHub/GVP_Science/data/output/pteropus_estimates_by_viral_family.csv')
pt.total = 55
pt.raw <- kevin[kevin$hHostNameFinal=='Pteropus_giganteus',]$vVirusNameCorrected

pt.est <- sum(pt.rare$Estimator) + (pt.total-sum(pt.rare$Observed))
pt.lci <- sum(pt.rare$LCL) + (pt.total-sum(pt.rare$Observed))
pt.uci <- sum(pt.rare$UCL) + (pt.total-sum(pt.rare$Observed))

pt.rate <- length(pt.raw)/c(pt.est,pt.lci,pt.uci)

rates <- (pt.rate+mac.rate)/2
names(rates) <- c('est','lower','upper')

corrected <- raw.ests[[1]]/rates
corrected

# 50% check

raw.ests.50 <- copredict.ci(virus,5291,100)
raw.ests.50[,3] <- dmm::unfactor(raw.ests.50[,3])
raw.ests.50[,4] <- dmm::unfactor(raw.ests.50[,4])

raw.ests.50

corrected.50 <- raw.ests.50[3,2:4]/rates
corrected.50








### PART B. Separate analyses for DNA and RNA viruses ### 

virusmeta <- read.csv("../Github/brevity/Olival Nature 2017 Raw Data/viruses.csv")

dnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='DNA',]$vVirusNameCorrected)
rnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='RNA',]$vVirusNameCorrected)


# Compile sub-association lists
rnavirus <- virus[virus$Parasite %in% rnalist,]
dnavirus <- virus[virus$Parasite %in% dnalist,]

# Estimation
raw.ests.dna <- copredict(dnavirus,5291,1000)
raw.ests.rna <- copredict(rnavirus,5291,1000)

raw.ests.dna
raw.ests.rna

raw.ests.50.dna <- copredict.ci(dnavirus,5291,1000)
raw.ests.50.rna <- copredict.ci(rnavirus,5291,1000)

raw.ests.50.dna
raw.ests.50.rna

corrected.dna <- raw.ests.dna[[1]]/rates
corrected.rna <- raw.ests.rna[[1]]/rates

corrected.dna
corrected.rna


raw.ests.50.dna[,3] <- dmm::unfactor(raw.ests.50.dna[,3])
raw.ests.50.dna[,4] <- dmm::unfactor(raw.ests.50.dna[,4])
corrected.50.dna <- raw.ests.50.dna[3,2:4]/rates
corrected.50.dna

raw.ests.50.rna[,3] <- dmm::unfactor(raw.ests.50.rna[,3])
raw.ests.50.rna[,4] <- dmm::unfactor(raw.ests.50.rna[,4])
corrected.50.rna <- raw.ests.50.rna[3,2:4]/rates
corrected.50.rna

corrected.dna+corrected.rna

corrected.50.dna+corrected.50.rna


corrected.dna*0.141
corrected.rna*0.417
corrected.dna*0.141 + corrected.rna*0.417


corrected.50.dna*0.141
corrected.50.rna*0.417
corrected.50.dna*0.141 + corrected.50.rna*0.417





#############
############# EVERYTHING BELOW THIS POINT IS SUPPLEMENT
############# WHERE WE TRIED DOING SEPARATE RATES OF SAMPLING
############# FOR DNA AND RNA VIRUSES 
#############










############## DNA NUMBERS



### Derive the rates of sampling ###

mac.rare <- read.csv('../GitHub/GVP_Science/data/output/macaca_estimates_by_viral_family.csv')
mac.rare <- mac.rare[mac.rare$Class=='DNA',]
mac.total = 23
mac.raw <- kevin[kevin$hHostNameFinal=='Macaca_mulatta',]$vVirusNameCorrected

mac.est <- sum(mac.rare$Estimator) + (mac.total-sum(mac.rare$Observed))
mac.lci <- sum(mac.rare$LCL) + (mac.total-sum(mac.rare$Observed))
mac.uci <- sum(mac.rare$UCL) + (mac.total-sum(mac.rare$Observed))

mac.rate <-  12/c(mac.est,mac.lci,mac.uci)


pt.rare <- read.csv('../GitHub/GVP_Science/data/output/pteropus_estimates_by_viral_family.csv')
pt.rare <- pt.rare[pt.rare$Class=='DNA',]
pt.total = 32
pt.raw <- kevin[kevin$hHostNameFinal=='Pteropus_giganteus',]$vVirusNameCorrected

pt.est <- sum(pt.rare$Estimator) + (pt.total-sum(pt.rare$Observed))
pt.lci <- sum(pt.rare$LCL) + (pt.total-sum(pt.rare$Observed))
pt.uci <- sum(pt.rare$UCL) + (pt.total-sum(pt.rare$Observed))

pt.rate <- 0/c(pt.est,pt.lci,pt.uci)

rates <- (pt.rate+mac.rate)/2
names(rates) <- c('est','upper','lower')

corrected <- raw.ests.dna[[1]]/rates
corrected


raw.ests.50.dna[,3] <- dmm::unfactor(raw.ests.50.dna[,3])
raw.ests.50.dna[,4] <- dmm::unfactor(raw.ests.50.dna[,4])
corrected.50.dna <- raw.ests.50.dna[3,2:4]/rates
corrected.50.dna


############## RNA NUMBERS



### Derive the rates of sampling ###

mac.rare <- read.csv('../GitHub/GVP_Science/data/output/macaca_estimates_by_viral_family.csv')
mac.rare <- mac.rare[mac.rare$Class=='RNA',]
mac.total = 161
mac.raw <- kevin[kevin$hHostNameFinal=='Macaca_mulatta',]$vVirusNameCorrected

mac.est <- sum(mac.rare$Estimator) + (mac.total-sum(mac.rare$Observed))
mac.lci <- sum(mac.rare$LCL) + (mac.total-sum(mac.rare$Observed))
mac.uci <- sum(mac.rare$UCL) + (mac.total-sum(mac.rare$Observed))

mac.rate <-  10/c(mac.est,mac.lci,mac.uci)


pt.rare <- read.csv('../GitHub/GVP_Science/data/output/pteropus_estimates_by_viral_family.csv')
pt.rare <- pt.rare[pt.rare$Class=='RNA',]
pt.total = 23
pt.raw <- kevin[kevin$hHostNameFinal=='Pteropus_giganteus',]$vVirusNameCorrected

pt.est <- sum(pt.rare$Estimator) + (pt.total-sum(pt.rare$Observed))
pt.lci <- sum(pt.rare$LCL) + (pt.total-sum(pt.rare$Observed))
pt.uci <- sum(pt.rare$UCL) + (pt.total-sum(pt.rare$Observed))

pt.rate <- 3/c(pt.est,pt.lci,pt.uci)

rates <- (pt.rate+mac.rate)/2
names(rates) <- c('est','upper','lower')

corrected <- raw.ests.rna[[1]]/rates
corrected


raw.ests.50.rna[,3] <- dmm::unfactor(raw.ests.50.rna[,3])
raw.ests.50.rna[,4] <- dmm::unfactor(raw.ests.50.rna[,4])
corrected.50.rna <- raw.ests.50.rna[3,2:4]/rates
corrected.50.rna