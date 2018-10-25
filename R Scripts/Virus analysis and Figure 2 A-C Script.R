library(ggplot2)
library(cowplot)
library(codependent)

### PART A. TOTAL VIRAL RICHNESS ###

# Read in data 
virus <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/associations.csv")
virus <- virus[,c(2:1)]
names(virus) <- c('Host','Parasite')
# Exclude humans
virus <- virus[!(virus$Host %in% c('Homo_sapiens')),]

### Predict total richness ### 
copredict(5291,virus,200,1)

# Number of hosts
length(unique(virus$Host))
# Number of viruses
length(unique(virus$Parasite))
# These get used to scale out estimates later

# Average host range of a virus
mean(data.frame(table(virus$Host))$Freq)









### PART B. Separate analyses for DNA and RNA viruses ### 

virusmeta <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/viruses.csv")

dnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='DNA',]$vVirusNameCorrected)
rnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='RNA',]$vVirusNameCorrected)

# Compile sub-association lists
rnavirus <- virus[virus$Parasite %in% rnalist,]
dnavirus <- virus[virus$Parasite %in% dnalist,]

# Estimation
copredict(5291,rnavirus,200,1)
copredict(5291,dnavirus,200,1)

# Derive total proportion of zoonotic RNA and DNA viruses 
table(virusmeta$vDNAoRNA,virusmeta$IsZoonotic)
29/(29+176)
159/(159+222)

# Assign zoonotic data to the virus dataset to derive by-group estimates
virus$zoonotic=0
for (i in 1:nrow(virus)) {
virus$zoonotic[i] <- virusmeta[virusmeta$vVirusNameCorrected==virus$vVirusNameCorrected[i],]$IsZoonotic
}

# Read in host information
hostmeta <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/hosts.csv")
virus$hostgroup = ''

# Pull host names into virus data
for (i in 1:nrow(virus)) {
  virus$hostgroup[i] <- as.character(hostmeta[hostmeta$hHostNameFinal==virus$hHostNameFinal[i],]$hOrder)
}

# Zoonotic rates by group
virus <- virus[,c(1,7,8)]
virus <- unique(virus)
df <- table(virus$hostgroup,virus$zoonotic)










### PART C. UPPER BOUNDS ###

# Upper bound on total viral richness
copredict.ci(5291,virus,200,1)

# Re-partition by DNA/RNA viruses
virusmeta <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/viruses.csv")

dnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='DNA',]$vVirusNameCorrected)
rnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='RNA',]$vVirusNameCorrected)

rnavirus <- virus[virus$Parasite %in% rnalist,]
dnavirus <- virus[virus$Parasite %in% dnalist,]

# Analysis for RNA viruses only
copredict.ci(5291,rnavirus,200,1)

# [1] "True number of species in entire network is 345"
# [1] "Estimated number of species in entire network is 393.590434663421"
# [1] "The lower 95% CI is 389.705103384448"
# [1] "The upper 95% CI is 397.514502410087"
# [1] "Extrapolated estimated number of species is 1156.54367524801"

14.92537*1157

17269*159/(159+222)

# [1] "The lower 95% CI is 1142.10898172487"

14.92537*1142 

17045*159/(159+222)
# [1] "The upper 95% CI is 1171.16080353039"

14.92537*1171

17478*159/(159+222)

# Analysis for DNA viruses only 
copredict.ci(5291,dnavirus,200,1)

# [1] "Extrapolated estimated number of species is 2092.91700248616"

14.92537*2093

31239*159/(159+222)
# [1] "The lower 95% CI is 2035.41012410777"

14.92537*2035

30373*159/(159+222)

# [1] "The upper 95% CI is 2152.0486350218"

14.92537*2152

32119*159/(159+222)













### PART D. DEVELOPING CODE FOR FIGURES ###

# Assemble sub-samples to plot
df.virus <- curve.df(virus,100)
df.rna <- curve.df(rnavirus,100)
df.dna <- curve.df(dnavirus,100)

# Fit curves for 2A
model.all <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.virus)
predicted.all <- data.frame(pred = predict(model.all), host = df.virus$n.host)

#Plot 2A
g1 <- ggplot(df.virus, aes(n.host, n.par)) + xlim(0,775) + ylim(0,650) + xlab('Mammals') + ylab('All viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('mediumpurple1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.all, aes(x=host, y=pred))
g1
model.all

# Fit curves for 2B
model.rna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.rna)
predicted.rna <- data.frame(pred = predict(model.rna), host = df.rna$n.host)

# Plot 2B
g2 <- ggplot(df.rna, aes(n.host, n.par)) + xlim(0,740) + ylim(0,370) + xlab('Mammals') + ylab('RNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('mediumorchid1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.rna, aes(x=host, y=pred))
g2
model.rna

# Fit curves for 2C
model.dna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.dna)
predicted.dna <- data.frame(pred = predict(model.dna), host = df.dna$n.host)

# Plot 2C
g3 <- ggplot(df.dna, aes(n.host, n.par)) + xlim(0,180) + ylim(0,170) + xlab('Mammals') + ylab('DNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('slateblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.dna, aes(x=host, y=pred))
g3
model.dna

# Plot everything together, for final figure assembly
plot1 <- plot_grid(g1, g2, g3,
                      labels=c("A", "B", "C"), nrow=1, ncol = 3)

plot2 <- plot_grid(NULL, NULL,
                      labels=c("D", "E"), nrow=1, ncol = 2)

plot3 <- plot_grid(blank, blank, nrow=1, ncol = 2)

plot_grid(plot1,plot2,plot2,nrow=3)













### PART E. SUPPLEMENTAL FIGURE 1. ###


virus <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/associations.csv")
virus <- virus[,c(2:1)]
names(virus) <- c('Host','Parasite')

df.virus <- binera(virus,100,plots=true)
