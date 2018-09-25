
library(reshape)

# 1. Robertson 1929: plant-pollinator interactions

rob1929_raw <- read.csv("~/Github/brevity/other networks/robertson1929.csv")
rob1929 <- rob1929_raw[,c('plant','poll')]
colnames(rob1929) <- c("Plant","Pollinator")

# 2. Schleuning 2010: seed-disperser interactions

sch2010_raw <- read.csv("~/Github/brevity/other networks/schleuning2010.csv")
sch2010_raw[sch2010_raw>1] <- 1
sch2010 <- melt(sch2010_raw, id=c("Plant.species"))
sch2010 <- sch2010[sch2010$value==1,]
sch2010 <- sch2010[,c(2,1)]
colnames(sch2010) <- c("Plant","Disperser")

# 3. Toju 2018: plant-arbuscular mycorrhizae

tojuraw.1 <- read.csv('~/Github/brevity/other networks/toju descriptors.csv')
tojuraw.2 <- read.csv('~/Github/brevity/other networks/toju.csv')
otu.list <- unique(tojuraw.1[tojuraw.1$Category.In.This.Study=='Arbuscular_Mycorrhizal',]$OTU.code)

tojuraw.2$Plant <- gsub('_TES', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_TOMA', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YKS', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_SGD', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YSD', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YAKU', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YONA', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_IRI', '', tojuraw.2$Plant)

toju <- aggregate(. ~ Plant, data=tojuraw.2, FUN=sum)
toju <- melt(toju, id=c("Plant"))
toju <- toju[toju$value>1,]
colnames(toju)[1:2] <- c("Plant","Microbe")
toju$Microbe <- gsub('X', '', toju$Microbe)
toju$Microbe <- gsub('\\.', ':', toju$Microbe)
toju <- toju[toju$Microbe %in% as.character(otu.list),]
toju <- toju[,c(1:2)]

# 4. Host-helminth relationships
