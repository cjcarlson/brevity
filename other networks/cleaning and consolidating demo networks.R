
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
colnames(sch2010) <- c('Plant','Disperser')
