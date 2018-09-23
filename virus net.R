
virus <- read.csv("~/Github/brevity/olival nature 2017/associations.csv")
virus <- virus[,c(1:2)]
names(virus) <- c('Host','Parasite')

copredict(5128,virus,20,1)

length(unique(virus$Host))
length(unique(virus$Parasite))
