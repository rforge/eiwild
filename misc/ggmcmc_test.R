

library(ggmcmc)
library(eiwild)
library(dplyr)

data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100,verbose=100)


s <- res$draws$cellCounts[,1:2]
S <- ggs(s)

S2 <- filter(S, grepl("CSU_1", Parameter))

str(S2)

ggs_histogram(S)
ggs_density(S)


data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
out1 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                 sample=1000, thinning=2, burnin=100, verbose=100)
out2 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                 sample=1000, thinning=2, burnin=100, verbose=100)
out3 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                 sample=1000, thinning=2, burnin=100, verbose=100)
out4 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                 sample=1000, thinning=2, burnin=100, verbose=100)


s1 <- out1$draws$cellCounts[,1:2]
s2 <- out2$draws$cellCounts[,1:2]
s3 <- out3$draws$cellCounts[,1:2]
s4 <- out4$draws$cellCounts[,1:2]


S1 <- ggs(s1)
S2 <- ggs(s2)
S2$Chain <- 2

S3 <- ggs(s3)
S3$Chain <- 3

S4 <- ggs(s4)
S4$Chain <- 4


S <- rbind(S1,S2,S3,S4)
attr(S, "nChains") <- 4

ggs_histogram(S)
ggs_density(S)
ggs_traceplot(S)
ggs_running(S)

ggplot(dm.rm, aes(x = Iteration, y = m))

ggplot(dm.rm, aes(x = Iteration, y = rm, colour = as.factor(Chain))) + 
  geom_line() + ylab("Running Mean") + facet_wrap(~Parameter)

ggs_autocorrelation(S)



