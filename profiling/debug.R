dataset.N <- 130
data.source.pollen <- tibble_Europe2$filtered.counts[[dataset.N]]
data.source.age <- tibble_Europe2$list_ages[[dataset.N]]
rand = 10
interest.treshold = 8000
intrapolate = T
BIN = 250
standardise = T
S.value = 150
sm.type = "grim" 
N.points = 5
range.age.max = 300
grim.N.max = 9
DC = "chisq"
Debug = F

data.source.extrap <- data.bin

data.source.pollen.extract <- data.source.pollen
data.source.age.extract <- data.source.age


data.source.pollen.check <- data.source.pollen
data.source.age.check <- data.source.age
proportion = F
Debug=F


tibble_Europe_Roc$ROC[[1]]


which(tibble_Europe2$dataset.id %in%  22936 )

dataset.N <- 130


test <- fc_ratepol( data.source.pollen =  tibble_Europe2$filtered.counts[[dataset.N]],
                    data.source.age = tibble_Europe2$list_ages[[dataset.N]],
                    rand = 10,
                    interest.treshold = 8000,
                    intrapolate = F,
                    BIN = 250,
                    standardise = T, 
                    S.value = 150, 
                    sm.type = "grim", 
                    N.points = 5, 
                    range.age.max = 300, 
                    grim.N.max = 9,
                    DC = "chisq",
                    Debug = F)

test %>% ggplot(aes( y=DF.RoC, 
                     x= DF.Newage))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=DF.RoC.05q, ymax=DF.RoC.95q), color="gray", alpha=1/5)+
  #geom_point(alpha=1/5)+
  geom_line(size=1)+
  geom_point(data = test[test$Peak==T,], color="blue", size=3)+
  #geom_hline(yintercept = mean(test$RoC.median), color="green")+
  geom_hline(yintercept = median(test$DF.RoC), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5))


data.sd <- fc_standar(data.work, S.value, Debug=Debug)
data.sd.check <- fc_check(data.sd, proportion = T, Debug=Debug)

data.sd.t <- fc_standar(data.work.t, S.value, Debug=Debug)
data.sd.check.t <- fc_check(data.sd.t, proportion = T, Debug=Debug)


data.sd.check$Age$newage
data.sd.check.t$Age$newage

data.sd.check$Pollen %>% colSums()
data.sd.check.t$Pollen %>% colSums()


r.m.full %>% ggplot(aes( y=DF.RoC, 
                      x= DF.Newage))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=DF.RoC.05q, ymax=DF.RoC.95q), color="gray", alpha=1/5)+
  #geom_point(alpha=1/5)+
  geom_line(size=1)+
  geom_point(data = r.m.full[r.m.full$Peak==T,], color="blue", size=3)+
  #geom_hline(yintercept = mean(test$RoC.median), color="green")+
  geom_hline(yintercept = median(r.m.full$DF.RoC), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5)) 





data.sub<-tibble_Europe2[c(1:10),]


data.frame(POLLEN=reshape2::melt(tibble_Europe2$filtered.counts[[dataset.N]]), 
           AGE=rep(tibble_Europe2$list_ages[[dataset.N]]$ages$age, ncol(tibble_Europe2$filtered.counts[[dataset.N]]))) %>%
  ggplot(aes( y=POLLEN.value, 
                  x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  #geom_point(alpha=1/5)+
  geom_line(aes(group=POLLEN.variable), alpha=1)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  #geom_density_2d()+
  xlab("Age")+ylab("Pollen")+
  facet_wrap(~POLLEN.variable)+
  coord_flip(ylim = c(0,8000))



n.sets <- 1:length(test)

names(n.sets) <- n.sets

for(i in n.sets[c(1:65,67:70,72:74,76:125,127:131,133:length(n.sets))])
{
  vec<- test[[i]]$age
  for(j in 2:length(vec))
    {
    z<- abs(vec[j]-vec[j-1])
    if(z==0) stop(names(n.sets)[i])
    }
}

# 66, 71, 75, 126, 132

tibble_Europe2$list_ages[[66]]$ages
tibble_Europe2$list_ages[[75]]$ages
tibble_Europe2$list_ages[[126]]$ages

tibble_Europe2$list_ages[[71]]$ages
tibble_Europe2$list_ages[[132]]$ages




r.m.full %>% ggplot(aes( y=RoC, 
                     x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), color="gray")+
  #geom_point(alpha=1/5)+
  geom_line()+
  geom_point(data = r.m.full[r.m.full$Peak==T,], color="red", size=3)+
  geom_hline(yintercept = median(r.m.full$RoC), color="blue")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip()

microbenchmark (
  .subset2(as.data.frame(t(data.source)),i),
  data.source[i,]
)

microbenchmark(
  data.source[,i],
  .subset2(data.source,i)
)


test <- as.data.frame(matrix(1:16,ncol=4,nrow = 4))

library(data.table)

test2 <- as.data.table(matrix(1:16,ncol=4,nrow = 4))

test2[,..i]  
test2[i,]
test2

microbenchmark(
  test[i,],
  slice(test,i),
  test2[i,]
)


microbenchmark(
  test[,i],
  .subset2(test,i),
  pull(test,i),
  select(test,i),
  test2[,..i]
)

microbenchmark(
  test[i,],
  slice(test,i),
  test2[i,],
  test[,i],
  .subset2(test,i),
  pull(test,i),
  select(test,i),
  test2[,..i]
)

microbenchmark(
  .subset2(test,2)[2],
  .subset(.subset2(test,2),2)
)


test3 <- 1:20

microbenchmark(
  .subset(test3,15),
  test3[15]
)




DF.size <- data.frame(matrix(nrow=nrow(data.sub),ncol=3))
names(DF.size) <- c("order","ID","size")
DF.size$order <- 1:nrow(DF.size)

for(i in 1:nrow(DF.size))
{
  DF.size$ID[i] <- data.sub$dataset.id[[i]]
  DF.size$size[i] <- nrow(data.sub$filtered.counts[[i]])
}

DF.size[order(DF.size$size),]


test<-  tibble_Europe2 %>%
  mutate( 
    MAT = map(list_ages, .f = function(x) {
       pluck(x,"age_position")
    }),
    MAT.ed = map(MAT, .f = function(x) {
      x[,8:ncol(x)]
    }),
    HasZero = map(MAT.ed, .f = function(x) {
      ifelse(any(x==0),T,F)
    })
  ) %>%
  select(c("dataset.id","HasZero"))

rm(test)


test2 <- test$HasZero %>% unlist()

c(1:length(test2))[test2]

N <- 238

tibble_Europe2$list_ages[[N]]$age_position %>%
  apply(., 2, FUN= function(x) {any(x==0)})


tibble_Europe2$list_ages[[N]]$age_quantile %>%
  t()

tibble_Europe2$list_ages[[N]]$age_position[,10]
  
tibble_Europe2$list_ages[[N]]$age_position %>%
  as.data.frame() %>%
  select(.,"Pos14")


tibble_Europe2[214,]$dataset.id



s.time <- Sys.time()

tibble_Europe_Roc <-  tibble_Europe2[1:5,] %>%
  mutate(., ROC = map2(filtered.counts,list_ages,
                       .f = function(.x,.y)
                       {res <- fc_ratepol(
                         data.source.pollen = .x,
                         data.source.age = .y,
                         interest.treshold = 8000,
                         rand = 1000,
                         standardise = T, 
                         S.value = 150, 
                         sm.type = "grim", 
                         N.points = 5, 
                         range.age.max = 500, 
                         grim.N.max = 9,
                         DC = "chisq",
                         Debug = F
                       )} ))

f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time


tibble_Europe2$list_ages[[1]]$age_position %>% apply(.,1, FUN= is.unsorted) 

tibble_Europe2$list_ages[[1]]$age_position[1,]


