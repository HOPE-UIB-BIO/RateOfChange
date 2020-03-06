which(tibble_Europe2$dataset.id %in%  45117 )


dataset.N <- 122
data.source.pollen <- tibble_Europe2$filtered.counts[[dataset.N]]
data.source.age <- tibble_Europe2$list_ages[[dataset.N]]
sm.type = "grim" 
N.points = 5
range.age.max = 500 
grim.N.max = 9
BIN = F
BIN.size = 500
Shiftbin = F
shift.value = 100
N.shifts = 5
rand = 100
standardise = T 
S.value = 150 
DC = "chisq"
interest.treshold = 8000
Debug = F

data.source.extrap <- data.bin


data.source <- data.sd

data.subset<- data.work
BINS <- SELECTED.BINS

data.source.bin<- data.smooth.check

data.source.check <- list.res


data.source.pollen.extract <- data.source.pollen
data.source.age.extract <- data.source.age

data.source.pollen.check <- data.source.pollen
data.source.age.check <- data.source.age


dataset.N <- 122
rm(test)

test <- fc_ratepol( data.source.pollen =  tibble_Europe2$filtered.counts[[dataset.N]],
                    data.source.age = tibble_Europe2$list_ages[[dataset.N]],
                    sm.type = "grim" ,
                    N.points = 5,
                    range.age.max = 500, 
                    grim.N.max = 9,
                    BIN = T,
                    BIN.size = 500,
                    Shiftbin = T,
                    N.shifts = 10,
                    rand = 100,
                    standardise = T, 
                    S.value = 150 ,
                    DC = "chisq",
                    interest.treshold = 8000,
                    Debug = F)


test %>% ggplot(aes( y=RUN.RoC, 
                     x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  geom_line(size=1)+
  geom_line(data=data.frame(RUN.RoC = predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = test)),
                            RUN.Age.Pos = test$RUN.Age.Pos),
            color="blue", size=1)+
  geom_point(data = test[test$Peak.gam==T,], color="red", size=3)+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5))


test %>% ggplot(aes(x= RUN.Age.Pos))+
  theme_classic()+
  coord_flip(xlim=c(0,8000))+
  scale_x_continuous(trans = "reverse")+
  geom_density(fill="gray80")



test %>% ggplot(aes( y=RUN.RoC, 
                     x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  #geom_point(alpha=1/5)+
  geom_line(size=1)+
  #geom_point(data = test[test$soft.Peak==T,], color="orange", size=3)+
  #geom_point(data = test[test$Peak==T,], color="blue", size=3)+
  #geom_hline(yintercept = mean(test$RoC.median), color="green")+
  #geom_hline(yintercept = median(test$RUN.RoC), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  geom_smooth(method = "loess")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5))

example <- fc_extract(tibble_Europe2$filtered.counts[[dataset.N]],
           tibble_Europe2$list_ages[[dataset.N]],
           Debug=F) %>%
  fc_smooth(., 
            sm.type = sm.type, 
            N.points = N.points,
            range.age.max = range.age.max,
            grim.N.max = grim.N.max,
            Debug=F) %>%
  fc_check(.,proportion = T)


N.taxa <- 10

Common.list <- example$Pollen %>%
  colSums() %>% 
  sort(decreasing = T) %>%
  .subset(.,1:N.taxa) %>%
  names()


data.frame(POLLEN=example$Pollen %>%
             select(.,Common.list) %>%
             reshape2::melt() , 
           AGE=rep(example$Age$newage, N.taxa))  %>%
  ggplot(aes( y=POLLEN.value, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=rep(0,length(POLLEN.value)),ymax=POLLEN.value, fill=POLLEN.variable), 
              color="gray20", alpha=1/5)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  #facet_wrap(~POLLEN.variable, ncol = N.taxa)+
  xlab("Age")+ylab("Pollen")+
  coord_flip(xlim=c(0,8000), ylim = c(0,1))


Common.list <- tibble_Europe2$filtered.counts[[dataset.N]] %>%
  colSums() %>% 
  sort(decreasing = T) %>%
  .subset(.,1:N.taxa) %>%
  names()


data.frame(POLLEN=tibble_Europe2$filtered.counts[[dataset.N]] %>%
                 select(.,Common.list) %>%
                 reshape2::melt() , 
               AGE=rep(tibble_Europe2$list_ages[[dataset.N]]$ages$age, N.taxa)) %>%
  ggplot(aes( y=POLLEN.value, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=rep(0,length(POLLEN.value)),ymax=POLLEN.value), color="gray20", fill="gray80")+
  facet_wrap(~POLLEN.variable)+
  xlab("Age")+ylab("Pollen")+
  coord_flip(xlim=c(0,8000), ylim = c(0,1000))




s.time <- Sys.time()

tibble_Europe_Roc <-  tibble_Europe2[1:20,] %>%
  mutate(., ROC = map2(filtered.counts,list_ages,
                       .f = function(.x,.y)
                       {res <- fc_ratepol(
                         data.source.pollen = .x,
                         data.source.age = .y,
                         interest.treshold = 8000,
                         rand = 10,
                         intrapolate = F,
                         BIN = 250,
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


r.m.full %>% ggplot(aes( y=RUN.RoC, 
                      x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  #geom_point(alpha=1/5)+
  geom_line(size=1)+
  #geom_point(data = r.m.full[r.m.full$soft.Peak==T,], color="blue", size=3)+
  #geom_hline(yintercept = mean(test$RoC.median), color="green")+
  geom_hline(yintercept = median(r.m.full$RUN.RoC), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,3)) 


r.m.full %>% ggplot(aes( x=RUN.Age.Pos))+
  geom_density(fill="gray80")+
  theme_classic()+
  coord_flip()+
  scale_x_continuous(trans = "reverse")
  
  

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

abe <- c(1:10)

t <- "abe"

get(noquote(t))

noquote(t)
