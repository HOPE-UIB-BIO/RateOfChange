dataset.N <- 2
sm.type = "age.w" 
N.points = 5
range.age.max = 500 
grim.N.max = 9
BIN = T
BIN.size = 500
Shiftbin = T
shift.value = 100
N.shifts = 5
rand = 1000
standardise = T 
S.value = 150 
DC = "chord"
interest.treshold = 8000
Debug = F


test <- fc_ratepol( data.source.pollen =  tibble_Europe2$filtered.counts[[dataset.N]],
                    data.source.age = tibble_Europe2$list_ages[[dataset.N]],
                    sm.type = sm.type,
                    N.points = N.points,
                    range.age.max = range.age.max, 
                    grim.N.max = grim.N.max,
                    BIN = BIN,
                    BIN.size = BIN.size,
                    Shiftbin = Shiftbin,
                    N.shifts = N.shifts,
                    rand = rand,
                    standardise = standardise, 
                    S.value = S.value ,
                    DC = DC,
                    interest.treshold = interest.treshold,
                    Debug = Debug)


test %>% ggplot(aes( y=RUN.RoC, 
                     x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  geom_line(size=1)+
  geom_line(data=data.frame(RUN.RoC = predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = test)),
                            RUN.Age.Pos = test$RUN.Age.Pos),
            color="blue", size=1)+
  geom_point(color="black", size=1)+
  geom_point(data = test[test$soft.Peak==T,], color="yellow", size=1)+
  geom_point(data = test[test$Peak==T,], color="orange", size=2)+
  geom_point(data = test[test$Peak.gam==T,], color="red", size=3)+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5))



ggsave(paste0("RUN_",dataset.N,"_",sm.type,"_BIN_",BIN,"_Shift_",Shiftbin,"_DC_",DC,"_R_",rand), device = "png")

test %>% ggplot(aes(x= RUN.Age.Pos))+
  theme_classic()+
  coord_flip(xlim=c(0,8000))+
  scale_x_continuous(trans = "reverse")+
  geom_histogram(fill="gray80", binwidth = 500,color="gray50")

