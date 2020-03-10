dataset.N <- 2
sm.type = "age.w" 
N.points = 5
range.age.max = 500 
grim.N.max = 9
interest.treshold = 8000
Debug = F


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
              color="gray20", alpha=1/5, size=0.1)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  #facet_wrap(~POLLEN.variable, ncol = N.taxa)+
  xlab("Age")+ylab("Pollen")+
  coord_flip(xlim=c(0,interest.treshold), ylim = c(0,1))+
  ggtitle(sm.type)

ggsave(paste0("RUN_",dataset.N,"_",sm.type,"_Pollen"), device = "png")
