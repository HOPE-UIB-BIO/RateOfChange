# res.df.plot <- read.csv("results20202014.csv", row.names = 1)
# load("C:/Users/ondre/Dropbox/HOPE_data/tibble_Europe_filtered13.02.20.RData")

# ----------------------------------------------
#             RESULT VISUALISATION 
# ----------------------------------------------


RoC_summary_p1 <- res.df.plot %>%
  ggplot(aes( y=RoC.median, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,30000), ylim = c(0,5))+
  geom_line(aes(group=as.factor(ID)),alpha=1/10, size=1)+
  geom_hline(yintercept = 0, color="red")+
  geom_smooth(color="green", method = "loess", se=F)+
  xlab("Age")+ylab("Rate of Change")
RoC_summary_p1


RoC_summary_p1b <-res.df.plot %>%
  ggplot(aes(x=age))+
  geom_density(fill="gray")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,30000))+
  xlab("")+ylab("Density of the samples")
RoC_summary_p1b


res.df.plot %>% 
  ggplot(aes(x=age))+
  geom_histogram(fill="gray",color="gray80", binwidth = 200)+
  theme_classic()+
  scale_x_continuous()+ #trans = "reverse"
  coord_cartesian (xlim=c(0,30000))+
  xlab("age")+ylab("Number of the samples")+
  ggtitle("binwidth = 200")


RoC_summary_p2 <- res.df.plot %>%
  filter(Peak==T) %>%
  ggplot(aes( y=RoC.median, 
              x= age))+
  theme_classic()+
  scale_color_gradient2(low="white",mid="darkblue",high="black", midpoint = 4)+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,30000), ylim = c(0,5))+
  geom_point(aes(color=RoC.median), alpha=1/10, size=3)+
  geom_hline(yintercept = 0, color="red")+
  geom_smooth(color="orange", method = "loess", se=F)+
  xlab("")+ylab("Rate of Change")+
  theme(legend.position = "none")
RoC_summary_p2

RoC_summary_p2b <- res.df.plot %>%
  filter (Peak==T) %>%
  ggplot(aes( x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,30000))+
  geom_density(fill="gray")+
  xlab("")+ylab("Density of Peak-points")
RoC_summary_p2b


RoC_summary<- ggarrange(RoC_summary_p1,RoC_summary_p1b,
                        RoC_summary_p2,RoC_summary_p2b,
                        ncol = 4, labels = c("A","B","C","D"))
RoC_summary

ggsave("RoC_summary.pdf")

which(tibble_Europe2$dataset.id %in%  22988)

dataset.N <- 130
max.age <- max(tibble_Europe2$list_ages[[dataset.N]]$ages$age)

N.taxa <- 10

Common.list <- tibble_Europe2$filtered.counts[[dataset.N]] %>%
  colSums() %>% 
  sort(decreasing = T) %>%
  .subset(.,1:N.taxa) %>%
  names()
  
p0<-data.frame(POLLEN=tibble_Europe2$filtered.counts[[dataset.N]] %>%
                 select(.,Common.list) %>%
                 reshape2::melt() , 
               AGE=rep(tibble_Europe2$list_ages[[dataset.N]]$ages$age, N.taxa)) %>%
  ggplot(aes( y=POLLEN.value, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=rep(0,length(POLLEN.value)),ymax=POLLEN.value), alpha=1)+
  geom_line(aes(group=POLLEN.variable), alpha=1/10)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  facet_wrap(~POLLEN.variable)+
  geom_vline(xintercept = 520, color="red")+
  xlab("Age")+ylab("Pollen")+
  coord_flip(xlim=c(0,max.age))
p0

ggsave("ExamplePlot01.pdf",plot= p0, width = 50, height = 30, units= "cm", dpi= 600)

p1a <-   res.df.plot %>%
  filter(ID==unique(res.df.plot$ID)[dataset.N]) %>%
  ggplot(aes( y=RoC.median, 
              x= age)) +
  theme_classic() +
  scale_x_continuous(trans = "reverse") +
  coord_flip(xlim=c(0,max.age), ylim = c(0,5)) +
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/5) +
  geom_line(aes(group=as.factor(ID)),alpha=1, size=2) +
  geom_point(data = . %>% filter(Peak==T),color="blue", alpha=1, size=3) +
  geom_hline(yintercept = c( 
    res.df.plot %>%
      filter(ID==unique(res.df.plot$ID)[dataset.N]) %>%
      select("RoC.median") %>%
      apply(.,2, FUN = median)
  ), color="blue") +
  geom_vline(xintercept = 520, color="red")+
  geom_hline(yintercept = 0, color="red") +
  xlab("Age")+ylab("Rate of Change")
p1a


p1b <- res.df.plot %>%
  filter(ID==unique(res.df.plot$ID)[dataset.N]) %>%
  filter(Peak==T)%>%
  ggplot(aes(x=age))+
  theme_classic()+
  geom_density(fill="gray")+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,max.age))+
  xlab("Age")+ylab("Density of Peak-Points")
p1b


p1<- ggarrange(p1a,p1b, ncol = 2, labels = c("A","B"))
p1

ggsave("ExamplePlot02.pdf",p1, width = 20, height = 15, units= "cm", dpi= 600)

p2<- res.df.plot %>%
  ggplot(aes( y=RoC.median, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,30000), ylim = c(0,5))+
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/5)+
  geom_line(aes(group=as.factor(ID)),alpha=1, size=1)+
  geom_point(data = res.df.plot[res.df.plot$Peak==T,],color="blue", alpha=1, size=1)+
  geom_hline(yintercept = 0, color="red")+
  facet_wrap(~ID)+
  xlab("Age")+ylab("Rate of Change")
p2

ggsave("PerSitePlot.pdf",p2, width = 50, height = 30, units= "cm", dpi= 600)



p1a <-   res.df.plot %>%
  filter(ID==unique(res.df.plot$ID)[dataset.N]) %>%
  ggplot(aes( y=RoC.median, 
              x= age)) +
  theme_classic() +
  scale_x_continuous(trans = "reverse") +
  coord_flip(xlim=c(0,max.age), ylim = c(0,5)) +
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/5) +
  geom_line(aes(group=as.factor(ID)),alpha=1, size=2) +
  geom_point(data = . %>% filter(Peak==T),color="blue", alpha=1, size=3) +
  geom_hline(yintercept = c( 
    res.df.plot %>%
      filter(ID==unique(res.df.plot$ID)[dataset.N]) %>%
      select("RoC.median") %>%
      apply(.,2, FUN = median)
  ), color="blue") +
  geom_vline(xintercept = 520, color="red")+
  geom_hline(yintercept = 0, color="red") +
  xlab("Age")+ylab("Rate of Change")
p1a


tibble_Europe2$list_ages[[dataset.N]]$age_position %>%
  scale(center=T, scale = T) %>%
  t() %>%
  as.data.frame() %>%
  mutate(
    q_95 = apply(.,1,FUN=function(x) quantile(x,c(0.975))),
    q_05 = apply(.,1,FUN=function(x) quantile(x,c(0.025))),
    q_83 = apply(.,1,FUN=function(x) quantile(x,c(0.8337))),
    q_17 = apply(.,1,FUN=function(x) quantile(x,c(0.1663)))
  ) %>%
  select(.,c("q_05","q_95","q_83","q_17")) %>%
  bind_cols(.,tibble_Europe2$list_ages[[dataset.N]]$ages) %>%
  ggplot(.,aes(x= age))+
    geom_ribbon(aes(ymax=q_95,ymin=q_05), alpha=1/5)+
  geom_ribbon(aes(ymax=q_83,ymin=q_17), alpha=1/5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip()+
  xlab("Age")+ylab("Age diference (SD)")
  


tibble_Europe2$list_ages[[dataset.N]]$age_quantile %>%
  t() %>%
  as.data.frame() %>%
  bind_cols(.,tibble_Europe2$list_ages[[dataset.N]]$ages ) %>%
  dplyr::rename(.,q_95='95%',q_05="5%") %>%
  ggplot(aes(
    x=depth,
    y=age
  ))+
  geom_ribbon(aes(ymax=q_95,ymin=q_05), alpha=1/5)+
  geom_line()+
  theme_classic()



# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())
