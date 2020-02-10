# res.df.plot <- read.csv("results20201002.csv", row.names = 1)
# load("~/HOPE/Data/tibble_Europe_filtered29.01.20.Rdata")

# ----------------------------------------------
#             RESULT VISUALISATION 
# ----------------------------------------------


RoC_summary_p1 <- res.df.plot[,] %>%
  ggplot(aes( y=RoC.mean, 
              x= DF.Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,20000), ylim = c(0,5))+
  geom_line(aes(group=as.factor(ID)),alpha=1/10, size=1)+
  geom_hline(yintercept = 0, color="red")+
  geom_smooth(color="green", method = "loess", se=F)+
  xlab("Age")+ylab("Rate of Change")
RoC_summary_p1

RoC_summary_p2 <- res.df.plot[,] %>%
  ggplot(aes( y=RoC.mean, 
              x= DF.Age))+
  theme_classic()+
  scale_color_gradient2(low="white",mid="darkblue",high="black", midpoint = 4)+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,20000), ylim = c(0,5))+
  geom_point(data = res.df.plot[res.df.plot$Peak==T,],aes(color=RoC.mean), alpha=1/10, size=3)+
  geom_hline(yintercept = 0, color="red")+
  geom_smooth(color="green", method = "loess", se=F)+
  xlab("Age")+ylab("Rate of Change")+
  theme(legend.position = "none")
RoC_summary_p2

RoC_summary_p3 <- res.df.plot[res.df.plot$Peak==T,] %>%
  ggplot(aes( x= DF.Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,20000))+
  geom_density(fill="gray")+
  xlab("Age")+ylab("Density of Peak-points")
RoC_summary_p3


RoC_summary<- ggarrange(RoC_summary_p1,RoC_summary_p2,RoC_summary_p3, ncol = 3, labels = c("A","B","C"))
RoC_summary

ggsave("RoC_summary.pdf")


p0<-data.frame(POLLEN=reshape2::melt(tibble_Europe2[2,]$filtered.counts[[1]]), 
           AGE=rep(tibble_Europe2[2,]$list_ages[[1]]$ages$age, ncol(tibble_Europe2[2,]$filtered.counts[[1]]))) %>%
  ggplot(aes( y=POLLEN.value, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=rep(0,12000),ymax=POLLEN.value), alpha=1)+
  geom_line(aes(group=POLLEN.variable), alpha=1/10)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  facet_wrap(~POLLEN.variable)+
  xlab("Age")+ylab("Pollen")+
  coord_flip()
p0

ggsave("ExamplePlot01.pdf",plot= p0, width = 50, height = 30, units= "cm", dpi= 600)

p1a <- res.df.plot[res.df.plot$ID=="17334",] %>%
  ggplot(aes( y=RoC.mean, 
              x= DF.Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,15000), ylim = c(0,1.5))+
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/2)+
  geom_line(aes(group=as.factor(ID)),alpha=1, size=1)+
  geom_point(data = res.df.plot[res.df.plot$Peak==T & res.df.plot$ID=="17334",],color="blue", alpha=1, size=3)+
  geom_hline(yintercept = median(res.df.plot[res.df.plot$ID=="17334",]$RoC.mean), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")
p1a

p1b <- res.df.plot[res.df.plot$Peak==T & res.df.plot$ID=="17334",] %>%
  ggplot(aes(x=DF.Age))+
  theme_classic()+
  geom_density(fill="gray")+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,15000))+
  xlab("Age")+ylab("Density of Peak-Points")
p1b


p1<- ggarrange(p1a,p1b, ncol = 2, labels = c("A","B"))
p1

ggsave("ExamplePlot02.pdf",p1, width = 20, height = 15, units= "cm", dpi= 600)

p2<- res.df.plot %>%
  ggplot(aes( y=RoC.mean, 
              x= DF.Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(0,20000), ylim = c(0,5))+
  geom_line(aes(group=as.factor(ID)),alpha=1, size=1)+
  geom_point(data = res.df.plot[res.df.plot$Peak==T,],color="blue", alpha=1, size=1)+
  geom_hline(yintercept = 0, color="red")+
  facet_wrap(~ID)+
  xlab("Age")+ylab("Rate of Change")
p2

ggsave("PerSitePlot.pdf",p2, width = 50, height = 30, units= "cm", dpi= 600)


# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())
