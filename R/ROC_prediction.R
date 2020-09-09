##################################################
### -------- Predicting ROC values  ---------- ###
##################################################

# cGAM modeling of ROC values

# ------------------------------------------------------------------------------

# ----------------------------------------------
#                     SETUP
# ----------------------------------------------
library(tidyverse)
library(mgcv)
library(gratia)
library(sp)
library(rgdal)
library(maps)
library(RColorBrewer)
theme_set(theme_classic())


countries <- readOGR(dsn = "C:/Users/omo084/OneDrive - University of Bergen/PRIVATE/GITHUB/NEOTOMA-SelectedPollenSites/data/spatial/Countries",
                     layer = "TM_WORLD_BORDERS-0.3")

# ----------------------------------------------
#                     DATA 
# ----------------------------------------------

# load HOPE pollen data
tibble_HOPE <- readRDS("C:/Users/omo084/OneDrive - University of Bergen/HOPE_data/R_Data/_HOPE_DATA/HOPE_data_smooth20200730.RDS")

ROC_HOPE <- readRDS("C:/Users/omo084/OneDrive - University of Bergen/HOPE_data/R_Data/_HOPE_DATA/_ROC/HOPE_Roc20200804.RDS")

HOPE_data = tibble_HOPE %>%
  left_join(ROC_HOPE, by="dataset.id")

HOPE_data_work  <- HOPE_data %>%
  unnest(cols = c(ROC)) %>%
  mutate(BIN = ceiling(AGE/500)) %>%
  mutate(BIN = (BIN-1)*500) %>%
  dplyr::select(REGION, dataset.id,Schulz_Bio, long, lat, BIN, AGE, ROC, PEAK)

HOPE_data_work_EU = HOPE_data_work %>%
  filter(REGION  == "Europe") %>%
  dplyr::select(-REGION)

# Site age distribution
HOPE_data_work_EU %>%
  group_by(dataset.id) %>%
  summarise(age_min=min(BIN),
            age_max = max(BIN)) %>%
  arrange(age_min, age_max) %>%
  mutate(ORDER = 1:nrow(.)) %>%
  ggplot(aes(x=age_min,xend=age_max, y=-ORDER,yend=-ORDER))+
  geom_segment()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(trans = "reverse")+
  labs(x="Age (cal yr BP)",
       y="")


EU_data_BIN = HOPE_data_work_EU %>%
  mutate(LONG = ceiling(long),
         LAT = ceiling(lat)) %>%
  group_by(LONG,LAT, BIN) %>%
  dplyr::summarise(N.samples=n(),
                   ROC.mean =  mean(ROC),
                   ROC.median = median(ROC),
                   ROC.upq = quantile(ROC, 0.95),
                   PEAK.m = mean(PEAK)) %>%
  mutate(PEAK.m = replace(PEAK.m, is.na(PEAK.m),0)) %>%
  filter(BIN < 9e3) 

# Site spatial distribution with ROC
EU_data_BIN %>%
  ggplot(aes(x=LONG, y=LAT))+
  borders(fill="gray90", colour="gray90")+
  geom_point(aes(color=ROC.upq, size=N.samples))+
  facet_wrap(~BIN)+
  coord_quickmap(xlim = c(-12,40), ylim=c(38,70))+
  scale_color_gradient(low="purple",high = "yellow")

# Site spatial distribution with peak points
EU_data_BIN %>%
  ggplot(aes(x=LONG, y=LAT))+
  borders(fill="gray90", colour="gray90")+
  geom_point(aes(color=PEAK.m, size=N.samples))+
  facet_wrap(~BIN)+
  coord_quickmap(xlim = c(-12,40), ylim=c(38,70))+
  scale_color_gradient(low="purple",high = "yellow")

summary(EU_data_BIN)

EU_data_BIN$BIN %>%
  unique() %>% length()


# ----------------------------------------------
#                 MODELING 
# ----------------------------------------------


active_cores = parallel::detectCores()-1

# test BAM without space
EU_bam_test = bam(ROC.upq~s(BIN, k=15),
             data = EU_data_BIN,
             weights = N.samples/mean(N.samples),
             family = Gamma(), #betar(eps=.Machine$double.eps*1e8)
             method = "fREML", 
             nthreads = active_cores,
             discrete = T,
             select = F)

gam.check(EU_bam_test)

summary(EU_bam_test)

appraise(EU_bam_test, method = "simulate")

draw(EU_bam_test, scales = "fixed")


# BAM with spatial interacion
EU_bam = bam(ROC.upq~s(BIN, k=15)+
               s(LONG, LAT, k = 200, bs = "ds")+
               ti(LONG, LAT, BIN, d=c(2,1), bs=c("ds","tp"),k=c(30,10)),
             data = EU_data_BIN,
             weights = N.samples/mean(N.samples),
             family = Gamma(), #betar(eps=.Machine$double.eps*1e8)
             method = "fREML", 
             nthreads = active_cores,
             discrete = T,
             select = F)

gam.check(EU_bam)

summary(EU_bam)

appraise(EU_bam, method = "simulate")

draw(EU_bam, scales = "fixed")

# ----------------------------------------------
#                   PLOT 
# ----------------------------------------------

EU_pred_data <- with(EU_data_BIN,
                     expand_grid(BIN = seq(0,max(EU_data_BIN$BIN),500),
                                 LONG = seq(-12,42,length.out = 100),
                                 LAT = seq(36,72,length.out = 100)))

EU_pred <- predict.bam(EU_bam, EU_pred_data, type = "response")

ind_EU <- exclude.too.far(EU_pred_data$LONG, EU_pred_data$LAT,
                          EU_data_BIN$LONG, EU_data_BIN$LAT, dist = 0.1)
table(ind_EU)

EU_pred[ind_EU] <- NA

EU_pred_data_fin <- cbind(EU_pred_data,fit= EU_pred)



spdf <- SpatialPointsDataFrame(coords = EU_pred_data_fin[, c("LONG", "LAT")], data = EU_pred_data_fin,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs "))

whatever <- spdf[!is.na(over(spdf, as(countries, "SpatialPolygons"))), ]
whatever <- as.data.frame(whatever)

EU_pred_data_fin = whatever %>% as_tibble()

quantile(EU_pred_data_fin$fit, seq(0.2,0.8,0.2), na.rm=T)

EU_pred_data_fin = EU_pred_data_fin %>%
  mutate(Fit.cat = ifelse(fit >0.4983, 5,
                          ifelse(fit>0.4172,4,
                                 ifelse(fit>0.3537,3,
                                        ifelse(fit>0.3038,2,1)))))
summary(EU_pred_data_fin)

EU_pred_data_fin %>%
  #filter(fit >= 0 & fit < 2.6) %>%
  ggplot(aes(x=LONG, y=LAT))+
  borders(fill="gray90", colour = "gray90")+
  geom_raster(aes(fill=as.factor(Fit.cat)))+
  #borders(fill=NA, colour = "gray30")+
  facet_wrap(~BIN)+
  coord_quickmap(xlim = c(-12,42), ylim=c(36,72))+
  scale_fill_manual(values=brewer.pal(5, "YlOrRd"), na.value=NA)+
  labs(y="latitude",
       x="longitude",
       fill="RoC category",
       title="Predicted values for EU")


EU_pred_data_fin %>%
  #filter(fit >= 0 & fit < 2.6) %>%
  ggplot(aes(x=LONG, y=LAT))+
  borders(fill="gray90", colour = "gray90")+
  geom_raster(aes(fill=log(fit+1)))+
  #borders(fill=NA, colour = "gray30")+
  facet_wrap(~BIN)+
  coord_quickmap(xlim = c(-12,42), ylim=c(36,72))+
  scale_fill_gradient(low="yellow",high = "red", na.value=NA)+
  labs(y="latitude",
       x="longitude",
       fill="RoC score",
       title="Predicted values for EU")


ggsave("fig/EU_pred_500.pdf",
       width = 40,height = 30, units = "cm")



