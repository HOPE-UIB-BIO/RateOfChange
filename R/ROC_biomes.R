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
library(ggpubr)
library(raster)
library(svMisc)
theme_set(theme_classic())


WWF <- readOGR(dsn = "C:/Users/omo084/Documents/GITHUB/NEOTOMA-SelectedPollenSites/data/spatial/WWF",
                     layer = "wwf_terr_biomes")


tibble_HOPE <- readRDS("C:/Users/omo084/OneDrive - University of Bergen/HOPE_data/R_Data/_HOPE_DATA/HOPE_data_smooth20200730.RDS")
ROC_HOPE <- readRDS("C:/Users/omo084/OneDrive - University of Bergen/HOPE_data/R_Data/_HOPE_DATA/_ROC/HOPE_Roc20200804.RDS")

FINAL_data_smooth_READY = tibble_HOPE %>%
  inner_join(ROC_HOPE, by ="dataset.id")

add_region <- function(x, region = regions){
  x2 <- x # you need a copy to be able to keep the data and cbind in the last line
  coordinates(x2) <- ~long+lat
  proj4string(x2) <- (CRS(proj4string(region))) # assumes projection are the same
  cbind(x, over(x = x2, y = region)) 
  
}

FINAL_data_smooth_READY = add_region(FINAL_data_smooth_READY, region = WWF) %>%
  #dplyr::select(-Id) %>%
  as_tibble()

FINAL_data_smooth_READY = FINAL_data_smooth_READY %>%
  filter(is.na(BIOM_NAME)==F)

# remove weird islands
sites_to_remove <- c(41532,1697)

FINAL_data_smooth_READY <- FINAL_data_smooth_READY %>%
  filter(!(dataset.id %in% sites_to_remove))

age.treshold = 12e3 
Roc.treshold = 2 

# save the result from the k.check fc
f <- function(b, k.sample = 10e3, k.rep = 1e3) {
  mgcv:::k.check(b, subsample = k.sample, n.rep = k.rep)
}

select_model <- function(var_y,var_x, family, data, weights=NA) {
  for(i in 3:(min(c(nrow(data)-1,100)))){
    print(paste("trying k=",i))
    formula.w <- paste(var_y,"~ s(",var_x,", k=",i,")")
    
    if (is.na(weights)==F){
      data = data %>%
        mutate(W= with(data,get(weights))) 
      data$W = data$W / mean(data$W)
      
    } else {
      data = data %>%
        mutate(W=rep(1,nrow(data)))
    }
    
    suppressWarnings(gam.w <- gam(formula = as.formula(formula.w), data = data,family =  noquote(family),
                                  weights = W,
                                  method = "REML", niterPQL=50))
    
    
    suppressWarnings(basis <- f(gam.w))
    
    if(is.na(basis[4])==F){
      if(basis[4] > 0.05) {break}  
    }
    
  }
  return(gam.w)
} 

predict_gam_values <- function(gam_model, var_x, deriv = F){
  
  new_data  <-  tibble(VAR = seq(gam_model$var.summary[[1]][1],gam_model$var.summary[[1]][3], length.out =100))
  names(new_data) = var_x
  
  new_data_pred <- cbind(new_data,
                         data.frame(predict(gam_model, se.fit = TRUE,newdata = new_data, type="response"))) %>%
    as_tibble()
  
  crit.t <- qt(0.975, df = df.residual(gam_model))
  
  new_data_pred <- new_data_pred %>%
    mutate(upper = fit + (crit.t * se.fit),
           lower = fit - (crit.t * se.fit))
  
  k_value = gam_model$formula %>%
    as.character(.) %>%
    .[3] %>%
    sub(".*=","",.) %>%
    sub(")","",.) %>%
    as.double(.)
  
  gam_model_summary = gam_model %>%
    summary()
  
  p_value = gam_model_summary$s.table[4]
  
  if (deriv == T & p_value <0.05){
    gam_model_deriv <- fderiv(gam_model, newdata = new_data, n=1000)
    
    gam_model_deriv_sint <- with(new_data, cbind(confint(gam_model_deriv, nsim = 1000,type="simultaneous", transform	= T),
                                                 BIN=BIN)) %>% as_tibble()
    Deriv = gam_model_deriv_sint %>%
      mutate(significante_change = lower > 0 | upper < 0 ) %>% 
      dplyr::select(BIN,significante_change)
  } else {Deriv = NA}
  
  return(list(data=new_data_pred,k=k_value,p=p_value, FirstDeriv = Deriv))
}

draw_gam_full <- function(data, pred_gam_up, pred_gam_pk, region, siluete=F , blue_samples = F, palette_x = pallete.1, deriv = F, y_cut = 1.8){
  
  p0 <- ggplot()+
    scale_x_continuous(trans = "reverse", breaks = seq(0,20e3, 2e3))+
    coord_flip(xlim=c(age.treshold,0))+
    theme_classic()+
    theme(legend.position = "none")+
    labs(x="Age (yr BP)")
  
  if (deriv == T){
    
    # ROC 
    
    p3 <- p0 +
      geom_vline(xintercept = seq(0,20e3, 2e3), color= "gray90", size = 0.1)+
      geom_ribbon(data = pred_gam_up$data, aes(x=BIN, y=fit, ymin=lower, ymax=upper),fill="gray80", alpha=ifelse(pred_gam_up$p<0.05,0.5,0))+
      geom_line(data = pred_gam_up$data, aes(x=BIN, y=fit), lty=ifelse(pred_gam_up$p<0.05,1,2))+
      scale_fill_manual(values = palette_x)+
      scale_color_manual(values = palette_x)+
      labs(y="95% quantile ROC score / proportion of peak points")+
      scale_y_continuous(limits = c(-0.1,y_cut), breaks = seq(0,2,0.2),  sec.axis =  sec_axis(name =  ,~ ./2, breaks=seq(0,2,0.2)) )
    
    
    if (is.na(pred_gam_up$FirstDeriv)==F){
      rect.df_upq <-  left_join(pred_gam_up$FirstDeriv,pred_gam_up$data, by = "BIN") %>%
        filter(significante_change==T)
      
      if(nrow(rect.df_upq) > 0){
        p3 <- p3 + geom_point(data= rect.df_upq, aes(x=BIN, y=fit), size=2, color="gray30", shape=8)
      }  
    }
    
    p3 <- p3 + geom_point(data=data,aes(x=BIN, y=ROC.upq,color=region), shape=15, size = 2)
    
    # peak point
    
    p3 <- p3 +geom_ribbon(data = pred_gam_pk$data, aes(x=BIN, y=fit*2, ymin=lower*2, ymax=upper*2), fill="gray80",alpha=ifelse(pred_gam_pk$p<0.05,0.5,0))+
      geom_line(data = pred_gam_pk$data, aes(x=BIN, y=fit*2),lty=ifelse(pred_gam_pk$p<0.05,3,2))
    
    
    if (is.na(pred_gam_pk$FirstDeriv)==F){
      rect.df_pk <- left_join(pred_gam_pk$FirstDeriv,pred_gam_pk$data, by = "BIN") %>%
        filter(significante_change==T)
      
      if(nrow(rect.df_pk) > 0){
        p3 <- p3 + geom_point(data= rect.df_pk, aes(x=BIN, y=fit*2), size=2, color="gray30", shape=8)
      }
    }
    
    p3 <- p3 +  geom_point(data=data,aes(x=BIN, y=VALUE.mean*2,color=region), shape=0, size = 2)
    
    
  } else {
    p3 <- p0 +
      geom_ribbon(data = pred_gam_up$data, aes(x=BIN, y=fit, ymin=lower, ymax=upper),fill="gray80", alpha=ifelse(pred_gam_up$p<0.05,0.5,0))+
      geom_point(data=data,aes(x=BIN, y=ROC.upq,color=region), shape=15)+
      geom_line(data = pred_gam_up$data, aes(x=BIN, y=fit), lty=ifelse(pred_gam_up$p<0.05,1,2))+
      geom_ribbon(data = pred_gam_pk$data, aes(x=BIN, y=fit*2, ymin=lower*2, ymax=upper*2), fill="gray80",alpha=ifelse(pred_gam_pk$p<0.05,0.5,0))+
      geom_point(data=data,aes(x=BIN, y=VALUE.mean*2,color=region), shape=0)+
      geom_line(data = pred_gam_pk$data, aes(x=BIN, y=fit*2),lty=ifelse(pred_gam_pk$p<0.05,3,2))+
      scale_color_manual(values = palette_x)+
      labs(y="95% quantile ROC score / proportion of peak points")+
      scale_y_continuous(limits = c(-0.1,y_cut),breaks = seq(0,2,0.2),  sec.axis =  sec_axis(name =  ,~ ./2, breaks=seq(0,2,0.2))  )  
  }
  
  if (siluete == T){
    region_coord_w <- region_coord %>%
      filter(REGION == region)
    
    p2<- ggplot()+
      borders(fill=pallete.1[names(pallete.1)==region], colour=NA, alpha=.3)+
      coord_fixed(xlim=c(region_coord_w$long_min,region_coord_w$long_max),
                  ylim = c(region_coord_w$lat_min, region_coord_w$lat_max))+
      theme_transparent()
    
    p2_g <- ggplotGrob(p2) 
    
    p3_a <- p3+rremove("xylab") + annotation_custom(grob = p2_g, xmin=-0, xmax=-5e3, ymin=0.9, ymax=y_cut) 
  } else {
    p3_a <- p3+rremove("xylab")
  }
  
  if (blue_samples == T){
    Max_samples = max(data$N.samples)
    
    p1<- ggplot()+
      coord_flip(xlim=c(age.treshold,0))+
      theme_classic()+
      labs(x="Age (yr BP)")+
      geom_ribbon(data=data, aes(x=BIN,ymax=N.samples), ymin=0, fill="lightblue1")+
      scale_x_continuous(trans = "reverse", position = "top")+
      scale_y_continuous(trans = "reverse", position = "right", breaks = c(0,Max_samples))+
      theme(axis.title.x = element_blank())+
      labs(y="Number of samples in each time bin")  
    
    p4 <- ggarrange(p3_a,p1+rremove("y.text")+rremove("y.title")+rremove("y.ticks"), align = "h", nrow = 1, widths = c(1,0.2))
    
  } else {
    p4 <- p3_a + theme(text = element_text(size = 16))
  }
  
  return(p4)
}

draw_region <- function(region, plot_data ,BIN_data){
  
  data.w <- plot_data %>%
    filter(REGION == region) %>%
    ungroup()
  
  data.w.sum <- BIN_data %>%
    filter(REGION == region) %>%
    ungroup()
  
  # MEDIAN ROC GAM
  #gam.median <- select_model(var_y ="ROC.median",var_x="BIN", family="Gamma", data = data.w.sum, weights = "N.samples")
  #pred_gam.median<-predict_gam_values(gam.median, "BIN")
  
  # 95% quantile ROC GAM
  gam.upq <- select_model(var_y ="ROC.upq",var_x="BIN",family="Gamma()", data = data.w.sum, weights="N.samples")
  
  pred_gam.upq<-predict_gam_values(gam_model = gam.upq, var_x= "BIN", deriv = T)
  
  # PEAK POINTS GAM
  gam.PEAK <- select_model(var_y ="VALUE.mean",var_x="BIN",family="betar(eps=.Machine$double.eps*1e8)", data = data.w.sum,weights="N.samples")
  
  pred_gam.PEAK <-predict_gam_values(gam_model = gam.PEAK,  var_x= "BIN", deriv = T)
  
  p_fin = draw_gam_full(data = data.w.sum, 
                        pred_gam_up = pred_gam.upq, 
                        pred_gam_pk = pred_gam.PEAK, 
                        region = region, 
                        blue_samples = F,
                        siluete=T, 
                        palette_x = pallete.1, 
                        deriv=T,
                        y_cut = 1.3)
  
  return(p_fin)
}

draw_cluster <- function(region){
  
  data.w <- FINAL_data_smooth_READY %>%
    #filter(REGION == region) %>%
    arrange(dataset.id)
  
  data.w$BIOM_NAME  <- with(data.w, reorder(BIOM_NAME ,-lat, mean,))
  
  getPalette = colorRampPalette(brewer.pal(min(8, length(levels(data.w$BIOM_NAME))), "Set2" ))
  Palette.w<- getPalette(length(levels(data.w$BIOM_NAME)))
  names(Palette.w)<- levels(data.w$BIOM_NAME)
  
  long_lim  = c(min(data.w$long), max(data.w$long))
  lat_lim  = c(min(data.w$lat), max(data.w$lat))
  
  p1 <- data.w %>%
    ggplot(aes(x=long, y=lat))+
    borders(fill = "gray90", colour = "gray90") +
    coord_quickmap(xlim=long_lim, ylim=lat_lim)+
    geom_point(aes(color=BIOM_NAME), size=3)+
    #geom_rug(aes(x=long), sides = "b")+
    #geom_rug(aes(x=lat), sides = "l")+
    theme_classic()+
    scale_color_manual(values = Palette.w)+
    labs(y= "latitude",
         x= "longitude")+
    theme(legend.position = "bottom")  
  p1
  
  
  data.w.ROC <- data.w %>%
    unnest(c(ROC)) %>%
    mutate(BIN =  ceiling(AGE/500)*500) %>%
    dplyr::select(REGION,BIOM_NAME, BIN, ROC) %>%
    group_by(REGION,BIOM_NAME, BIN) %>%
    summarise(.groups="drop",
              N= n(),
              ROC.m = median(ROC),
              ROC.sd = sd(ROC),
              ROC.up = quantile(ROC,0.975),
              ROC.dw = quantile(ROC,0.025),
              ROC.SE = ROC.sd / sqrt(N)) %>%
    left_join(.,data.w %>%
                arrange(dataset.id) %>%
                unnest(c(ROC)) %>%
                mutate(BIN =  ceiling(AGE/500)*500 ) %>%
                dplyr::select(REGION,BIOM_NAME,dataset.id, BIN, PEAK) %>%
                ungroup() %>%
                group_by(REGION,BIOM_NAME,dataset.id, BIN) %>%
                dplyr::summarise(.groups="keep", 
                                 PEAK.m = mean(PEAK)) %>%
                group_by(REGION,BIOM_NAME, BIN) %>%
                dplyr::summarise(.groups="keep",
                                 N = n(),
                                 VALUE.mean = mean(PEAK.m),
                                 VALUE.sd = sd(PEAK.m),
                                 VALUE.SE = VALUE.sd / sqrt(N)),
              by=c("REGION","BIOM_NAME","BIN")) %>%
    mutate(VALUE.mean = replace(VALUE.mean, is.na(VALUE.mean),0),
           VALUE.sd = replace(VALUE.sd, is.na(VALUE.sd),0))
  
  cluster_result = data.w %>%
    group_by(BIOM_NAME) %>%
    summarise(.groups = "drop",
              N=n()) %>%
    #filter(N >= 10) %>%
    mutate(C_char = as.character(BIOM_NAME)) %>%
    mutate(data =  purrr::map(C_char, .f=function(x){
      data.w.ROC.c <- data.w.ROC  %>%
        filter(BIOM_NAME == as.character(x)) %>%
        mutate(N.samples = N.x,
               ROC.upq = ROC.up)
      return(data.w.ROC.c)
    })) %>%
    mutate(gam_ROC_upq = purrr::map(data, .f=function(x){
      gam_ROC_upq <- select_model(var_y = "ROC.up",var_x= "BIN", family = "Gamma()", data = x,weights = "N.samples")
      return(gam_ROC_upq)
    })) %>%
    mutate(pred_gam_ROC_upq = purrr::map(gam_ROC_upq, .f=function(x){
      pred_gam_ROC_upq <-predict_gam_values(x, "BIN", deriv = T)
      return(pred_gam_ROC_upq)
    })) %>%
    mutate(gam_ROC_PEAK = purrr::map(data, .f=function(x){
      gam_ROC_PEAK <- select_model(var_y = "VALUE.mean",var_x= "BIN", family = "betar(eps=.Machine$double.eps*1e8)", data = x,weights = "N.samples")
      return(gam_ROC_PEAK)
    })) %>%
    mutate(pred_gam_ROC_PEAK = purrr::map(gam_ROC_PEAK, .f=function(x){
      pred_gam_ROC_PEAK <-predict_gam_values(x, "BIN", deriv = T)    
      return(pred_gam_ROC_PEAK)
    })) %>%
    mutate(Plot = purrr::pmap(list(data, pred_gam_ROC_upq,pred_gam_ROC_PEAK,C_char), .f=function(x,y,z,u){
      plot.res = draw_gam_full(x,y,z,region =as.character(u), siluete = F, palette_x = Palette.w, deriv = T, y_cut=1.75, blue_samples = T)
      return(plot.res)
    }))
  
  p2 <- ggarrange(plotlist = cluster_result$Plot) %>%
    annotate_figure(., bottom = "95% quantile of RoC", top="percentage of peak points")
  
  p.fin <- ggarrange(p1,p2,nrow=1, widths = c(1,.75))
  
  return(list(data_sites = data.w, data_c = cluster_result, plot = annotate_figure(p.fin, top=paste(region,",N =",nrow(data.w)))))
}

re_Asia = draw_cluster("Asia")

re_LA = draw_cluster("Latin America")
