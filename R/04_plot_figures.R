#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                    Graphical outcomes
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 1. Load model and  data  -----
#----------------------------------------------------------#

# datasets
data_correct <- 
  read_rds("data/output/datasets/subsets/data_correct.rds")

data_false <- 
  read_rds("data/output/datasets/subsets/data_false.rds")

data_detail_correct <- 
  read_rds("data/output/datasets/subsets/data_detail_correct.rds")

data_detail_false <- 
  read_rds("data/output/datasets/subsets/data_detail_false.rds")

# models
mod_correct_select <-
  read_rds("data/output/models/mod_correct_select.rds")

mod_false_select <-
  read_rds("data/output/models/mod_false_select.rds")

mod_detail_correct_select <-
  read_rds("data/output/models/mod_detail_correct_select.rds")

mod_detail_false_select <-
  read_rds("data/output/models/mod_detail_false_select.rds")


#----------------------------------------------------------#
# 2. calculate EMM  -----
#----------------------------------------------------------#

# 2.1 for all WU and Peaks ----- 
emmeans_correct_full <-
  emmeans(
    mod_correct_select,
    pairwise ~ WU * Peak,
    type = "response") 

emmeans_correct_full_tibble <- 
  tibble(
    emmeans_correct_full$emmeans %>% 
      as_tibble(),
    segment = "correct detection") 

emmeans_false_full <-
  emmeans(
    mod_false_select,
    pairwise ~ WU * Peak,
    type = "response") 

emmeans_false_full_tibble <- 
  tibble(
    emmeans_false_full$emmeans %>% 
      as_tibble(),
    segment = "false positives") 

emmeans_correct_WU <- 
  emmeans(
    mod_correct_select,
    ~ WU,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "correct detection")

emmeans_false_WU <- 
  emmeans(
    mod_false_select,
    ~ WU,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "false positives")

emmeans_correct_peak <- 
  emmeans(
    mod_correct_select,
    ~ Peak,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "correct detection")
  
emmeans_false_peak <- 
  emmeans(
    mod_false_select,
    ~ Peak,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "false positives")


# 2.2 detail datset type and RoC setting ----- 

emmeans_detail_correct_full <-
  emmeans(
    mod_detail_correct_select,
    pairwise ~  DC + position +  smooth + 
      DC:position + position:smooth,
    type = "response") 

emmeans_detail_false_full <-
  emmeans(
    mod_detail_false_select,
    pairwise ~ DC + position +  smooth + 
      DC:position + DC:smooth + position:smooth +
      DC:position:smooth,
    type = "response") 

emmeans_detail_correct_full_tibble <- 
  tibble(
    emmeans_detail_correct_full$emmeans %>% 
      as_tibble(),
    segment = "correct detection") 

emmeans_detail_false_full_tibble <- 
  tibble(
    emmeans_detail_false_full$emmeans %>% 
      as_tibble(),
    segment = "false positives")

#----------------------------------------------------------#
# 3. Plot figures  -----
#----------------------------------------------------------#

#---------------------------------------------#
# 3.1 Peak points & WU (Fig. 2) -----
#---------------------------------------------#

(plot_full_WU <-
  bind_rows(
    emmeans_correct_WU,
    emmeans_false_WU) %>% 
  ggplot(
    aes(
      y = response,
      x = WU,
      col = segment,
      fill = segment)) + 
  
  geom_bar(
    col = "gray30",
    stat = "identity",
    width = 0.75,
    position = position_dodge(width = 0.75)) + 
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL),
    col = "gray30",
    position = position_dodge(width = 0.75),
    width = 0.2) +
  
  geom_point(
    shape = 0,
    col = "gray30",
    position = position_dodge(width = 0.75),
    size = 1) +
  
  #scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_segment) +
  scale_fill_manual(values = color_legen_segment) +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "right",
    text = element_text(size = text_size))+
   labs(
     y = "proportion of detected Working Units",
     x = "Working Unit type",
     color = "",
     fill = ""
   ))

(plot_full_peak <-
    bind_rows(
      emmeans_correct_peak,
      emmeans_false_peak) %>% 
    ggplot(
      aes(
        y = response,
        x = Peak,
        col = segment,
        fill = segment)) + 
    
    geom_bar(
      col = "gray30",
      stat = "identity",
      width = 0.75,
      position = position_dodge(width = 0.75)) + 
    
    geom_errorbar(
      aes(
        ymin = lower.CL,
        ymax = upper.CL),
      col = "gray30",
      position = position_dodge(width = 0.75),
      width = 0.2) +
    
    geom_point(
      shape = 0,
      col = "gray30",
      position = position_dodge(width = 0.75),
      size = 1) +
 
    #scale_y_continuous(limits = c(0,1))+
    scale_color_manual(values = color_legen_segment) +
    scale_fill_manual(values = color_legen_segment) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      legend.position = "right",
      text = element_text(size = text_size))+
    labs(
      y = "proportion of detected Working Units",
      x = "Peak point detection method",
      color = "",
      fill = ""))



bind_rows(
  emmeans_correct_full_tibble,
  emmeans_false_full_tibble
  ) %>%  
  ggplot(
    aes(
      y = response,
      x = Peak,
      col = segment,
      fill = segment)) + 
  
  geom_bar(
    col = "gray30",
    stat = "identity",
    width = 0.75,
    position = position_dodge(width = 0.75)) + 
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL),
    col = "gray30",
    position = position_dodge(width = 0.75),
    width = 0.2) +
  
  geom_point(
    shape = 0,
    col = "gray30",
    position = position_dodge(width = 0.75),
    size = 1) +
  
  facet_wrap( ~ WU, nrow = 1) +
  
  #scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_segment) +
  scale_fill_manual(values = color_legen_segment) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -90),
    text = element_text(size = text_size))


bind_rows(
  emmeans_detail_correct_full_tibble,
  emmeans_detail_false_full_tibble) %>%  
  ggplot(
    aes(
      y = response,
      x = position,
      col = smooth,
      fill = smooth)) + 
  
  geom_bar(
    col = "gray30",
    stat = "identity",
    width = 0.75,
    position = position_dodge(width = 0.75)) + 
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL),
    col = "gray30",
    position = position_dodge(width = 0.75),
    width = 0.2) +
  
  geom_point(
    shape = 0,
    col = "gray30",
    position = position_dodge(width = 0.75),
    size = 1) +
  
  facet_grid(segment ~ DC,
             scales = "free_y") +
  
  #scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_smooth) +
  scale_fill_manual(values = color_legen_smooth) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -90),
    text = element_text(size = text_size))





