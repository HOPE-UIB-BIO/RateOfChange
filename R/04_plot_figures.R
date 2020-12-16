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
data_success_focus <- 
  read_rds("data/output/datasets/data_success_focus.rds")

data_success_empty <- 
  read_rds("data/output/datasets/data_success_empty.rds")

data_detail_focus <- 
  read_rds("data/output/datasets/data_detail_focus.rds")

data_detail_empty <- 
  read_rds("data/output/datasets/data_detail_empty.rds")

# models
mod_success_focus_final <-
  read_rds("data/output/models/temp_model_success_focus_select.rds")

mod_success_empty_final <-
  read_rds("data/output/models/temp_model_success_empty_select.rds")

mod_detail_focus_final <-
  read_rds("data/output/models/temp_model_detail_focus_select.rds")

mod_detail_empty_final <-
  read_rds("data/output/models/temp_model_detail_empty_select.rds")


#----------------------------------------------------------#
# 2. calculate EMM  -----
#----------------------------------------------------------#

emmeans_success_focus_full <-
  emmeans(
    mod_success_focus_final,
    pairwise ~ WU * Peak  * dataset_type,
    type = "response") 

emmeans_success_focus_full_tibble <- 
  tibble(
    emmeans_success_focus_full$emmeans %>% 
      as_tibble(),
    segment = "correct detection") 

emmeans_success_empty_full <-
  emmeans(
    mod_success_empty_final,
    pairwise ~ WU * Peak  * dataset_type,
    type = "response") 

emmeans_success_empty_full_tibble <- 
  tibble(
    emmeans_success_empty_full$emmeans %>% 
      as_tibble(),
    segment = "false positives") 

emmeans_success_focus_WU <- 
  emmeans(
    mod_success_focus_final,
    ~ WU,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "correct detection")

emmeans_success_empty_WU <- 
  emmeans(
    mod_success_empty_final,
    ~ WU,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "false positives")

emmeans_success_focus_peak <- 
  emmeans(
    mod_success_focus_final,
    ~ Peak,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "correct detection")
  
emmeans_success_empty_peak <- 
  emmeans(
    mod_success_empty_final,
    ~ Peak,
    type = "response") %>% 
  as_tibble() %>% 
  mutate(segment = "false positives")


emmeans_detail_focus_full <-
  emmeans(
    mod_detail_focus_final,
    pairwise ~ position +  smooth + position:smooth,
    type = "response") 

emmeans_detail_empty_full <-
  emmeans(
    mod_detail_empty_final,
    pairwise ~ position +  smooth + position:smooth,
    type = "response") 

emmeans_detail_focus_full_tibble <- 
  tibble(
    emmeans_detail_focus_full$emmeans %>% 
      as_tibble(),
    segment = "correct detection") 

emmeans_detail_empty_full_tibble <- 
  tibble(
    emmeans_detail_empty_full$emmeans %>% 
      as_tibble(),
    segment = "false positives")

#----------------------------------------------------------#
# 3. Plot figures  -----
#----------------------------------------------------------#

#---------------------------------------------#
# 3.1 Peak points & WU (Fig. 2) -----
#---------------------------------------------#

bind_rows(
  emmeans_success_focus_full_tibble,
  emmeans_success_empty_full_tibble) %>%  
  ggplot(
    aes(
      y = response,
      x = dataset_type,
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
  
  facet_grid(WU ~ Peak) +
  
  scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_segment) +
  scale_fill_manual(values = color_legen_segment) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -90),
    text = element_text(size = text_size))


bind_rows(
  emmeans_detail_focus_full_tibble,
  emmeans_detail_empty_full_tibble) %>%  
  ggplot(
    aes(
      y = response,
      x = position,
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
  
  facet_grid( ~ smooth) +
  
  scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_segment) +
  scale_fill_manual(values = color_legen_segment) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -90),
    text = element_text(size = text_size))



bind_rows(
  emmeans_success_focus_WU,
  emmeans_success_empty_WU) %>% 
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
  
  scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_segment) +
  scale_fill_manual(values = color_legen_segment) +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "right",
    text = element_text(size = text_size)) 


bind_rows(
  emmeans_success_focus_peak,
  emmeans_success_empty_peak) %>% 
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
  
  scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = color_legen_segment) +
  scale_fill_manual(values = color_legen_segment) +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "right",
    text = element_text(size = text_size)) 



