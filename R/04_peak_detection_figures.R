#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                    Graphical outcomes
#                  of peak point detection
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
    pairwise ~  DC + position + smooth + DC:position + DC:smooth + 
      position:smooth + DC:position:smooth ,
    type = "response") 

emmeans_detail_false_full <-
  emmeans(
    mod_detail_false_select,
    pairwise ~ DC + position + smooth + DC:position + DC:smooth + 
      position:smooth + DC:position:smooth,
    type = "response") 

emmeans_detail_correct_full_tibble <- 
  tibble(
    emmeans_detail_correct_full$emmeans %>% 
      as_tibble(),
    segment = "correct detection") %>% 
  mutate(smooth = fct_relevel(smooth,"None", "Shep", "M_avg", "Age_w", "Grimm"))

emmeans_detail_false_full_tibble <- 
  tibble(
    emmeans_detail_false_full$emmeans %>% 
      as_tibble(),
    segment = "false positives") %>% 
  mutate(smooth = fct_relevel(smooth,"None", "Shep", "M_avg", "Age_w", "Grimm"))



#----------------------------------------------------------#
# 3. Plot figures  -----
#----------------------------------------------------------#

#---------------------------------------------#
# 3.1 (Fig. 2) Peak points & WU  -----
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
     col = gray_dark,
     stat = "identity",
     width = 0.75,
     size = line_size,
     position = position_dodge(width = 0.75)) + 
   
   geom_errorbar(
     aes(
       ymin = lower.CL,
       ymax = upper.CL),
     col = gray_dark,
     size = line_size,
     position = position_dodge(width = 0.75),
     width = 0.2) +
   
   geom_point(
     shape = 0,
     col = gray_dark,
     position = position_dodge(width = 0.75),
     size = 1) +
   
   scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.05)) +
   scale_color_manual(values = color_legen_segment) +
   scale_fill_manual(values = color_legen_segment) +
   theme(
     axis.text.x = element_text(angle = 90, 
                                hjust = 0.95, 
                                vjust = 0.35),
     line = element_line(size = line_size),
     legend.position = "right",
     text = element_text(size = text_size))+
   labs(
     y = "Proportion of detected Working Units",
     x = "Working Unit type",
     color = "",
     fill = ""
   ))

emmeans_correct_WU

emmeans_false_WU


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
      col = gray_dark,
      stat = "identity",
      width = 0.75,
      size = line_size,
      position = position_dodge(width = 0.75)) + 
    
    geom_errorbar(
      aes(
        ymin = lower.CL,
        ymax = upper.CL),
      col = gray_dark,
      size = line_size,
      position = position_dodge(width = 0.75),
      width = 0.2) +
    
    geom_point(
      shape = 0,
      col = gray_dark,
      position = position_dodge(width = 0.75),
      size = 1) +
    
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.05)) +
    scale_color_manual(values = color_legen_segment) +
    scale_fill_manual(values = color_legen_segment) +
    theme(
      axis.text.x = element_text(angle = 90, 
                                 hjust = 0.95, 
                                 vjust = 0.35),
      line = element_line(size = line_size),
      legend.position = "right",
      text = element_text(size = text_size))+
    labs(
      y = "Proportion of detected Working Units",
      x = "Peak-point detection method",
      color = "",
      fill = ""))

emmeans_correct_peak

(figure_02 <- 
    ggarrange(
      plot_full_WU,
      plot_full_peak + rremove("ylab"),
      nrow = 1,
      align = "hv",
      labels = LETTERS[1:2],
      common.legend = TRUE,
      legend = "right"
    ))


ggsave("data/output/figures/fig_2_raw.pdf",
       figure_02,
       width = pdf_width, height = pdf_height, units = pdf_units)


#---------------------------------------------#
# 3.2 (Fig. 3) Properties of RoC and dataset type  -----
#---------------------------------------------#

(figure_03 <- 
   bind_rows(
     emmeans_detail_correct_full_tibble,
     emmeans_detail_false_full_tibble) %>% 
   mutate(
     position = stringr::str_replace(position, " density level", "") ) %>% 
   mutate(
     position = factor(position, levels = c("low","high") )) %>% 
   ggplot(
     aes(
       y = response,
       x = position,
       col = segment,
       fill = segment)) + 
   
   geom_bar(
     col = gray_dark,
     stat = "identity",
     width = 0.75,
     size = line_size,
     position = position_dodge(width = 0.75)) + 
   
   geom_errorbar(
     aes(
       ymin = lower.CL,
       ymax = upper.CL),
     col = gray_dark,
     size = line_size, 
     position = position_dodge(width = 0.75),
     width = 0.2) +
   
   geom_point(
     shape = 0,
     col = gray_dark,
     position = position_dodge(width = 0.75),
     size = 1) +
   
   facet_grid(DC ~ smooth) +
   
   #scale_y_continuous(limits = c(0,1))+
   scale_color_manual(values = color_legen_segment) +
   scale_fill_manual(values = color_legen_segment) +
   theme(
     legend.position = "right",
     axis.text.x = element_text(angle = 90, 
                                hjust = 0.95, 
                                vjust = 0.35),
     line = element_line(size = line_size),
     text = element_text(size = text_size))+
   labs(
     y = "Proportion of detected Working Units",
     x = "Density of levels in position of change",
     color = "",
     fill = ""))


ggsave("data/output/figures/fig_3_raw.pdf",
       figure_03,
       width = pdf_width, height = pdf_height, units = pdf_units)

