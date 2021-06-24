#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                   Plotting examples
#         of data generation and peak point detection
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 0. load data and tibble creation  -----
#----------------------------------------------------------#

# ROC results
ROC_all <- read_rds("data/output/datasets/RoC/ROC_all.rds")

# site A
data_site_A <- 
  read_rds("data/output/datasets/pollen_sites/data_site_A.rds")

data_site_A_merged <- 
  inner_join(
    data_site_A$age_data,
    data_site_A$community_data,
    by = "sample.id")

bin_size_result <- 
  read_rds("data/output/datasets/RoC/bin_size_result.rds")

mod_detail_correct_select <-
  read_rds("data/output/models/mod_detail_correct_select.rds")

mod_detail_false_select <-
  read_rds("data/output/models/mod_detail_false_select.rds")

data_detail_correct <- 
  read_rds("data/output/datasets/subsets/data_detail_correct.rds")

data_detail_false <- 
  read_rds("data/output/datasets/subsets/data_detail_false.rds")


emmeans_detail_correct_smooth <-
  emmeans(
    mod_detail_correct_select,
    pairwise ~  smooth ,
    type = "response") 

emmeans_detail_false_smooth <-
  emmeans(
    mod_detail_false_select,
    pairwise ~  smooth ,
    type = "response") 

emmeans_detail_correct_smooth_tibble <- 
  tibble(
    emmeans_detail_correct_smooth$emmeans %>% 
      as_tibble(),
    segment = "correct detection") %>% 
  mutate(smooth = fct_relevel(smooth,"None", "Shep", "M_avg", "Age_w", "Grimm"))


emmeans_detail_false_smooth_tibble <- 
  tibble(
    emmeans_detail_false_smooth$emmeans %>% 
      as_tibble(),
    segment = "false positives") %>% 
  mutate(smooth = fct_relevel(smooth,"None", "Shep", "M_avg", "Age_w", "Grimm"))

emmeans_detail_correct_DC <-
  emmeans(
    mod_detail_correct_select,
    pairwise ~  DC ,
    type = "response") 

emmeans_detail_false_DC <-
  emmeans(
    mod_detail_false_select,
    pairwise ~  DC ,
    type = "response") 

emmeans_detail_correct_DC_tibble <- 
  tibble(
    emmeans_detail_correct_DC$emmeans %>% 
      as_tibble(),
    segment = "correct detection") 


emmeans_detail_false_DC_tibble <- 
  tibble(
    emmeans_detail_false_DC$emmeans %>% 
      as_tibble(),
    segment = "false positives") 


emmeans_detail_correct_position  <-
  emmeans(
    mod_detail_correct_select,
    pairwise ~  position  ,
    type = "response") 

emmeans_detail_false_position  <-
  emmeans(
    mod_detail_false_select,
    pairwise ~  position  ,
    type = "response") 

emmeans_detail_correct_position_tibble <- 
  tibble(
    emmeans_detail_correct_position$emmeans %>% 
      as_tibble(),
    segment = "correct detection") 


emmeans_detail_false_position_tibble <- 
  tibble(
    emmeans_detail_false_position$emmeans %>% 
      as_tibble(),
    segment = "false positives") 



#----------------------------------------------------------#
# 1. (Fig S1) Data generation figure  -----
#----------------------------------------------------------#

# visual definition for charts
env_lmits <- c(95, 120)
grain_limit <- c(0, 150)


# 1.1. data generation ----

# 1.1.1. late -----

random_data_late <- 
  .simulate.pollen.data(
    time = time_seq,
    nforc = N_env, 
    mean = 100, 
    sdev = .15,
    nprox = high_diversity,
    var = 20,
    range = 20,
    manual_edit = TRUE,
    breaks = breaks_late,
    jitter = TRUE,
    rarity = TRUE,
    transform_to_counts = T,
    N_pollen_grains = 300)

random_data_late_env_var <-
  random_data_late$env_var %>% 
  as.data.frame()

names(random_data_late_env_var) <- paste0("V",1:N_env)

plot_random_data_late_env <-
  bind_cols(
    random_data_late_env_var,
    random_data_late$list_ages) %>% 
  dplyr::select(-sample.id) %>%
  pivot_longer(., cols = -c(age)) %>%
  arrange(age, value) %>%
  ggplot(
    aes(
      x = age,
      y = value)) +
  
  geom_vline(
    xintercept = breaks_late,
    color = gray_light,
    size = line_size) +
  
  geom_line(
    aes(color = name),
    size = line_size) +
  
  coord_flip(
    xlim = c(age_lim, 0)) +
  
  scale_x_continuous(trans = "reverse") +
  scale_y_continuous(limits = env_lmits) +
  theme(
    line = element_line(size = line_size),
    text = element_text(size = text_size),
    legend.position = "none") +
  
  labs(
    y = "",
    x = "")

plot_random_data_late_pollen <-
  bind_cols(
    random_data_late$community_data,
    random_data_late$list_ages) %>% 
  dplyr::select(-sample.id) %>%
  pivot_longer(., cols = -c(age)) %>%
  arrange(age, value) %>%
  ggplot(
    aes(
      x = age,
      y = value)) +
  
  geom_ribbon(
    aes(
      ymin = rep(0,length(value)),
      ymax = value,
      fill = name), 
    color = gray_dark,
    alpha=1/5,
    size = line_size)+
  
  geom_vline(
    xintercept = breaks_late,
    color = gray_light,
    size = line_size) +
  
  coord_flip(
    xlim = c(age_lim, 0)) +
  
  scale_x_continuous(trans = "reverse") +
  scale_y_continuous(limits = grain_limit) +
  theme(
    line = element_line(size = line_size),
    text = element_text(size = text_size),
    legend.position = "none") +
  
  labs(
    y = "",
    x = "")

plot_random_data_late_full <-
  ggarrange(
    plot_random_data_late_env + 
      rremove("xy.text") + rremove("ticks"),
    plot_random_data_late_pollen + 
      rremove("xy.text") + rremove("ticks"),
    nrow = 1,
    align = "h"
  )

plot_random_data_late_full


# 1.1.2. early ----

random_data_early <- 
  .simulate.pollen.data(
    time = time_seq,
    nforc = N_env, 
    mean = 100, 
    sdev = .15,
    nprox = high_diversity,
    var = 20,
    range = 20,
    manual_edit = TRUE,
    breaks = breaks_early,
    jitter = TRUE,
    rarity = TRUE,
    transform_to_counts = T,
    N_pollen_grains = 300)

random_data_early_env_var <-
  random_data_early$env_var %>% 
  as.data.frame()

names(random_data_early_env_var) <- paste0("V",1:N_env)

plot_random_data_early_env <-
  bind_cols(
    random_data_early_env_var,
    random_data_early$list_ages) %>% 
  dplyr::select(-sample.id) %>%
  pivot_longer(., cols = -c(age)) %>%
  arrange(age, value) %>%
  ggplot(
    aes(
      x = age,
      y = value)) +
  
  geom_vline(
    xintercept = breaks_early,
    color = gray_light,
    size = line_size) +
  
  geom_line(
    aes(color = name),
    size = line_size) +
  
  coord_flip(
    xlim = c(age_lim, 0)) +
  
  scale_x_continuous(trans = "reverse") +
  scale_y_continuous(limits = env_lmits) +
  theme(
    line = element_line(size = line_size),
    text = element_text(size = text_size),
    legend.position = "none") +
  
  labs(
    y = "Value of env. variable",
    x = "")

plot_random_data_early_pollen <-
  bind_cols(
    random_data_early$community_data,
    random_data_early$list_ages) %>% 
  dplyr::select(-sample.id) %>%
  pivot_longer(., cols = -c(age)) %>%
  arrange(age, value) %>%
  ggplot(
    aes(
      x = age,
      y = value)) +
  
  geom_ribbon(
    aes(
      ymin = rep(0,length(value)),
      ymax = value,
      fill = name), 
    color = gray_dark,
    alpha=1/5,
    size = line_size)+
  
  geom_vline(
    xintercept = breaks_early,
    color = gray_light,
    size = line_size) +
  
  coord_flip(
    xlim = c(age_lim, 0)) +
  
  scale_x_continuous(trans = "reverse") +
  scale_y_continuous(limits = grain_limit) +
  theme(
    line = element_line(size = line_size),
    text = element_text(size = text_size),
    legend.position = "none") +
  
  labs(
    y = "Number of pollen grains",
    x = "")

plot_random_data_early_full <-
  ggarrange(
    plot_random_data_early_env + 
      rremove("y.text") + rremove("y.ticks"),
    plot_random_data_early_pollen + 
      rremove("y.text") + rremove("y.ticks"),
    nrow = 1,
    align = "h"
  )

plot_random_data_early_full

# 1.2. Plot figure  -----

(figure_S1 <- 
   ggarrange(
     plot_random_data_late_full,
     plot_random_data_early_full,
     nrow = 2,
     labels = LETTERS[1:2]) %>% 
   annotate_figure(.,
                   left = text_grob(
                     "Age (cal yr BP)",
                     size = text_size,
                     rot = 90)))

ggsave(
  "data/output/figures/fig_S1_raw.pdf",
  figure_S1,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

# save the data

write_rds(
  random_data_late,
  "data/output/datasets/simulated/fig_S1_random_data_late.rds")

write_rds(
  random_data_early,
  "data/output/datasets/simulated/fig_S1_random_data_early.rds")



#----------------------------------------------------------#
# 2. (Fig S2) Peak point detection figure  -----
#----------------------------------------------------------#


# findount one example of calculation
calculation_number_select <- 
  ROC_all %>% 
  filter(
    smooth == "none",
    DC == "chisq",
    diversity == "high_diversity",
    position == "breaks_late",
    WU == "MW") %>% 
  arrange(calculation_ID) %>% 
  distinct(calculation_ID) %>% 
  slice(4) %>% 
  pluck(1)

(figure_S2 <- 
    ROC_all %>% 
    filter(calculation_ID == calculation_number_select) %>% 
    dplyr::select(Age,ROC,ROC_dw,ROC_up, contains("Peak")) %>% 
    pivot_longer(
      cols = -c(Age, ROC, ROC_dw, ROC_up),
      names_to = "peak_type",
      values_to = "Peak") %>% 
    mutate(peak_type = str_replace(peak_type, "Peak_", "")) %>% 
    mutate(
      peak_type = factor(
        peak_type,
        levels = c(
          "threshold",
          "trend_linear",
          "trend_non_linear",
          "trend_GAM_deriv",
          "SNI"))) %>% 
    .plot.roc.curve(., roc_max = 1)+
    facet_wrap(~peak_type, nrow = 1) + 
    geom_rect(
      data = tibble(
        Age = 0,
        ROC = 0,
        Age_min = breaks_late[1] - 500,
        Age_max = breaks_late[2] + 500,
        RoC_min = -Inf,
        RoC_max = Inf),
      aes(
        xmin = Age_min,
        xmax = Age_max,
        ymin = RoC_min,
        ymax = RoC_max),
      alpha = 1/3,
      fill = gray_light)+
    geom_vline(xintercept = c(breaks_late[1], breaks_late[2]),
               color = gray_dark))

ggsave(
  "data/output/figures/fig_S2_raw.pdf",
  figure_S2,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

#----------------------------------------------------------#
# 3. (Fig S3) Comparison of success detection   -----
#----------------------------------------------------------#

dataset_id_select_late <-  
  ROC_all %>% 
  filter(
    smooth == "none",
    DC == "chisq",
    diversity == "high_diversity",
    position == "breaks_late") %>% 
  arrange(dataset_ID) %>% 
  distinct(dataset_ID) %>% 
  slice(1) %>% 
  pluck(1)

dataset_id_select_early <-  
  ROC_all %>% 
  filter(
    smooth == "none",
    DC == "chisq",
    diversity == "high_diversity",
    position == "breaks_early") %>% 
  arrange(dataset_ID) %>% 
  distinct(dataset_ID) %>% 
  slice(1) %>% 
  pluck(1)



comparison_data_select <-
  ROC_all %>% 
  dplyr::filter(dataset_ID  %in% c(dataset_id_select_late,
                                   dataset_id_select_early) ) %>% 
  dplyr::mutate(
    WU = factor(WU, levels = c("levels", "bins", "MW")),
    DC = factor(DC, levels = c("chord", "chisq"))) %>% 
  dplyr::mutate(
    DC = fct_recode(
      DC,
      "Chord" = "chord",
      "Chisq" = "chisq"))



(figure_S3_A <- 
    comparison_data_select %>% 
    dplyr::filter(dataset_ID == dataset_id_select_late) %>% 
    ggplot(
      aes(
        y = ROC,
        x = Age))+
    geom_rect(
      data = tibble(
        Age = 0,
        ROC = 0,
        Age_min = breaks_late[1] - 500,
        Age_max = breaks_late[2] + 500,
        RoC_min = -Inf,
        RoC_max = Inf),
      aes(
        xmin = Age_min,
        xmax = Age_max,
        ymin = RoC_min,
        ymax = RoC_max),
      alpha = 1/3,
      fill = gray_light)+
    geom_vline(
      xintercept = c(breaks_late[1], breaks_late[2]),
      color = gray_dark)+
    geom_rug(
      data = data_site_A$age_data,
      aes(x = age, y = NULL)) + 
    geom_line(
      aes(group = calculation_ID),
      alpha = 1/3,
      color = gray_dark)+
    geom_point(alpha = 1/3, colour = gray_dark)+
    geom_point( data = comparison_data_select %>% 
                  dplyr::filter(dataset_ID == dataset_id_select_late) %>% 
                  dplyr::filter(Peak_trend_non_linear  == TRUE),
                alpha = 1/3, colour = "green")+
    scale_x_continuous(breaks = seq(0,8e3,1e3))+
    facet_wrap( ~ WU, nrow = 3, scales = "free_y" ))


(figure_S3_B <- 
    comparison_data_select %>% 
    dplyr::filter(dataset_ID == dataset_id_select_early) %>% 
    ggplot(
      aes(
        y = ROC,
        x = Age))+
    geom_rect(
      data = tibble(
        Age = 0,
        ROC = 0,
        Age_min = breaks_early[1] - 500,
        Age_max = breaks_early[2] + 500,
        RoC_min = -Inf,
        RoC_max = Inf),
      aes(
        xmin = Age_min,
        xmax = Age_max,
        ymin = RoC_min,
        ymax = RoC_max),
      alpha = 1/3,
      fill = gray_light)+
    geom_vline(
      xintercept = c(breaks_early[1], breaks_early[2]),
      color = gray_dark)+
    geom_rug(
      data = data_site_A$age_data,
      aes(x = age, y = NULL)) + 
    geom_line(
      aes(group = calculation_ID),
      alpha = 1/3,
      color = gray_dark)+
    geom_point(alpha = 1/3, colour = gray_dark)+
    geom_point( data = comparison_data_select %>% 
                  dplyr::filter(dataset_ID == dataset_id_select_early) %>% 
                  dplyr::filter(Peak_trend_non_linear  == TRUE),
                alpha = 1/3, colour = "green")+
    scale_x_continuous(breaks = seq(0,8e3,1e3))+
    facet_wrap( ~ WU, nrow = 3, scales = "free_y" ))


(figure_S3 <- 
    ggarrange(
      figure_S3_A + ggpubr::rremove("xylab"),
      figure_S3_B + ggpubr::rremove("xylab"),
      nrow = 1,
      labels = LETTERS[1:2]) %>% 
    annotate_figure(.,
                    bottom = text_grob(
                      "Age (cal yr BP)",
                      size = text_size),
                    left = text_grob(
                      "Rate−of−Change score",
                      size = text_size,
                      rot = 90)))

ggsave(
  "data/output/figures/fig_S3_raw.pdf",
  figure_S3,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)


#----------------------------------------------------------#
# 4. (Fig S4) resolution sensitivity  -----
#----------------------------------------------------------#

bin_size <- 1000
res <- c(20, 50, 70, 100, 150, 200, 250, 300, 400, 500)

resolution_levels <- 
  .thin.data(WU = "levels")

resolution_bins <- 
  .thin.data(WU = "bins", size_of_bin = bin_size)

resolution_MW <- 
  .thin.data(WU = "MW", size_of_bin = bin_size)

resolution_merge <-
  bind_rows(
    resolution_levels %>% 
      mutate(
        WU = "individual levels"),
    resolution_bins %>% 
      mutate(
        WU = "binning"),
    resolution_MW %>% 
      mutate(
        WU = "mowing window")) %>% 
  dplyr::mutate(
    working_unit = factor(
      WU, 
      levels = c("individual levels", "binning", "mowing window")))

write_rds(
  resolution_merge,
  "data/output/datasets/sensitivity/resolution_merge.rds")

(figure_S4 <-
    resolution_merge %>% 
    ggplot() + 
    geom_rug(
      data = data_site_A$age_data,
      aes(x = age)) +
    geom_line(
      aes(x = Age, 
          y = ROC, 
          colour = resolution)) +
    geom_point(
      data = resolution_merge %>% 
        dplyr::filter(Peak  == TRUE),
      aes(x = Age, 
          y = ROC, 
          colour = resolution))+
    scale_x_continuous(breaks = seq(0,8e3,1e3))+
    viridis::scale_colour_viridis(direction = -1, discrete= TRUE)+
    labs(
      y = "Rate−of−Change score",
      x = "Age (cal yr BP)", 
      colour = "resolution (yr)") +
    facet_wrap(~working_unit, nrow = 3, scales = "free_y")+
    theme(
      line = element_line(size = line_size),
      text = element_text(size = text_size)))

ggsave(
  "data/output/figures/fig_S4_raw.pdf",
  figure_S4,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

#----------------------------------------------------------#
# 5. (Fig S5) hiatus effect  -----
#----------------------------------------------------------#

hiatus_length_vector <- c(0, 100, 200, 500, 1000, 2000)
hiatus_age <- 3e3

hiatus_levels <- 
  .create.hiatus(WU = "levels", hiatus_age = hiatus_age)

hiatus_bins <- 
  .create.hiatus(WU = "bins", size_of_bin = 500, hiatus_age = hiatus_age)

hiatus_MW <- 
  .create.hiatus(WU = "MW", size_of_bin = 500, hiatus_age = hiatus_age)

hiatus_merge <-
  bind_rows(
    hiatus_levels %>% 
      mutate(
        WU = "individual levels"),
    hiatus_bins %>% 
      mutate(
        WU = "binning"),
    hiatus_MW %>% 
      mutate(
        WU = "mowing window")) %>% 
  dplyr::mutate(
    working_unit = factor(
      WU, 
      levels = c("individual levels", "binning", "mowing window")),
    hiatus_length = factor(
      hiatus_length,
      levels = hiatus_length_vector) )

(figure_S5 <- 
    hiatus_merge %>% 
    ggplot() + 
    geom_rug(
      data = data_site_A$age_data,
      aes(x = age)) +
    geom_vline(xintercept = hiatus_age, color = gray_light)+
    geom_line(
      aes(x = Age, 
          y = ROC, 
          colour = hiatus_length)) +
    geom_point(
      data = hiatus_merge %>% 
        dplyr::filter(Peak  == TRUE),
      aes(x = Age, 
          y = ROC, 
          colour = hiatus_length))+
    scale_x_continuous(breaks = seq(0,8e3,1e3))+
    viridis::scale_colour_viridis(direction = -1, discrete= TRUE)+
    labs(
      y = "Rate−of−Change score",
      x = "Age (cal yr BP)", 
      colour = "hiatus length (yr)") +
    facet_wrap(~working_unit, nrow = 3, scales = "free_y")+
    theme(
      line = element_line(size = line_size),
      text = element_text(size = text_size)))

ggsave(
  "data/output/figures/fig_S5_raw.pdf",
  figure_S5,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

#----------------------------------------------------------#
# 6. (Fig S6) missing top of the core  -----
#----------------------------------------------------------#

top_core_vector <- c(0, 100, 200, 500, 1000, 2000)

top_core_levels <- 
  .remove.top.core(WU = "levels")

top_core_bins <- 
  .remove.top.core(WU = "bins", size_of_bin = 500)

top_core_MW <- 
  .remove.top.core(WU = "MW", size_of_bin = 500)

top_core_merge <-
  bind_rows(
    top_core_levels %>% 
      mutate(
        WU = "individual levels"),
    top_core_bins %>% 
      mutate(
        WU = "binning"),
    top_core_MW %>% 
      mutate(
        WU = "mowing window")) %>% 
  dplyr::mutate(
    working_unit = factor(
      WU, 
      levels = c("individual levels", "binning", "mowing window")),
    top_core_age = factor(
      top_core_age ,
      levels = top_core_vector) )


(figure_S6 <- 
    top_core_merge %>% 
    ggplot() + 
    geom_rug(
      data = data_site_A$age_data,
      aes(x = age)) +
    geom_line(
      aes(x = Age, 
          y = ROC, 
          colour = top_core_age)) +
    geom_point(
      data = top_core_merge %>% 
        dplyr::filter(Peak  == TRUE),
      aes(x = Age, 
          y = ROC, 
          colour = top_core_age))+
    scale_x_continuous(breaks = seq(0,8e3,1e3))+
    viridis::scale_colour_viridis(direction = -1, discrete= TRUE)+
    labs(
      y = "Rate−of−Change score",
      x = "Age (cal yr BP)", 
      colour = "removed top core (yr)") +
    facet_wrap(~working_unit, nrow = 3, scales = "free_y")+
    theme(
      line = element_line(size = line_size),
      text = element_text(size = text_size)))

ggsave(
  "data/output/figures/fig_S6_raw.pdf",
  figure_S6,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

#----------------------------------------------------------#
# 7. (Fig S7) binning affect on RoC -----
#----------------------------------------------------------#


bin_size_vector <- c(100,200,500,750,1e3,1500,3e3)

(figure_S7 <- 
    bin_size_result %>%
    ggplot() +
    
    geom_rug(
      data = data_site_A$age_data,
      aes(x = age),
      sides = "b") +
    
    geom_line(
      aes(
        x = Age,
        y= ROC,
        group = bin,
        color = factor(bin, levels = bin_size_vector))) +
    
    viridis::scale_colour_viridis(direction = -1, discrete= TRUE) +
    
    theme(
      line = element_line(size = line_size),
      text = element_text(size = text_size)) +
    
    labs(
      x ="Age (cal yr BP)",
      y = "Rate-of-Change score \n [Chisq-distance per 500 yr]",
      color = "time bin size"))


ggsave(
  "data/output/figures/fig_S7_raw.pdf",
  figure_S7,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

#----------------------------------------------------------#
# 8 (Fig. S8) Comparion of methods  -----
#----------------------------------------------------------#

(plot_detail_smooth <- 
   bind_rows(
     emmeans_detail_correct_smooth_tibble,
     emmeans_detail_false_smooth_tibble) %>% 
   ggplot(
     aes(
       y = response,
       x = smooth,
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
   
   #scale_y_continuous(limits = c(0,1))+
   scale_color_manual(values = color_legen_segment) +
   scale_fill_manual(values = color_legen_segment) +
   theme(
     legend.position = "right",
     line = element_line(size = line_size),
     text = element_text(size = text_size))+
   labs(
     y = "Proportion of detected Working Units",
     x = "Smoothing method",
     color = "",
     fill = ""))


(plot_detail_DC <- 
    bind_rows(
      emmeans_detail_correct_DC_tibble,
      emmeans_detail_false_DC_tibble) %>% 
    ggplot(
      aes(
        y = response,
        x = DC,
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
    
    #scale_y_continuous(limits = c(0,1))+
    scale_color_manual(values = color_legen_segment) +
    scale_fill_manual(values = color_legen_segment) +
    theme(
      legend.position = "right",
      line = element_line(size = line_size),
      text = element_text(size = text_size))+
    labs(
      y = "Proportion of detected Working Units",
      x = "Dissimilarity coeficient",
      color = "",
      fill = ""))

(plot_detail_position <- 
    bind_rows(
      emmeans_detail_correct_position_tibble,
      emmeans_detail_false_position_tibble) %>% 
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
    
    #scale_y_continuous(limits = c(0,1))+
    scale_color_manual(values = color_legen_segment) +
    scale_fill_manual(values = color_legen_segment) +
    theme(
      legend.position = "right",
      line = element_line(size = line_size),
      text = element_text(size = text_size))+
    labs(
      y = "Proportion of detected Working Units",
      x = "Position of enviromental change",
      color = "",
      fill = ""))

(figure_S8 <- 
    ggarrange(
      plot_detail_DC +  rremove("ylab"),
      plot_detail_smooth + rremove("ylab"),
      plot_detail_position + rremove("ylab"),
      nrow = 3,
      labels = LETTERS[1:3],
      common.legend = TRUE,
      legend = "right"
    ) %>% 
    annotate_figure(
      left = text_grob(
        "Proportion of detected Working Units",
        size = text_size,
        rot = 90
      )
    )) 


ggsave(
  "data/output/figures/fig_S8_raw.pdf",
  figure_S8,
  width = pdf_width,
  height = pdf_height * 2,
  units = pdf_units)


