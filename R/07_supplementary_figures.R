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
# 1. (Fig S1) Data generation figure  -----
#----------------------------------------------------------#

# visual definition for charts
env_lmits <- c(95, 120)
grain_limit <- c(0, 150)


# 1.1. data generation ----

# 1.1.1. recent -----

random_data_recent <- 
  .simulate.pollen.data(
    time = time_seq,
    nforc = N_env, 
    mean = 100, 
    sdev = .15,
    nprox = high_diversity,
    var = 20,
    range = 20,
    manual_edit = TRUE,
    breaks = breaks_recent,
    jitter = TRUE,
    rarity = TRUE,
    transform_to_counts = T,
    N_pollen_grains = 300)

random_data_recent_env_var <-
  random_data_recent$env_var %>% 
  as.data.frame()

names(random_data_recent_env_var) <- paste0("V",1:N_env)

plot_random_data_recent_env <-
  bind_cols(
    random_data_recent_env_var,
    random_data_recent$list_ages) %>% 
  dplyr::select(-sample.id) %>%
  pivot_longer(., cols = -c(age)) %>%
  arrange(age, value) %>%
  ggplot(
    aes(
      x = age,
      y = value)) +
  
  geom_vline(
    xintercept = breaks_recent,
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

plot_random_data_recent_pollen <-
  bind_cols(
    random_data_recent$community_data,
    random_data_recent$list_ages) %>% 
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
    xintercept = breaks_recent,
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

plot_random_data_recent_full <-
  ggarrange(
    plot_random_data_recent_env + 
      rremove("xy.text") + rremove("ticks"),
    plot_random_data_recent_pollen + 
      rremove("xy.text") + rremove("ticks"),
    nrow = 1,
    align = "h"
  )

plot_random_data_recent_full


# 1.1.2 late ----

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
    y = "Value of env. variable",
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
    y = "Number of pollen grains",
    x = "")

plot_random_data_late_full <-
  ggarrange(
    plot_random_data_late_env + 
      rremove("y.text") + rremove("y.ticks"),
    plot_random_data_late_pollen + 
      rremove("y.text") + rremove("y.ticks"),
    nrow = 1,
    align = "h"
  )

plot_random_data_late_full

# 1.2. Plot figure  -----

figure_S1 <- 
  ggarrange(
    plot_random_data_recent_full,
    plot_random_data_late_full,
    nrow = 2,
    labels = LETTERS[1:2]) %>% 
  annotate_figure(.,
                  left = text_grob(
                    "Age (cal yr BP)",
                    size = text_size,
                    rot = 90))

ggsave(
  "data/output/figures/fig_S1_raw.pdf",
  figure_S1,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)


#----------------------------------------------------------#
# 2. (Fig S2) Peak point detection figure  -----
#----------------------------------------------------------#

# load data
ROC_all <- read_rds("data/output/datasets/RoC/ROC_all.rds")

# findount one example of calculation
calculation_number_select <- 
  ROC_all %>% 
  filter(
    smooth == "none",
    DC == "chisq",
    diversity == "high_diversity",
    position == "breaks_recent",
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
  .plot.roc.curve(., roc_max = 1) +
  facet_wrap(~peak_type, nrow = 1) + 
  geom_rect(
    data = tibble(
      Age = 0,
      ROC = 0,
      Age_min = breaks_recent[1] - 500,
      Age_max = breaks_recent[2] + 500,
      RoC_min = -Inf,
      RoC_max = Inf),
    aes(
      xmin = Age_min,
      xmax = Age_max,
      ymin = RoC_min,
      ymax = RoC_max),
    alpha = 1/3,
    fill = gray_light))

ggsave(
  "data/output/figures/fig_S2_raw.pdf",
  figure_S2,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

