#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                 Plotting RoC figures
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 1. Load data  -----
#----------------------------------------------------------#

data_site_A <- read_rds("data/output/datasets/pollen_sites/data_site_A.rds")
data_site_B <- read_rds("data/output/datasets/pollen_sites/data_site_B.rds")
data_site_C <- read_rds("data/output/datasets/pollen_sites/data_site_C.rds")
data_site_D <- read_rds("data/output/datasets/pollen_sites/data_site_D.rds")

data_site_A_RoC_levels <- 
  read_rds("data/output/example_data_roc/data_site_A_RoC_levels.rds")
data_site_A_RoC_bins <-
  read_rds("data/output/example_data_roc/data_site_A_RoC_bins.rds")
data_site_A_RoC_MW <-
  read_rds("data/output/example_data_roc/data_site_A_RoC_MW.rds")

data_site_B_RoC_levels <- 
  read_rds("data/output/example_data_roc/data_site_B_RoC_levels.rds")
data_site_B_RoC_bins <-
  read_rds("data/output/example_data_roc/data_site_B_RoC_bins.rds")
data_site_B_RoC_MW <-
  read_rds("data/output/example_data_roc/data_site_B_RoC_MW.rds")

data_site_C_RoC_levels <- 
  read_rds("data/output/example_data_roc/data_site_C_RoC_levels.rds")
data_site_C_RoC_bins <-
  read_rds("data/output/example_data_roc/data_site_C_RoC_bins.rds")
data_site_C_RoC_MW <-
  read_rds("data/output/example_data_roc/data_site_C_RoC_MW.rds")

data_site_D_RoC_levels <- 
  read_rds("data/output/example_data_roc/data_site_D_RoC_levels.rds")
data_site_D_RoC_bins <-
  read_rds("data/output/example_data_roc/data_site_D_RoC_bins.rds")
data_site_D_RoC_MW <-
  read_rds("data/output/example_data_roc/data_site_D_RoC_MW.rds")

ROC_site_A_all_settings <-
  read_rds("data/output/datasets/RoC/ROC_site_A_all_settings.rds")

bin_size_result <- 
  read_rds("data/output/datasets/RoC/bin_size_result.rds")


#----------------------------------------------------------#
# 2. Extract pollen data  -----
#----------------------------------------------------------#

# common taxa in each dataset
data_site_A_dom <- .extract.dominant.pollen.taxa(data_site_A)
data_site_B_dom <- .extract.dominant.pollen.taxa(data_site_B)
data_site_C_dom <- .extract.dominant.pollen.taxa(data_site_C)
data_site_D_dom <- .extract.dominant.pollen.taxa(data_site_D)

# Common pollen taxa in all datasets
common_taxa <- 
  c(data_site_A_dom,
    data_site_B_dom,
    data_site_C_dom,
    data_site_D_dom) %>%
  unique() %>%
  sub("/",".",.) %>%
  sub("-",".",.) %>%
  sub(")",".",.) %>%
  sub(".\\(","..",.) 

# prepare pollen data for the common taxa
data_site_A_pollen <- .get.pollen.data(data_site_A, common_taxa)
data_site_B_pollen <- .get.pollen.data(data_site_B, common_taxa)
data_site_C_pollen <- .get.pollen.data(data_site_C, common_taxa)
data_site_D_pollen <- .get.pollen.data(data_site_D, common_taxa)

#----------------------------------------------------------#
# 4. Create individual plots  -----
#----------------------------------------------------------#

# 4.1. Sample density -----
plot_site_A_density <- .plot.density.of.samples(data_site_A)
plot_site_B_density <- .plot.density.of.samples(data_site_B)
plot_site_C_density <- .plot.density.of.samples(data_site_C)
plot_site_D_density <- .plot.density.of.samples(data_site_D)

# 4.2. pollen distribution -----
plot_site_A_pollen <- 
  .plot.pollen.data(data_site_A_pollen, common_taxa, strip_names = T)

plot_site_B_pollen <- 
  .plot.pollen.data(data_site_B_pollen, common_taxa)

plot_site_C_pollen <- 
  .plot.pollen.data(data_site_C_pollen, common_taxa)

plot_site_D_pollen <- 
  .plot.pollen.data(data_site_D_pollen, common_taxa)


# 4.3. RoC curve -----
max_roc <- 1.5

plot_site_A_roc_levels <- .plot.roc.curve(data_site_A_RoC_levels, roc_max = max_roc)
plot_site_B_roc_levels <- .plot.roc.curve(data_site_B_RoC_levels, roc_max = max_roc)
plot_site_C_roc_levels <- .plot.roc.curve(data_site_C_RoC_levels, roc_max = max_roc)
plot_site_D_roc_levels <- .plot.roc.curve(data_site_D_RoC_levels, roc_max = max_roc)

plot_site_A_roc_bins <- .plot.roc.curve(data_site_A_RoC_bins, roc_max = max_roc)
plot_site_B_roc_bins <- .plot.roc.curve(data_site_B_RoC_bins, roc_max = max_roc)
plot_site_C_roc_bins <- .plot.roc.curve(data_site_C_RoC_bins, roc_max = max_roc)
plot_site_D_roc_bins <- .plot.roc.curve(data_site_D_RoC_bins, roc_max = max_roc)

plot_site_A_roc_MW <- .plot.roc.curve(data_site_A_RoC_MW, roc_max = max_roc)
plot_site_B_roc_MW <- .plot.roc.curve(data_site_B_RoC_MW, roc_max = max_roc)
plot_site_C_roc_MW <- .plot.roc.curve(data_site_C_RoC_MW, roc_max = max_roc)
plot_site_D_roc_MW <- .plot.roc.curve(data_site_D_RoC_MW, roc_max = max_roc)


#----------------------------------------------------------#
# 5. (Fig 4) Building the figure 4 -----
#----------------------------------------------------------#

rel_w_density <- 1
rel_w_pollen <- 3
rel_w_roc <- 0.8

plot_site_A_full <-
  plot_grid(
    plot_site_A_density + 
      rremove("x.title") + rremove("x.text") + rremove("x.ticks"),
    plot_site_A_pollen + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_A_roc_levels + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_A_roc_bins + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_A_roc_MW + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    ncol = 5,
    align = "h",
    axis = "bt",
    rel_widths = c(rel_w_density, rel_w_pollen, rel_w_roc, rel_w_roc, rel_w_roc))

plot_site_B_full <-
  plot_grid(
    plot_site_B_density + 
      rremove("x.title") + rremove("x.text") + rremove("x.ticks"),
    plot_site_B_pollen + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_B_roc_levels + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_B_roc_bins + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_B_roc_MW + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    ncol = 5,
    align = "h",
    axis = "bt",
    rel_widths = c(rel_w_density, rel_w_pollen, rel_w_roc, rel_w_roc, rel_w_roc))

plot_site_C_full <-
  plot_grid(
    plot_site_C_density + 
      rremove("x.title") + rremove("x.text") + rremove("x.ticks"),
    plot_site_C_pollen + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_C_roc_levels + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_C_roc_bins + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    plot_site_C_roc_MW + 
      rremove("xy.title") + rremove("xy.text") + rremove("ticks"),
    ncol = 5,
    align = "h",
    axis = "bt",
    rel_widths = c(rel_w_density, rel_w_pollen, rel_w_roc, rel_w_roc, rel_w_roc))

plot_site_D_full <-
  plot_grid(
    plot_site_D_density, 
    plot_site_D_pollen + 
      rremove("y.title") + rremove("y.text") + rremove("y.ticks"),
    plot_site_D_roc_levels + 
      rremove("xy.title") + rremove("y.text") + rremove("y.ticks"),
    plot_site_D_roc_bins + 
      rremove("y.title") + rremove("y.text") + rremove("y.ticks"),
    plot_site_D_roc_MW + 
      rremove("xy.title") + rremove("y.text") + rremove("y.ticks"),
    ncol = 5,
    align = "h",
    axis = "bt",
    rel_widths = c(rel_w_density, rel_w_pollen, rel_w_roc, rel_w_roc, rel_w_roc))

figure_4 <-
  ggarrange(
    plot_site_A_full,
    plot_site_B_full,
    plot_site_C_full,
    plot_site_D_full,
    nrow = 4,
    labels = LETTERS[1:4],
    heights = c(1.5, 1, 1, 1.1),
    legend = "none")

ggsave(
  "data/output/figures/fig_4_raw.pdf",
  figure_4,
  height = pdf_height * 2,
  width = pdf_width,
  units = pdf_units)


#----------------------------------------------------------#
# 6. (Fig 5) Roc in multiple settings -----
#----------------------------------------------------------#

(figure_5 <- 
  ROC_site_A_all_settings %>% 
  unnest(ROC) %>% 
  mutate(
    smooth_type = fct_relevel(
      smooth_type,
      "none", "shep", "m.avg", "age.w", "grim"),
    DC = fct_relevel(
      DC,
      "chord","chisq"
    )) %>%
  mutate(
    smooth_type = fct_recode(
      smooth_type,
      "None" = "none",
      "Shep" = "shep",
      "M_avg" = "m.avg",
      "Age_w" = "age.w",
      "Grimm" = "grim"),
    DC = fct_recode(
      DC,
      "Chord" = "chord",
      "Chisq" = "chisq")) %>% 
  .plot.roc.curve(., roc_max = 1) +
  facet_grid(DC~smooth_type))

ggsave(
  "data/output/figures/fig_5_raw.pdf",
  figure_5,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)

#----------------------------------------------------------#
# 6. (Fig S3) binning affect on RoC -----
#----------------------------------------------------------#

(figure_S3 <- 
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
       color = bin)) +
   
   viridis::scale_colour_viridis(direction = -1) +
   
   theme(
     line = element_line(size = line_size),
     text = element_text(size = text_size)) +
   
   labs(
     x ="Age (cal yr BP)",
     y = "Rate-of-Change score",
     color = "Bin size"))

ggsave(
  "data/output/figures/fig_S3_raw.pdf",
  figure_S3,
  height = pdf_height,
  width = pdf_width,
  units = pdf_units)
