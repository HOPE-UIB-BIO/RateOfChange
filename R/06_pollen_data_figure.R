#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                     Calculation RoC 
#                     in example data
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
plot_site_A_roc_levels <- .plot.roc.curve(data_site_A_RoC_levels)
plot_site_B_roc_levels <- .plot.roc.curve(data_site_B_RoC_levels)
plot_site_C_roc_levels <- .plot.roc.curve(data_site_C_RoC_levels)
plot_site_D_roc_levels <- .plot.roc.curve(data_site_D_RoC_levels)

plot_site_A_roc_bins <- .plot.roc.curve(data_site_A_RoC_bins)
plot_site_B_roc_bins <- .plot.roc.curve(data_site_B_RoC_bins)
plot_site_C_roc_bins <- .plot.roc.curve(data_site_C_RoC_bins)
plot_site_D_roc_bins <- .plot.roc.curve(data_site_D_RoC_bins)

plot_site_A_roc_MW <- .plot.roc.curve(data_site_A_RoC_MW)
plot_site_B_roc_MW <- .plot.roc.curve(data_site_B_RoC_MW)
plot_site_C_roc_MW <- .plot.roc.curve(data_site_C_RoC_MW)
plot_site_D_roc_MW <- .plot.roc.curve(data_site_D_RoC_MW)


#----------------------------------------------------------#
# 5. Building the figure  -----
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

