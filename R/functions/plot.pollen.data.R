.plot.pollen.data <- function(data, common_list, strip_names = F){
  
  # Create a custom pallete
  .get.palette <- colorRampPalette(brewer.pal(8, "Set2"))
  color_legen_pollen_taxa <- .get.palette(length(common_list))
  names(color_legen_pollen_taxa)<- sort(common_list)  
  
  # plot
  pollen_plot <-  
    data %>%
    bind_rows(
      .,tibble(
        sample.id = NA,
        name = common_list,
        value = 0,
        age = 10e3,
        newage = NA)) %>%
    ggplot(
      aes( 
        y = value, 
        x = age))+
    
    geom_hline(
      yintercept = c(0,1),
      color = gray_light,
      size = line_size) +
    
    geom_vline(
      xintercept = seq(from = 0, to = age_lim, by=2000),
      color = gray_light,
      size = line_size) +
    
    geom_ribbon(
      aes(
        ymin = rep(0, length(value)),
        ymax = value, 
        fill = name), 
      color = gray_dark,
      alpha = 1,
      size = line_size)+
    
    facet_wrap( ~ name, 
                ncol = length(common_list))+
    
    scale_x_continuous(trans = "reverse")+
    scale_y_continuous(breaks = c(0,1))+
    
    coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
    
    scale_fill_manual("pollen taxa",values = color_legen_pollen_taxa)+
    
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(angle = 45),
          line = element_line(size = line_size),
          text = element_text(size = text_size))+
    
    labs(
      x="Age (cal yr BP)",
      y="Proportion of total pollen")
  
  if(strip_names == F){
    pollen_plot <-
      pollen_plot + 
      theme(
        strip.background = element_blank(),
        strip.text = element_blank()
      )
  }
  
  return(pollen_plot)
  
}


