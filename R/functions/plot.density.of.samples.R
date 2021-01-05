.plot.density.of.samples <- function(data){
  
  figure_density <- 
    data$age_data %>%
    filter(age < age_lim+1e3) %>%
    ggplot(
      aes(x=age))+
    
    geom_hline(
      yintercept = c(0,3e-4),
      color = gray_light,
      size = line_size)+
    
    geom_vline(
      xintercept = seq(from = 0,to = age_lim, by=2e3),
      color = gray_light,
      size = line_size)+
    
    geom_density(
      color = gray_dark,
      fill = gray_light)+
    
    geom_rug(sides = "b")+
    
    scale_x_continuous(trans = "reverse")+
    scale_y_continuous(breaks = c(0,3e-4))+
    
    coord_flip(
      xlim = c(age_lim, 0),
      ylim = c(0, 3e-4))+  
  
  theme(
    line = element_line(size = line_size),
    text = element_text(size = text_size))+
    
    labs(x="Age (cal yr BP)",
         y="Density of samples")
  
  return(figure_density)
      
}

