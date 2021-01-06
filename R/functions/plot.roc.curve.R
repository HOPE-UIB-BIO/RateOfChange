.plot.roc.curve <- function(data, roc_max = 2){
  
  data_peak <-
    data %>% 
    filter(Peak == T)
  
  roc_plot <- 
    data %>%
    ggplot(
      aes(
        y =ROC, 
        x = Age))+
    
    geom_vline(
      xintercept = seq(from = 0,to = age_lim, by=2000),
      color = gray_light,
      size = line_size) +
    
    geom_ribbon(
      aes(
        ymin = ROC_dw, 
        ymax = ROC_up),
      alpha=1/2,
      color = gray_light,
      fill = gray_light) +
    
    geom_line(
      alpha = 1,
      size = line_size) +
    
    geom_point(
      data = data_peak,
      color = "green",
      size = 2,
      shape=16,
      alpha=2/3) +
    
    geom_hline(
      yintercept = 0,
      color = gray_dark,
      size = line_size)+
    
    scale_x_continuous(trans = "reverse")+
    coord_flip(
      xlim = c(age_lim, 0),
      ylim = c(0, roc_max)) +
    
    theme(
      legend.position = "none",
      text = element_text(size = text_size),
      line = element_line(size = line_size)) +
    
    labs(
      x = "Age (cal yr BP)",
      y = "Rate-of-Change score")
    
  return(roc_plot)
}


