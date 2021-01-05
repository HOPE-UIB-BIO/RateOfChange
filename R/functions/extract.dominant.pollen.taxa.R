.extract.dominant.pollen.taxa <- function(data, N_species = 5){
  
  cmn_df <- data$community_data[ ,-1]
  cmn_df_proportion <- cmn_df/rowSums(cmn_df)*100
  species_list <- 
    colSums(cmn_df_proportion) %>%
    sort(decreasing =T)
  sel_taxa_names <- names(species_list)[1:N_species]
  
  return(sel_taxa_names)
}
