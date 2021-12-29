library(dplyr)

#' Determine the upper limit of the y axis zoom facet 
#' 
#' @param data Preprocessed data
determine_zoom_y_limit <- function(data) {

  proportion_limit <- data %>%
    filter(variant=="B.1.1.529 (Omikron)") %>%
    summarise(proportion = max(proportion))
  
  ceiling(min(1, proportion_limit$proportion*1.5)*100)/100
}
