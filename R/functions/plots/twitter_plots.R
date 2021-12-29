#' Export a plot with dimensions and DPI optimized for usage on Twitter
#' 
#' @param data Preprocssed data
#' @param plot Ggplot object to export
#' @param filename Name of the file to use for storage in the output/ folder
#' @param subtitle Subtitle to add to the plot using [ggplot2::ggtitle()]
#' 
#' @return Relative output path of the stored file
export_plot_with_title <- function(data, plot, filename="austria-variants.png", subtitle="") {
  p_enriched <- plot +
    ggtitle("SARS-CoV-2 Variants - Austria\nSequencing and targeted PCR combined", subtitle = subtitle) +
    labs(caption = glue::glue("Data up to {max(data$Kalenderwoche)}
         Source: AGES, https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/"))
  
  output_path <- paste0("output/", filename)
  ggsave(plot = p_enriched, output_path, width = 10, height = 5)
  
  output_path
}


export_plot_quasibinomial_fit_ci_zoom <- function(data, plot) {
  p_enriched <- plot +
    ggtitle("SARS-CoV-2 Variants - Austria", subtitle = "Sequencing and targeted PCR combined\nOverdispersed binomial GLM fit") +
    labs(x="Date, week ending",
         y="Proportion of cases with known/specified variant",
         color=NULL, 
         fill=NULL, 
         caption = glue::glue("Data up to {max(data$Kalenderwoche)}
         Source: AGES, https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/
       Shaded area showing 95% CI of quasibinomial model"))
  
  output_path <- "output/austria-variants-binomial-glm.png"
  ggsave(plot = p_enriched, output_path, width = 10, height = 5)
  
  output_path
}

export_plot_log_odds <- function(data, plot) {
  p_enriched <- plot +
    ggtitle("SARS-CoV-2 Variants - Austria", subtitle = "Sequencing and targeted PCR combined\nOverdispersed binomial GLM fit")
  
  output_path <- "output/austria-variants-binomial-glm-log-odds.png"
  ggsave(plot = p_enriched, output_path, width = 10, height = 5)
  
  output_path
}
