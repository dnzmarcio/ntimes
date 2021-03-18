#'Normality tests

#'
#'@description Perform Anderson-Darling, Shapiro-Francia,
#'Kolmogorov-Smirnov, Cramer-vonMises, and Pearson normality tests
#'for several variables. In addition, it also plots a Quantile-Quantile plot and a histogram.
#'
#'@param data a data frame with the variables.
#'@param group an optional character indicating the group variable.
#'@param norm.test a function with a numeric vector as input and a list as output containing an object named \code{p.value} similar to \link[ntimes]{helper_sf_test}.
#'
#'@return \code{tab} a table of p-values for each of the five normality tests.
#'@return \code{plot} a Quantile-Quantile plot.
#'
#'@examples
#'library(magrittr)
#'data(iris)
#'
#'iris %>% select(-Species) %>% nt_norm_test()
#'iris %>% nt_norm_test(group = Species)
#'
#'@importFrom purrr map2
#'@importFrom tidyr pivot_longer
#'@importFrom dplyr mutate bind_cols
#'
#'@export
nt_norm_test <- function(data, group = NULL, norm.test = helper_sf_test){

  group <- enquo(group)

  if (!quo_is_null(group)){
    vars <- select(.data = data, -{{group}})
    group.var <- select(.data = data, {{group}})
    group.name <- names(group.var)
    group.label <- map2(group.var, group.name, extract_label)
  } else {
    vars <-  data
    group.var <- NULL
    group.name <- NULL
    group.label <- NULL
  }

  vars.name <- names(vars)
  vars.label <- map2(vars, vars.name, extract_label)

  if (!is.null(group.var)){
    tab <- data %>%
      pivot_longer(cols = -{{group}}, names_to = "Variable", values_to = "value") %>%
      group_by(.data$Variable, {{group}}) %>%
      summarise(p.value = norm.test(.data$value)$p.value) %>%
      rename(Group = {{group}}, 'p value' = .data$p.value)

  } else {
    tab <- vars %>%
      pivot_longer(cols = everything(), names_to = "Variable", values_to = "value") %>%
      group_by(.data$Variable) %>%
      summarise(p.value = norm.test(.data$value)$p.value) %>%
      rename('p value' = .data$p.value)
  }

  qq_plot <- map2(vars, vars.label, qq_plot, group = group.var[[1]],
                  group.label = group.label[[1]])

  hist <- map2(vars, vars.label, hist_plot, group = group.var[[1]],
               group.label = group.label[[1]])

  out <- list(tab = tab, qq_plot = qq_plot, hist = hist)
  return(out)
}




#'@import ggplot2
#'@importFrom stats qnorm
#'@importFrom dplyr summarise
qq_plot <-  function(var, var.label,
                     group = NULL, group.label = NULL){

  qqline_slope <- function(var){
    y <- quantile(var, c(0.25, 0.75), na.rm = TRUE)
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
  }

  qqline_int <- function(var){
    y <- quantile(var, c(0.25, 0.75), na.rm = TRUE)
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1] - slope * x[1]
  }


  if (is.null(group)){
    data_plot <- data.frame(var = var)
    slope <- qqline_slope(var)
    int <- qqline_int(var)

    out <- ggplot(data_plot, aes_string(sample = 'var')) + stat_qq() +
      geom_abline(slope = slope, intercept = int) + theme_bw()

  } else {
    if (is.null(group.label))
      group.label <- "group"

    group <- paste(group.label, group, sep = ": ")
    group <- as.factor(group)
    data_plot <- data.frame(var, group)
    data_qqline <- data_plot %>% group_by(group) %>%
      summarise(slope = qqline_slope(var),
                int = qqline_int(var))


    out <- ggplot(data_plot, aes_string(sample = 'var')) + stat_qq() +
      geom_abline(data = data_qqline, aes_string(intercept = "int", slope = "slope")) +
      theme_bw() +
      ggtitle(paste0(var.label)) + facet_grid(. ~ group)

  }
  return(out)
}

hist_plot <-  function(var, var.label,
                       group = NULL, group.label = NULL){

  if (is.null(group)){
    data_plot <- data.frame(var = var)

    out <- ggplot(data_plot, aes_string(x = 'var')) +
      geom_histogram() + theme_bw()

  } else {
    if (is.null(group.label))
      group.label <- "group"

    group <- paste(group.label, group, sep = ": ")
    group <- as.factor(group)
    data_plot <- data.frame(var, group)

    out <- ggplot(data_plot, aes_string(x = 'var')) +
      geom_histogram() + theme_bw() +
      facet_grid(. ~ group) +
      ggtitle(paste0(var.label))

  }
  return(out)
}

