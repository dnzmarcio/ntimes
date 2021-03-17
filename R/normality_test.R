#'Normality tests
#'
#'@importFrom purrr map2
#'@importFrom tidyr gather nest unnest
#'@importFrom dplyr mutate
#'@importFrom tibble data_frame as_data_frame
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
#'@export
nt_norm_test <- function(data, group = NULL, test = nt_norm_test(test = "sf"),
                         digits = 3, pvalue.plot = TRUE){

  data <- as_data_frame(data)
  group <- enquo(group)

  if (!quo_is_null(group)){
    vars <- select(.data = data, -!!group)
    group <- select(.data = data, !!group)
    group.name <- names(group)
    group.label <- extract_label(group, group.name)
  } else {
    vars <-  data
    group <- NULL
    group.name <- NULL
    group.label <- NULL
  }

  vars.name <- names(vars)
  vars.label <- map2(vars, vars.name, extract_label)

  plot <- list()

  if (!is.null(group)){
    tab <- vars %>%
      gather(key = "Variable", value = "value") %>%
      nest(-.data$Variable) %>%
      mutate(Variable = unlist(vars.label)) %>%
      mutate(p.value = map(.data$data,
                       ~ norm_test(var = .$value, group = group[[1]],
                                   test = test, digits = digits))) %>%
      unnest(.data$p.value, .drop = TRUE)
  } else {
    tab <- vars %>%
      gather(key = "Variable", value = "value") %>%
      nest(-.data$Variable) %>%
      mutate(Variable = unlist(vars.label)) %>%
      mutate(p.value = map(.data$data,
                       ~ norm_test(var = .$value, group = NULL,
                                   test = test, digits = digits))) %>%
      unnest(.data$p.value, .drop = TRUE)
  }

  qq_plot <- map2(vars, vars.label, qq_plot, group = group[[1]],
                  test = test,
                  group.label = group.label[[1]],
                  digits = digits,
                  pvalue.plot  = pvalue.plot)

  hist <- map2(vars, vars.label, hist_plot, group = group[[1]],
               test = test,
               group.label = group.label[[1]],
               digits = digits,
               pvalue.plot  = pvalue.plot)

  out <- list(pvalue = tab, qq_plot = qq_plot, hist = hist)
  return(out)
}




#'@import ggplot2
#'@importFrom stats qnorm
#'@importFrom dplyr summarise
qq_plot <-  function(var, var.label,
                     group = NULL, group.label = NULL,
                     test = nt_norm_test(test = "sf"),
                     digits = 3, pvalue.plot = TRUE){

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

    p <- norm_test(var = var, test = test)
    p <- ifelse(round(p[[1]], digits) != 0,
                        paste0("= ", round(p[[1]], digits)), "< 0.001")
    pvalue <- paste("p-value", p)

    out <- ggplot(data_plot, aes_string(sample = 'var')) + stat_qq() +
      geom_abline(slope = slope, intercept = int) + theme_bw()
    if (pvalue.plot)
      out <- out +
      annotate(geom = "text", label = pvalue,
               x = -Inf, y = Inf,
               hjust = -0.5, vjust = 2, size = 4)

  } else {
    if (is.null(group.label))
      group.label <- "group"

    group <- paste(group.label, group, sep = ": ")
    group <- as.factor(group)
    data_test <- data.frame(var, group)
    data_qqline <- data.frame(var, group) %>% group_by(group) %>%
      summarise(slope = qqline_slope(var),
                int = qqline_int(var))

    p <- norm_test(var = var, group = group, test = test) %>%
      select(-.data$Group)
    p <- ifelse(round(p[[1]], digits) != 0,
                paste("=", round(p[[1]], digits)), "< 0.001")
    pvalue <- paste("p-value", p)

    out <- ggplot(data_test, aes_string(sample = 'var')) + stat_qq() +
      geom_abline(data = data_qqline, aes(intercept = int, slope = slope)) + theme_bw() +
      facet_grid(. ~ group) +
      ggtitle(paste0(var.label))
    if (pvalue.plot)
      out <- out +
        annotate(geom = "text", label = pvalue,
                 x = -Inf, y = Inf,
                 hjust=-0.5, vjust=2, size = 4)

  }
  return(out)
}

hist_plot <-  function(var, var.label,
                       group = NULL, group.label = NULL,
                       test = "sf", digits = 3, pvalue.plot = TRUE){

  if (is.null(group)){
    data_plot <- data.frame(var = var)

    p <- norm_test(var = var, test = test)
    p <- ifelse(round(p[[1]], digits) != 0,
                paste0("= ", round(p[[1]], digits)), "< 0.001")
    pvalue <- paste("p-value", p)

    out <- ggplot(data_plot, aes_string(x = 'var')) +
      geom_histogram() + theme_bw()
    if (pvalue.plot)
      out <- out +
      annotate(geom = "text", label = pvalue,
               x = -Inf, y = Inf,
               hjust = -0.5, vjust = 2, size = 4)

  } else {
    if (is.null(group.label))
      group.label <- "group"

    group <- paste(group.label, group, sep = ": ")
    group <- as.factor(group)
    data_test <- data.frame(var, group)

    p <- norm_test(var = var, group = group, test = test) %>%
      select(-.data$Group)
    p <- ifelse(round(p[[1]], digits) != 0,
                paste("=", round(p[[1]], digits)), "< 0.001")
    pvalue <- paste("p-value", p)

    out <- ggplot(data_test, aes_string(x = 'var')) +
      geom_histogram() + theme_bw() +
      facet_grid(. ~ group) +
      ggtitle(paste0(var.label))
    if (pvalue.plot)
      out <- out +
      annotate(geom = "text", label = pvalue,
               x = -Inf, y = Inf,
               hjust=-0.5, vjust=2, size = 4)

  }
  return(out)
}

