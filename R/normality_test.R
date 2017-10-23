#'Normality tests
#'
#'@importFrom purrr map2
#'@importFrom tidyr gather nest unnest
#'@importFrom dplyr mutate
#'@importFrom tibble data_frame as_data_frame
#'
#'@description Perform Anderson-Darling, Shapiro-Francia,
#'Kolmogorov-Smirnov, Cramer-vonMises, and Pearson normality tests
#'for several variables. In addition, it also plots a Quantile-Quantile plot.
#'
#'@param data a data frame with the variables.
#'@param group an optional character indicating the group variable.
#'@param test a character value specifying the normality test to be performed.
#'The options are Anderson-Darling (\code{ad}), Shapiro-Francia (\code{"sf"}),
#'Kolmogorov-Smirnov (\code{ks}),  Cramer-vonMises (\code{cvm}) and
#'Pearson (\code{ps}). The default is Shapiro-Francia (\code{"sf"}).
#'@param digits a integer value indicating the numer of decimal places.
#'@param pvalue.plot a logical value indicating if the p-value
#'should be presented in the Quantile-Quantile plot.
#'
#'@details The function is a wrapper of \link[nortest]{ad.test},
#'\link[nortest]{sf.test}, \link[nortest]{lillie.test},
#'\link[nortest]{cvm.test}, and \link[nortest]{pearson.test}.
#'
#'@return \code{tab} a table of p-values for each of the five normality tests.
#'@return \code{plot} a Quantile-Quantile plot.
#'
#'@examples
#'nt_norm_test(iris, group = Species)
#'
#'@export
nt_norm_test <- function(data, group = NULL, test = "sf",
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
  vars.label <- unlist(map2(vars, vars.name, extract_label))

  plot <- list()

  if (!is.null(group)){
    tab <- vars %>%
      gather(key = "Variable", value = "value") %>%
      nest(-.data$Variable) %>% mutate(Variable = vars.label) %>%
      mutate(tab = map(.data$data,
                       ~ norm_test(var = .$value, group = group[[1]],
                                   test = test, digits = digits))) %>%
      unnest(tab, .drop = TRUE)
  } else {
    tab <- vars %>%
      gather(key = "Variable", value = "value") %>%
      nest(-.data$Variable) %>% mutate(Variable = vars.label) %>%
      mutate(tab = map(.data$data,
                       ~ norm_test(var = .$value, group = NULL,
                                   test = test, digits = digits))) %>%
      unnest(tab, .drop = TRUE)
  }

  # for (i in 1:nv){
  #
  #   if (!is.null(group)){
  #     plot[[paste(vars.name[i])]] <-
  #       norm_plot(var = vars[[i]], group = group[[1]],
  #                 test = test,
  #                 var.label = vars.label[[i]],
  #                 group.label = group.label[[1]],
  #                 digits = digits,
  #                 pvalue.plot  = pvalue.plot)
  #   } else {
  #     plot[[paste(vars.name[i])]] <-
  #       norm_plot(var = vars[[i]], group = NULL,
  #                 test = test,
  #                 var.label = vars.label[[i]],
  #                 group.label = group.label[[1]],
  #                 digits = digits,
  #                 pvalue.plot  = pvalue.plot)
  #   }
  # }

  out <- list(pvalue = tab, plot = plot)
  return(out)
}

#'@importFrom nortest ad.test sf.test lillie.test cvm.test pearson.test
#'@importFrom dplyr select
#'@importFrom purrr map
#'@importFrom tidyr nest unnest
#'@importFrom tibble data_frame
norm_test <-  function(var, group = NULL,
                       test,
                       digits = 3){

  if (!is.numeric(var))
    stop("'var' contains non-numeric variable.")

  aux_norm_test <- function(var, test){

    if (test == "ad")
      result <- try(ad.test(var), silent = TRUE)
    if (test == "sf")
      result <- try(sf.test(var), silent = TRUE)
    if (test == "ks")
      result <- try(lillie.test(var), silent = TRUE)
    if (test == "cvm")
      result <- try(cvm.test(var), silent = TRUE)
    if (test == "ps")
      result <- try(pearson.test(var), silent = TRUE)

    out <- result$p.value

    return(out)
  }

  if (is.null(group)){
    out <- aux_norm_test(var, test = test)
  } else {
    data.test <- data_frame(var, g = group)
    out <- data.test %>% nest(-.data$g) %>%
      mutate(p = map(.data$data, ~ aux_norm_test(.$var, test = test))) %>%
      unnest(.data$p, drop = TRUE) %>% select(Group = .data$g, `p value` = .data$p)
  }

  return(out)
}

#'@import ggplot2
norm_plot <-  function(var,
                       group = NULL,
                       test = "sf",
                       var.label = NULL,
                       group.label = NULL,
                       digits = 3,
                       pvalue.plot = TRUE){

  if (is.null(group)){
    data_plot <- data.frame(var = (var - mean(var))/sd(var))

    p <- norm_test(var = var, test = test)
    p <- ifelse(round(p[[1]], digits) != 0,
                        paste0("= ", round(p[[1]], digits)), "< 0.001")
    pvalue <- paste("p-value", p)

    out <- ggplot(data_plot, aes_string(sample = 'var')) + stat_qq() +
      geom_abline(intercept = 0, slope = 1) + theme_bw()
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

    p <- norm_test(var = var, test = test, group = group) %>% dplyr::select(-group)
    p <- ifelse(round(p[[1]], digits) != 0,
                paste("=", round(p[[1]], digits)), "< 0.001")
    pvalue <- paste("p-value", p)

    out <- ggplot(data_test, aes_string(sample = 'var')) + stat_qq() +
      geom_abline(intercept = 0, slope = 1) + theme_bw() +
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


