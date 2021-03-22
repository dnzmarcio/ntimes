#'Descriptive measures
#'
#'@description It calculates a descriptive table for quantitative and
#'qualitative variables in a publishable format.
#'
#'@param data a data frame with the variables.
#'@param group an optional  with the group variable.
#'@param labels a list of labels with components given by their variable names.
#'@param measures_qt a list of functions to summarise quantitative variables. See more in details.
#'@param measures_ql a list of functions to summarise qualitative variables. See more in details.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'
#'@details For quantitative variables, mean +/- sd, median (quantile 0.25 - quantile 0.75),
#'median (minimum - maximum) and number of missing observations are calculated using the
#'functions \code{\link[ntimes]{helper_mean_sd}}, \code{\link[ntimes]{helper_median_iqr}},
#'\code{\link[ntimes]{helper_median_range}} and \code{\link[ntimes]{helper_missing}}.
#'For qualitative variables, percentage (frequency) is calculated using
#'\code{\link[ntimes]{helper_perc_count}}.
#'
#'@return a data frame with summary for all variables by group.
#'
#'@examples
#'data(iris)
#'
#'iris %>% nt_describe(group = Species,
#'                     labels = list(Sepal.Length = "Sepal Length",
#'                                   Sepal.Width = "Sepal Width",
#'                                   Petal.Length = "Petal Length",
#'                                   Petal.Width = "Petal Width"))
#'
#'@importFrom purrr pmap
#'@importFrom dplyr filter select
#'@importFrom utils write.csv
#'@importFrom magrittr %>%
#'@importFrom rlang := .data quo_is_null enquo
#'@importFrom tibble tibble
#'
#'@export
nt_describe <- function(data,
                        group = NULL,
                        labels = NULL,
                        measures_qt = list(helper_mean_sd,
                                           helper_median_iqr,
                                           helper_median_range,
                                           helper_missing),
                        measures_ql = list(helper_perc_count),
                        digits = 2,
                        save = FALSE,
                        file = "descriptive_analysis"){

  group <- enquo(group)

  if (!quo_is_null(group)){
    vars <- select(.data = data, -!!group)
    group <- select(.data = data, !!group)
    group.name <- names(group)
  } else {
    vars <- data
    group <- NULL
    group.name <- NULL
  }

  vars.name <- names(vars)
  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)
    if (!is.null(group)){
      group <- data_labeller(group, labels)
      group.label <- extract_label(group[[1]], group.name)
    }
  } else {
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)
    if (!is.null(group))
      group.label <- extract_label(group[[1]], group.name)
  }

  temp <- pmap(.l = list(vars, vars.name, vars.label),
               .f = aux_describe,
               group = group[[1]],
               group.name = group.name,
               group.label = group.label,
               digits = digits,
               measures_qt = measures_qt,
               measures_ql = measures_ql)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  attr(out, "ntimes") <- "descriptive"
  return(out)
}

aux_describe <- function(var, var.name, var.label,
                         group, group.name, group.label,
                         digits, measures_qt, measures_ql){

  if (!is.null(group)){
    if (!is.factor(group)){
      group <- as.factor(group)
      warning(paste0("Group variable was transformed into a factor."))
    }
  }

  if (is.numeric(var)){
    out <- describe_quantitative(var = var,
                                 group = group,
                                 digits = digits,
                                 var.label = var.label,
                                 group.label = group.label,
                                 measures_qt = measures_qt)
  } else {
    out <- describe_qualitative(var = var,
                                group = group,
                                digits = digits,
                                var.label = var.label,
                                group.label = group.label,
                                measures_ql = measures_ql)


  }
  return(out)
}

describe_quantitative <- function(var, group,
                                  digits,
                                  var.label, group.label,
                                  measures_qt){

  if (is.null(group)) {
    desc <- quantitative_measures(var,
                                  digits = digits,
                                  measures_qt = measures_qt)
    out <- format_quantitative(desc = desc, group = NULL,
                               var.label = var.label,
                               group.label = group.label)
  } else {

    aux <- function(x, g = NULL){
      out <- format_quantitative(x, group = g,
                                 var.label = var.label,
                                 group.label = group.label)
      return(out)
    }

    desc <- tapply(var, group, quantitative_measures,
                   digits = digits,
                   measures_qt = measures_qt)
    group.lv <- setNames(as.list(levels(group)), levels(group))
    temp <- mapply(aux, desc, group.lv, SIMPLIFY = FALSE)

    out <- Reduce(function(x, y)
      merge(x, y, by = "Variable", all = TRUE, sort = F),
      temp)
  }

  out <- as.data.frame(out)
  return(out)
}

quantitative_measures <- function(x, digits, measures_qt){

  out <- lapply(measures_qt, function(f) f(x, digits = digits))
  out$n <- list(name = " n", value = length(x))

  return(out)
}

format_quantitative <- function(desc,
                                group,
                                var.label,
                                group.label = NULL){

  if (length(desc[-length(desc)]) > 1){
    aux_variable <- c(var.label,
                      paste0("  \t ",
                             Reduce(c, lapply(desc[-length(desc)],
                                              function(x)
                                                x$name))))
    aux_measures <- c("", Reduce(c, lapply(desc[-length(desc)],
                                           function(x)
                                             x$value)))
  } else {
    aux_variable <- var.label
    aux_measures <- desc[[1]]$value
  }

  out <- data.frame(Variable = aux_variable, Measure = aux_measures)

  if (is.null(desc$n$value))
    desc$n$value <- 0

  if (!is.null(group)){
    colnames(out)[2] <- paste0(group.label, ": ", group,
                               " (n = ", desc$n$value, ")")
  } else {
    colnames(out)[2] <- paste0("All (n = ", desc$n$value, ")")
  }
  return(out)
}

describe_qualitative <- function(var, group = NULL,
                                 digits = 2,
                                 var.label, group.label,
                                 measures_ql){

  if (!is.factor(var)){
    var <- as.factor(var)
    warning(paste0(var.label, " was transformed into a factor."))
  }

  lv <- levels(var)

  if (is.null(group)) {
    desc <- qualitative_measures(h = var, digits = digits,
                                 measures_ql = measures_ql)
    out <- format_qualitative(desc = desc, group = NULL,
                              var.label = var.label)

  } else {

    aux <- function(x, g = NULL){
      out <- format_qualitative(x, group = g,
                                var.label = var.label,
                                group.label = group.label)
      return(out)
    }

    desc <- tapply(var, group, qualitative_measures,
                   digits = digits, measures_ql = measures_ql)
    group.lv <- setNames(as.list(levels(group)), levels(group))
    temp <- mapply(aux, desc, group.lv, SIMPLIFY = FALSE)

    out <- Reduce(function(x, y)
      merge(x, y, by = "Variable", all = TRUE, sort = F), temp)
  }

  return(out)
}

qualitative_measures <- function(h, digits, measures_ql){

  aux <- function(f) {
    out <- f(h, digits = digits)
  }
  out <- lapply(measures_ql, aux)
  out$n <- list(name = " n", value = length(h))

  return(out)

}

format_qualitative <- function(desc, group,
                               var.label, group.label){

  if (is.null(desc[[1]]$name)){
    aux_variable <- var.label
  } else {
    aux_variable <- c(var.label, paste0("  \t ",
                                        as.character(desc[[1]]$name)))
  }

  aux_measure <- c("", desc[[1]]$value)

  out <- data.frame(Variable = aux_variable, Measure = aux_measure)


  if (is.null(desc$n))
    desc$n$value <- 0

  if (!is.null(group)){
    colnames(out)[2] <- paste0(group.label, ": ", group,
                               " (n = ", desc$n$value, ")")
  } else {
    colnames(out)[2] <- paste0("All (n = ", desc$n$value, ")")
  }
  return(out)

}



