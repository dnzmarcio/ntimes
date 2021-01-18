#'Descriptive measures
#'
#'@description Calculating a descriptive table for quantitative and
#'qualitative variables in the usual format to present in scientific article.
#'
#'@importFrom purrr map2
#'@importFrom dplyr filter select
#'@importFrom utils write.csv
#'@importFrom magrittr %>%
#'@importFrom rlang := .data quo_is_null enquo
#'@importFrom tibble tibble
#'
#'@param data a data frame with the variables.
#'@param group an optional  with the group variable.
#'@param measures a list of functions to summarise quantitative variables. See more in details.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@return a data frame with summary for all variables by group.
#'
#'@examples
#'data(iris)
#'
#'iris %>% nt_describe(group = Species)
#'
#'@export
nt_describe <- function(data,
                        group = NULL,
                        measures_qt = list(nt_mean_sd,
                                        nt_median_iqr,
                                        nt_median_range,
                                        nt_missing),
                        measures_ql = list(nt_perc_count),
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

  temp <- map2(.x = vars, .y = vars.name, .f = aux_describe,
               group = group[[1]], group.name = group.name, digits = digits,
               measures_qt = measures_qt,
               measures_ql = measures_ql)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  attr(out, "ntimes") <- "descriptive"
  return(out)
}

#'@importFrom stats sd
#'@export
nt_mean_sd <- function(var, digits){
  mean <- format(round(mean(var, na.rm = TRUE), digits), nsmall = digits)
  sd <- format(round(sd(var, na.rm = TRUE), digits), nsmall = digits)
  name <- paste0("Mean", " \U00b1 ", "SD")
  measure <- paste0(mean, " \U00b1 ", sd)
  out <- list(name = name, measure = measure)
  return(out)
}

#'@importFrom stats median quantile
#'@export
nt_median_iqr <- function(var, digits){
  median <- format(round(median(var, na.rm = TRUE), digits),
                   nsmall = digits)
  q25 <- format(round(quantile(var, probs = 0.25, na.rm = TRUE), digits),
                nsmall = digits)
  q75 <- format(round(quantile(var, probs = 0.75, na.rm = TRUE), digits),
                nsmall = digits)
  name <- "Median (Q25% ; Q75%)"
  measure <- paste0(median," (", q25, " ; ", q75, ")")
  out <- list(name = name, measure = measure)
  return(out)
}

#'@importFrom stats median
#'@export
nt_median_range <- function(var, digits){
  median <- format(round(median(var, na.rm = TRUE), digits),
                   nsmall = digits)
  min <- format(round(min(var, na.rm = TRUE), digits),
                nsmall = digits)
  max <- format(round(max(var, na.rm = TRUE), digits),
                nsmall = digits)
  name <- "Median (Min ; Max)"
  measure <- paste0(median," (", min, " ; ", max, ")")
  out <- list(name = name, measure = measure)
  return(out)
}

#'@export
nt_missing <- function(var, digits){
  out <- list(name = "Missing", measure = sum(is.na(var)))
  return(out)
}

#'@importFrom forcats fct_explicit_na
#'@export
nt_perc_count <- function(var, digits){
  h <- fct_explicit_na(var, na_level = "Missing")
  lh <- levels(var)

  count <- tapply(h, h, length)
  count <- ifelse(is.na(count), 0, count)
  n <- length(h)
  perc <- 100*prop.table(count)
  perc <- ifelse(!is.finite(perc), NA, format(round(perc, digits), nsmall = digits))

  perc_count <- paste0(perc, " (", count, ")")

  if (!("Missing" %in% lh)){
    lh <- c(lh, "Missing")
    perc_count <- c(perc_count, "0 (0)")
  }

  out <- list(name= lh, measure = perc_count)
}

aux_describe <- function(var, var.name, group, group.name,
                         digits, measures_qt, measures_ql){

  var.label <- extract_label(var, var.name)
  unit.label <- extract_unit(var)

  if (!is.null(group)){
    group.label <- extract_label(group, group.name)

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
                                 unit.label = unit.label,
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
                                  var.label, unit.label, group.label,
                                  measures_qt){

  if (is.null(group)) {
    desc <- quantitative_measures(var,
                                  digits = digits,
                                  measures_qt = measures_qt)
    out <- format_quantitative(desc = desc, group = NULL,
                               var.label = var.label,
                               unit.label = unit.label,
                               group.label = group.label)
  } else {

    aux <- function(x, g = NULL){
      out <- format_quantitative(x, group = g,
                                 var.label = var.label,
                                 unit.label = unit.label,
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
  out$n <- list(name = " n", measure = length(x))

  return(out)
}

format_quantitative <- function(desc,
                                group,
                                var.label,
                                unit.label,
                                group.label = NULL){

  var.label <- ifelse(unit.label == "", var.label,
                      paste0(var.label, " (", unit.label, ")"))

  if (length(desc[-length(desc)]) > 1){
    aux_variable <- c(var.label,
                      paste0("  \t ", Reduce(c, lapply(desc[-length(desc)],
                                       function(x) x$name))))
    aux_measures <- c("", Reduce(c, lapply(desc[-length(desc)],
                                           function(x) x$measure)))
  } else {
    aux_variable <- var.label
    aux_measures <- desc[[1]]$measure
  }

  out <- data.frame(Variable = aux_variable, Measure = aux_measures)

  if (is.null(desc$n$measure))
    desc$n$measure <- 0

  if (!is.null(group)){
    colnames(out)[2] <- paste0(group.label, ": ", group,
                               " (n = ", desc$n$measure, ")")
  } else {
    colnames(out)[2] <- paste0("All (n = ", desc$n$measure, ")")
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

  out <- lapply(measures_ql, function(f) f(h, digits = digits))
  out$n <- list(name = " n", measure = length(h))

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

  aux_measure <- c("", desc[[1]]$measure)

  out <- data.frame(Variable = aux_variable, Measure = aux_measure)


  if (is.null(desc$n))
    desc$n$measure <- 0

  if (!is.null(group)){
    colnames(out)[2] <- paste0(group.label, ": ", group,
                               " (n = ", desc$n$measure, ")")
  } else {
    colnames(out)[2] <- paste0("All (n = ", desc$n$measure, ")")
  }
  return(out)

}



