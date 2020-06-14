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
#'@importFrom tibble data_frame as_data_frame
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
                        measures = list(nt_mean_sd,
                                        nt_median_iqr,
                                        nt_median_range,
                                        nt_missing),
                        digits = 2,
                        save = FALSE,
                        file = "descriptive_analysis"){

  data <- as_data_frame(data)
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
               measures = measures)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  return(out)
}

#'@importFrom stats sd
#'@export
nt_mean_sd <- function(var, digits){
  mean <- round(mean(var, na.rm = TRUE), digits)
  sd <- round(sd(var, na.rm = TRUE), digits)
  name <- paste0("  Mean", " \U00b1 ", "SD")
  measure <- paste0(mean, " \U00b1 ", sd)
  out <- list(name = name, measure = measure)
  return(out)
}

#'@importFrom stats median quantile
#'@export
nt_median_iqr <- function(var, digits){
  median <- round(median(var, na.rm = TRUE), digits)
  q25 <- round(quantile(var, probs = 0.25, na.rm = TRUE), digits)
  q75 <- round(quantile(var, probs = 0.75, na.rm = TRUE), digits)
  name <- "  Median (Q25% ; Q75%)"
  measure <- paste0(median," (", q25, " ; ", q75, ")")
  out <- list(name = name, measure = measure)
  return(out)
}

#'@export
nt_missing <- function(var, digits){
  out <- list(name = " Missing", measure = sum(is.na(var)))
  return(out)
}

#'@importFrom stats median
#'@export
nt_median_range <- function(var, digits){
  median <- round(median(var, na.rm = TRUE), digits)
  min <- round(min(var, na.rm = TRUE), digits)
  max <- round(max(var, na.rm = TRUE), digits)
  name <- "  Median (Min ; Max)"
  measure <- paste0(median," (", min, " ; ", max, ")")
  out <- list(name = name, measure = measure)
  return(out)
}

aux_describe <- function(var, var.name, group, group.name, digits, measures){

  var.label <- extract_label(var, var.name)
  unit.label <- extract_unit(var)

  if (!is.null(group))
    group.label <- extract_label(group, group.name)

  if (is.numeric(var)){
    out <- describe_quantitative(var = var,
                                 group = group,
                                 digits = digits,
                                 var.label = var.label,
                                 unit.label = unit.label,
                                 group.label = group.label,
                                 measures = measures)
  } else {
    out <- describe_qualitative(var = var,
                                group = group,
                                digits = digits,
                                var.label = var.label,
                                group.label = group.label)


  }
  return(out)
}

describe_quantitative <- function(var, group,
                                  digits,
                                  var.label, unit.label, group.label,
                                  measures){

  if (is.null(group)) {
    desc <- quantitative_measures(var,
                                  digits = digits,
                                  measures = measures)
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
                   measures = measures)
    group.lv <- setNames(as.list(levels(group)), levels(group))
    temp <- mapply(aux, desc, group.lv, SIMPLIFY = FALSE)

    out <- Reduce(function(x, y)
      merge(x, y, by = "Variable", all = TRUE, sort = F),
      temp)
  }

  out <- as.data.frame(out)
  return(out)
}

quantitative_measures <- function(x, digits, measures){

  out <- lapply(measures, function(f) f(x, digits = digits))
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
                      Reduce(c, lapply(desc[-length(desc)],
                                       function(x) x$name)))
    aux_measures <- c("", Reduce(c, lapply(desc[-length(desc)],
                                           function(x) x$measure)))
  } else {
    aux_variable <- var.label
    aux_measures <- desc[[1]]$measure
  }

  out <- data.frame(Variable = aux_variable, Measure = aux_measures)

  if (!is.null(group)){
    colnames(out)[2] <- paste0(group.label, ":", group,
                               " (n = ", desc$n$measure, ")")
  } else {
    colnames(out)[2] <- paste0("(n = ", desc$n$measure, ")")
  }
  return(out)
}

describe_qualitative <- function(var, group = NULL,
                                 digits = 2,
                                 var.label, group.label){

  if (!is.factor(var))
    stop(paste0(var.label, " is not a factor."))

  lv <- levels(var)

  if (is.null(group)) {
    desc <- qualitative_measures(h = var, digits = digits)
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
                   digits = digits)
    group.lv <- setNames(as.list(levels(group)), levels(group))
    temp <- mapply(aux, desc, group.lv, SIMPLIFY = FALSE)

    out <- Reduce(function(x, y)
      merge(x, y, by = "Variable", all = TRUE, sort = F), temp)
  }

  return(out)
}

#'@importFrom forcats fct_explicit_na
qualitative_measures <- function(h, digits){
  h <- fct_explicit_na(h, na_level = "Missing")
  levels <- levels(h)
  count <- tapply(h, h, length)
  n <- length(h)
  missing <- count[length(count)]
  perc <- round(100*prop.table(count), digits)
  perc <- ifelse(!is.finite(perc), NA, perc)

  perc_count <- paste0(perc, " (", count, ")")

  out <- list(levels = levels, perc_count = perc_count, n = n)

  return(out)
}

format_qualitative <- function(desc, group,
                               var.label, group.label = NULL){

  aux_variable <- c(var.label, paste(" ", as.character(desc$levels)))
  aux_measure <- c("", desc$perc_count)

  out <- data.frame(Variable = aux_variable, Measure = aux_measure)

  if (!is.null(group)){
    colnames(out)[2] <- paste0(group.label, ":", group,
                               " (n = ", desc$n, ")")
  } else {
    colnames(out)[2] <- paste0("(n = ", desc$n, ")")
  }
  return(out)

}



