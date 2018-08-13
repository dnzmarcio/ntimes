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
#'@param measures a character value indicating which measures should be presented: mean.sd, mean.iqr, median.range, and missing.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@return A data frame with mean, standard deviation, median, quantile 25\%,
#'quantile 75\%, minimum, maximum, sample size and missing data for
#'quantitative variables, and frequency, percentage, sample size and
#'missing data for qualitative variables.
#'
#'@examples
#'library(dplyr)
#'library(magrittr)
#'data(iris)
#'
#'iris_nt <- iris %>%
#'  mutate(Species = ql_var(Species,
#'                          from = c("setosa", "versicolor", "virginica"),
#'                          to = c("Setosa", "Versicolor", "Virginica"),
#'                          order = c("Virginica", "Setosa", "Versicolor")))
#'iris_nt %>% nt_describe(group = Species)
#'
#'@export
nt_describe <- function(data,
                        group = NULL,
                        measures = c("mean.sd", "median.iqr",
                                     "median.range", "missing"),
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
               group = group, group.name = group.name, digits = digits,
               measures = measures)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  return(out)
}

aux_describe <- function(var, var.name, group, group.name, digits, measures){

  var.label <- extract_label(var, var.name)
  unit.label <- extract_unit(var)

  if (!is.null(group))
    group.label <- extract_label(group[[1]], group.name)

  if (is.numeric(var)){
    out <- describe_quantitative(var = var,
                                 group = group,
                                 digits = digits,
                                 var.name = var.name,
                                 var.label = var.label,
                                 unit.label = unit.label,
                                 group.label = group.label,
                                 measures = measures)
  } else {
    out <- describe_qualitative(var = var,
                                group = group,
                                digits = digits,
                                var.name = var.name,
                                var.label = var.label,
                                group.label = group.label)
  }
  return(out)
}

#'@importFrom purrr map
#'@importFrom dplyr mutate filter select transmute everything
#'@importFrom tidyr gather nest unnest
#'@importFrom magrittr %>%
describe_quantitative <- function(var, group,
                                  digits,
                                  var.name, var.label,
                                  unit.label, group.label, measures){

  if (is.null(group)) {

    temp <- quantitative_measures(var, digits = digits)
    temp <- temp %>% mutate(variable = var.name) %>%
      select(.data$variable, everything()) %>%
      transmute(variable = .data$variable,
                mean.sd =
                  paste0(.data$mean, " \U00b1 ", .data$sd),
                median.q25.q75 =
                  paste0(.data$median,
                         " (", .data$q25, " ; ", .data$q75, ")"),
                median.min.max =
                  paste0(.data$median,
                         " (", .data$min, " ; ", .data$max, ")"),
                n = .data$n, missing = .data$missing)

    out <- format_quantitative(temp, group = FALSE,
                               var.label = var.label,
                               unit.label = unit.label,
                               measures = measures)

  } else {
    temp <- data_frame(!!var.name := var, g = fct_drop(group[[1]])) %>%
      gather(key = "key", value = "value", -.data$g) %>%
      filter(!is.na(.data$g)) %>% nest(.data$value) %>%
      mutate(desc = map(.data$data,
                        ~ quantitative_measures(x = .$value,
                                                digits = digits))) %>%
      unnest(.data$desc, .drop = TRUE) %>%
      mutate(variable = .data$key, group = .data$g) %>%
      select(.data$variable, .data$group, everything()) %>%
      select(-.data$key, -.data$g) %>%
      transmute(variable = .data$variable,
                group = .data$group,
                mean.sd =
                  paste0(.data$mean, " \U00b1 ", .data$sd),
                median.q25.q75 =
                  paste0(.data$median,
                         " (", .data$q25, " ; ", .data$q75, ")"),
                median.min.max =
                  paste0(.data$median,
                         " (", .data$min, " ; ", .data$max, ")"),
                n = .data$n, missing = .data$missing)
    temp <- split(temp, temp$group)
    temp <- map(.x = temp, .f = format_quantitative, group = TRUE,
                var.label = var.label,
                unit.label = unit.label,
                group.label = group.label,
                measures = measures)

    out <- Reduce(function(x, y)
      merge(x, y, by = "Variable", all = TRUE, sort = F),
      temp)
  }

  out <- as.data.frame(out)
  return(out)
}

#'@importFrom stats quantile sd median
#'@importFrom dplyr summarise_all
quantitative_measures <- function(x, digits){

  if (all(is.na(x))){
    x <- data_frame(x = as.numeric(x))
  } else {
    x <- data_frame(x)
  }

  q25 <- function(x, na.rm) quantile(x, probs = 0.25, na.rm = na.rm)
  q75 <- function(x, na.rm) quantile(x, probs = 0.75, na.rm = na.rm)
  sample_size <- function(x, na.rm) length(x)
  missing <- function(x, na.rm) sum(is.na(x))

  funs <- c("mean", "median", "sd", "q25", "q75", "min",
            "max", "sample_size", "missing")

  out <- round(summarise_all(.tbl = x, .funs = funs, na.rm = TRUE), digits)
  rownames(out) <- NULL
  colnames(out)[length(funs) - 1] <- "n"

  return(out)
}

format_quantitative <- function(var, group,
                                var.label, unit.label, group.label = NULL,
                                measures){

  var.label <- ifelse(unit.label == "", var.label,
                      paste0(var.label, " (", unit.label, ")"))

  if (length(measures) > 1){
    aux_variable <- c(var.label, NA)
    aux_measures <- c("", NA)
    index <- 2

    if (any(measures == "mean.sd")) {
      aux_variable[index] <- paste("  Mean", "SD", sep = " \U00b1 ")
      aux_measures[index] <- var$mean.sd
      index <- index + 1
    }
    if (any(measures == "median.iqr")) {
      aux_variable[index] <- "  Median (Q25% ; Q75%)"
      aux_measures[index] <- var$median.q25.q75
      index <- index + 1
    }
    if (any(measures == "median.range")) {
      aux_variable[index] <- "  Median (Min ; Max)"
      aux_measures[index] <- var$median.min.max
      index <- index + 1
    }
    if (any(measures == "missing")) {
      aux_variable[index] <- "  Missing"
      aux_measures[index] <- unique(var$missing)
    }
  } else {
    aux_variable <- var.label
    aux_measures <-  switch(measures,
                            "mean.sd" = var$mean.sd,
                            "median.iqr" = var$median.q25.q75,
                            "median.range" = var$median.min.max,
                            "missing" = unique(var$missing))
  }

  out <- data_frame(Variable = aux_variable, Measure = aux_measures)

  if (group){
    colnames(out)[2] <- paste0(group.label, ":", var$group,
                               " (n = ", var$n, ")")
  } else {
    colnames(out)[2] <- paste0("(n = ", var$n, ")")
  }
  return(out)
}

#'@importFrom purrr map
#'@importFrom dplyr mutate filter select transmute everything
#'@importFrom tidyr gather nest unnest
#'@importFrom magrittr %>%
describe_qualitative <- function(var, group = NULL,
                                 digits = 2,
                                 var.name, var.label, group.label){

  lv <- levels(as.factor(var))

  if (is.null(group)) {
    temp <- qualitative_measures(h = var, digits = digits, levels = lv)
    temp <- temp %>% mutate(variable = var.name) %>%
      select(.data$variable, everything()) %>%
      select(category = .data$category, percent = .data$perc,
             frequency = .data$freq, missing = .data$missing,
             n = .data$n)

    out <- format_qualitative(temp, group = FALSE, var.label = var.label)

  } else {
    temp <- data_frame(!!var.name := var, g = fct_drop(group[[1]])) %>%
      gather(key = "key", value = "value", -.data$g) %>%
      filter(!is.na(.data$g)) %>% nest(.data$value) %>%
      mutate(desc = map(.data$data,
                        ~ qualitative_measures(h = .$value,
                                               digits = digits,
                                               levels = lv))) %>%
      unnest(.data$desc, .drop = TRUE) %>%
      select(group = .data$g, variable = .data$key, category = .data$category,
             percent = .data$perc, frequency = .data$freq,
             .data$n, missing = .data$missing)
    temp <- split(temp, temp$group)
    temp <- map(.x = temp, .f = format_qualitative, group = TRUE,
                var.label = var.label,
                group.label = group.label)

    out <- Reduce(function(x, y)
      merge(x, y, by = "Variable", all = TRUE, sort = F), temp)
  }

  out <- as.data.frame(out)
  return(out)
}

#'@importFrom forcats fct_explicit_na fct_count
#'@importFrom dplyr mutate filter select
qualitative_measures <- function(h, digits, levels){
  h <- factor(h, levels = levels)
  h <- fct_explicit_na(h, na_level = "Missing")

  out <- fct_count(h) %>%
      mutate(freq = .data$n, n = sum(.data$freq),
             missing = sum(.data$freq[.data$f == "Missing"])) %>%
      filter(.data$f != "Missing") %>%
      mutate(perc = round(100*prop.table(.data$freq), digits)) %>%
      select(category = .data$f, .data$perc, .data$freq, .data$n, .data$missing)

  return(out)
}


format_qualitative <- function(var, group,
                               var.label, group.label = NULL){

  aux_variable <- c(var.label, paste(" ", as.character(var$category)), "Missing")
  aux_measure <- c("", paste0(var$frequency, " (", var$percent, ")"),
                   unique(var$missing))

  out <- data_frame(Variable = aux_variable, Measure = aux_measure)

  if (group){
    colnames(out)[2] <- paste0(group.label, ":", unique(var$group),
                               " (n = ", unique(var$n), ")")
  } else {
    colnames(out)[2] <- paste0("(n = ", unique(var$n), ")")
  }
  return(out)

}



