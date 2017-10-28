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
#'@param measures a character vector indicating which measures should be presented.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@return A data frame with mean, standard deviation, median, quantile 25\%,
#'quantile 75\%, minimum, maximum, sample size and missing data for
#'quantitative variables, and frequency, percentage, sample size and
#'missing data for qualitative variables.
#'
#'@examples
#'library(magrittr)
#'iris %>% nt_describe(group = Species)
#'
#'@export
nt_describe <- function(data,
                        group = NULL,
                        measures = c("mean.sd", "median.quantiles", "median.range"),
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
               group = group, group.name = group.name,
               format = format, digits = digits)

  out <- Reduce(rbind, temp)

  all.measures <- c("mean.sd", "median.quantiles", "median.range")
  if (!all(measures != measures)){
    if (!any(measures == "mean.sd"))
      out <- out %>% filter(.data$Variable != " Mean \U00b1 SD")
    if (!any(measures == "median.quantiles"))
      out <- out %>% filter(.data$Variable != " Median (Q25% ; Q75%)")
    if (!any(measures == "median.range"))
      out <- out %>% filter(.data$Variable != " Median (Min ; Max)")
  }

  if (save)
    write.csv(out, file = paste0(file, ".csv"))
  return(out)
}

aux_describe <- function(var, var.name, group, group.name, format, digits){

  var.label <- extract_label(var, var.name)
  unit.label <- extract_unit(var)

  if (!is.null(group))
    group.label <- extract_label(group, group.name)

  if (is.numeric(var)){
    out <- describe_quantitative(var = var,
                                 group = group,
                                 format = format,
                                 digits = digits,
                                 var.name = var.name,
                                 var.label = var.label,
                                 unit.label = unit.label,
                                 group.label = group.label)
  } else {
    out <- describe_qualitative(var = var,
                                group = group,
                                format = format,
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
                                  format, digits,
                                  var.name, var.label,
                                  unit.label, group.label){

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
                               unit.label = unit.label)

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
                group.label = group.label)

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

  x <- data_frame(x)

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
                                var.label, unit.label, group.label = NULL){

  var.label <- ifelse(unit.label == "", var.label,
                      paste0(" (", unit.label, ")"))

  aux_variable <- c(var.label,
                    paste(" Mean", "SD", sep = " \U00b1 "),
                    " Median (Q25% ; Q75%)",
                    " Median (Min ; Max)", " Missing")
  aux_measures <- c("", var$mean.sd, var$median.q25.q75,
                    var$median.min.max, unique(var$missing))

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
                                 format = TRUE,
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



