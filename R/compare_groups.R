#'Compare two groups
#'
#'@description Performing comparisons between two groups.
#'
#'@importFrom forcats fct_drop
#'@importFrom purrr map2
#'@importFrom dplyr select
#'@importFrom utils write.csv
#'@importFrom rlang := .data quo_is_null enquo
#'@importFrom tibble data_frame as_data_frame
#'
#'@param data a data frame with the variables.
#'@param group a data frame with the group variable.
#'@param alternative a character value indicating the alternative hypothesis,
#'must be one of "two.sided", "greater" or "less".
#'@param test a character value indicating the tests to be performed.
#'The options are "auto", "par" and "npar". See more in details.
#'@param conf.level a character value specifying the confidence level of the confidence interval for
#'the difference between the two groups.
#'@param paired a logical value indicating whether a paired test should be used.
#'@param norm.test a character value specifying the normality test to be performed.
#'The options are Anderson-Darling (\code{ad}), Shapiro-Francia (\code{"sf"}),
#'Kolmogorov-Smirnov (\code{ks}),  Cramer-vonMises (\code{cvm}) and
#'Pearson (\code{p}). The default is Shapiro-Francia (\code{"sf"}). It is only used if
#'\code{test = "automatic"}.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits.ci the number of digits to present the confidence intervals.
#'@param digits.p the number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'
#'@details If \code{test = "automatic"}, the normality assumption will be verified by
#'a normality test (Anderson-Daling (\link[nortest]{ad.test}),
#'Shapiro-Francia (\link[nortest]{sf.test}),
#''Kolmogorov-Smirnov (\link[nortest]{lillie.test}),
#'Cramer-vonMises (\link[nortest]{cvm.test}),
#'and Pearson (\link[nortest]{pearson.test})) and
#'Levene test (\link[car]{leveneTest}) will evaluate the assumption of
#'homocedasticity at a significance level of 0.05.
#'If the data satisfies both assumptions, then t-test is chosen;
#'if only normality is satisfied, then Welch t-test; if only homoscedasticity, then
#'Mann-Whitney; if neither assumptions, then Brunner-Munzel t test.
#'
#'@examples
#'library(dplyr)
#'library(forcats)
#'data(iris)
#'
#'iris_nt <- iris %>% filter(Species != "setosa") %>% mutate(Species = fct_drop(Species))
#'iris_nt %>% nt_compare_tg(group = Species)
#'
#'@export
nt_compare_tg <- function(data, group,
                          alternative = "two.sided",
                          test = "auto",
                          conf.level = 0.95,
                          paired = FALSE,
                          norm.test = "sf",
                          format = TRUE,
                          digits.ci = 3,
                          digits.p = 5,
                          save = FALSE,
                          file = "nt_compare_tg"){

  group <- enquo(group)

  vars <- select(.data = data, -!!group)
  group <- select(.data = data, !!group)

  vars.name <- names(vars)
  group.name <- names(group)

  if (nlevels(fct_drop(group[[1]])) != 2)
    stop("'group' should have only two levels.")

  temp <- map2(.x = vars, .y = vars.name, .f = aux_compare_tg,
               group = group, group.name = group.name,
               test = test,
               alternative = alternative,
               conf.level = conf.level,
               paired = paired,
               norm.test = norm.test,
               format = format,
               digits.p = digits.p,
               digits.ci = digits.ci)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))
  return(out)
}

aux_compare_tg <- function(var, var.name, group, group.name = group.name,
                           test, alternative, conf.level,
                           paired = paired, norm.test,
                           format, digits.p, digits.ci){

  var.label <- extract_label(var, var.name)
  group.label <- extract_label(group, group.name)
  if (is.numeric(var)){
      out <- nt_dist_qt_tg(var = var,
                           group = group,
                           test = test,
                           alternative = alternative,
                           conf.level = conf.level,
                           paired = paired,
                           norm.test = norm.test,
                           format = format,
                           digits.p = digits.p,
                           digits.ci = digits.ci,
                           var.name = var.name,
                           var.label = var.label,
                           group.label = group.label)

  } else {
    out <- nt_dist_ql_tg(var = var,
                         group = group,
                         alternative = alternative,
                         conf.level = conf.level,
                         paired = paired,
                         format = format,
                         digits.p = digits.p,
                         digits.ci = digits.ci,
                         var.name = var.name,
                         var.label = var.label,
                         group.label = group.label)

  }

  return(out)
}




#'Compare more than two groups
#'
#'@description Performing comparisons among three or more groups.
#'
#'@importFrom forcats fct_drop
#'@importFrom purrr map2
#'@importFrom dplyr select
#'@importFrom utils write.csv
#'@importFrom rlang := .data quo_is_null enquo
#'@importFrom tibble data_frame as_data_frame
#'
#'@param data a data frame with the variables.
#'@param group a data frame with the group variable.
#'@param test a character value indicating the tests to be performed.
#'The options are "auto", "par" and "npar". See more in details.
#'@param norm.test a character value specifying the normality test to be performed.
#'The options are Anderson-Darling (\code{ad}), Shapiro-Francia (\code{"sf"}),
#'Kolmogorov-Smirnov (\code{ks}),  Cramer-vonMises (\code{cvm}) and
#'Pearson (\code{p}). The default is Shapiro-Francia (\code{"sf"}). It is only used if
#'\code{test = "automatic"}.
#'@param digits.p the number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@param keep.data a logical value indicating if
#'
#'@details If \code{test = "automatic"}, the normality assumption will be verified by
#'a normality test (Anderson-Daling (\link[nortest]{ad.test}),
#'Shapiro-Francia (\link[nortest]{sf.test}),
#''Kolmogorov-Smirnov (\link[nortest]{lillie.test}),
#'Cramer-vonMises (\link[nortest]{cvm.test}),
#'and Pearson (\link[nortest]{pearson.test})) and
#'Levene test (\link[car]{leveneTest}) will evaluate the assumption of
#'homocedasticity at a significance level of 0.05.
#'If the data satisfies both assumptions, then ANOVA is chosen;
#'if only normality is satisfied, then Welch ANOVA; if only homoscedasticity
#'or neither assumptions, then Kruskal-Wallis.
#'
#'@examples
#'library(magrittr)
#'data(iris)
#'
#'iris %>% nt_compare_mg(group = Species)
#'
#'@export
nt_compare_mg <- function(data, group,
                          test = "auto",
                          norm.test = "sf",
                          digits.p = 5,
                          save = FALSE,
                          file = "nt_compare_mg",
                          keep.data = FALSE){

  group <- enquo(group)

  vars <- select(.data = data, -!!group)
  group <- select(.data = data, !!group)

  vars.name <- names(vars)
  group.name <- names(group)

  if (nlevels(fct_drop(group[[1]])) == 2)
    stop("'group' should have more than two levels.")

  temp <- map2(.x = vars, .y = vars.name, .f = aux_compare_mg,
               group = group, group.name = group.name,
               test = test,
               norm.test = norm.test, digits.p = digits.p)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  if (keep.data){
    result <- out
    out <- list()
    out$result <- result
    out$data <- data
  }

  return(out)
}


aux_compare_mg <- function(var, var.name, group, group.name = group.name,
                           test, norm.test,
                           format, digits.p){

  var.label <- extract_label(var, var.name)
  group.label <- extract_label(group, group.name)

  if (is.numeric(var)){
    out <- nt_dist_qt_mg(var = var,
                         group = group,
                         test = test,
                         norm.test = norm.test,
                         digits.p = digits.p,
                         var.name = var.name,
                         var.label = var.label,
                         group.label = group.label)

  } else {
    out <- nt_dist_ql_mg(var = var,
                         group = group,
                         digits.p = digits.p,
                         var.name = var.name,
                         var.label = var.label,
                         group.label = group.label)

  }

  return(out)
}
