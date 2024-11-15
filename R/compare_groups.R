#'Compare two groups
#'
#'@description Performing comparisons between two groups.
#'
#'@importFrom forcats fct_drop
#'@importFrom purrr map2
#'@importFrom dplyr select
#'@importFrom utils write.csv
#'@importFrom rlang := .data quo_is_null enquo
#'@importFrom tibble tibble
#'
#'@param data a data frame with the variables.
#'@param group a data frame with the group variable.
#'@param labels a list of labels with components given by their variable names.
#'@param alternative a character value indicating the alternative hypothesis,
#'must be one of "two.sided", "greater" or "less".
#'@param norm_test a function with a numeric vector as input and a list as output containing an object named \code{p_value} similar to \link[ntimes]{helper_sf_test}.
#'@param var_test a function with a numeric vector, group vector and paired logical variable as input and a list as output containing an object named \code{p_value} similar to \link[ntimes]{helper_levene_test}.
#'@param qt_test a list of functions for four possible cases: (1) normality and homoscedasticity,
#'(2) normality and heteroscedasticity, (3) non-normality and homoscedasticity and (4) normality and heteroscedasticity.
#'@param conf_level a character value specifying the confidence level of the confidence interval for
#'the difference between the two groups.
#'@param paired a logical value indicating whether a paired test should be used.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits_ci the number of digits to present the confidence intervals.
#'@param digits_p the number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@param ... a list with additional arguments to be passed to the helper functions.
#'
#'@examples
#'data(iris)
#'library(dplyr)
#'
#'iris |> filter(Species != "setosa") |>
#'  mutate(Species = droplevels(Species)) |>
#'  nt_compare_tg(group = Species,
#'                labels = list(Sepal.Length = "Sepal Length",
#'                              Sepal.Width = "Sepal Width",
#'                              Petal.Length = "Petal Length",
#'                              Petal.Width = "Petal Width"))
#'@export
nt_compare_tg <- function(data, group, labels = NULL,
                          alternative = "two.sided",
                          norm_test = helper_sf_test,
                          var_test = helper_levene_test,
                          qt_test =
                            list(helper_student_t, helper_welch_t,
                                 helper_mann_whitney, helper_brunner_munzel),
                          paired = FALSE,
                          conf_level = 0.95,
                          format = TRUE,
                          digits_ci = 3,
                          digits_p = 5,
                          save = FALSE,
                          file = "nt_compare_tg",
                          ...){

  group <- enquo(group)

  vars <- select(.data = data, -!!group)
  group <- select(.data = data, !!group)

  vars_name <- names(vars)
  group_name <- names(group)

  vars_name <- names(vars)
  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars_label <- map2(.x = vars, .y = as.list(vars_name),
                       .f = extract_label)
    if (!is.null(group)){
      group <- data_labeller(group, labels)
      group_label <- extract_label(group[[1]], group_name)
    }
  } else {
    vars_label <- map2(.x = vars, .y = as.list(vars_name),
                       .f = extract_label)
    if (!is.null(group))
      group_label <- extract_label(group[[1]], group_name)
  }

  if (!is.factor(group[[1]])){
    group[[1]] <- as.factor(group[[1]])
    warning(paste(group_label, "was transformed into a factor."))
  }

  if (nlevels(group[[1]]) != 2)
    stop("'group' should have only two levels.")
  temp <- pmap(.l = list(vars, vars_name, vars_label),
               .f = aux_compare_tg,
               group = group[[1]],
               group_name = group_name,
               group_label = group_label,
               norm_test = norm_test,
               var_test = var_test,
               qt_test = qt_test,
               paired = paired,
               alternative = alternative,
               conf_level = conf_level,
               format = format,
               digits_p = digits_p,
               digits_ci = digits_ci,
               ...)

  out <- Reduce(rbind, temp)

  if (format){
    out <- out |> mutate('95% CI' =
                            paste0("(", .data$Lower, " ; ",
                                   .data$Upper, ")")) |>
      select(.data$Variable, .data$Group, .data$Hypothesis,
             .data$Test, .data$`95% CI`, `p value` = .data$`p_value`)
  }


  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  attr(out, "ntimes") <- "two_groups"

  return(out)
}

aux_compare_tg <- function(var, var_name, var_label,
                           group, group_name, group_label,
                           norm_test, var_test, qt_test,
                           paired = paired, alternative, conf_level,
                           format, digits_p, digits_ci){

  if (is.numeric(var)){
    out <- dist_qt_tg(var = var,
                      group = group,
                      var_label = var_label,
                      group_label = group_label,
                      norm_test = norm_test,
                      var_test = var_test,
                      qt_test = qt_test,
                      paired = paired,
                      alternative = alternative,
                      conf_level = conf_level,
                      digits_p = digits_p,
                      digits_ci = digits_ci)

  } else {
    if (!is.factor(var)){
      var <- as.factor(var)
      warning(paste(var_label, "was transformed into a factor."))
    }

    out <- dist_ql_tg(var = var,
                      group = group,
                      var_label = var_label,
                      group_label = group_label,
                      paired = paired,
                      alternative = alternative,
                      conf_level = conf_level,
                      digits_p = digits_p,
                      digits_ci = digits_ci)

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
#'@importFrom tibble tibble
#'
#'@param data a data frame with the variables.
#'@param group a data frame with the group variable.
#'@param labels a list of labels with components given by their variable names.
#'@param norm_test a function with a numeric vector as input and a list as output containing an object named \code{p_value} similar to \link[ntimes]{helper_sf_test}.
#'@param var_test a function with a numeric vector, group vector and paired logical variable as input and a list as output containing an object named \code{p_value} similar to \link[ntimes]{helper_levene_test}.
#'@param qt_test a list of functions for three possible cases: (1) normality and homoscedasticity,
#'(2) normality and heteroscedasticity, (3) non-normality and homoscedasticity/heteroscedasticity.
#'@param contrast a matrix of contrasts. See more details in \code{\link[multcomp]{glht}}.
#'@param alternative a character value indicating the alternative hypothesis,
#'must be one of "two.sided", "greater" or "less".
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits_ci the number of digits to present the confidence intervals.
#'@param digits_p the number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@param multiple_comparisons a logical value indicating if pairwise comparisons should be performed.
#'
#'@details If \code{test = "automatic"}, the normality assumption will be verified by
#'\code{norm_test} and homoscedasticity assumption will evaluate the assumption of
#'\code{var_test} at a significance level of 0.05.
#'If the data satisfies both assumptions, then \code{qt_test[[1]]} is chosen;
#'if only normality is satisfied, then \code{qt_test[[2]]}; if only homoscedasticity
#'or neither assumptions, then \code{qt_test[[3]]}.
#'
#'@examples
#'data(iris)
#'
#'iris |> nt_compare_mg(group = Species)
#'
#'@export
nt_compare_mg <- function(data, group, labels = NULL,
                          norm_test = helper_sf_test,
                          var_test = helper_levene_test,
                          qt_test =
                            list(helper_anova, helper_welch_anova,
                                 helper_kruskal_wallis),
                          contrast = "Tukey",
                          alternative = "two.sided",
                          format = TRUE, digits_p = 3, digits_ci = 2,
                          save = FALSE, file = "nt_compare_mg",
                          multiple_comparisons = FALSE){

  group <- enquo(group)

  vars <- select(.data = data, -!!group)
  group <- select(.data = data, !!group)

  vars_name <- names(vars)
  group_name <- names(group)

  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars_label <- map2(.x = vars, .y = as.list(vars_name),
                       .f = extract_label)
    if (!is.null(group)){
      group <- data_labeller(group, labels)
      group_label <- extract_label(group[[1]], group_name)
    }
  } else {
    vars_label <- map2(.x = vars, .y = as.list(vars_name),
                       .f = extract_label)
    if (!is.null(group))
      group_label <- extract_label(group[[1]], group_name)
  }

  if (nlevels(fct_drop(group[[1]])) == 2)
    stop("'group' should have more than two levels.")

  temp <- pmap(.l = list(vars, vars_name, vars_label),
               .f = aux_compare_mg,
               group = group[[1]],
               group_name = group_name,
               group_label = group_label,
               norm_test = norm_test, var_test = var_test,
               qt_test = qt_test,
               digits_p = digits_p,
               mc = multiple_comparisons)

  omnibus_test <- Reduce(rbind, temp)

  if (multiple_comparisons){
    aux <- omnibus_test |> filter(.data$p_value < 0.05)
    vars_name <- unlist(aux$Variable)
    vars <- data |> select(all_of(vars_name))
    vars_label.mc <- vars_label[vars_name]
    group <- data |> select(all_of(group_name))
    test <- omnibus_test |> pull(.data$Test)

    temp <- pmap(list(vars, vars_name, vars_label.mc, test),
                 .f = aux_compare_mc,
                 group = group,
                 group_name = group_name, group_label = group_label,
                 alternative = alternative, contrast = contrast,
                 digits_p = digits_p, digits_ci = digits_ci)
    mc.test <- Reduce(rbind, temp)

    if (format){
      mc.test <- mc.test |>
        mutate('95% CI' =
                 paste0("(", .data$Lower, " ; ",
                        .data$Upper, ")"),) |>
        select(.data$Variable, .data$Group, .data$Hypothesis,
               .data$Test, .data$`95% CI`, `p value` = .data$`p_value`)

      if (save)
        write.csv(mc.test, file = paste0(file, "_mc_test.csv"))
    }

    if (!is.null(labels))
      omnibus_test <- omnibus_test |>
        mutate(Variable = str_replace_all(.data$Variable, unlist(vars_label)))
  }

  if (format){
    omnibus_test <- omnibus_test |>
      select(.data$Variable, .data$Group, .data$Test, .data$Hypothesis,
             `p value` = .data$`p_value`)
  }

  if (save)
    write.csv(omnibus_test, file = paste0(file, "_omnibus_test.csv"))

  if (!multiple_comparisons){
    out <- list(omnibus_test = omnibus_test)
    attr(out, "ntimes") <- "multiple_groups"
  } else {
    out <- list(omnibus_test = omnibus_test, mc.test = mc.test)
    attr(out, "ntimes") <- "multiple_comparisons"
  }



  return(out)
}


aux_compare_mg <- function(var, var_name, var_label,
                           group, group_name, group_label,
                           norm_test, var_test,
                           qt_test,
                           format, digits_p, mc){
  if (mc){
    var_label <- var_name
    group_label <- group_name
  }

  if (is.numeric(var)){
    out <- dist_qt_mg(var = var,
                      group = group,
                      qt_test = qt_test,
                      norm_test = norm_test,
                      digits_p = digits_p,
                      var_label = var_label,
                      group_label = group_label)

  } else {
    out <- dist_ql_mg(var = var,
                      group = group,
                      digits_p = digits_p,
                      var_label = var_label,
                      group_label = group_label)

  }

  return(out)
}

aux_compare_mc <- function(var, var_name, var_label,
                           omnibus_test,
                           group, group_name, group_label,
                           alternative, contrast, digits_p, digits_ci){


  if (is.numeric(var)){
    out <- dist_qt_mc(var = var,
                         omnibus_test = omnibus_test,
                         group = group[[1]],
                         alternative = alternative,
                         contrast = contrast,
                         digits_p = digits_p,
                         digits_ci = digits_ci,
                         var_label = var_label,
                         group_label = group_label)

  } else {
    out <- dist_ql_mc(var = var,
                         group = group[[1]],
                         alternative = alternative,
                         contrast = contrast,
                         digits_p = digits_p,
                         digits_ci = digits_ci,
                         var_label = var_label,
                         group_label = group_label)
  }

  return(out)
}
