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
#'@param norm.test a function with a numeric vector as input and a list as output containing an object named \code{p.value} similar to \link[ntimes]{helper_sf_test}.
#'@param var.test a function with a numeric vector, group vector and paired logical variable as input and a list as output containing an object named \code{p.value} similar to \link[ntimes]{helper_levene_test}.
#'@param qt.test a list of functions for four possible cases: (1) normality and homoscedasticity,
#'(2) normality and heteroscedasticity, (3) non-normality and homoscedasticity and (4) normality and heteroscedasticity.
#'@param conf.level a character value specifying the confidence level of the confidence interval for
#'the difference between the two groups.
#'@param paired a logical value indicating whether a paired test should be used.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits.ci the number of digits to present the confidence intervals.
#'@param digits.p the number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'
#'@examples
#'library(maggritr)
#'data(iris)
#'
#'iris %>% filter(Species != "setosa") %>%
#'  mutate(Species = as.factor(Species)) %>%
#'  nt_compare_tg(group = Species,
#'                labels = list(Sepal.Length = "Sepal Length",
#'                              Sepal.Width = "Sepal Width",
#'                              Petal.Length = "Petal Length",
#'                              Petal.Width = "Petal Width"))
#'@export
nt_compare_tg <- function(data, group, labels = NULL,
                          alternative = "two.sided",
                          norm.test = helper_sf_test,
                          var.test = helper_levene_test,
                          qt.test =
                            list(helper_student_t, helper_welch_t,
                                 helper_mann_whitney, helper_brunner_munzel),
                          paired = FALSE,
                          conf.level = 0.95,
                          format = TRUE,
                          digits.ci = 3,
                          digits.p = 5,
                          save = FALSE,
                          file = "nt_compare_tg",
                          ...){

  group <- enquo(group)

  vars <- select(.data = data, -!!group)
  group <- select(.data = data, !!group)

  vars.name <- names(vars)
  group.name <- names(group)

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

  if (nlevels(fct_drop(group[[1]])) != 2)
    stop("'group' should have only two levels.")
  temp <- pmap(.l = list(vars, vars.name, vars.label),
               .f = aux_compare_tg,
               group = group[[1]],
               group.name = group.name,
               group.label = group.label,
               norm.test = norm.test,
               var.test = var.test,
               qt.test = qt.test,
               paired = paired,
               alternative = alternative,
               conf.level = conf.level,
               format = format,
               digits.p = digits.p,
               digits.ci = digits.ci,
               ...)

  out <- Reduce(rbind, temp)

  if (format){
    out <- out %>% mutate('95% CI' =
                            paste0("(", .data$Lower, " ; ",
                                   .data$Upper, ")")) %>%
      select(.data$Variable, .data$Group, .data$Hypothesis,
             .data$Test, .data$`95% CI`, `p value` = .data$`p.value`)
  }


  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  attr(out, "ntimes") <- "two_groups"

  return(out)
}

aux_compare_tg <- function(var, var.name, var.label,
                           group, group.name, group.label,
                           norm.test, var.test, qt.test,
                           paired = paired, alternative, conf.level,
                           format, digits.p, digits.ci){

  if (is.numeric(var)){
    out <- dist_qt_tg(var = var,
                      group = group,
                      var.label = var.label,
                      group.label = group.label,
                      norm.test = norm.test,
                      var.test = var.test,
                      qt.test = qt.test,
                      paired = paired,
                      alternative = alternative,
                      conf.level = conf.level,
                      digits.p = digits.p,
                      digits.ci = digits.ci)

  } else {
    if (!is.factor(var)){
      var <- as.factor(var)
      warning(paste(var.label, "was transformed into a factor."))
    }

    out <- dist_ql_tg(var = var,
                      group = group,
                      var.label = var.label,
                      group.label = group.label,
                      paired = paired,
                      alternative = alternative,
                      conf.level = conf.level,
                      digits.p = digits.p,
                      digits.ci = digits.ci)

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
#'@param norm.test a function with a numeric vector as input and a list as output containing an object named \code{p.value} similar to \link[ntimes]{helper_sf_test}.
#'@param var.test a function with a numeric vector, group vector and paired logical variable as input and a list as output containing an object named \code{p.value} similar to \link[ntimes]{helper_levene_test}.
#'@param qt.test a list of functions for three possible cases: (1) normality and homoscedasticity,
#'(2) normality and heteroscedasticity, (3) non-normality and homoscedasticity/heteroscedasticity.
#'@param contrast a matrix of contrasts. See more details in \code{\link[multcomp]{glht}}.
#'@param alternative a character value indicating the alternative hypothesis,
#'must be one of "two.sided", "greater" or "less".
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits.ci the number of digits to present the confidence intervals.
#'@param digits.p the number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@param mc a logical value indicating if pairwise comparisons should be performed.
#'
#'@details If \code{test = "automatic"}, the normality assumption will be verified by
#'\code{norm.test} and homoscedasticity assumption will evaluate the assumption of
#'\code{var.test} at a significance level of 0.05.
#'If the data satisfies both assumptions, then \code{qt.test[[1]]} is chosen;
#'if only normality is satisfied, then \code{qt.test[[2]]}; if only homoscedasticity
#'or neither assumptions, then \code{qt.test[[3]]}.
#'
#'@examples
#'library(magrittr)
#'data(iris)
#'
#'iris %>% nt_compare_mg(group = Species)
#'
#'@export
nt_compare_mg <- function(data, group,
                          norm.test = helper_sf_test,
                          var.test = helper_levene_test,
                          qt.test =
                            list(helper_anova, helper_welch_anova,
                                 helper_kruskal_wallis),
                          contrast = "Tukey",
                          alternative = "two.sided",
                          format = TRUE, digits.p = 3, digits.ci = 2,
                          save = FALSE, file = "nt_compare_mg", mc = FALSE){

  group <- enquo(group)

  vars <- select(.data = data, -!!group)
  group <- select(.data = data, !!group)

  vars.name <- names(vars)
  group.name <- names(group)

  if (nlevels(fct_drop(group[[1]])) == 2)
    stop("'group' should have more than two levels.")

  temp <- map2(.x = vars, .y = vars.name, .f = aux_compare_mg,
               group = group, group.name = group.name,
               norm.test = norm.test, var.test = var.test,
               qt.test = qt.test,
               digits.p = digits.p, mc = mc)

  omnibus.test <- Reduce(rbind, temp)

  if (mc){
    aux <- omnibus.test[which(omnibus.test$p.value < 0.05), ]
    vars.name <- aux$Variable
    vars <- data[vars.name]
    group.name <- unique(aux$Group)
    group <- data[group.name]
    test <- omnibus.test$Test

    temp <- pmap(list(vars, vars.name, test), .f = aux_compare_mc,
                 group = group, group.name = group.name,
                 alternative = alternative, contrast = contrast,
                 digits.p = digits.p, digits.ci = digits.ci)
    mc.test <- Reduce(rbind, temp)

    if (format){
      mc.test <- mc.test %>%
        mutate('95% CI' =
                 paste0("(", .data$Lower, " ; ",
                        .data$Upper, ")"),) %>%
        select(.data$Variable, .data$Group, .data$Hypothesis,
               .data$Test, .data$`95% CI`, `p value` = .data$`p.value`)

      if (save)
        write.csv(mc.test, file = paste0(file, "_mc_test.csv"))
    }
  }

  if (format){
    omnibus.test <- omnibus.test %>%
      select(.data$Variable, .data$Group, .data$Test, .data$Hypothesis,
             `p value` = .data$`p.value`)
  }

  if (save)
    write.csv(omnibus.test, file = paste0(file, "_omnibus_test.csv"))

  if (!mc){
    out <- list(omnibus.test = omnibus.test)
    attr(out, "ntimes") <- "multiple_groups"
  } else {
    out <- list(omnibus.test = omnibus.test, mc.test = mc.test)
    attr(out, "ntimes") <- "multiple_comparisons"
  }



  return(out)
}


aux_compare_mg <- function(var, var.name, group, group.name,
                           norm.test, var.test,
                           qt.test,
                           format, digits.p, mc){
  if (mc){
    var.label <- var.name
    group.label <- group.name
  } else {
    unit.label <- extract_unit(var)
    var.label <- extract_label(var, var.name)
    group.label <- extract_label(group, group.name)
  }

  if (is.numeric(var)){
    out <- dist_qt_mg(var = var,
                         group = group[[1]],
                         qt.test = qt.test,
                         norm.test = norm.test,
                         digits.p = digits.p,
                         var.label = var.label,
                         group.label = group.label)

  } else {
    out <- dist_ql_mg(var = var,
                         group = group[[1]],
                         digits.p = digits.p,
                         var.label = var.label,
                         group.label = group.label)

  }

  return(out)
}

aux_compare_mc <- function(var, var.name, omnibus.test, group, group.name,
                           alternative, contrast, digits.p, digits.ci){

  var.label <- extract_label(var, var.name)
  group.label <- extract_label(group, group.name)

  if (is.numeric(var)){
    out <- dist_qt_mc(var = var,
                         omnibus.test = omnibus.test,
                         group = group[[1]],
                         alternative = alternative,
                         contrast = contrast,
                         digits.p = digits.p,
                         digits.ci = digits.ci,
                         var.label = var.label,
                         group.label = group.label)

  } else {
    out <- dist_ql_mc(var = var,
                         group = group[[1]],
                         alternative = alternative,
                         contrast = contrast,
                         digits.p = digits.p,
                         digits.ci = digits.ci,
                         var.label = var.label,
                         group.label = group.label)
  }

  return(out)
}
