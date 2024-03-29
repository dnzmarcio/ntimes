#'Mean and standard deviation
#'
#'@description It calculates mean and standard deviation, concatenating them to
#'present on a table.
#'
#'@param var a numeric vector.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with the first element \code{name} as the measure name and the
#'second element as the \code{value} for a given variable.
#'
#'@importFrom stats sd
#'@export
helper_mean_sd <- function(var, digits, ...){

  ldots <- list(...)

  mean <- format(round(mean(var, na.rm = TRUE), digits), nsmall = digits)
  sd <- format(round(sd(var, na.rm = TRUE), digits), nsmall = digits)
  name <- paste0("\t Mean", " \U00b1 ", "SD")
  value <- paste0(mean, " \U00b1 ", sd)
  out <- list(name = name, value = value)
  return(out)
}

#'Median with first and third quantiles
#'
#'@description It calculates median with quantiles 25% and 75%, concatenating them
#'to present on a table.
#'
#'@param var a numeric vector.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with the first element \code{name} as the measure name and the
#'second element as the \code{value} for a given variable.
#'
#'@importFrom stats median quantile
#'@export
helper_median_iqr <- function(var, digits, ...){

  ldots <- list(...)

  median <- format(round(median(var, na.rm = TRUE), digits),
                   nsmall = digits)
  q25 <- format(round(quantile(var, probs = 0.25, na.rm = TRUE), digits),
                nsmall = digits)
  q75 <- format(round(quantile(var, probs = 0.75, na.rm = TRUE), digits),
                nsmall = digits)
  name <- "\t Median (Q25% ; Q75%)"
  value <- paste0(median," (", q25, " ; ", q75, ")")
  out <- list(name = name, value = value)
  return(out)
}


#'Median with minimum and maximum
#'
#'@description It calculates median with minimum and maximum, concatenating them
#'to present on a table.
#'
#'@param var a numeric vector.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with the first element \code{name} as the measure name and the
#'second element as the \code{value} for a given variable.
#'
#'@importFrom stats median
#'@export
helper_median_range <- function(var, digits, ...){

  ldots <- list(...)

  median <- format(round(median(var, na.rm = TRUE), digits),
                   nsmall = digits)
  min <- format(round(min(var, na.rm = TRUE), digits),
                nsmall = digits)
  max <- format(round(max(var, na.rm = TRUE), digits),
                nsmall = digits)
  name <- "\t Median (Min ; Max)"
  value <- paste0(median," (", min, " ; ", max, ")")
  out <- list(name = name, value = value)
  return(out)
}

#'Number of missing observations
#'
#'@description It calculates the number of missing observations.
#'
#'@param var a numeric vector.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with the first element \code{name} as the measure name and the
#'second element as the \code{value} for a given variable.
#'
#'@export
helper_missing <- function(var, ...){

  ldots <- list(...)

  out <- list(name = "\t Missing",
              value = sum(is.na(var)))

  return(out)
}

#'Percentages and frequencies
#'
#'@description It calculates percentages and frequencies, concatenating them to
#'present on a table.
#'
#'@param var a numeric vector.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with the first element \code{name} as the measure name and the
#'second element as the \code{value} for a given variable.
#'
#'@importFrom forcats fct_na_value_to_level
#'@export
helper_perc_count <- function(var, digits, ...){

  ldots <- list(...)

  h <- fct_na_value_to_level(var, level = "Missing")
  lh <- levels(h)

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

  lh <- paste0("\t ", lh)

  out <- list(name = lh, value = perc_count)
}

#'Frequencies and percentages
#'
#'@description It calculates percentages and frequencies, concatenating them to
#'present on a table.
#'
#'@param var a numeric vector.
#'@param digits a numeric value specifying the number of digits to present the results.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with the first element \code{name} as the measure name and the
#'second element as the \code{value} for a given variable.
#'
#'@importFrom forcats fct_na_value_to_level
#'@export
helper_count_perc <- function(var, digits, ...){

  ldots <- list(...)

  h <- fct_na_value_to_level(var, level = "Missing")
  lh <- levels(h)

  count <- tapply(h, h, length)
  count <- ifelse(is.na(count), 0, count)
  n <- length(h)
  perc <- 100*prop.table(count)
  perc <- ifelse(!is.finite(perc), NA, format(round(perc, digits), nsmall = digits))

  perc_count <- paste0(count, " (", perc, ")")

  if (!("Missing" %in% lh)){
    lh <- c(lh, "Missing")
    perc_count <- c(perc_count, "0 (0)")
  }

  lh <- paste0("\t ", lh)

  out <- list(name = lh, value = perc_count)
}

#'Helper for Shapiro Francia test for Normality
#'
#'@description It performs the Shapiro-Francia to test the assumption of normality.
#'
#'@param x a numeric vector.
#'
#'@return a list indicating the test and its p-value.
#'
#'@importFrom nortest sf.test
#'@export
helper_sf_test <-  function(x){

  p.value <- nortest::sf.test(x)$p.value
  test <- "Shapiro-Francia"

  out <- list(test = test, p.value = p.value)
}

#'Helper for Levene test
#'
#'@description It performs the Levene test to test homoscedasticity.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'@param paired a logical value.
#'
#'@return a list indicating the test and its p-value.
#'
#'@importFrom tidyr drop_na
#'@importFrom lawstat levene.test
#'@importFrom PairedData levene.Var.test.paired
#'@export
helper_levene_test <- function(x, g, paired){

  temp <- na.exclude(data.frame(x, g))
  if (!paired){
    p.value <- lawstat::levene.test(temp$x, temp$g)$p.value
    test <- "Levene"
  } else {
    lg <- levels(g)
    p.value <- PairedData::levene.Var.test.paired(temp$x[temp$g == lg[1]],
                                                  temp$x[temp$g == lg[2]])$p.value
  }

  out <- list(test = test, p.value = p.value)
}

#'Helper for Analysis of Variance
#'
#'@description It performs the Analysis of Variance.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed and its respective \code{p.value}.
#'
#'@export
helper_anova <- function(x, g){

  data.test <- data.frame(x, g)
  result <- stats::oneway.test(x ~ g, data = data.test, var.equal = TRUE)
  test <- "ANOVA"
  p.value <- result$p.value

  out <- list(test = test, p.value = p.value)
}

#'Helper for Welch's Analysis of Variance
#'
#'@description It performs the Welch's Analysis of Variance.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed and its respective \code{p.value}.
#'
#'@export
helper_welch_anova <- function(x, g){

  data.test <- data.frame(x, g)
  result <- stats::oneway.test(x ~ g, data = data.test, var.equal = FALSE)
  test <- "Welch's ANOVA"
  p.value <- result$p.value

  out <- list(test = test, p.value = p.value)
}

#'Helper for Kruskal-Wallis test
#'
#'@description It performs the Kruskal-Wallis test.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed and its respective \code{p.value}.
#'
#'@export
helper_kruskal_wallis <- function(x, g){

  data.test <- data.frame(x, g)
  result <- stats::kruskal.test(x ~ g, data = data.test)
  test <- "Kruskal-Wallis"
  p.value <- result$p.value

  out <- list(test = test, p.value = p.value)
}

#'Helper for Student's t-test
#'
#'@description It performs Student's t-test.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'@param paired a logical value.
#'@param alternative a character value.
#'@param conf.level a numeric value.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed, its respective \code{p.value} and \code{lower} and \code{upper}
#'limits of the confidence interval for the effect size considered in the test.
#'
#'@export
helper_student_t <- function(x, g, paired, alternative, conf.level){

  data.test <- data.frame(x, g)
  result <- stats::t.test(data.test$x ~ data.test$g,
                          var.equal = TRUE,
                          alternative = alternative,
                          paired = paired,
                          conf.level = conf.level)

  p.value <- result$p.value
  lower <- result$conf.int[1]
  upper <- result$conf.int[2]
  test <- ifelse(paired, "Paired Student's t-test", "Student's t-test")

  out <- list(test = test, p.value = p.value,
              lower = lower, upper = upper)
}

#'Helper for Welch's t-test
#'
#'@description It performs the Welch's t-test.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'@param paired a logical value.
#'@param alternative a character value.
#'@param conf.level a numeric value.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed, its respective \code{p.value} and \code{lower} and \code{upper}
#'limits of the confidence interval for the effect size considered in the test.
#'
#'@export
helper_welch_t <- function(x, g, paired, alternative, conf.level){

  data.test <- data.frame(x, g)
  result <- stats::t.test(data.test$x ~ data.test$g,
                          var.equal = FALSE,
                          alternative = alternative,
                          paired = paired,
                          conf.level = conf.level)

  p.value <- result$p.value
  lower <- result$conf.int[1]
  upper <- result$conf.int[2]
  test <- ifelse(paired, "Paired Welch's t-test", "Welch's t-test")

  out <- list(test = test, p.value = p.value,
              lower = lower, upper = upper)
}

#'Helper for Mann-Whitney test
#'
#'@description It performs the Mann-Whitney test.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'@param paired a logical value.
#'@param alternative a character value.
#'@param conf.level a numeric value.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed, its respective \code{p.value} and \code{lower} and \code{upper}
#'limits of the confidence interval for the effect size considered in the test.
#'
#'@export
helper_mann_whitney <- function(x, g, paired, alternative, conf.level){

  data.test <- data.frame(x, g)
  result <- stats::wilcox.test(data.test$x ~ data.test$g,
                               alternative = alternative,
                               paired = paired, conf.int = TRUE,
                               conf.level = conf.level)

  p.value <- result$p.value
  lower <- result$conf.int[1]
  upper <- result$conf.int[2]
  test <- ifelse(paired, "Paired Mann-Whitney", "Mann-Whitney")

  out <- list(test = test, p.value = p.value,
              lower = lower, upper = upper)

}

#'Helper for Brunner-Munzel t-test
#'
#'@description It performs the Brunner-Munzel t-test.
#'
#'@param x a numeric vector.
#'@param g a factor.
#'@param paired a logical value.
#'@param alternative a character value.
#'@param conf.level a numeric value.
#'
#'@details This function can be modified by the user,
#'but input and output should be kept the same.
#'
#'@return a list with \code{test} indicating the test that was
#'performed, its respective \code{p.value} and \code{lower} and \code{upper}
#'limits of the confidence interval for the effect size considered in the test.
#'
#'@export
helper_brunner_munzel <- function(x, g, paired, alternative, conf.level){

  data.test <- data.frame(x, g)

  if (!paired){
    result <- nparcomp::npar.t.test(x ~ g,
                                    data = data.test,
                                    method = "t.app",
                                    conf.level = conf.level,
                                    alternative = alternative,
                                    info = FALSE)

    p.value <- result$Analysis[, 6]
    lower <- ifelse(result$Analysis[, 3] < 0, 0, result$Analysis[, 3])
    upper <- ifelse(result$Analysis[, 4] > 1, 1, result$Analysis[, 4])
    test <- "Brunner-Munzel t-test"
  } else {
    result <- nparcomp::npar.t.test.paired(x ~ g,
                                           data = data.test,
                                           alternative = alternative,
                                           conf.level = conf.level,
                                           info = FALSE)

    p.value <- result$Analysis[1, 5]
    lower <- ifelse(result$Analysis[1, 1] < 0, 0, result$Analysis[1, 1])
    upper <- ifelse(result$Analysis[1, 3] > 1, 1, result$Analysis[1, 3])
    test <- "Paired Brunner-Munzel t-test"
  }

  out <- list(test = test, p.value = p.value,
              lower = lower, upper = upper)

}
