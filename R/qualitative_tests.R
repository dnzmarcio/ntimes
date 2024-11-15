#'@importFrom forcats fct_drop
#'@importFrom dplyr bind_cols
#'@importFrom stats prop.test fisher.test chisq.test mcnemar.test mantelhaen.test
#'@importFrom methods is
dist_ql_tg <-  function(var, group,
                        alternative,
                        conf_level, paired,
                        digits_p, digits_ci,
                        var_name, var_label, group_label){

  data.test <- data.frame(x = var, g = fct_drop(group))
  lg <- levels(data.test$g)

  if (!paired){
    tab <- table(data.test$g, data.test$x)

    if (!is.finite(max(dim(tab))) | min(dim(tab)) < 2) {
      test <- NA
      p_value <- NA
      hypothesis <- NA
      lower <- NA
      upper <- NA
    } else {
      result <- try(stats::fisher.test(tab), silent = TRUE)

      if (any(is(result) != "try-error")){
        p_value <- result$p.value

        if (max(dim(tab)) == 2){
          pt <- try(stats::prop.test(tab, conf_level = conf_level), silent = TRUE)

          if (any(is(pt) != "try-error")){
            lower <- pt$conf.int[[1]]
            upper <- pt$conf.int[[2]]
            test <- "Fisher's exact test"
          } else {
            lower <- NA
            upper <- NA
            test <- "Fisher's exact test"
          }

        } else {
          lower <- NA
          upper <- NA
          test <- "Fisher-Freeman-Halton's exact test"
        }
      } else {
        result <- try(stats::chisq.test(tab), silent = TRUE)
        lower <- NA
        upper <- NA
        if (any(is(result) != "try-error")){
          p_value <- result$p.value
          test <- "Chi-Square test"
        } else {
          p_value <- NA
          test <- NA
        }
      }

      alt <- switch(alternative,
                    "greater" = ">",
                    "less" = "<",
                    "two.sided" = "!=")
      hypothesis <- paste(lg[2], "=", lg[1])
    }

  } else {
    tab <- table(data.test$x, data.test$g)
    result <- mcnemar.test(tab)
    p_value <- result$p.value
    test <- ifelse(max(dim(tab)) == 2, "McNemar test",
                   "McNemar-Bowker test")
    lower <- NA
    upper <- NA
    hypothesis <- "Marginal Homogeneity"
  }

  out <- data.frame(Variable = var_label,
                    Group = group_label,
                    Hypothesis = hypothesis,
                    Lower = round(lower, digits_ci),
                    Upper = round(upper, digits_ci),
                    Test = test,
                    p_value = round(p_value, digits_p))


  return(out)
}


#'@importFrom forcats fct_drop
#'@importFrom dplyr bind_cols
#'@importFrom stats fisher.test chisq.test
dist_ql_mg <-  function(var, group,
                        digits_p, var_name, var_label, group_label){

  data.test <- data.frame(x = var, g = fct_drop(group))
  lg <- levels(data.test$g)

  tab <- table(data.test$g, data.test$x)

  if (!is.finite(max(dim(tab))) | min(dim(tab)) < 2) {
    test <- NA
    p_value <- NA
    hypothesis <- NA
    lower <- NA
    upper <- NA
  } else {
    result <- try(stats::fisher.test(tab), silent = TRUE)

    if (any(is(result) != "try-error")){
      p_value <- result$p.value
      test <- "Fisher's exact test"
    } else {
      result <- try(stats::chisq.test(tab), silent = TRUE)

      if (any(is(result) != "try-error")){
        p_value <- result$p.value
        test <- "Chi-square test"
      } else {
        p_value <- NA
        test <- NA
      }
    }
  }

  hypothesis <- "Association"

  out <- data.frame(Variable = var_label[[1]], Group = group_label[[1]],
                    Hypothesis = hypothesis,  Test = test,
                    `p value` =  round(p_value, digits_p))

  return(out)
}

#'@importFrom multcomp glht mcp
#'@importFrom stats confint
dist_ql_mc <-  function(var, group, alternative, contrast,
                        digits_p, digits_ci, var_label, group_label) {

  data.test <- data.frame(x = var, g = group)

  av <- glm(x ~ g, data = data.test, family = "binomial")
  mc <- glht(av, linfct = mcp(g = contrast), alternative = alternative)
  sm <- summary(mc)
  p_value <- round(sm$test$pvalues, digits_p)

  test <- paste(contrast)
  dif <- as.character(rownames(sm[2]$linfct))
  alt <- switch(alternative,
                "two.sided" = " != ",
                "greater" = " < ",
                "less" = " > ")
  hypothesis <- gsub(" - ", alt, dif)

  test <- paste0("Parametric ", contrast)

  lower <- round(exp(confint(mc)$confint[, 2]), digits_ci)
  upper <- round(exp(confint(mc)$confint[, 3]), digits_ci)

  out <- data.frame(Variable = var_label, Group = group_label,
                    Hypothesis = hypothesis, Lower = lower, Upper = upper,
                    Test = test, `p value` = p_value)

  return(out)
}
