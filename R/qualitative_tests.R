#'@importFrom forcats fct_drop
#'@importFrom dplyr bind_cols
#'@importFrom stats prop.test fisher.test chisq.test mcnemar.test mantelhaen.test
dist_ql_tg <-  function(var, group,
                        alternative,
                        conf.level, paired,
                        digits.p, digits.ci,
                        var.name, var.label, group.label){

  data.test <- data.frame(x = var, g = fct_drop(group))
  lg <- levels(data.test$g)

  if (!paired){
    tab <- table(data.test$g, data.test$x)

    if (!is.finite(max(dim(tab))) | min(dim(tab)) < 2) {
      test <- NA
      p.value <- NA
      hypothesis <- NA
      lower <- NA
      upper <- NA
    } else {
      result <- try(stats::fisher.test(tab), silent = TRUE)

      if (class(result) != "try-error"){
        p.value <- result$p.value

        if (max(dim(tab)) == 2){
          pt <- try(stats::prop.test(tab, conf.level = conf.level), silent = TRUE)

          if (class(pt) != "try-error"){
            lower <- pt$conf.int[[1]]
            upper <- pt$conf.int[[2]]
            test <- "Fisher's Exact"
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
        if (class(result) != "try-error"){
          p.value <- result$p.value
          test <- "Chi-Square"
        } else {
          p.value <- NA
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
    p.value <- result$p.value
    test <- ifelse(max(dim(tab)) == 2, "McNemar test", "McNemar-Bowker test")
    lower <- NA
    upper <- NA
    hypothesis <- "Marginal Homogeneity"
  }

  out <- data.frame(Variable = var.label,
                    Group = group.label,
                    Hypothesis = hypothesis,
                    Lower = round(lower, digits.ci),
                    Upper = round(upper, digits.ci),
                    Test = test,
                    p.value = round(p.value, digits.p))


  return(out)
}


#'@importFrom forcats fct_drop
#'@importFrom dplyr bind_cols
#'@importFrom stats fisher.test chisq.test
dist_ql_mg <-  function(var, group,
                        digits.p, var.name, var.label, group.label){

  data.test <- data.frame(x = var, g = fct_drop(group))
  lg <- levels(data.test$g)

  tab <- table(data.test$g, data.test$x)

  if (!is.finite(max(dim(tab))) | min(dim(tab)) < 2) {
    test <- NA
    p.value <- NA
    hypothesis <- NA
    lower <- NA
    upper <- NA
  } else {
    result <- try(stats::fisher.test(tab), silent = TRUE)

    if (class(result) != "try-error"){
      p.value <- result$p.value
      test <- "Fisher's exact test"
    } else {
      result <- try(stats::chisq.test(tab), silent = TRUE)

      if (class(result) != "try-error"){
        p.value <- result$p.value
        test <- "Chi-square test"
      } else {
        p.value <- NA
        test <- NA
      }
    }
  }

  hypothesis <- "Association"

  out <- data.frame(Variable = var.label[[1]], Group = group.label[[1]],
                    Hypothesis = hypothesis,  Test = test,
                    `p value` =  round(p.value, digits.p))

  return(out)
}

#'@importFrom multcomp glht mcp
#'@importFrom stats confint
dist_ql_mc <-  function(var, group, alternative, contrast,
                        digits.p, digits.ci, var.label, group.label) {

  data.test <- data.frame(x = var, g = group)

  av <- glm(x ~ g, data = data.test, family = "binomial")
  mc <- glht(av, linfct = mcp(g = contrast), alternative = alternative)
  sm <- summary(mc)
  p.value <- round(sm$test$pvalues, digits.p)

  test <- paste(contrast)
  dif <- as.character(rownames(sm[2]$linfct))
  alt <- switch(alternative,
                "two.sided" = " != ",
                "greater" = " < ",
                "less" = " > ")
  hypothesis <- gsub(" - ", alt, dif)

  test <- paste0("Parametric ", contrast)

  lower <- round(exp(confint(mc)$confint[, 2]), digits.ci)
  upper <- round(exp(confint(mc)$confint[, 3]), digits.ci)

  out <- data.frame(Variable = var.label, Group = group.label,
                    Hypothesis = hypothesis, Lower = lower, Upper = upper,
                    Test = test, `p value` = p.value)

  return(out)
}
