#'@importFrom forcats fct_drop
#'@importFrom dplyr bind_cols
#'@importFrom stats prop.test fisher.test chisq.test mcnemar.test mantelhaen.test
nt_dist_ql_tg <-  function(var, group,
                           alternative,
                           conf.level, paired,
                           format, digits.p, digits.ci,
                           var.name, var.label, group.label){

  data.test <- data_frame(x = var, g = fct_drop(group[[1]]))
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
      result <- try(fisher.test(tab), silent = TRUE)

      if (class(result) != "try-error"){
        p.value <- result$p.value

        if (max(dim(tab)) == 2){
          pt <- prop.test(tab, conf.level = conf.level)
          lower <- pt$conf.int[[1]]
          upper <- pt$conf.int[[2]]
          test <- "Fisher's Exact"
        } else {
          lower <- NA
          upper <- NA
          test <- "Fisher-Freeman-Halton's Exact"
        }
      } else {
        result <- chisq.test(tab)
        lower <- NA
        upper <- NA
        p.value <- result$p.value
        test <- "Chi-Square"
      }

      alt <- switch(alternative,
                    "greater" = ">",
                    "less" = "<",
                    "two.sided" = "!=")
      hypothesis <- paste(lg[2], "=", lg[1])
    }

  } else {
    id <- data_frame(id = rep(1:(length(var)/2), 2))
    data.test <- bind_cols(id, data.test)

    tab <- table(data.test$x, data.test$g, data.test$id)
    result <- mantelhaen.test(tab)
    p.value <- result$p.value
    test <- ifelse(dim(tab) == 2, "McNemar", "Cochran-Mantel-Haenszel")
    lower <- NA
    upper <- NA
    hypothesis <- "Marginal Heterogeneity"
  }

  out <- data_frame(Variable = var.label[[1]], Group = group.label[[1]],
                    Hypothesis = hypothesis, Lower = lower, Upper = upper,
                    Test = test, p.value) %>%
    mutate(Lower = round(.data$Lower, digits.ci),
           Upper = round(.data$Upper, digits.ci),
           'p value' = round(.data$p.value, digits.p))

  if (format){
    out <- out %>% mutate(`95% CI` = paste0("(",.data$Lower,
                                            " ; ",
                                            .data$Upper, ")")) %>%
      select(-.data$Lower, -.data$Upper, -.data$p.value)
  }

  return(out)
}


#'@importFrom forcats fct_drop
#'@importFrom dplyr bind_cols
#'@importFrom stats fisher.test chisq.test
nt_dist_ql_mg <-  function(var, group,
                           format, digits.p, var.name, var.label, group.label){

  data.test <- data_frame(x = var, g = fct_drop(group[[1]]))
  lg <- levels(data.test$g)

  tab <- table(data.test$g, data.test$x)

  if (!is.finite(max(dim(tab))) | min(dim(tab)) < 2) {
    test <- NA
    p.value <- NA
    hypothesis <- NA
    lower <- NA
    upper <- NA
  } else {
    result <- try(fisher.test(tab), silent = TRUE)

    if (class(result) != "try-error"){
      p.value <- result$p.value
      test <- "Fisher's Exact test"
    } else {
      result <- chisq.test(tab)

      p.value <- result$p.value
      test <- "Chi-Square"
    }
  }

  hypothesis <- "Association"

  out <- data_frame(Variable = var.label[[1]], Group = group.label[[1]],
                    Hypothesis = hypothesis,  Test = test,
                    `p value` =  round(p.value, digits.p))

  return(out)
}

#'@importFrom multcomp glht mcp
#'@export
nt_dist_ql_mc <-  function(var, group, alternative, contrast,
                           format, digits.p, digits.ci, var.label, group.label) {

  data.test <- data_frame(x = var, g = group[[1]])

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

  out <- data_frame(Variable = var.label, Group = group.label,
                    Hypothesis = hypothesis, Lower = lower, Upper = upper,
                    Test = test, `p value` = p.value)

  if (format){
    out <- out %>% mutate('95% CI' = paste0("(", .data$Lower, " ; ",
                                            .data$Upper, ")")) %>%
      select(.data$Variable, .data$Group, .data$Hypothesis,
             .data$Test, .data$`95% CI`, .data$`p value`)
  }

  return(out)
}
