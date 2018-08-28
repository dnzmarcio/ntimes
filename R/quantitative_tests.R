#'@importFrom forcats fct_drop
#'@importFrom dplyr mutate select
#'@importFrom tidyr nest unnest
#'@importFrom purrr map
#'@importFrom magrittr %>%
#'@importFrom car leveneTest
#'@importFrom stats t.test wilcox.test
#'@importFrom nparcomp npar.t.test npar.t.test.paired
nt_dist_qt_tg <-  function(var, group, test,
                           alternative, conf.level, paired,
                           norm.test,
                           format, digits.p, digits.ci,
                           var.label, group.label) {

  data.test <- data_frame(x = var, g = fct_drop(group[[1]]))
  nlg <- nlevels(fct_drop(data.test$g[!is.na(data.test$x)]))
  lg <- levels(data.test$g)

  if (nlg == 2)  {

    if (test == "auto"){
      p.norm <- data.test %>% nest(-.data$g) %>%
        mutate(p = map(.data$data, ~ norm_test(var = .x$x, test = norm.test))) %>%
        select(.data$p) %>% unnest(.data$p)
      p.norm <- ifelse(is.na(p.norm[[1]]), 0, p.norm[[1]])

      p.var <- leveneTest(data.test$x ~ data.test$g)
      p.var <- p.var$'Pr(>F)'[[1]]

      if (all(p.norm > 0.05)) {
        if (p.var > 0.05) {
          result <- t.test(data.test$x ~ data.test$g,
                           var.equal = TRUE,
                           alternative = alternative,
                           paired = paired,
                           conf.level = conf.level)
          test <- ifelse(paired, "Paired t", "t")
        } else {
          result <- t.test(data.test$x ~ data.test$g,
                           var.equal = TRUE,
                           alternative = alternative,
                           paired = paired,
                           conf.level = conf.level)
          test <- ifelse(paired, "Paired Welch t", "Welch t")
        }

        p.value <- result$p.value
        lower <- result$conf.int[[1]]
        upper <- result$conf.int[[2]]

      } else {
        if (p.var > 0.05) {
          result <- wilcox.test(data.test$x ~ data.test$g,
                                alternative = alternative,
                                paired = paired, conf.int = TRUE,
                                conf.level = conf.level)

          p.value <- result$p.value
          lower <- result$conf.int[[1]]
          upper <- result$conf.int[[2]]
          test <- ifelse(paired, "Paired Mann-Whitney", "Mann-Whitney")

        } else {
          if (!paired){
            result <- npar.t.test(x ~ g,
                                  data = data.test,
                                  method = "t.app",
                                  rounds = max(digits.p, digits.ci),
                                  alternative = alternative,
                                  info = FALSE)

            p.value <- result$Analysis[, 6]
            lower <- ifelse(result$Analysis[, 3] < 0, 0, result$Analysis[, 3])
            upper <- ifelse(result$Analysis[, 4] > 1, 1, result$Analysis[, 4])
            test <- "Brunner-Munzel t test"
          } else {
            result <- npar.t.test.paired(x ~ g,
                                         data = data.test,
                                         rounds = max(digits.p, digits.ci),
                                         alternative = alternative,
                                         info = FALSE,
                                         plot.simci = FALSE)

            p.value <- result$Analysis[1, 5]
            lower <- ifelse(result$Analysis[1, 1] < 0, 0, result$Analysis[1, 1])
            upper <- ifelse(result$Analysis[1, 3] > 1, 1, result$Analysis[1, 3])
            test <- "Paired Brunner-Munzel t test"
          }
        }
      }
    }

    if (test == "par"){
      p.var <- leveneTest(data.test$x ~ data.test$g)
      p.var <- p.var$'Pr(>F)'[[1]]

      if (p.var > 0.05) {
        result <- t.test(data.test$x ~ data.test$g,
                         var.equal = TRUE,
                         alternative = alternative,
                         paired = paired,
                         conf.level = conf.level)
        test <- ifelse(paired, "Paired t", "t")
      } else {
        result <- t.test(data.test$x ~ data.test$g,
                         var.equal = TRUE,
                         alternative = alternative,
                         paired = paired,
                         conf.level = conf.level)
        test <- ifelse(paired, "Paired Welch t", "Welch t")
      }

      p.value <- result$p.value
      lower <- result$conf.int[[1]]
      upper <- result$conf.int[[2]]
    }

    if (test == "npar"){
      p.var <- leveneTest(data.test$x ~ data.test$g)
      p.var <- p.var$'Pr(>F)'[[1]]

      if (p.var > 0.05) {
        result <- wilcox.test(data.test$x ~ data.test$g,
                              alternative = alternative,
                              paired = paired, conf.int = TRUE,
                              conf.level = conf.level)

        p.value <- result$p.value
        lower <- result$conf.int[[1]]
        upper <- result$conf.int[[2]]
        test <- ifelse(paired, "Paired Mann-Whitney", "Mann-Whitney")

      } else {
        if (!paired){
          result <- npar.t.test(x ~ g,
                                data = data.test,
                                method = "t.app",
                                rounds = max(digits.p, digits.ci),
                                alternative = alternative,
                                info = FALSE)

          p.value <- result$Analysis[1, 6]
          lower <- ifelse(result$Analysis[, 3] < 0, 0, result$Analysis[, 3])
          upper <- ifelse(result$Analysis[, 4] > 1, 1, result$Analysis[, 4])
          test <- "Brunner-Munzel t test"
        } else {
          result <- npar.t.test.paired(x ~ g,
                                       data = data.test,
                                       alternative = alternative,
                                       rounds = max(digits.p, digits.ci),
                                       info = FALSE, plot.simci = FALSE)

          p.value <- result$Analysis[1, 5]
          lower <- ifelse(result$Analysis[1, 1] < 0, 0, result$Analysis[1, 1])
          upper <- ifelse(result$Analysis[1, 3] > 1, 1, result$Analysis[1, 3])
          test <- "Paired Brunner-Munzel t test"
        }
      }
    }

    if (test == "t"){
      result <- t.test(data.test$x ~ data.test$g,
                       var.equal = TRUE,
                       alternative = alternative,
                       paired = paired,
                       conf.level = conf.level)
      test <- ifelse(paired, "Paired t", "t")

      p.value <- result$p.value
      lower <- result$conf.int[[1]]
      upper <- result$conf.int[[2]]
    }

    if (test == "wt"){
      result <- t.test(data.test$x ~ data.test$g,
                       var.equal = TRUE,
                       alternative = alternative,
                       paired = paired,
                       conf.level = conf.level)
      test <- ifelse(paired, "Paired Welch t", "Welch t")

      p.value <- result$p.value
      lower <- result$conf.int[[1]]
      upper <- result$conf.int[[2]]
    }

    if (test == "mw"){
      result <- wilcox.test(data.test$x ~ data.test$g,
                            alternative = alternative,
                            paired = paired, conf.int = TRUE,
                            conf.level = conf.level)

      p.value <- result$p.value
      lower <- result$conf.int[[1]]
      upper <- result$conf.int[[2]]
      test <- ifelse(paired, "Paired Mann-Whitney", "Mann-Whitney")
    }

    if (test == "bm"){
      if (!paired){
        result <- npar.t.test(x ~ g,
                              data = data.test,
                              method = "t.app",
                              rounds = max(digits.p, digits.ci),
                              alternative = alternative,
                              info = FALSE)

        p.value <- result$Analysis[1, 6]
        lower <- ifelse(result$Analysis[, 3] < 0, 0, result$Analysis[, 3])
        upper <- ifelse(result$Analysis[, 4] > 1, 1, result$Analysis[, 4])
        test <- "Brunner-Munzel t test"
      } else {
        result <- npar.t.test.paired(x ~ g,
                                     data = data.test,
                                     alternative = alternative,
                                     rounds = max(digits.p, digits.ci),
                                     info = FALSE, plot.simci = FALSE)

        p.value <- result$Analysis[1, 5]
        lower <- ifelse(result$Analysis[1, 1] < 0, 0, result$Analysis[1, 1])
        upper <- ifelse(result$Analysis[1, 3] > 1, 1, result$Analysis[1, 3])
        test <- "Paired Brunner-Munzel t test"
      }
    }

    alt <- switch(alternative,
                  "two.sided" = " != ",
                  "greater" = " < ",
                  "less" = " > ")
    hypothesis <- paste(lg[2], alt, lg[1])

  } else {
    p.value <- NA
    test <- NA
    lower <- NA
    upper <- NA
    hypothesis <- NA
  }


  out <- data_frame(Variable = var.label, Group = group.label,
                    Hypothesis = hypothesis, Lower = lower, Upper = upper,
                    Test = test, p.value) %>%
    mutate(Lower = round(.data$Lower, digits.ci),
           Upper = round(.data$Upper, digits.ci),
           'p value' = round(.data$p.value, digits.p), -.data$p.value)

  if (format){
    out <- out %>% mutate('95% CI' = paste0("(", .data$Lower, " ; ",
                                           .data$Upper, ")")) %>%
      select(.data$Variable, .data$Group, .data$Hypothesis,
             .data$Test, .data$`95% CI`, .data$`p value`)
  }

  return(out)
}


#'@importFrom forcats fct_drop
#'@importFrom dplyr mutate select
#'@importFrom tidyr nest unnest
#'@importFrom purrr map
#'@importFrom magrittr %>%
#'@importFrom car leveneTest
#'@importFrom stats oneway.test kruskal.test
nt_dist_qt_mg <-  function(var, group, test,
                           norm.test, digits.p,
                           var.label, group.label) {

  data.test <- data_frame(x = var, g = fct_drop(group[[1]]))
  nlg <- nlevels(fct_drop(data.test$g[!is.na(data.test$x)]))
  lg <- levels(data.test$g)

  if (nlg > 2)  {

    if (test == "auto"){
      p.norm <- data.test %>% nest(-.data$g) %>%
        mutate(p = map(.data$data, ~ norm_test(var = .x$x, test = norm.test))) %>%
        select(.data$p) %>% unnest(.data$p)
      p.norm <- ifelse(is.na(p.norm[[1]]), 0, p.norm[[1]])

      p.var <- leveneTest(data.test$x ~ data.test$g)
      p.var <- p.var$'Pr(>F)'[[1]]

      if (all(p.norm > 0.05)) {
        if (p.var > 0.05) {
          result <- oneway.test(x ~ g, data = data.test, var.equal = TRUE)
          test <- "ANOVA"
        } else {
          result <- oneway.test(x ~ g, data = data.test, var.equal = FALSE)
          test <- "Welch-ANOVA"
        }
      } else {
        result <- kruskal.test(x ~ g, data = data.test)
        test <- "Krukal-Wallis"
      }
    }

    if (test == "par"){
      p.var <- leveneTest(data.test$x ~ data.test$g)
      p.var <- p.var$'Pr(>F)'[[1]]

      if (p.var > 0.05) {
        result <- oneway.test(x ~ g, data = data.test, var.equal = TRUE)
        test <- "ANOVA"
      } else {
        result <- oneway.test(x ~ g, data = data.test, var.equal = FALSE)
        test <- "Welch-ANOVA"
      }
    }

    if (test == "npar"){
      result <- kruskal.test(x ~ g, data = data.test)
      test <- "Krukal-Wallis"
    }

    p.value <- result$p.value
    hypothesis <- "At least one group is different"

  } else {
    p.value <- NA
    test <- NA
    hypothesis <- NA
  }

  out <- data_frame(Variable = var.label, Group = group.label,
                    Hypothesis = hypothesis, Test = test,
                    `p value` = round(p.value, digits.p))

    return(out)
}

#'@importFrom multcomp glht mcp
#'@importFrom nparcomp nparcomp
#'@importFrom stats confint aov
nt_dist_qt_mc <-  function(var, otest, group, alternative, contrast,
                           format, digits.p, digits.ci, var.label, group.label) {

  data.test <- data_frame(x = var, g = group[[1]])

  if(otest == "ANOVA"){
    av <- aov(x ~ g, data = data.test)
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

    lower <- round(confint(mc)$confint[, 2], digits.ci)
    upper <- round(confint(mc)$confint[, 3], digits.ci)

    out <- data_frame(Variable = var.label, Group = group.label,
                      Hypothesis = hypothesis, Lower = lower, Upper = upper,
                      Test = test, `p value` = p.value)
  } else{

    alternative <- switch(alternative,
                          "less" = "greater",
                          "greater" = "less",
                          "two.sided" = "two.sided")

    npc <- try(nparcomp(x ~ g, data = data.test,  type = contrast,
                        alternative = alternative, rounds = digits.p,
                        asy.method = "mult.t",
                        plot.simci = FALSE, info = FALSE),
               silent = TRUE)

    if (class(npc) != "try.error") {
      p.value <- round(npc$Analysis[,6], 3)
      lower <- ifelse(npc$Analysis[, 3] < 0, 0, round(npc$Analysis[, 3], 2))
      upper <- ifelse(npc$Analysis[, 4] > 1, 1, round(npc$Analysis[, 4], 2))
    } else {
      p.value <- NA
    }

    test <- paste("Non-parametric", contrast)
    dif <- rownames(npc$Contrast)
    alt <- switch(alternative,
                  "two.sided" = " != ",
                  "greater" = " > ",
                  "less" = " < ")
    hypothesis <- gsub(" - ", alt, dif)

    out <- data_frame(Variable = var.label, Group = group.label,
                      Hypothesis = hypothesis, Lower = lower, Upper = upper,
                      Test = test, `p value` = p.value)

  }

  if (format){
    out <- out %>% mutate('95% CI' = paste0("(", .data$Lower, " ; ",
                                            .data$Upper, ")")) %>%
      select(.data$Variable, .data$Group, .data$Hypothesis,
             .data$Test, .data$`95% CI`, .data$`p value`)
  }

  return(out)
}


