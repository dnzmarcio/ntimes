dist_qt_tg <-  function(var, var_label, group, group_label,
                        norm_test, var_test, qt_test,
                        alternative, conf_level, paired,
                        digits_p, digits_ci,
                        ...) {

  data_test <- data.frame(x = var, g = droplevels(group))
  nlg <- nlevels(droplevels(data_test$g[!is.na(data_test$x)]))
  lg <- levels(data_test$g)

  if (nlg == 2)  {

    # Checking normality
    aux <- function(x){
      result <- try(norm_test(x), silent = TRUE)
      out <- ifelse(any(is(result) == "try-error"), 0, result$p_value)
    }
    p.norm <- tapply(data_test$x, data_test$g, aux)

    # Checking homoscedasticity
    result <- try(var_test(x = data_test$x, g = data_test$g, paired = paired), silent = TRUE)
    p.var <- ifelse(any(is(result) == "try-error"), 0, result$p_value)

    # Possible cases
    if (all(p.norm > 0.05)) {
      if (p.var > 0.05) {
        result <- try(qt_test[[1]](data_test$x, data_test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf_level = conf_level), silent = TRUE)

      } else {
        result <- try(qt_test[[2]](data_test$x, data_test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf_level = conf_level), silent = TRUE)
      }

    } else {
      if (p.var > 0.05) {
        result <- try(qt_test[[3]](data_test$x, data_test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf_level = conf_level), silent = TRUE)

      } else {
        result <- try(qt_test[[4]](x = data_test$x,
                                   g = data_test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf_level = conf_level), silent = TRUE)
      }
    }

    if (any(is(result) != "try-error")){
      test <- result$test
      p_value <- result$p_value
      lower <- result$lower
      upper <- result$upper

      alt <- switch(alternative,
                    "two.sided" = " \u2260 ",
                    "greater" = " < ",
                    "less" = " > ")
      hypothesis <- paste(lg[2], alt, lg[1])

    } else {
      test <- NA
      p_value <- NA
      lower <- NA
      upper <- NA
      hypothesis <- NA
    }

  } else {
    test <- NA
    p_value <- NA
    lower <- NA
    upper <- NA
    hypothesis <- NA
  }

  out <- data.frame(Variable = var_label,
                    Group = group_label,
                    Hypothesis = hypothesis,
                    Lower = round(lower, digits_ci),
                    Upper = round(upper, digits_ci),
                    Test = test,
                    p_value =round(p_value, digits_p))

  return(out)
}

dist_qt_mg <-  function(var, group,
                        norm_test, var_test, qt_test,
                        digits_p,
                        var_label, group_label) {

  data_test <- data.frame(x = var, g = droplevels(group))
  nlg <- nlevels(droplevels(data_test$g[!is.na(data_test$x)]))
  lg <- levels(data_test$g)

  if (nlg > 2)  {

    # Checking normality
    aux <- function(x){
      result <- try(norm_test(x), silent = TRUE)
      out <- ifelse(any(is(result) == "try-error"), 0, result$p_value)
    }
    p.norm <- tapply(data_test$x, data_test$g, aux)

    # Checking homocedasticity
    result <- try(var_test(x = data_test$x, g = data_test$g, paired = FALSE), silent = TRUE)
    p.var <- ifelse(any(is(result) == "try-error"), 0, result$p_value)

    if (all(p.norm > 0.05)) {
      if (p.var > 0.05) {
        result <- try(qt_test[[1]](data_test$x, data_test$g), silent = TRUE)
      } else {
        result <- try(qt_test[[2]](data_test$x, data_test$g), silent = TRUE)
      }
    } else {
      result <- try(qt_test[[3]](data_test$x, data_test$g), silent = TRUE)
    }

    if (any(is(result) != "try-error")){
      test <- result$test
      p_value <- result$p_value
      hypothesis <- "At least one group is different"

    } else {
      test <- NA
      p_value <- NA
      hypothesis <- NA
    }

  } else {
    p_value <- NA
    test <- NA
    hypothesis <- NA
  }

  out <- data.frame(Variable = var_label, Group = group_label,
                    Hypothesis = hypothesis, Test = test,
                    p_value =round(p_value, digits_p))

  return(out)
}

#'@importFrom multcomp glht mcp
#'@importFrom nparcomp nparcomp
#'@importFrom stats confint aov
#'@importFrom methods is
dist_qt_mc <-  function(var, omnibus_test, group,
                        alternative, contrast,
                        digits_p, digits_ci,
                        var_label, group_label) {

  data_test <- data.frame(x = var, g = group)

  if(omnibus_test == "ANOVA" | omnibus_test == "Welch's ANOVA"){
    av <- aov(x ~ g, data = data_test)
    mc <- try(glht(av, linfct = mcp(g = contrast), alternative = alternative),
              silent = TRUE)

    if (any(is(mc) != "try.error")) {
      sm <- summary(mc)
      test <- paste(contrast)
      dif <- as.character(rownames(sm[2]$linfct))
      alt <- switch(alternative,
                    "two.sided" = " \u2260 ",
                    "greater" = " > ",
                    "less" = " < ")
      hypothesis <- gsub(" - ", alt, dif)
      test <- paste0("Parametric ", contrast)

      p_value <- round(sm$test$pvalues, digits_p)
      lower <- round(confint(mc)$confint[, 2], digits_ci)
      upper <- round(confint(mc)$confint[, 3], digits_ci)

    } else {
      hypothesis <- NA
      p_value <- NA
      lower <- NA
      upper <- NA
    }

    out <- data.frame(Variable = var_label, Group = group_label,
                      Hypothesis = hypothesis, Lower = lower, Upper = upper,
                      Test = test, p_value)
  } else {

    mc <- try(nparcomp(x ~ g, data = data_test,  type = contrast,
                        alternative = alternative, rounds = digits_p,
                        asy.method = "mult.t",
                        plot.simci = FALSE, info = FALSE),
               silent = TRUE)

    if (any(is(mc) != "try.error")) {
      p_value <- round(mc$Analysis[,6], 3)
      lower <- ifelse(mc$Analysis[, 3] < 0, 0, round(mc$Analysis[, 3], 2))
      upper <- ifelse(mc$Analysis[, 4] > 1, 1, round(mc$Analysis[, 4], 2))

      test <- paste("Non-parametric", contrast)
      dif <- rownames(mc$Contrast)
      alt <- switch(alternative,
                    "two.sided" = " \u2260 ",
                    "greater" = " > ",
                    "less" = " < ")
      hypothesis <- gsub(" - ", alt, dif)
    } else {
      test <- NA
      hypothesis <- NA
      p_value <- NA
      lower <- NA
      upper <- NA
    }

    out <- data.frame(Variable = var_label, Group = group_label,
                      Hypothesis = hypothesis, Lower = lower, Upper = upper,
                      Test = test, p_value)

  }
  rownames(out) <- NULL

  return(out)
}

