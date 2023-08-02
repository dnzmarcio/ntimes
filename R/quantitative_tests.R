dist_qt_tg <-  function(var, var.label, group, group.label,
                        norm.test, var.test, qt.test,
                        alternative, conf.level, paired,
                        digits.p, digits.ci,
                        ...) {

  data.test <- data.frame(x = var, g = droplevels(group))
  nlg <- nlevels(droplevels(data.test$g[!is.na(data.test$x)]))
  lg <- levels(data.test$g)

  if (nlg == 2)  {

    # Checking normality
    aux <- function(x){
      result <- try(norm.test(x), silent = TRUE)
      out <- ifelse(any(is(result) == "try-error"), 0, result$p.value)
    }
    p.norm <- tapply(data.test$x, data.test$g, aux)

    # Checking homoscedasticity
    result <- try(var.test(x = data.test$x, g = data.test$g, paired = paired), silent = TRUE)
    p.var <- ifelse(any(is(result) == "try-error"), 0, result$p.value)

    # Possible cases
    if (all(p.norm > 0.05)) {
      if (p.var > 0.05) {
        result <- try(qt.test[[1]](data.test$x, data.test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf.level = conf.level), silent = TRUE)

      } else {
        result <- try(qt.test[[2]](data.test$x, data.test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf.level = conf.level), silent = TRUE)
      }

    } else {
      if (p.var > 0.05) {
        result <- try(qt.test[[3]](data.test$x, data.test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf.level = conf.level), silent = TRUE)

      } else {
        result <- try(qt.test[[4]](x = data.test$x,
                                   g = data.test$g,
                                   alternative = alternative,
                                   paired = paired,
                                   conf.level = conf.level), silent = TRUE)
      }
    }

    if (any(is(result) != "try-error")){
      test <- result$test
      p.value <- result$p.value
      lower <- result$lower
      upper <- result$upper

      alt <- switch(alternative,
                    "two.sided" = " \u2260 ",
                    "greater" = " < ",
                    "less" = " > ")
      hypothesis <- paste(lg[2], alt, lg[1])

    } else {
      test <- NA
      p.value <- NA
      lower <- NA
      upper <- NA
      hypothesis <- NA
    }

  } else {
    test <- NA
    p.value <- NA
    lower <- NA
    upper <- NA
    hypothesis <- NA
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

dist_qt_mg <-  function(var, group,
                        norm.test, var.test, qt.test,
                        digits.p,
                        var.label, group.label) {

  data.test <- data.frame(x = var, g = droplevels(group))
  nlg <- nlevels(droplevels(data.test$g[!is.na(data.test$x)]))
  lg <- levels(data.test$g)

  if (nlg > 2)  {

    # Checking normality
    aux <- function(x){
      result <- try(norm.test(x), silent = TRUE)
      out <- ifelse(any(is(result) == "try-error"), 0, result$p.value)
    }
    p.norm <- tapply(data.test$x, data.test$g, aux)

    # Checking homocedasticity
    result <- try(var.test(x = data.test$x, g = data.test$g, paired = FALSE), silent = TRUE)
    p.var <- ifelse(any(is(result) == "try-error"), 0, result$p.value)

    if (all(p.norm > 0.05)) {
      if (p.var > 0.05) {
        result <- try(qt.test[[1]](data.test$x, data.test$g), silent = TRUE)
      } else {
        result <- try(qt.test[[2]](data.test$x, data.test$g), silent = TRUE)
      }
    } else {
      result <- try(qt.test[[3]](data.test$x, data.test$g), silent = TRUE)
    }

    if (any(is(result) != "try-error")){
      test <- result$test
      p.value <- result$p.value
      hypothesis <- "At least one group is different"

    } else {
      test <- NA
      p.value <- NA
      hypothesis <- NA
    }

  } else {
    p.value <- NA
    test <- NA
    hypothesis <- NA
  }

  out <- data.frame(Variable = var.label, Group = group.label,
                    Hypothesis = hypothesis, Test = test,
                    p.value = round(p.value, digits.p))

  return(out)
}

#'@importFrom multcomp glht mcp
#'@importFrom nparcomp nparcomp
#'@importFrom stats confint aov
#'@importFrom methods is
dist_qt_mc <-  function(var, omnibus.test, group,
                        alternative, contrast,
                        digits.p, digits.ci,
                        var.label, group.label) {

  data.test <- data.frame(x = var, g = group)

  if(omnibus.test == "ANOVA" | omnibus.test == "Welch's ANOVA"){
    av <- aov(x ~ g, data = data.test)
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

      p.value <- round(sm$test$pvalues, digits.p)
      lower <- round(confint(mc)$confint[, 2], digits.ci)
      upper <- round(confint(mc)$confint[, 3], digits.ci)

    } else {
      hypothesis <- NA
      p.value <- NA
      lower <- NA
      upper <- NA
    }

    out <- data.frame(Variable = var.label, Group = group.label,
                      Hypothesis = hypothesis, Lower = lower, Upper = upper,
                      Test = test, p.value)
  } else {

    mc <- try(nparcomp(x ~ g, data = data.test,  type = contrast,
                        alternative = alternative, rounds = digits.p,
                        asy.method = "mult.t",
                        plot.simci = FALSE, info = FALSE),
               silent = TRUE)

    if (any(is(mc) != "try.error")) {
      p.value <- round(mc$Analysis[,6], 3)
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
      p.value <- NA
      lower <- NA
      upper <- NA
    }

    out <- data.frame(Variable = var.label, Group = group.label,
                      Hypothesis = hypothesis, Lower = lower, Upper = upper,
                      Test = test, p.value)

  }
  rownames(out) <- NULL

  return(out)
}

