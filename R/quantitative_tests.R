dist_qt_tg <-  function(var, group,
                           norm.test, var.test, distr.test,
                           alternative, conf.level, paired,
                           digits.p, digits.ci,
                           var.label, group.label) {

  data.test <- data.frame(x = var, g = droplevels(group))
  nlg <- nlevels(droplevels(data.test$g[!is.na(data.test$x)]))
  lg <- levels(data.test$g)

  if (nlg == 2)  {

    # Checking normality
    aux <- function(x){
      result <- try(norm.test(x), silent = TRUE)
      out <- ifelse(class(result) == "try-error", 0, result$p.value)
    }
    p.norm <- tapply(data.test$x, data.test$g, aux)

    # Checking homocedasticity
    result <- try(var.test(data.test$x, data.test$g), silent = TRUE)
    p.var <- ifelse(class(result) == "try-error", 0, result$p.value)

    # Possible cases
    if (all(p.norm > 0.05)) {
      if (p.var > 0.05) {
        result <- try(distr.test[[1]](data.test$x, data.test$g,
                                  alternative = alternative,
                                  paired = paired,
                                  conf.level = conf.level), silent = TRUE)

      } else {
        result <- try(distr.test[[2]](data.test$x, data.test$g,
                                  alternative = alternative,
                                  paired = paired,
                                  conf.level = conf.level), silent = TRUE)
      }

    } else {
      if (p.var > 0.05) {
        result <- try(distr.test[[3]](data.test$x, data.test$g,
                                  alternative = alternative,
                                  paired = paired,
                                  conf.level = conf.level), silent = TRUE)

      } else {
        result <- try(distr.test[[4]](data.test$x, data.test$g,
                                  alternative = alternative,
                                  paired = paired,
                                  conf.level = conf.level), silent = TRUE)
      }
    }

    if (class(result) != "try-error"){
      test <- result$test
      p.value <- result$p.value
      lower <- result$lower
      upper <- result$upper

      alt <- switch(alternative,
                    "two.sided" = " != ",
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
                        norm.test, var.test, distr.test,
                        digits.p,
                        var.label, group.label) {

  data.test <- data.frame(x = var, g = droplevels(group))
  nlg <- nlevels(droplevels(data.test$g[!is.na(data.test$x)]))
  lg <- levels(data.test$g)

  if (nlg > 2)  {

    # Checking normality
    aux <- function(x){
      result <- try(norm.test(x), silent = TRUE)
      out <- ifelse(class(result) == "try-error", 0, result$p.value)
    }
    p.norm <- tapply(data.test$x, data.test$g, aux)

    # Checking homocedasticity
    result <- try(var.test(data.test$x, data.test$g), silent = TRUE)
    p.var <- ifelse(class(result) == "try-error", 0, result$p.value)

    if (all(p.norm > 0.05)) {
      if (p.var > 0.05) {
        result <- try(distr.test[[1]](data.test$x, data.test$g), silent = TRUE)
      } else {
        result <- try(distr.test[[2]](data.test$x, data.test$g), silent = TRUE)
      }
    } else {
      result <- try(distr.test[[3]](data.test$x, data.test$g), silent = TRUE)
    }

    if (class(result) != "try-error"){
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
dist_qt_mc <-  function(var, omnibus.test, group,
                        alternative, contrast,
                        digits.p, digits.ci,
                        var.label, group.label) {

  data.test <- data.frame(x = var, g = group)

  if(omnibus.test == "ANOVA" | omnibus.test == "Welch's ANOVA"){
    av <- aov(x ~ g, data = data.test)
    mc <- try(glht(av, linfct = mcp(g = contrast), alternative = alternative),
              silent = TRUE)

    if (class(mc) != "try.error") {
      sm <- summary(mc)
      test <- paste(contrast)
      dif <- as.character(rownames(sm[2]$linfct))
      alt <- switch(alternative,
                    "two.sided" = " != ",
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

    if (class(mc) != "try.error") {
      p.value <- round(mc$Analysis[,6], 3)
      lower <- ifelse(mc$Analysis[, 3] < 0, 0, round(mc$Analysis[, 3], 2))
      upper <- ifelse(mc$Analysis[, 4] > 1, 1, round(mc$Analysis[, 4], 2))

      test <- paste("Non-parametric", contrast)
      dif <- rownames(mc$Contrast)
      alt <- switch(alternative,
                    "two.sided" = " != ",
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
  return(out)
}

#'@export
nt_norm_test <-  function(test){

  if (test == "ad")
    out <- function(x) {
      out <- nortest::ad.test(x)
      return(out)
    }
  if (test == "sf")
    out <- function(x) {
      out <- nortest::sf.test(x)
      return(out)
    }
  if (test == "ks")
    out <- function(x) {
      out <- nortest::lillie.test(x)
      return(out)
    }
  if (test == "cvm")
    out <- function(x) {
      out <- nortest::cvm.test(x)
      return(out)
    }
  if (test == "ps")
    out <- function(x) {
      out <- nortest::pearson.test(x)
      return(out)
    }

  return(out)
}

#'@export
nt_var_test <- function(test, paired = FALSE){

  if (test == "levene"){
    if (!paired){
      out <- function(x, g) {
        out <- lawstat::levene.test(x, g)
        return(out)
      }
    } else {
      out <- function(x, g) {
        lg <- levels(g)
        out <- Paired::levene.Var.test(x[g == lg[1]], x[g == lg[2]])
        return(out)
      }
    }
  }

  return(out)
}

#'@export
nt_anova <- function(x, g){

  data.test <- data.frame(x, g)
  result <- stats::oneway.test(x ~ g, data = data.test, var.equal = TRUE)
  test <- "ANOVA"
  p.value <- result$p.value

  out <- list(test = test, p.value = p.value)
}

#'@export
nt_welch_anova <- function(x, g){

  data.test <- data.frame(x, g)
  result <- stats::oneway.test(x ~ g, data = data.test, var.equal = FALSE)
  test <- "Welch's ANOVA"
  p.value <- result$p.value

  out <- list(test = test, p.value = p.value)
}

#'@export
nt_kruskal_wallis <- function(x, g){
  data.test <- data.frame(x, g)
  result <- stats::kruskal.test(x ~ g, data = data.test)
  test <- "Kruskal-Wallis"
  p.value <- result$p.value

  out <- list(test = test, p.value = p.value)
}

#'@export
nt_student_t <- function(x, g, paired, alternative, conf.level){

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

#'@export
nt_welch_t <- function(x, g, paired, alternative, conf.level){

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


#'@export
nt_mann_whitney <- function(x, g, paired, alternative, conf.level){

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

#'@export
nt_brunner_munzel <- function(x, g, alternative, paired, conf.level){

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
