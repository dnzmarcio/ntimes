#'Simple Generalized Linear Models
#'
#'@description Fit simple GLM.
#'
#'@param data a data frame with the variables.
#'@param response a character value indicating the response variable.
#'@param ...  character values indicating confounding variables.
#'@param family a character indicating family distribution. See more \code{\link[stats]{family}}.
#'@param robust_variance a function yielding a covariance matrix or a covariance matrix. See more \code{\link[lmtest]{coeftest}}.
#'@param increment a named list indicating the magnitude of increments to calculate odds ratio for continuous covariates.
#'@param ci_type a character value indicating the procedure to calculate confidence intervals: likelihood ratio (\code{profile}) or wald (\code{Wald})
#'@param conf_level a numerical value indicating the confidence level for parameters of interest.
#'@param exponentiate a logical value indicating whether coefficients should be exponentiated.
#'@param format a logical value indicating whether the output should be formatted.
#'@param labels a list of labels with components given by their variable names.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits_p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@examples
#'library(titanic)
#'library(dplyr)
#'
#'data(titanic_train)
#'titanic_nt <- titanic_train |>
#'  mutate(Sex = factor(Sex,
#'                      levels = c("male", "female"),
#'                      labels = c("Male", "Female")),
#'         Pclass = factor(Pclass,
#'                         levels = 1:3,
#'                         labels = c("I", "II", "III")),
#'         Embarked = factor(Embarked,
#'                           levels = c("C", "Q", "S"),
#'                           labels = c("Cherbourg", "Queenstown", "Southampton")))
#'
#'titanic_nt |> select(Survived, Sex, Age, Pclass, Embarked) |>
#'  nt_simple_glm(response = Survived, Age,
#'                     family = binomial(link = "logit"),
#'                     exponentiate = TRUE,
#'                     labels = list(Pclass = "Passenger class"))
#'
#'@import titanic
#'@importFrom purrr pmap map2
#'@importFrom dplyr select
#'@importFrom rlang enquo quos .data
#'@importFrom tibble as_data_frame
#'@importFrom utils write.csv
#'
#'@export
nt_simple_glm <- function(data, response, ...,
                          family, robust_variance = NULL,
                          increment = NULL,
                          ci_type = "Wald", conf_level = 0.95,
                          exponentiate = FALSE,
                          format = TRUE, labels = NULL,
                          digits = 2, digits_p = 3,
                          save = FALSE, file = "simple_logistic"){

  response <- enquo(response)
  aux <- quos(...)

  vars <- select(.data = data, -!!response)
  response <- select(.data = data, !!response)

  if (length(aux) > 0){
    for (i in 1:length(aux)){
      vars <- select(.data = vars, -!!aux[[i]])
    }
    add <- select(.data = data, !!!aux)
    add.name <- names(add)
    add.label <- map2(add, add.name, extract_label)
  } else {
    add <- NULL
    add.name <- NULL
    add.label <- NULL
  }

  vars.name <- names(vars)
  response.name <- names(response)

  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)

    response <- data_labeller(response, labels)
    response.label <- extract_label(response, response.name)
  } else {
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)
  }

  temp <- pmap(.l = list(vars, vars.name, vars.label),
               .f = aux_simple_glm,
               response = response[[1]], response.label = response.label,
               add = add, add.name = add.name, add.label = add.label,
               family = family, robust_variance = robust_variance,
               increment = increment, exponentiate = exponentiate,
               conf_level = conf_level, ci_type = ci_type,
               format = format)

  out <- Reduce(rbind, temp)

  if (format){
    out <- out |> mutate(null.deviance = round(.data$null.deviance, digits),
                          logLik = round(.data$logLik, digits),
                          AIC = round(.data$AIC, digits),
                          BIC = round(.data$BIC, digits),
                          deviance = round(.data$deviance, digits)) |>
      transmute(Variable = .data$term, Group = .data$group,
                `Estimate (95% CI)` = ifelse(.data$estimate == 1 &
                                              is.na(.data$conf.low) &
                                              is.na(.data$conf.high),
                                             "Reference",
                                 paste0(round(.data$estimate, digits), " (",
                                 round(.data$conf.low, digits), " ; ",
                                 round(.data$conf.high, digits), ")")),
                `Wald p value` = ifelse(round(.data$p.value, digits_p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits_p))),
                `LR p value` = ifelse(round(.data$p.value.lr, digits_p) == 0, "< 0.001",
                                        as.character(round(.data$p.value.lr, digits_p))),
                n = as.character(.data$n),
                null.deviance = as.character(.data$null.deviance),
                logLik = as.character(.data$logLik),
                AIC = as.character(.data$AIC),
                BIC = as.character(.data$BIC),
                deviance =as.character(.data$deviance)) |>
      mutate(`Estimate (95% CI)` =
               recode(.data$`Estimate (95% CI)`,
                      `NA (NA ; NA)` = "Reference")) |>
      replace_na(list(`Wald p value` = "", `LR p value` = "",
                      n = "", null.deviance = "", df.null = "",
                      logLik = "", AIC = "", BIC = "",
                      deviance = "", df.residual = ""))
  }

  if (save){
    write.csv(out, file = paste0(file, "_regression.csv"))
  }

  return(out)
}

aux_simple_glm <- function(var, var.name, var.label,
                           response, response.label,
                           add, add.name, add.label,
                           family, robust_variance,
                           increment, exponentiate,
                           conf_level, ci_type,
                           format){

  if (is.factor(var)){
    var <- droplevels(var)
  }

  if (!is.null(add)){
    aux <- data.frame(add, var)
    aux.labels <- c(add.label, var = var.label)
    aux.names <- c(add.name, var.name)
    data.model <- data.frame(response = response, add = add, var = var)
  } else {
    aux <- data.frame(var)
    aux.labels <- c(var = var.label)
    aux.names <- var.name
    data.model <- data.frame(response = response, var = var)
  }

  tab.labels <- list()
  tab.levels <- list()

  for (i in 1:ncol(aux)){
    if (is.numeric(aux[, i])){
      tab.levels[[names(aux.labels)[i]]] <- ifelse(is.null(increment[[aux.names[[i]]]]),
                              "every 1 unit of change",
                              paste0("every ",
                                     increment[[aux.names[i]]],
                                     " unit of change"))
    } else {
      if (!is.factor(aux[, i]))
        aux[, i] <- as.factor(aux[, i])

      lv <- levels(aux[, i])
      tab.levels[[names(aux.labels)[i]]] <- lv[1:length(lv)] #paste0(lv[2:length(lv)], "/", lv[1])
    }

    tab.labels[[names(aux.labels)[i]]] <- rep(aux.labels[[i]], length(tab.levels[[names(aux.labels)[i]]]))
  }

  out <- fit_simple_glm(data.model, family,
                        tab.labels, tab.levels, var.label,
                        robust_variance, increment[[var.name]],
                        exponentiate, conf_level, ci_type)

  if (format){
    out$p.value.lr = ifelse(duplicated(out$term), NA, out$p.value.lr)
    out$term = ifelse(duplicated(out$term), "", out$term)
  }

  return(out)
}

#'@importFrom broom tidy glance
#'@importFrom stats na.exclude glm anova
#'@importFrom lmtest coeftest
#'@importFrom stringr str_which
fit_simple_glm <- function(data, family,
                           tab.labels, tab.levels, var.label,
                           robust_variance, increment,
                           exponentiate, conf_level, ci_type){

  data <- na.exclude(data)

  if (!is.null(increment))
      data[["var"]] <- data[["var"]]/increment

  fit <- glm(response ~ ., data = data, family = family)

  if (is.null(robust_variance)){
    temp <- tidy(fit, exponentiate = exponentiate,
                 conf_level = conf_level,
                 conf.type = ci_type,
                 conf.int = TRUE)

  } else {
    step <- coeftest(fit, vcov. = robust_variance)
    temp <- tidy(step, exponentiate = exponentiate,
                 conf_level = conf_level,
                 conf.int = TRUE)
  }

  nlv <- lapply(tab.levels, length)

  for (i in 1:length(nlv)){
    if (nlv[[i]] > 1){
      add.row <- str_which(temp$term, names(tab.levels)[i])[1]
      temp <- rbind(temp[1:(add.row-1), ], NA,
                    temp[add.row:nrow(temp), ])
      temp[add.row, 1] <- paste0(names(tab.levels)[i], tab.levels[[i]][1])
      temp[add.row, 2] <- 1
    }

  }

  temp$term <- c("(Intercept)", unlist(tab.labels))
  temp$group <- c("", unlist(tab.levels))
  temp <- temp[, c(1, 8, 2, 6, 7, 5)]

  if (exponentiate)
    temp <- temp[-1, ]

  p.value.lr <- rep(NA, length(unique(unlist(tab.labels))))

  for (i in 2:ncol(data)){
    if (ncol(data) > 2){
      fit0 <- glm(response ~ ., data = data[, -i], family = family)
      p.value.lr[(i-1)] <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]
    } else {
      fit0 <- glm(response ~ 1, data = data, family = family)
      p.value.lr[(i-1)] <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]
    }
  }

  aux <- data.frame(p.value.lr = p.value.lr, n = nrow(data), glance(fit))
  out <- merge(temp, aux, by = 0, all = TRUE, sort = TRUE)[, -1]

  return(out)
}



#'Multivariable Generalized Linear models
#'
#'@description Tabulating results from multivariable GLMs.
#'
#'@param fit a glm object.
#'@param exponentiate a logical value indicating whether coefficients should be exponentiated.
#'@param ci_type a character value indicating the procedure to calculate confidence intervals: likelihood ratio (\code{profile}) or wald (\code{wald}).
#'@param contrast.qt a character indicating whether the contrast for quantitative covariates. Options are every one-unit of change (\code{one-unit}), quartiles (\code{quartiles}) or provided by the user (\code{user}).
#'@param user.contrast a variable named list of numerical vectors indicating contrast for a covariate.
#'@param user.contrast.interaction a variable named list of numerical vectors indicating a contrast for interaction.
#'@param table.reference a logical value indicating whether the output should be presented with a line indicating the reference category.
#'@param format a logical value indicating whether the output should be formatted.
#'@param labels a list of labels with components given by their variable names.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits_p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'
#'@examples
#'library(titanic)
#'library(dplyr)
#'
#'data(titanic_train)
#'dt <- titanic_train |> mutate(Sex = factor(Sex,
#'                                            levels = c("male", "female"),
#'                                            labels = c("Male", "Female")),
#'                               Pclass = factor(Pclass,
#'                                               levels = 1:3,
#'                                               labels = c("I", "II", "III"))
#'                               )
#'
#'fit <- glm(Survived ~ Age + Sex + Pclass, data = dt, family = "binomial")
#'
#'nt_multiple_glm(fit, exponentiate = TRUE)
#'
#'@importFrom purrr map
#'@importFrom utils write.csv
#'@export
nt_multiple_glm <- function(fit, exponentiate = FALSE,
                            ci_type = "Wald", contrast.qt = "one-unit",
                            user.contrast = NULL, user.contrast.interaction = NULL,
                            table.reference = TRUE,
                            format = TRUE, labels = NULL,
                            digits = 2, digits_p = 3,
                            save = FALSE, file = "nt_multiple_glm"){

  out <- aux_multiple_glm(fit = fit,
                               exponentiate = exponentiate,
                               ci_type = ci_type,
                               contrast.qt = contrast.qt,
                               user.contrast = user.contrast,
                               user.contrast.interaction = user.contrast.interaction,
                               format = format, table.reference = table.reference)
  ref <- reference_df(fit)$ref

  if (format)
    out$effect <- out$effect |>
    transmute(Variable = .data$variable, Group = .data$group,
              `Estimate (95% CI)` = ifelse(is.na(.data$estimate),
                                           "Reference",
                                           paste0(round(.data$estimate, digits), " (",
                                                  round(.data$conf.low, digits), " ; ",
                                                  round(.data$conf.high, digits), ")")),
              `p value` = ifelse(is.na(.data$p.value), "",
                                 ifelse(round(.data$p.value, digits_p) == 0, "< 0.001",
                                        as.character(round(.data$p.value, digits_p))))) |>
    replace_na(list(`p value` = ""))

  if (!is.null(labels)){
    aux_labels <- labels
    names(aux_labels) <- paste0("^", names(aux_labels), "$")

    out$effect <- out$effect |>
      mutate(Variable =
               str_replace_all(.data$Variable, unlist(aux_labels)))
  }

  if (save)
    write.csv(out$effect, file = paste0(file, ".csv"))

  out <- list(effect = out$effect, coef = out$coef, ref = ref)

  return(out)
}

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
aux_multiple_glm <- function(fit, exponentiate, robust_variance,
                             ci_type, contrast.qt,
                             user.contrast, user.contrast.interaction,
                             format, table.reference){

  aux <- extract_data(fit)

  effect <- fit_multiple_glm(fit, fit.vars = aux,
                             exponentiate = exponentiate,
                             robust_variance = robust_variance,
                             type = ci_type,
                             contrast.qt = contrast.qt,
                             user.contrast = user.contrast,
                             user.contrast.interaction = user.contrast.interaction,
                             table.reference = table.reference)
  effect <- effect |>
    separate(.data$term, into = c("variable", "group"), sep = ":")

  if (format)
    effect <- effect |> group_by(.data$variable) |>
    mutate(aux_variable = ifelse(duplicated(.data$variable), "", .data$variable)) |>
    ungroup(.data$variable) |> select(-.data$variable) |>
    rename(variable = .data$aux_variable)

  # temp <- unlist(aux$var.labels)
  # labels <- paste0(temp, " ")
  # names(labels) <- names(temp)
  # coef <- tidy(fit) |>
  #   mutate(term = str_replace_all(.data$term, labels),
  #          term = sub(" $", "", x = .data$term))

  coef <- summary(fit)

  out <- list(effect = effect, coef = coef)

  return(out)
}

#'@importFrom stats model.matrix formula setNames anova vcov glm update.formula
fit_multiple_glm <- function(fit, fit.vars, exponentiate, robust_variance,
                             type, contrast.qt,
                             user.contrast, user.contrast.interaction,
                             table.reference){

  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta_var <- as.matrix(vcov(fit))
  term_labels <- attr(fit$terms, "term.labels")

  interaction <- colnames(attr(fit$terms, "factors"))[attr(fit$terms, "order") > 1]

  for (i in 1:length(fit.vars$var)){

    if (length(interaction) > 0){
      cond.interaction <- grepl(fit.vars$var[i], x = interaction, fixed = TRUE)
    } else {
      cond.interaction <- FALSE
    }

    if (all(!cond.interaction)){
      temp <- contrast_df(data = fit.vars$data,
                          var = fit.vars$var[i],
                          ref = ref,
                          contrast.qt = contrast.qt,
                          user.contrast = user.contrast,
                          table.reference = table.reference)
      design_matrix <- model.matrix(formula(fit), data = temp$new.data)

      drop <- which(grepl(fit.vars$var[i], x = as.character(term_labels), fixed = TRUE))
      fit0 <- glm(update.formula(fit$formula,
                                 paste0(" ~ . - ", paste(term_labels[drop],
                                                         collapse = " - "))),
                  data = na.exclude(fit.vars$data), family = fit$family)

      contrast <- contrast_calc(fit = fit, fit0 = fit0,
                                design_matrix = design_matrix,
                                beta = beta, beta_var = beta_var,
                                type = type)

      aux <- grep("one-unit change", temp$label)
      if (table.reference & length(aux) == 0)
        contrast <- rbind(NA, contrast)

      temp <- data.frame(term = temp$label, contrast)

      if (type == "lr")
        temp[1:2, 5] <- temp[2:1, 5]

      if (i > 1)
        temp <- rbind(out, temp)

      out <- temp

    } else {
      for (k in which(cond.interaction)){

        interaction.vars <- fit.vars$var[sapply(fit.vars$var, grepl,
                                                x = as.character(interaction[k]),
                                                fixed = TRUE)]
        others <- interaction.vars[interaction.vars != fit.vars$var[i]]
        temp <- contrast_df(data = fit.vars$data, var = fit.vars$var[i],
                            ref = ref, user.contrast = user.contrast,
                            interaction = others,
                            user.contrast.interaction = user.contrast.interaction)

        design_matrix <- sapply(temp$new.data, function(x) model.matrix(fit, x), simplify = FALSE)

        drop <- which(grepl(fit.vars$var[i], x = as.character(term_labels), fixed = TRUE))
        fit0 <- glm(update.formula(fit$formula, paste0(" ~ . - ", paste(term_labels[drop], collapse = " - "))),
                    data = fit.vars$data, family = fit$family)

        contrast <- contrast_calc(fit = fit, fit0 = fit0, design_matrix = design_matrix,
                                  beta = beta, beta_var = beta_var,
                                  type = type)

        if (table.reference){
          aux.contrast <- list()
          index <- 1
          for (j in 1:length(temp$label)){

            if (j %in% temp$seq_nl){
              aux.contrast[[j]] <- matrix(NA, ncol = ncol(contrast))
              colnames(aux.contrast[[j]]) <- c("estimate", "lower", "upper", "p.value")
            } else {
              aux.contrast[[j]] <- contrast[index, ]
              index <- index + 1
            }
          }
          contrast <- Reduce(rbind, aux.contrast)
        }

        contrast <- Reduce(rbind, aux.contrast)

        temp <- data.frame(term = temp$label, contrast)

        if (i > 1)
          temp <- rbind(out, temp)

        out <- temp

      }
    }
  }

  if (exponentiate)
    out[, 2:4] <- exp(out[, 2:4])
  rownames(out) <- NULL
  colnames(out) <- c("term", "estimate", "conf.low", "conf.high", "p.value")
  return(out)

}


