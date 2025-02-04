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
                          exponentiate = TRUE,
                          type = "wald", contrast_qt = "one-unit",
                          user_contrast = NULL,
                          table_reference = TRUE,
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
    add_name <- names(add)
    add_label <- map2(add, add_name, extract_label)
  } else {
    add <- NULL
    add_name <- NULL
    add_label <- NULL
  }

  vars_name <- names(vars)
  response_name <- names(response)

  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars_label <- map2(.x = vars, .y = as.list(vars_name),
                       .f = extract_label)


    response <- data_labeller(response, labels)
    response_label <- extract_label(response, response_name)

  } else {
    vars_label <- map2(.x = vars, .y = as.list(vars_name),
                       .f = extract_label)
  }

  temp <- pmap(.l = list(vars, vars_name, vars_label),
               .f = aux_simple_glm,
               response = response[[1]],
               response_label = response_label,
               add = add,
               add_name = add_name,
               add_label = add_label,
               family = family,
               exponentiate = exponentiate,
               robust_variance = robust_variance,
               type = type,
               contrast_qt = contrast_qt,
               user_contrast = user_contrast,
               table_reference = table_reference,
               format = format)

  out <- Reduce(rbind, temp)

  if (format){
    out <- out |>
      # mutate(null.deviance = round(.data$null.deviance, digits),
      #                     logLik = round(.data$logLik, digits),
      #                     AIC = round(.data$AIC, digits),
      #                     BIC = round(.data$BIC, digits),
      #                     deviance = round(.data$deviance, digits)) |>
      transmute(Variable = .data$variable, Contrast = .data$group,
                `Estimate (95% CI)` = ifelse(is.na(.data$estimate) &
                                               is.na(.data$conf_low) &
                                               is.na(.data$conf_high),
                                             "Reference",
                                             paste0(round(.data$estimate, digits), " (",
                                                    round(.data$conf_low, digits), " ; ",
                                                    round(.data$conf_high, digits), ")")),
                `Wald p value` = ifelse(is.na(.data$p_value_wald), "",
                                   ifelse(round(.data$p_value_wald, digits_p) == 0, "< 0.001",
                                 as.character(round(.data$p_value_wald, digits_p)))),
                `LR p value` = ifelse(is.na(.data$p_value_lr), "",
                                        ifelse(round(.data$p_value_lr, digits_p) == 0, "< 0.001",
                                               as.character(round(.data$p_value_lr, digits_p))))
                # `LR p value` = ifelse(round(.data$p_value.lr, digits_p) == 0, "< 0.001",
                #                         as.character(round(.data$p_value.lr, digits_p))),
                # n = as.character(.data$n),
                # null.deviance = as.character(.data$null.deviance),
                # logLik = as.character(.data$logLik),
                # AIC = as.character(.data$AIC),
                # BIC = as.character(.data$BIC),
                # deviance =as.character(.data$deviance)
                ) #|>
      # replace_na(list(`Wald p value` = "", `LR p value` = "",
      #                 n = "", null.deviance = "", df.null = "",
      #                 logLik = "", AIC = "", BIC = "",
      #                 deviance = "", df.residual = ""))
  }

  if (save){
    write.csv(out, file = paste0(file, "_regression.csv"))
  }

  return(out)
}

aux_simple_glm <- function(var, var_name, var_label,
                           response, response_label,
                           add, add_name, add_label,
                           family,
                           exponentiate,
                           robust_variance,
                           type,
                           contrast_qt,
                           user_contrast,
                           table_reference,
                           format){

  if (is.factor(var)){
    var <- droplevels(var)
  }

  if (!is.null(add)){
    aux <- data.frame(add, var)
    aux_labels <- c(add_label, var = var_label)
    aux_names <- c(add_name, var_name)
    data_model <- data.frame(response = response, add = add, var = var)
  } else {
    aux <- data.frame(var)
    aux_labels <- c(var = var_label)
    aux_names <- var_name
    data_model <- data.frame(response = response, var = var)
  }


  out <- fit_simple_glm(data_model, family,
                        exponentiate, robust_variance,
                        type, contrast_qt,
                        user_contrast,
                        table_reference,
                        var_label)

  out <- out |>
    separate(.data$term, into = c("variable", "group"), sep = ":")

  if (format)
    out <- out |>
    group_by(.data$variable) |>
    mutate(aux_variable = ifelse(duplicated(.data$variable), "", .data$variable)) |>
    ungroup(.data$variable) |>
    select(-.data$variable) |>
    rename(variable = .data$aux_variable)

  return(out)
}

#'@importFrom broom tidy glance
#'@importFrom stats na.exclude glm anova
#'@importFrom lmtest coeftest
#'@importFrom stringr str_which
fit_simple_glm <- function(data_model, family,
                           exponentiate, robust_variance,
                           type, contrast_qt,
                           user_contrast,
                           table_reference,
                           var_label){

  data_model <- na.exclude(data_model)

  fit <- glm(response ~ ., data = data_model, family = family)
  aux <- extract_data(fit, data = data_model)

  ref <- reference_df(fit, data = data_model)$df
  beta <- as.numeric(fit$coefficients)
  beta_var <- as.matrix(vcov(fit))
  term_labels <- attr(fit$terms, "term.labels")

  for (i in 1:length(aux$var)){

    temp <- contrast_df(data = aux$data,
                        var = aux$var[i],
                        ref = ref,
                        contrast_qt = contrast_qt,
                        user_contrast = user_contrast,
                        table_reference = table_reference)
    temp$label <- gsub("var:", paste0(var_label, ":"), temp$label)

    design_matrix <- model.matrix(formula(fit), data = temp$new.data)

    drop <- which(grepl(aux$var[i], x = as.character(term_labels), fixed = TRUE))
    fit0 <- update(fit, paste("~ . -", term_labels[drop]))

    contrast <- contrast_calc(fit = fit, fit0 = fit0,
                              design_matrix = design_matrix,
                              beta = beta, beta_var = beta_var,
                              type = type)

    label_index <- grep("every 1 unit of change", temp$label)
    if (table_reference & length(label_index) == 0)
      contrast <- rbind(NA, contrast)

    temp <- data.frame(term = temp$label, contrast)

    if (table_reference & length(label_index) == 0)
      temp[1:2, 6] <- temp[2:1, 6]

    if (i > 1)
      temp <- rbind(out, temp)

    out <- temp

  }

  if (exponentiate)
    out[, 2:4] <- exp(out[, 2:4])

  rownames(out) <- NULL
  colnames(out) <- c("term", "estimate", "conf_low", "conf_high", "p_value_wald", "p_value_lr")

  return(out)
}



#'Multivariable Generalized Linear models
#'
#'@description Tabulating results from multivariable GLMs.
#'
#'@param fit a glm object.
#'@param exponentiate a logical value indicating whether coefficients should be exponentiated.
#'@param ci_type a character value indicating the procedure to calculate confidence intervals: likelihood ratio (\code{profile}) or wald (\code{wald}).
#'@param contrast_qt a character indicating whether the contrast for quantitative covariates. Options are every one-unit of change (\code{one-unit}), quartiles (\code{quartiles}) or provided by the user (\code{user}).
#'@param user_contrast a variable named list of numerical vectors indicating contrast for a covariate.
#'@param user_contrast_interaction a variable named list of numerical vectors indicating a contrast for interaction.
#'@param table_reference a logical value indicating whether the output should be presented with a line indicating the reference category.
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
                            ci_type = "Wald", contrast_qt = "one-unit",
                            user_contrast = NULL, user_contrast_interaction = NULL,
                            table_reference = TRUE,
                            format = TRUE, labels = NULL,
                            digits = 2, digits_p = 3,
                            save = FALSE, file = "nt_multiple_glm"){

  out <- aux_multiple_glm(fit = fit,
                               exponentiate = exponentiate,
                               ci_type = ci_type,
                               contrast_qt = contrast_qt,
                               user_contrast = user_contrast,
                               user_contrast_interaction = user_contrast_interaction,
                               format = format, table_reference = table_reference)
  ref <- reference_df(fit)$ref

  if (format)
    out$effect <- out$effect |>
    transmute(Variable = .data$variable, Group = .data$group,
              `Estimate (95% CI)` = ifelse(is.na(.data$estimate),
                                           "Reference",
                                           paste0(round(.data$estimate, digits), " (",
                                                  round(.data$conf_low, digits), " ; ",
                                                  round(.data$conf_high, digits), ")")),
              `Wald p value` = ifelse(is.na(.data$p_value_wald), "",
                                 ifelse(round(.data$p_value_wald, digits_p) == 0, "< 0.001",
                                        as.character(round(.data$p_value_wald, digits_p)))),
              `LR p value` = ifelse(is.na(.data$p_value_lr), "",
                                      ifelse(round(.data$p_value_lr, digits_p) == 0, "< 0.001",
                                             as.character(round(.data$p_value_lr, digits_p))))) |>
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

  out <- list(effect = out$effect, output = out$output, ref = ref)

  return(out)
}

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
aux_multiple_glm <- function(fit, exponentiate, robust_variance,
                             ci_type, contrast_qt,
                             user_contrast, user_contrast_interaction,
                             format, table_reference){

  aux <- extract_data(fit)

  effect <- fit_multiple_glm(fit, fit_vars = aux,
                             exponentiate = exponentiate,
                             robust_variance = robust_variance,
                             type = ci_type,
                             contrast_qt = contrast_qt,
                             user_contrast = user_contrast,
                             user_contrast_interaction = user_contrast_interaction,
                             table_reference = table_reference)
  effect <- effect |>
    separate(.data$term, into = c("variable", "group"), sep = ":")

  if (format)
    effect <- effect |>
    group_by(.data$variable) |>
    mutate(aux_variable = ifelse(duplicated(.data$variable), "", .data$variable)) |>
    ungroup(.data$variable) |>
    select(-.data$variable) |>
    rename(variable = .data$aux_variable)

  output <- summary(fit)
  statistics <- glance(fit)

  out <- list(effect = effect, output = output)

  return(out)
}

#'@importFrom stats model.matrix formula setNames anova vcov glm update.formula
fit_multiple_glm <- function(fit, fit_vars, exponentiate, robust_variance,
                             type, contrast_qt,
                             user_contrast, user_contrast_interaction,
                             table_reference){

  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta_var <- as.matrix(vcov(fit))
  term_labels <- attr(fit$terms, "term.labels")

  interaction <- colnames(attr(fit$terms, "factors"))[attr(fit$terms, "order") > 1]

  for (i in 1:length(fit_vars$var)){

    if (length(interaction) > 0){
      cond.interaction <- grepl(fit_vars$var[i], x = interaction, fixed = TRUE)
    } else {
      cond.interaction <- FALSE
    }

    if (all(!cond.interaction)){
      temp <- contrast_df(data = fit_vars$data,
                          var = fit_vars$var[i],
                          ref = ref,
                          contrast_qt = contrast_qt,
                          user_contrast = user_contrast,
                          table_reference = table_reference)
      design_matrix <- model.matrix(formula(fit), data = temp$new.data)

      drop <- which(grepl(fit_vars$var[i], x = as.character(term_labels), fixed = TRUE))
      fit0 <- glm(update.formula(fit$formula,
                                 paste0(" ~ . - ", paste(term_labels[drop],
                                                         collapse = " - "))),
                  data = na.exclude(fit_vars$data), family = fit$family)

      contrast <- contrast_calc(fit = fit, fit0 = fit0,
                                design_matrix = design_matrix,
                                beta = beta, beta_var = beta_var,
                                type = type)

      label_index <- grep("every 1 unit of change", temp$label)
      if (table_reference & length(label_index) == 0)
        contrast <- rbind(NA, contrast)

      temp <- data.frame(term = temp$label, contrast)

      if (length(label_index) == 0)
        temp[1:2, 6] <- temp[2:1, 6]

      if (i > 1)
        temp <- rbind(out, temp)

      out <- temp

    } else {
      for (k in which(cond.interaction)){

        interaction_vars <- fit_vars$var[sapply(fit_vars$var, grepl,
                                                x = as.character(interaction[k]),
                                                fixed = TRUE)]
        others <- interaction_vars[interaction_vars != fit_vars$var[i]]
        temp <- contrast_df(data = fit_vars$data,
                            var = fit_vars$var[i],
                            ref = ref,
                            contrast_qt = contrast_qt,
                            user_contrast = user_contrast,
                            interaction = others,
                            user_contrast_interaction = user_contrast_interaction,
                            table_reference = table_reference)

        design_matrix <- sapply(temp$new.data, function(x) model.matrix(fit, x), simplify = FALSE)

        drop <- which(grepl(fit_vars$var[i], x = as.character(term_labels), fixed = TRUE))
        fit0 <- glm(update.formula(fit$formula, paste0(" ~ . - ", paste(term_labels[drop], collapse = " - "))),
                    data = fit_vars$data, family = fit$family)

        contrast <- contrast_calc(fit = fit, fit0 = fit0, design_matrix = design_matrix,
                                  beta = beta, beta_var = beta_var,
                                  type = type)

        if (table_reference){
          aux.contrast <- list()
          index <- 1
          for (j in 1:length(temp$label)){

            if (j %in% temp$seq_nl){
              aux.contrast[[j]] <- matrix(NA, ncol = ncol(contrast))
              colnames(aux.contrast[[j]]) <- c("estimate", "lower", "upper", "p_value")
            } else {
              aux.contrast[[j]] <- contrast[index, ]
              index <- index + 1
            }
          }
          contrast <- Reduce(rbind, aux.contrast)
        }

        contrast <- Reduce(rbind, aux.contrast)
        label_index <- grep("every 1 unit of change", temp$label)
        temp <- data.frame(term = temp$label, contrast)

        if (length(label_index) == 0)
          temp[1:2, 6] <- temp[2:1, 6]

        if (i > 1)
          temp <- rbind(out, temp)

        out <- temp

      }
    }
  }

  if (exponentiate)
    out[, 2:4] <- exp(out[, 2:4])
  rownames(out) <- NULL
  colnames(out) <- c("term", "estimate", "conf_low", "conf_high", "p_value_wald", "p_value_lr")
  return(out)

}


