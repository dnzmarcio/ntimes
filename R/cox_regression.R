#'Simple Cox Regression
#'
#'@description Performing simple Cox regression.
#'
#'@param data a data frame with the variables.
#'@param time a numeric vector with the follow-up time.
#'@param status a numeric vector indicating status, 0 = censored, 1 = event at time.
#'@param ...  character values indicating confounding variables.
#'@param labels a list of labels with components given by their variable names.
#'@param increment a named list indicating the magnitude of increments to calculate odds ratio for continuous covariates.
#'@param cluster a character vector containing the cluster variable.
#'@param strata a character vector containing the strata variable.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits_p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'
#'@examples
#'library(survival)
#'library(dplyr)
#'data(ovarian)
#'
#'ovarian_nt <- ovarian |> mutate(resid.ds = factor(resid.ds,
#'                                                   levels = 1:2,
#'                                                   labels = c("no", "yes")),
#'                              ecog.ps = factor(ecog.ps,
#'                                               levels = 1:2,
#'                                               labels = c("I", "II")),
#'                              rx =factor(rx,
#'                                          levels = 1:2,
#'                                          labels = c("t1", "t2")))
#'ovarian_nt |> nt_simple_cox(time = futime, status = fustat,
#'                             labels = list(resid.ds = "Residual Disease",
#'                                           ecog.ps = "ECOG-PS",
#'                                           rx = "Treatment",
#'                                           age = "Age"))
#'
#'@importFrom rlang enquo quos quo_get_expr .data
#'@importFrom dplyr select mutate transmute rename
#'@importFrom purrr pmap map2
#'@importFrom survival survfit
#'@importFrom broom tidy
#'@importFrom tidyr replace_na
#'@importFrom utils write.csv
#'
#'@export
nt_simple_cox <- function(data, time, status, ...,
                          cluster = FALSE, strata = NULL,
                          exponentiate = TRUE,
                          type = "wald",
                          contrast_qt = "one-unit",
                          user_contrast = NULL,
                          table_reference = TRUE,
                          format = TRUE, labels = NULL,
                          digits = 2, digits_p = 3,
                          save = FALSE, file = "simple_cox"){

  time <- enquo(time)
  status <- enquo(status)
  strata <- enquo(strata)
  aux <- quos(...)

  vars <- select(.data = data, -!!time)
  vars <- select(.data = vars, -!!status)

  if (!is.null(quo_get_expr(strata))){
    vars <- select(.data = vars, -!!strata)
    strata_var <- select(.data = data, !!strata)
    strata_var <- strata_var[[1]]
  } else {
    strata_var <- NULL
  }

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

  vars.name <- names(vars)

  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)
  } else {
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)
  }

  time <- select(.data = data, !!time)
  time <- time[[1]]
  status <- select(.data = data, !!status)
  status <- as.numeric(as.factor(status[[1]])) - 1

  fit <- survfit(Surv(time, status) ~ 1)
  survival <- tidy(fit) |> select(-.data$std.error) |>
    mutate(estimate = round(100*.data$estimate, digits),
           conf_low = round(100*.data$conf.low, digits),
           conf_high = round(100*.data$conf.high, digits)) |>
    transmute(Time = .data$time, .data$n.risk, .data$n.event, .data$n.censor,
              `Survival (95% CI)` = paste0(.data$estimate, " (",
                                           .data$conf_low, " ; ",
                                           .data$conf_high, ")"))

  temp <- pmap(.l = list(vars, vars.name, vars.label),
               .f = aux_simple_cox,
               time = time,
               status = status,
               add = add,
               add_name = add_name,
               add_label = add_label,
               strata_var = strata_var,
               exponentiate = exponentiate,
               type = type,
               contrast_qt = contrast_qt,
               user_contrast = user_contrast,
               table_reference = table_reference,
               format = format)

  out <- bind_rows(temp)

  if (format) {
    out <- out |>
      # mutate(concordance = round(.data$concordance, digits),
      #        r.squared = round(.data$r.squared, digits),
      #        AIC = round(.data$AIC, digits),
      #        ph.assumption = round(.data$ph.assumption, digits_p)) |>
      transmute(Variable = .data$variable,
                Contrast = .data$group,
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
                # ,
                # n = .data$n, n.event = .data$n.event,
                # concordance = .data$concordance, r.squared = .data$r.squared,
                # AIC = .data$AIC, ph.assumption  = .data$ph.assumption
                )
  }

  if (save){
    write.csv(out, file = paste0(file, "_regression.csv"))
    write.csv(survival, file = paste0(file, "_survival.csv"))
  }

  out <- list(survival = survival, effect = out)

  return(out)
}

#'@importFrom purrr map
#'@importFrom stats setNames
#'@importFrom dplyr mutate bind_cols
#'@importFrom rlang .data
aux_simple_cox <- function(var, var_name, var_label,
                           time, status,
                           add, add_name, add_label,
                           strata_var,
                           exponentiate,
                           type,
                           contrast_qt,
                           user_contrast,
                           table_reference,
                           format){

  var_label <- extract_label(var, var_name)

  if (!is.null(add)){
    aux <- data.frame(add, var)
    aux.labels <- c(add_label, var = var_label)
    aux.names <- c(add_name, var_name)
  } else {
    aux <- data.frame(var)
    aux.labels <- c(var = var_label)
    aux.names <- var_name
  }

  data_model <- bind_cols(time = time, status = status, var = var, add = add)

  out <- fit_simple_cox(data_model, strata_var,
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

#'@importFrom survival coxph Surv cox.zph
#'@importFrom broom tidy glance
#'@importFrom tidyr separate replace_na
#'@importFrom dplyr select mutate
#'@importFrom stats na.exclude update.formula anova
#'@importFrom stringr str_replace_all
#'@importFrom methods is
fit_simple_cox <- function(data_model, strata_var,
                           exponentiate, robust_variance,
                           type, contrast_qt,
                           user_contrast,
                           table_reference,
                           var_label){

  if (any(is.na(data_model)))
    strata_var <- strata_var[-which(is.na(data), arr.ind = TRUE)[, 1]]

  data_model <- na.exclude(data_model)

  if (is.null(strata_var)){
    fit <- try(coxph(Surv(time, status) ~ ., data = data_model), silent = TRUE)
  } else {
    fit <- try(coxph(Surv(time, status) ~ ., data = data_model), silent = TRUE)
  }

  aux <- extract_data(fit, data = data_model)
  ref <- reference_df(fit, data = data_model)$df
  beta <- as.numeric(fit$coefficients)
  beta_var <- as.matrix(vcov(fit))
  term_labels <- attr(fit$terms, "term.labels")

  for (i in 1:length(aux$var)){

    temp <- contrast_df(data = aux$data, var = aux$var[i],
                        ref = ref,
                        contrast_qt = contrast_qt,
                        user_contrast = user_contrast,
                        table_reference = table_reference)
    temp$label <- gsub("var:", paste0(var_label, ":"), temp$label)

    design_matrix <- model.matrix(fit, temp$new.data)

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

    if (length(label_index) == 0)
      temp[1:2, 6] <- temp[2:1, 6]

    if (i > 1)
      temp <- rbind(out, temp)

    out <- temp
  }

  if (exponentiate)
    out[, 2:4] <- exp(out[, 2:4])

  rownames(out) <- NULL
  colnames(out) <- c("term", "estimate", "conf_low", "conf_high", "p_value_wald", "p_value_lr")

  # if (is(fit) != "try-error"){
  #   aux01 <- data.frame(p_value.lr = p_value.lr)
  #   zph.table <- cox.zph(fit)$table
  #
  #   aux02 <- glance(fit) |> select(.data$n, n.event = .data$nevent,
  #                                 .data$concordance, .data$r.squared, .data$AIC) |>
  #     mutate(ph.assumption = zph.table[nrow(zph.table), 3])
  #
  #   aux <- merge(aux01, aux02, by = 0, all = TRUE)[-1]
  #   out <- merge(temp, aux, by = 0, all = TRUE)[, -1]
  # } else {
  #   out <- data.frame(term = tab.labels$var, group = NA, estimate = NA,
  #                     p_value = NA, conf.high = NA, p_value = NA,
  #                     p_value.lr = NA, n = NA, n.event = NA,
  #                     concordance = NA, r.squared = NA, AIC = NA, ph.assumption = NA)
  # }

  return(out)
}



#'Proportional Hazards Cox regression table
#'
#'@description Tabulating results from fitted Proportional Hazards Cox models.
#'
#'@param fit a coxph object.
#'@param ci_type a character value indicating the procedure to calculate confidence intervals: likelihood ratio (\code{lr}) or wald (\code{wald}).
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
#'library(survival)
#'library(dplyr)
#'
#'data(ovarian)
#'dt <- ovarian |> mutate(resid.ds = factor(resid.ds,
#'                                                levels = 1:2,
#'                                                labels = c("no", "yes")),
#'                              ecog.ps = factor(ecog.ps,
#'                                               levels = 1:2,
#'                                               labels  = c("I", "II")),
#'                              rx = factor(rx,
#'                                          levels = 1:2,
#'                                          labels = c("t1", "t2")))
#'
#'
#'fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps*rx, data = dt)
#'
#'nt_multiple_cox(fit)
#'
#'@importFrom purrr map2 map
#'@importFrom utils write.csv
#'@importFrom dplyr transmute bind_rows
#'@importFrom tidyr replace_na
#'@importFrom methods is
#'@export
nt_multiple_cox <- function(fit,
                            ci_type = "wald",
                            contrast_qt = "one-unit",
                            user_contrast = NULL,
                            user_contrast_interaction = NULL,
                            table_reference = TRUE,
                            format = TRUE, labels = NULL,
                            digits = 2, digits_p = 3,
                            save = FALSE, file = "nt_multiple_cox"){

  if (any(is(fit) != "coxph"))
    stop("fit object is not a coxph class")

  out <- aux_multiple_cox(fit = fit, ci_type = ci_type,
                          contrast_qt = contrast_qt,
                          user_contrast = user_contrast,
                          user_contrast_interaction = user_contrast_interaction,
                          format = format,
                          table_reference = table_reference)
  ref <- reference_df(fit)$ref

  if (format){
    out$effect <-  out$effect |>
    transmute(Variable = .data$variable, HR = .data$hr,
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

  } else {
    if (!is.null(labels)){
      aux_labels <- labels
      names(aux_labels) <- paste0("^", names(aux_labels), "$")

      out$effect <- out$effect |>
        mutate(Variable =
                 str_replace_all(.data$variable, unlist(aux_labels)))
    }
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
aux_multiple_cox <- function(fit, ci_type, contrast_qt,
                             user_contrast, user_contrast_interaction,
                             format, table_reference){

  aux <- extract_data(fit)

  effect <- fit_multiple_cox(fit, fit_vars = aux, type = ci_type,
                             contrast_qt = contrast_qt,
                         user_contrast = user_contrast,
                         user_contrast_interaction = user_contrast_interaction,
                         table_reference = table_reference)

  effect <- effect |>
    separate(.data$term, into = c("variable", "hr"), sep = ":")

  if (format)
    effect <- effect |> group_by(.data$variable) |>
    mutate(aux_variable = ifelse(duplicated(.data$variable), "",
                                 .data$variable)) |>
    ungroup(.data$variable) |> select(-.data$variable) |>
    rename(variable = .data$aux_variable)

  # temp <- unlist(aux$var_labels)
  # labels <- paste0(temp, " ")
  # names(labels) <- names(temp)
  # coef <- tidy(fit) |>
  #   mutate(term = str_replace_all(.data$term, labels),
  #          term = sub(" $", "", x = .data$term))

  output <- summary(fit)

  out <- list(effect = effect, output = output)

  return(out)
}


#'@importFrom stats model.matrix formula setNames anova vcov update.formula
#'@importFrom survival coxph
#'@importFrom stringr str_split
fit_multiple_cox <- function(fit, fit_vars, type, contrast_qt,
                         user_contrast, user_contrast_interaction,
                         table_reference){

  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta_var <- as.matrix(vcov(fit))
  term_labels <- attr(fit$terms, "term.labels")

  interaction <- colnames(attr(fit$terms, "factors"))[attr(fit$terms, "order") > 1]

  if (any(attr(fit$terms, "order") > 2)){
    temp <- str_split(interaction, ":")
    temp <- sapply(temp, FUN =
                     function(x) sapply(temp, FUN = function(y) x %in% y))
    index <- unique(sapply(temp, FUN =
                             function(x) max(which(apply(x, FUN =
                                                           function(x) all(x), 2)))))
    interaction <- interaction[index]
  }

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
      design_matrix <- model.matrix(fit, temp$new.data)

      drop <- which(grepl(fit_vars$var[i], x = as.character(term_labels), fixed = TRUE))
      fit0 <- coxph(update.formula(fit$formula, paste0(" ~ . - ", paste(term_labels[drop], collapse = " - "))),
                    data = na.exclude(fit_vars$data))

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

    } else {
      for (k in which(cond.interaction)){
        interaction_vars <- fit_vars$var[sapply(fit_vars$var, grepl, x = as.character(interaction[k]), fixed = TRUE)]
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
        fit0 <- coxph(update.formula(fit$formula,
                                     paste0(" ~ . - ", paste(term_labels[drop], collapse = " - "))),
                      data = fit_vars$data)

        contrast <- contrast_calc(fit = fit, fit0 = fit0, design_matrix = design_matrix,
                                  beta = beta, beta_var = beta_var,
                                  type = type)

        if (table_reference){
          aux.contrast <- list()
          index <- 1
          for (j in 1:length(temp$label)){

            if (j %in% temp$seq_nl){
              aux.contrast[[j]] <- matrix(NA, ncol = ncol(contrast))
              colnames(aux.contrast[[j]]) <- c("estimate", "lower", "upper", "p_value_wald", "p_value_lr")
            } else {
              aux.contrast[[j]] <- contrast[index, ]
              index <- index + 1
            }
          }
          contrast <- Reduce(rbind, aux.contrast)
        }

        temp <- data.frame(term = temp$label, contrast)
        temp[1:2, 6] <- temp[2:1, 6]

        if (i > 1)
          temp <- rbind(out, temp)

        out <- temp

      }
    }
  }

  out[, 2:4] <- exp(out[, 2:4])
  rownames(out) <- NULL
  colnames(out) <- c("term", "estimate", "conf_low", "conf_high", "p_value_wald", "p_value_lr")

  return(out)
}


