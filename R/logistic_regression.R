#'Simple Logistic Regression
#'
#'@description Performing simple Logistic regression.
#'
#'@param data a data frame with the variables.
#'@param response a character value indicating the response variable.
#'@param ...  character values indicating confounding variables.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits.p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@examples
#'library(titanic)
#'library(magrittr)
#'library(dplyr)
#'
#'data(titanic_train)
#'titanic_nt <- titanic_train %>%
#'  mutate(Sex = factor(Sex,
#'                      levels = c("male", "female"),
#'                      labels = c("Male", "Female")),
#'         Pclass = factor(Pclass,
#'                         from = 1:3,
#'                         levels = c("I", "II", "III"),
#'                         labels = "Passenger Class"),
#'         Embarked = factor(Embarked,
#'                           levels = c("C", "Q", "S"),
#'                           labels = c("Cherbourg", "Queenstown", "Southampton")))
#'
#'titanic_nt %>% select(Survived, Sex, Age, Pclass, Embarked) %>%
#'  nt_simple_logistic(response = Survived, Age)
#'
#'@import titanic
#'@importFrom purrr map2
#'@importFrom dplyr select
#'@importFrom rlang enquo quos .data
#'@importFrom tibble as_data_frame
#'@importFrom utils write.csv
#'
#'@export
nt_simple_logistic <- function(data, response, ..., increment = NULL,
                               format = TRUE, digits = 2, digits.p = 3,
                               save = FALSE, file = "simple_logistic"){

  response <- enquo(response)
  aux <- quos(...)

  vars <- select(.data = data, -!!response)
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

  response <- select(.data = data, !!response)
  response.name <- names(response)
  response.label <- extract_label(response, response.name)
  response <- response[[1]]

  temp <- map2(.x = vars, .y = vars.name, .f = aux_simple_logistic,
               response = response, response.label = response.label,
               add = add, add.name = add.name, add.label = add.label,
               increment = increment,
               format = format)

  out <- Reduce(rbind, temp)

  if (format){
    out <- out %>% mutate(null.deviance = round(.data$null.deviance, digits),
                          logLik = round(.data$logLik, digits),
                          AIC = round(.data$AIC, digits),
                          BIC = round(.data$BIC, digits),
                          deviance = round(.data$deviance, digits)) %>%
      transmute(Variable = .data$term, OR = .data$group,
                `Estimate (95% CI)` = paste0(round(.data$estimate, digits), " (",
                                 round(.data$conf.low, digits), " ; ",
                                 round(.data$conf.high, digits), ")"),
                `Wald p value` = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits.p))),
                `LR p value` = ifelse(round(.data$p.value.lr, digits.p) == 0, "< 0.001",
                                        as.character(round(.data$p.value.lr, digits.p))),
                n = .data$n, null.deviance = .data$null.deviance,
                logLik = .data$logLik, AIC = .data$AIC, BIC = .data$BIC,
                deviance = .data$deviance) %>%
      mutate(`Estimate (95% CI)` =
               recode(`Estimate (95% CI)`, `NA (NA ; NA)` = "Reference")) %>%
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

aux_simple_logistic <- function(var, var.name, response, response.label,
                                add, add.name, add.label, increment,
                                format){

  var.label <- extract_label(var, var.name)

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
      tab.levels[[i]] <- ifelse(is.null(increment[[aux.names[[i]]]]),
                              "every 1 unit of change",
                              paste0("every ",
                                     increment[[aux.names[i]]],
                                     " unit of change"))
    } else {
      lv <- levels(var)
      tab.levels[[i]] <- paste0(lv[2:length(lv)], "/", lv[1])
    }

    tab.labels[[names(aux.labels)[i]]] <- rep(aux.labels[[i]], length(tab.levels[[i]]))
  }

  out <- fit_logistic(data.model, tab.labels, tab.levels, var.label, increment[[var.name]])

  # if (any(duplicated(out$term))){
  #   index <- which(duplicated(out$term, fromLast = TRUE))
  #   out[(index+1):(index + length(lv) - 1), ] <- out[index:(index + length(lv) - 2), ]
  #   out[index, 3:ncol(out)] <- NA
  #   out[index:(index+length(lv)-1), 2] <- lv
  # }

  if (format){
    out$p.value.lr = ifelse(duplicated(out$term), NA, out$p.value.lr)
    out$term = ifelse(duplicated(out$term), "", out$term)
  }

  return(out)
}

#'@importFrom broom tidy glance
#'@importFrom stats na.exclude glm anova
fit_logistic <- function(data, tab.labels, tab.levels, var.label, increment){

  data <- na.exclude(data)

  if (!is.null(increment))
    #for (i in 1:length(increment)){
      data[[tab.labels[[i]]]] <- data[[tab.labels[[i]]]]/increment
    #}

  fit <- glm(response ~ ., data = data, family = "binomial")

  temp <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  temp$term <- c("(Intercept)", unlist(tab.labels))
  temp$group <- c("", unlist(tab.levels))
  temp <- temp[-1, c(1, 8, 2, 6, 7, 5)]

  p.value.lr <- rep(NA, length(unique(temp$term)))
  for (i in 2:ncol(data)){
    if (ncol(data) > 2){
      fit0 <- glm(response ~ ., data = data[, -i], family = "binomial")
      p.value.lr[(i-1)] <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]
    } else {
      fit0 <- glm(response ~ 1, data = data, family = "binomial")
      p.value.lr[(i-1)] <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]
    }
  }

  aux01 <- data.frame(term = unique(temp$term), p.value.lr = p.value.lr)
  aux02 <- data.frame(n = nrow(data), glance(fit))
  aux <- merge(aux01, aux02, by = 0, all = TRUE)[-1]

  out <- merge(temp, aux, by = "term", all = TRUE)

  return(out)
}



#'Logistic regression table
#'
#'@description Tabulating results from Logistic models.
#'
#'@param fit a fitted model.
#'@param ci.type a character value indicating the procedure to calculate confidence intervals: likelihood ratio (\code{lr}) or wald (\code{wald}).
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits.p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'
#'@examples
#'library(titanic)
#'library(magrittr)
#'library(dplyr)
#'
#'data(titanic_train)
#'dt <- titanic_train %>% mutate(Sex = ql_var(Sex,
#'                                            from = c("male", "female"),
#'                                            to = c("Male", "Female")),
#'                               Pclass = ql_var(Pclass,
#'                                               from = 1:3,
#'                                               to = c("I", "II", "III"),
#'                                               label = "Passenger Class"))
#'
#'fit <- glm(Survived ~ Age + Sex + Pclass, data = dt, family = "binomial")
#'
#'nt_multiple_logistic(fit)
#'
#'@importFrom purrr map
#'@importFrom utils write.csv
#'@export
nt_multiple_logistic <- function(fit, ci.type = "lr", user.contrast = NULL, user.contrast.interaction = NULL,
                                 format = TRUE, digits = 2, digits.p = 3,
                                 save = FALSE, file = "nt_multiple_logistic"){

  out <- aux_multiple_logistic(fit = fit, ci.type = ci.type,
                               user.contrast = user.contrast,
                               user.contrast.interaction = user.contrast.interaction,
                               format = format)
  ref <- reference_df(fit)$ref

  if (format)
    out$effect <- out$effect %>%
    transmute(Variable = .data$variable, OR = .data$or,
              `Estimate (95% CI)` = paste0(round(.data$estimate, digits), " (",
                                           round(.data$conf.low, digits), " ; ",
                                           round(.data$conf.high, digits), ")"),
              `p value` = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits.p)))) %>%
    replace_na(list(`p value` = ""))

  if (save)
    write.csv(out$effect, file = paste0(file, ".csv"))

  out <- list(effect = out$effect, coef = out$coef, ref = ref)

  return(out)
}

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
aux_multiple_logistic <- function(fit, ci.type, user.contrast, user.contrast.interaction, format){

  aux <- extract_data(fit)

  effect <- effect.glm(fit, fit.vars = aux, type = ci.type,
                       user.contrast = user.contrast,
                       user.contrast.interaction = user.contrast.interaction) %>%
    mutate(term = str_replace_all(.data$term, unlist(aux$var.labels))) %>%
    separate(.data$term, into = c("variable", "or"), sep = ":")

  if (format)
    effect <- effect %>% group_by(.data$variable) %>%
    mutate(p.value. = ifelse(duplicated(.data$variable), NA, .data$p.value),
           aux_variable = ifelse(duplicated(.data$variable), "", .data$variable)) %>%
    ungroup(.data$variable) %>% select(-.data$variable) %>%
    rename(variable = .data$aux_variable)

  temp <- unlist(aux$var.labels)
  labels <- paste0(temp, " ")
  names(labels) <- names(temp)
  coef <- tidy(fit) %>%
    mutate(term = str_replace_all(.data$term, labels),
           term = sub(" $", "", x = term))

  out <- list(effect = effect, coef = coef)

  return(out)
}

#'@importFrom stats model.matrix formula setNames anova vcov glm update.formula
effect.glm <- function(fit, fit.vars, type, user.contrast, user.contrast.interaction){

  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta.var <- as.matrix(vcov(fit))
  term.labels <- attr(fit$terms, "term.labels")

  interaction <- colnames(attr(fit$terms, "factors"))[attr(fit$terms, "order") > 1]

  for (i in 1:length(fit.vars$var)){

    if (length(interaction) > 0){
      cond.interaction <- grepl(fit.vars$var[i], x = interaction, fixed = TRUE)
    } else {
      cond.interaction <- FALSE
    }

    if (all(!cond.interaction)){
      temp <- contrast_df(data = fit.vars$data, var = fit.vars$var[i],
                          ref = ref, user.contrast = user.contrast)
      design.matrix <- model.matrix(formula(fit), data = temp$new.data)

      drop <- which(grepl(fit.vars$var[i], x = as.character(term.labels), fixed = TRUE))
      fit0 <- glm(update.formula(fit$formula,
                                 paste0(" ~ . - ", paste(term.labels[drop],
                                                         collapse = " - "))),
                  data = na.exclude(fit.vars$data), family = "binomial")

      contrast <- contrast_calc(fit = fit, fit0 = fit0,
                                design.matrix = design.matrix,
                                beta = beta, beta.var = beta.var,
                                type = type)

      temp <- data.frame(term = temp$label, contrast)

      if (i > 1)
        temp <- rbind(out, temp)

      out <- temp

    } else {
      for (k in which(cond.interaction)){

        interaction.vars <- fit.vars$var[sapply(fit.vars$var, grepl, x = as.character(interaction[k]), fixed = TRUE)]
        others <- interaction.vars[interaction.vars != fit.vars$var[i]]
        temp <- contrast_df(data = fit.vars$data, var = fit.vars$var[i],
                            ref = ref, user.contrast = user.contrast,
                            interaction = others,
                            user.contrast.interaction = user.contrast.interaction)

        design.matrix <- sapply(temp$new.data, function(x) model.matrix(fit, x), simplify = FALSE)

        drop <- which(grepl(fit.vars$var[i], x = as.character(term.labels), fixed = TRUE))
        fit0 <- glm(update.formula(fit$formula, paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
                    data = fit.vars$data, family = "binomial")

        contrast <- contrast_calc(fit = fit, fit0 = fit0, design.matrix = design.matrix,
                                  beta = beta, beta.var = beta.var,
                                  type = type)

        temp <- data.frame(term = temp$label, contrast)

        if (i > 1)
          temp <- rbind(out, temp)

        out <- temp

      }
    }
  }

  out[, 2:4] <- exp(out[, 2:4])
  colnames(out) <- c("term", "estimate", "conf.low", "conf.high", "p.value")
  return(out)

}


