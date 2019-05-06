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
#'  mutate(Sex = ql_var(Sex,
#'                      from = c("male", "female"),
#'                      to = c("Male", "Female")),
#'         Pclass = ql_var(Pclass,
#'                         from = 1:3,
#'                         to = c("I", "II", "III"),
#'                         label = "Passenger Class"),
#'         Embarked = ql_var(Embarked,
#'                           from = c("C", "Q", "S"),
#'                           to = c("Cherbourg", "Queenstown", "Southampton")))
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
nt_simple_logistic <- function(data, response, ...,
                               format = TRUE, digits = 2, digits.p = 3,
                               save = FALSE, file = "simple_logistic"){

  data <- as_data_frame(data)
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
               format = format)

  out <- Reduce(rbind, temp)

  if (format){
    out <- out %>% mutate(null.deviance = round(.data$null.deviance, digits),
                          logLik = round(.data$logLik, digits),
                          AIC = round(.data$AIC, digits),
                          BIC = round(.data$BIC, digits),
                          deviance = round(.data$deviance, digits)) %>%
      replace_na(list(n = "", null.deviance = "", df.null = "",
                      logLik = "", AIC = "", BIC = "",
                      deviance = "", df.residual = "")) %>%
      transmute(Variable = .data$term, Group = .data$group,
                OR.95CI = paste0(round(.data$estimate, digits), " (",
                                 round(.data$conf.low, digits), " ; ",
                                 round(.data$conf.high, digits), ")"),
                p.value = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits.p))),
                n = .data$n, null.deviance = .data$null.deviance,
                logLik = .data$logLik, AIC = .data$AIC, BIC = .data$BIC,
                deviance = .data$deviance) %>%
      replace_na(list(p.value = "")) %>%
      rename(`OR (95% CI)` = .data$OR.95CI, `p value` = .data$p.value)
  }

  if (save){
    write.csv(out, file = paste0(file, "_regression.csv"))
  }

  return(out)
}

#'@importFrom dplyr bind_cols rename
aux_simple_logistic <- function(var, var.name, response, response.label,
                                add, add.name, add.label,
                                format){

  var.label <- extract_label(var, var.name)
  aux <- cbind(add, var)
  aux.label <- c(add.label, var.label)

  if(ncol(aux) > 1) {
    var.class <- unlist(map(cbind(add, var), is.numeric))
  } else {
    var.class <- list(var = is.numeric(var))
  }

  tab.labels <- list()
  tab.levels <- list()

  for (i in 1:length(var.class)){
    if (var.class[[i]]){
      tab.labels[[colnames(aux)[i]]] <- aux.label[[i]]
      tab.levels[[colnames(aux)[i]]] <- ""
    } else {
      tab.labels[[colnames(aux)[i]]] <- paste0(aux.label[[i]], ":")
      lv <- levels(var)
      tab.levels[[colnames(aux)[i]]] <- paste0(lv[2:length(lv)], "/", lv[1])
    }
  }

  if (!is.list(tab.labels))
    tab.labels <- setNames(as.list(tab.labels), "var")

  data.model <- bind_cols(response = response, add = add, var = var)
  out <- fit_logistic(data.model, tab.labels, tab.levels, var.label)

  if (format){
    out <- out %>% mutate(term = ifelse(duplicated(.data$term), "", .data$term))
    if (is.factor(var))
      if (length(levels(var)) > 2)
        out <- out %>% mutate(p.value = .data$p.value.lh)
  }

  return(out)
}

#'@importFrom broom tidy glance
#'@importFrom tidyr separate replace_na
#'@importFrom dplyr select transmute mutate bind_cols
#'@importFrom stringr str_replace_all
#'@importFrom stats na.exclude glm anova
fit_logistic <- function(data, tab.labels, tab.levels, var.label){

  data <- na.exclude(data)

  fit <- glm(response ~ ., data = data, family = "binomial")

  temp <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(term = str_replace_all(.data$term, unlist(tab.labels))) %>%
    select(-.data$std.error, -.data$statistic) %>%
    filter(.data$term != "(Intercept)") %>%
    separate(.data$term, into = c("term", "group"), sep = ":", fill = "right") %>%
    mutate(group = unlist(tab.levels))

  fit0 <- glm(response ~ . - var, data = data, family = "binomial")
  p.value.lh <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]

  aux <- bind_cols(p.value.lh = p.value.lh, n = nrow(na.exclude(data)), glance(fit))

  out <- merge(data.frame(temp, row.names=NULL),
               data.frame(aux, row.names=NULL),
               by = 0, all = TRUE)[-1]

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
#'fit <- glm(Survived ~ Age + Sex + Pclass, data = dt)
#'
#'nt_multiple_logistic(fit)
#'
#'@importFrom purrr map
#'@importFrom utils write.csv
#'@export
nt_multiple_logistic <- function(fit, ci.type = "lr",
                                 format = TRUE, digits = 2, digits.p = 3,
                                 save = FALSE, file = "nt_multiple_logistic"){

  out <- aux_multiple_logistic(fit, format, ci.type)
  ref <- reference_df(fit)$ref

  if (format)
    out$effect <- out$effect %>%
    transmute(Variable = .data$variable, OR = .data$or,
              'Estimate (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                           round(.data$conf.low, digits), " ; ",
                                           round(.data$conf.high, digits), ")"),
              'p value LR' = ifelse(round(.data$p.value.lr, digits.p) == 0, "< 0.001",
                                    as.character(round(.data$p.value.lr, digits.p)))) %>%
    replace_na(list('p value LR' = ""))

  if (save)
    write.csv(out$effect, file = paste0(file, ".csv"))

  out <- list(effect = out$effect, coef = out$coef, ref = ref)

  return(out)
}

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
aux_multiple_logistic <- function(fit, format, ci.type){

  aux <- extract_data(fit)

  effect <- effect.glm(fit, type = ci.type, exponentiate = TRUE) %>%
    mutate(term = str_replace_all(.data$term, unlist(aux$var.labels))) %>%
    separate(.data$term, into = c("variable", "or"), sep = ":")

  if (format)
    effect <- effect %>% group_by(.data$variable) %>%
    mutate(aux_variable = ifelse(duplicated(.data$variable), "", .data$variable),
           p.value = ifelse(duplicated(.data$p.value), NA, .data$p.value)) %>%
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
effect.glm <- function(fit, type, exponentiate){

  aux <- extract_data(fit)
  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta.var <- as.matrix(vcov(fit))
  term.labels <- attr(fit$terms, "term.labels")

  interaction <- colnames(attr(fit$terms, "factors"))[attr(fit$terms, "order") > 1]

  for (i in 1:length(aux$var)){

    if (length(interaction) > 0){
      cond.interaction <- grepl(aux$var[i], x = interaction, fixed = TRUE)
    } else {
      cond.interaction <- FALSE
    }

    if (all(!cond.interaction)){
      temp <- contrast_df(aux$data, aux$var[i], ref)
      design.matrix <- model.matrix(formula(fit), data = temp$new.data)

      drop <- which(grepl(aux$var[i], x = as.character(term.labels), fixed = TRUE))
      fit0 <- glm(update.formula(fit$formula, paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
                  data = aux$data, family = "binomial")
      p.value <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]

      contrast <- contrast_calc(fit = fit, design.matrix = design.matrix,
                                beta = beta, beta.var = beta.var,
                                p.value = p.value, type = type)

      temp <- data.frame(term = temp$label, contrast)

      if (i > 1)
        temp <- rbind(out, temp)

      out <- temp

    } else {
      for (k in which(cond.interaction)){
        interaction.vars <- aux$var[sapply(aux$var, grepl, x = as.character(interaction[k]), fixed = TRUE)]
        others <- interaction.vars[interaction.vars != aux$var[i]]
        temp <- contrast_df(aux$data, aux$var[i], ref, others)
        design.matrix <- sapply(temp$new.data, function(x) model.matrix(fit, x), simplify = FALSE)

        drop <- which(grepl(aux$var[i], x = as.character(term.labels), fixed = TRUE))
        fit0 <- glm(update.formula(fit$formula, paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
                    data = aux$data, family = "binomial")
        p.value <- anova(fit0, fit, test = "Chisq")$`Pr(>Chi)`[2]

        contrast <- contrast_calc(fit = fit, design.matrix = design.matrix,
                                  beta = beta, beta.var = beta.var,
                                  p.value = p.value, type = type)

        temp <- data.frame(term = temp$label, contrast)

        if (i > 1)
          temp <- rbind(out, temp)

        out <- temp

      }
    }
  }

  if (exponentiate)
    out[, apply(out, 2, is.numeric)] <- exp(out[, apply(out, 2, is.numeric)])

  colnames(out) <- c("term", "estimate", "conf.low", "conf.high", "p.value.lr")

  return(out)

}



