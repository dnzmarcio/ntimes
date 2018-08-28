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
#'titanic_train <- titanic_train %>%
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
#'titanic_train %>% select(Survived, Sex, Age, Pclass, Embarked) %>%
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
               digits = digits, digits.p = digits.p, format = format)

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
                                digits, digits.p, format){

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
  out <- fit_logistic(data.model, tab.labels, tab.levels, var.label) %>%
    mutate(term = ifelse(duplicated(.data$term), "", .data$term))
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

  if (length(tab.levels[["var"]]) > 1)
    temp[temp$term == var.label, ]$p.value <-
    c(p.value.lh, rep(NA, (length(tab.levels[["var"]]) - 1)))

  aux <- bind_cols(n = nrow(na.exclude(data)), glance(fit))

  out <- merge(data.frame(temp, row.names=NULL),
               data.frame(aux, row.names=NULL),
               by = 0, all = TRUE)[-1]

  return(out)
}



#'Logistic regression table
#'
#'@description Tabulating results from Logistic models.
#'
#'@param fit.list a list of fitted models.
#'@param data a data frame containing the variables used to fit the models listed in fit.list.
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
#'fit.list <- list()
#'
#'fit.list[[1]] <- glm(Survived ~ Sex, data = dt)
#'fit.list[[2]] <- glm(Survived ~ Age + Sex, data = dt)
#'fit.list[[3]] <- glm(Survived ~ Age + Sex + Pclass, data = dt)
#'
#'nt_multiple_logistic(fit.list)
#'
#'@importFrom purrr map
#'@importFrom utils write.csv
#'@export
nt_multiple_logistic <- function(fit.list, fit.labels = NULL, type = "or",
                                 format = TRUE, digits = 2, digits.p = 3,
                                 save = FALSE, file = "nt_multiple_logistic"){

  if (class(fit.list) != "list")
    fit.list <- list(fit.list)
  if (is.null(fit.labels))
    fit.labels <- 1:length(fit.list)

  if (type == "or"){
    temp <- map2(fit.list, fit.labels, aux_multiple_logistic,
                 format = format, type = "or")
    tab <- Reduce(rbind, temp)
  } else {
    temp <- map2(fit.list, fit.labels, aux_multiple_logistic,
                 format = format, type = "coef")
    tab <- Reduce(rbind, temp)
  }

  ref <- map(fit.list, ~ reference_df(.x)$ref)

  if (format)
    tab <-  tab %>%
    transmute(Model = .data$model,
              Variable = .data$variable, OR = .data$or,
              'Estimate (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                           round(.data$conf.low, digits), " ; ",
                                           round(.data$conf.high, digits), ")"),
              'p value' = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits.p)))) %>%
    replace_na(list('p value' = ""))

  if (save)
    write.csv(tab, file = paste0(file, ".csv"))

  out <- list(tab = tab, ref = ref)

  return(out)
}

#'@importFrom dplyr select mutate transmute bind_cols full_join
#'@importFrom purrr map map2
#'@importFrom broom tidy glance
#'@importFrom tidyr separate unite
#'@importFrom tibble data_frame
#'@importFrom gsubfn gsubfn
aux_multiple_logistic <- function(fit, model.label, format, type){

  aux <- extract_data(fit)

  if (type == "or"){
    temp <- table_fit(fit, exponentiate = TRUE)
    out <- temp %>%
      mutate(model = model.label,
             term = str_replace_all(.data$term, unlist(aux$var.labels))) %>%
      separate(.data$term, into = c("variable", "or"), sep = ":")

    if (format)
      out <- out %>% group_by(.data$variable) %>%
      mutate(aux_variable = ifelse(duplicated(.data$variable), "", .data$variable),
             p.value = ifelse(duplicated(.data$p.value), NA, .data$p.value)) %>%
      ungroup(.data$variable) %>% select(-.data$variable) %>%
      rename(variable = .data$aux_variable)

  } else {
    out <- tidy(fit)
  }
  return(out)
}


