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
#'titanic_train <- titanic_train %>% mutate(Sex = ql_var(Sex,
#'                                                       from = c("male", "female"),
#'                                                       to = c("Male", "Female")),
#'                                          Pclass = ql_var(Pclass,
#'                                                          from = 1:3,
#'                                                          to = c("I", "II", "III"),
#'                                                          label = "Passenger Class"),
#'                                          Embarked = ql_var(Embarked,
#'                                                            from = c("C", "Q", "S"),
#'                                                            to = c("Cherbourg", "Queenstown", "Southampton"))                )
#'
#'dt <- titanic_train %>% select(Survived, Sex, Age, Pclass, Embarked)
#'
#'dt %>% nt_simple_logistic(response = Survived, Age)
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
  if (!is.null(add)){
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
    out <- out %>%
      transmute(Variable = term,
                OR.95CI = paste0(round(.data$estimate, digits), " (",
                                 round(.data$conf.low, digits), " ; ",
                                 round(.data$conf.high, digits), ")"),
                p.value = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits.p))),
                null.deviance = round(.data$null.deviance, digits),
                logLik = round(.data$logLik, digits),
                AIC = round(.data$AIC, digits),
                BIC = round(.data$BIC, digits),
                deviance = round(deviance, digits)) %>%
      rename(`OR (95% CI)` = OR.95CI, `p value` = p.value) %>%
      mutate(`OR (95% CI)` = replace(`OR (95% CI)`, 1, "(Reference)"),
             `p value` = replace(`p value`, 1, ""))
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
  var.class <- unlist(map(cbind(add, var), is.numeric))
  fit.labels <- ifelse(var.class,
                       c(add.label, var.label),
                       paste0(c(add.label, var.label), ": "))

  data.model <- bind_cols(response = response, add = add, var = var)
  out <- fit_logistic(data.model, fit.labels)

  return(out)
}

#'@importFrom broom tidy glance
#'@importFrom tidyr replace_na
#'@importFrom dplyr select transmute mutate
#'@importFrom stringr str_replace_all
fit_logistic <- function(data, fit.labels){

  data.model <- na.exclude(data)

  fit <- glm(response ~ ., data = data.model, family = "binomial")

  null.fit <- glm(response ~ 1, data = data.model, family = "binomial")
  p.value.F <- anova(null.fit, fit)$"p.value"

  temp <- tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(term = str_replace_all(.data$term, unlist(fit.labels))) %>%
    select(-.data$std.error, -.data$statistic) %>%
    mutate(estimate = replace(.data$estimate, 1, 1),
           conf.low = replace(.data$conf.low, 1, NA),
           conf.high = replace(.data$conf.high, 1, NA),
           p.value = replace(.data$p.value, 1, NA))

  aux <- bind_cols(n = nrow(na.exclude(data.model)), glance(fit))

  out <- merge(data.frame(temp, row.names=NULL), data.frame(aux, row.names=NULL),
               by = 0, all = TRUE)[-1] %>%
    replace_na(list(null.deviance = "", df.null = "",
                    logLik = "", AIC = "", BIC = "",
                    deviance = "", df.residual = ""))

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
#'nt_table_coxph(fit.list, data = dt)
#'
#'@importFrom purrr map
#'@importFrom utils write.csv
#'@export
nt_multiple_logistic <- function(fit.list, data,
                                 format = TRUE, digits = 2, digits.p = 3,
                                 save = FALSE, file = "nt_multiple_logistic"){

  temp <- map(fit.list, aux_multiple_logistic,
              format, digits, digits.p, data = data)

  out <- Reduce(rbind, temp)

  if (save)
    write.csv(out, file = paste0(file, ".csv"))

  return(out)
}

#'@importFrom dplyr select mutate transmute bind_cols full_join
#'@importFrom purrr map map2
#'@importFrom broom tidy glance
#'@importFrom tidyr separate unite
#'@importFrom tibble data_frame
#'@importFrom gsubfn gsubfn
aux_multiple_logistic <- function(fit, format, digits, digits.p, data){

  var.name <- names(data)
  var.label <- extract_label(var, var.name)
  var.class <- unlist(map(cbind(add, var), is.numeric))
  fit.labels <- ifelse(var.class,
                       c(add.label, var.label),
                       paste0(c(add.label, var.label), ": "))

  temp <- tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(term = str_replace_all(term, unlist(fit.labels))) %>%
    select(-.data$std.error, -.data$statistic) %>%
    mutate(estimate = replace(.data$estimate, 1, 1),
           conf.low = replace(.data$conf.low, 1, NA),
           conf.high = replace(.data$conf.high, 1, NA),
           p.value = replace(.data$p.value, 1, NA))

  aux <- glance(fit)

  if (format){
    temp <- temp %>%
      transmute(Variable = .data$term,
                OR.95CI = paste0(round(.data$estimate, digits), " (",
                                 round(.data$conf.low, digits), " ; ",
                                 round(.data$conf.high, digits), ")"),
                p.value = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value, digits.p))))

    aux <- aux %>%
      mutate(null.deviance = round(.data$null.deviance, digits),
             logLik = round(.data$logLik, digits),
             AIC = round(.data$AIC, digits),
             BIC = round(.data$BIC, digits),
             deviance = round(.data$deviance, digits))
  }

  out <- merge(data.frame(temp, row.names=NULL), data.frame(aux, row.names=NULL),
               by = 0, all = TRUE)[-1] %>%
    replace_na(list(null.deviance = "", df.null = "",
                    logLik = "", AIC = "", BIC = "",
                    deviance = "", df.residual = ""))

  return(out)
}


