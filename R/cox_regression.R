#'Simple Cox Regression
#'
#'@description Performing simple Cox regression.
#'
#'@param data a data frame with the variables.
#'@param time a numeric vector with the follow-up time.
#'@param status a numeric vector indicating status, 0 = censored, 1 = event at time.
#'@param ...  character values indicating confounding variables.
#'@param cluster a character vector containing the cluster variable.
#'@param strata a character vector containing the strata variable.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits.p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@examples
#'library(survival)
#'library(dplyr)
#'library(magrittr)
#'data(ovarian)
#'
#'ovarian_nt <- ovarian %>% mutate(resid.ds = ql_var(resid.ds,
#'                                                from = 1:2,
#'                                                to = c("no", "yes"),
#'                                                label = "Residual Disease"),
#'                              ecog.ps = ql_var(ecog.ps,
#'                                               from = 1:2,
#'                                               to = c("I", "II"),
#'                                               label = "ECOG-PS"),
#'                              rx = ql_var(rx,
#'                                          from = 1:2,
#'                                          to = c("t1", "t2"),
#'                                          label = "Treatment"),
#'                              age = qt_var(age,
#'                                           label = "Age"))
#'ovarian_nt %>% nt_simple_cox(time = futime, status = fustat)
#'
#'@importFrom rlang enquo quos quo_get_expr .data
#'@importFrom dplyr select mutate transmute rename
#'@importFrom purrr map2
#'@importFrom survival survfit
#'@importFrom broom tidy
#'@importFrom tidyr replace_na
#'@importFrom utils write.csv
#'
#'@export
nt_simple_cox <- function(data, time, status, ...,
                          cluster = FALSE, strata = NULL, format = TRUE,
                          digits = 2, digits.p = 3, save = FALSE,
                          file = "simple_cox"){

  data <- as_data_frame(data)
  time <- enquo(time)
  status <- enquo(status)
  strata <- enquo(strata)
  aux <- quos(...)

  if (ncol(data) > 2){
    vars <- select(.data = data, -!!time)
    vars <- select(.data = vars, -!!status)
    if (!is.null(quo_get_expr(strata))){
      vars <- select(.data = vars, -!!strata)
      strata.var <- select(.data = data, !!strata)
      strata.var <- strata.var[[1]]
    } else {
      strata.var <- NULL
    }
    vars.name <- names(vars)

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
  }

  time <- select(.data = data, !!time)
  time <- time[[1]]
  status <- select(.data = data, !!status)
  status <- as.numeric(as.factor(status[[1]])) - 1

  fit <- survfit(Surv(time, status) ~ 1)
  survival <- tidy(fit) %>% select(-.data$std.error) %>%
    mutate(estimate = round(100*.data$estimate, digits),
           conf.low = round(100*.data$conf.low, digits),
           conf.high = round(100*.data$conf.high, digits)) %>%
    transmute(Time = .data$time, .data$n.risk, .data$n.event, .data$n.censor,
              `Survival (95% CI)` = paste0(.data$estimate, " (",
                                           .data$conf.low, " ; ",
                                           .data$conf.high, ")"))

  temp <- map2(.x = vars, .y = vars.name, .f = aux_simple_cox,
               time = time, status = status,
               add = add, add.name = add.name, add.label = add.label,
               strata.var = strata.var, format = format)

  cox <- bind_rows(temp)

  if (format) {

  cox <- cox %>% mutate(concordance = round(.data$concordance, digits),
                        r.squared = round(.data$r.squared, digits),
                        AIC = round(.data$AIC, digits),
                        ph.assumption = round(.data$ph.assumption, digits.p)) %>%
    transmute(Variable = .data$term,  Group = .data$group,
              HR.95CI = paste0(round(.data$estimate, digits), " (",
                               round(.data$conf.low, digits), " ; ",
                               round(.data$conf.high, digits), ")"),
              p.value = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                               as.character(round(.data$p.value, digits.p))),
              n = .data$n, n.event = .data$n.event,
              concordance = .data$concordance, r.squared = .data$r.squared,
              AIC = .data$AIC, ph.assumption  = .data$ph.assumption) %>%
    replace_na(list(p.value = "")) %>%
    rename(`HR (95% CI)` = .data$HR.95CI, `p value` = .data$p.value)
  }

  if (save){
    write.csv(cox, file = paste0(file, "_regression.csv"))
    write.csv(survival, file = paste0(file, "_survival.csv"))
  }

  out <- list(survival = survival, cox = cox)

  return(out)
}

#'@importFrom purrr map
#'@importFrom magrittr %>%
#'@importFrom stats setNames
#'@importFrom dplyr mutate bind_cols
#'@importFrom rlang .data
aux_simple_cox <- function(var, var.name, time, status,
                           add, add.name, add.label,
                           strata.var, format){

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

  data.model <- bind_cols(time = time, status = status, var = var, add = add)
  out <- fit_cox(data.model, tab.labels, tab.levels, strata.var)

  if (format){
    out <- out %>% mutate(term = ifelse(duplicated(.data$term), "", .data$term))
    if (is.factor(var))
      if (length(levels(var)) > 2)
        out <- out %>% mutate(p.value = .data$p.value.lh)
  }

  return(out)
}

#'@importFrom survival coxph Surv cox.zph
#'@importFrom broom tidy glance
#'@importFrom tidyr separate replace_na
#'@importFrom dplyr select mutate
#'@importFrom stats na.exclude update.formula anova
#'@importFrom stringr str_replace_all
fit_cox <- function(data, tab.labels, tab.levels, strata.var){

  if (any(is.na(data)))
    strata.var <- strata.var[-which(is.na(data), arr.ind = TRUE)[, 1]]
  data <- na.exclude(data)

  if (is.null(strata.var)){
    fit <- coxph(Surv(time, status) ~ ., data = data)
    temp <- tidy(fit, exponentiate = TRUE) %>%
      mutate(term = str_replace_all(.data$term, unlist(tab.labels))) %>%
      select(-.data$std.error, -.data$statistic) %>%
      separate(.data$term, into = c("term", "group"), sep = ":", fill = "right") %>%
      mutate(group = unlist(tab.levels))
  } else {
    fit <- coxph(Surv(time, status) ~ strata(strata.var) + ., data = data)
    temp <- tidy(fit, exponentiate = TRUE) %>%
      mutate(term = str_replace_all(.data$term, unlist(tab.labels))) %>%
      select(-.data$std.error, -.data$statistic) %>%
      separate(.data$term, into = c("term", "group"), sep = ":", fill = "right") %>%
      mutate(group = unlist(tab.levels))
  }

  fit0 <- coxph(update.formula(fit$formula, paste0(" ~ . - var")),
                data = data)

  p.value.lh <- anova(fit0, fit, test = "Chisq")$`P(>|Chi|)`[2]

  zph.table <- cox.zph(fit)$table

  aux <- glance(fit) %>% select(.data$n, n.event = .data$nevent,
                                .data$concordance, .data$r.squared, .data$AIC) %>%
    mutate(p.value.lh = p.value.lh, ph.assumption = zph.table[nrow(zph.table), 3])

  out <- merge(data.frame(temp, row.names=NULL), data.frame(aux, row.names=NULL),
               by = 0, all = TRUE)[-1]

  return(out)
}



#'Proportional Hazards Cox regression table
#'
#'@description Tabulating results from fitted Proportional Hazards Cox models.
#'
#'@param fit.list a list of fitted models.
#'@param fit.labels a character vector labelling the models in \code{fit.list}
#'@param type a character value indicating either hazard ratio (\code{hr}) or coefficients (\code{coef}) will be shown in the output.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param digits.p a numerical value defining number of digits to present the p-values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character indicating the name of output file in csv format to be saved.
#'
#'@examples library(survival)
#'library(magrittr)
#'library(dplyr)
#'
#'data(ovarian)
#'ovarian_nt <- ovarian %>% mutate(resid.ds = ql_var(resid.ds,
#'                                                from = 1:2,
#'                                                to = c("no", "yes"),
#'                                                label = "Residual Disease"),
#'                              ecog.ps = ql_var(ecog.ps,
#'                                               from = 1:2,
#'                                               to = c("I", "II"),
#'                                               label = "ECOG-PS"),
#'                              rx = ql_var(rx,
#'                                          from = 1:2,
#'                                          to = c("t1", "t2"),
#'                                          label = "Treatment"),
#'                              age = qt_var(age,
#'                                           label = "Age"))
#'
#'fit.list <- list()
#'
#'fit.list[[1]] <- coxph(Surv(futime, fustat) ~ ecog.ps, data = ovarian_nt)
#'fit.list[[2]] <- coxph(Surv(futime, fustat) ~ resid.ds + rx, data = ovarian_nt)
#'fit.list[[3]] <- coxph(Surv(futime, fustat) ~ age + ecog.ps*rx, data = ovarian_nt)
#'
#'nt_multiple_cox(fit.list)
#'
#'@importFrom purrr map2 map
#'@importFrom utils write.csv
#'@importFrom dplyr transmute bind_rows
#'@importFrom tidyr replace_na
#'@export
nt_multiple_cox <- function(fit.list, fit.labels = NULL, type = "hr",
                            format = TRUE, digits = 2, digits.p = 3,
                            save = FALSE, file = "nt_multiple_cox"){

  if (class(fit.list) != "list")
    fit.list <- list(fit.list)
  if (is.null(fit.labels))
    fit.labels <- 1:length(fit.list)

  if (type == "hr"){
    temp <- map2(fit.list, fit.labels, aux_multiple_cox,
                 format = format, type = "hr")
    tab <- bind_rows(temp)
  } else {
    temp <- map2(fit.list, fit.labels, aux_multiple_cox,
                 format = format, type = "coef")
    tab <- bind_rows(temp)
  }

  ref <- map(fit.list, ~ reference_df(.x)$ref)

  if (format)
    tab <-  tab %>%
    transmute(Model = .data$model,
              Variable = .data$variable, HR = .data$hr,
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

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
aux_multiple_cox <- function(fit, model.label, format, type){

  aux <- extract_data(fit)

  if (type == "hr"){
    temp <- table_fit(fit, exponentiate = TRUE)
    out <- temp %>%
      mutate(model = model.label,
             term = str_replace_all(.data$term, unlist(aux$var.labels))) %>%
      separate(.data$term, into = c("variable", "hr"), sep = ":")

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


