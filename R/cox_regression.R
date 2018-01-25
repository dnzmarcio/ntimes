#'Simple Cox Regression
#'
#'@description Performing simple Cox regression.
#'
#'@param data a data frame with the variables.
#'@param time a numeric vector with the follow-up time.
#'@param status a numeric vector indicating status, 0 = censored, 1 = event at time.
#'@param strata a logical vector indicating if strata should be considered.
#'@param id a character vector containing the strata.
#'
#'@examples
#'library(survival)
#'library(dplyr)
#'library(magrittr)
#'data(ovarian)
#'
#'ovarian_nt <- ovarian %>% mutate(age = qt_var(age, label = "Age"),
#'                                 resid.ds = ql_var(resid.ds,
#'                                                   from = 1:2,
#'                                                   to = c("no", "yes")),
#'                                 ecog.ps = ql_var(ecog.ps,
#'                                                  from = 1:2,
#'                                                  to = c("I", "II")),
#'                                 rx = ql_var(rx,
#'                                             from = 1:2,
#'                                             to = c("t1", "t2")))
#'ovarian_nt %>% nt_simple_cox(time = futime, status = fustat)
#'
#'@export
nt_simple_cox <- function(data, time, status,
                          cluster = FALSE, strata = FALSE, id = NULL,
                          digits = 2, digits.p = 3, save = FALSE){

  data <- as_data_frame(data)
  time <- enquo(time)
  status <- enquo(status)

  if (ncol(data) > 2){
    vars <- select(.data = data, -!!time)
    vars <- select(.data = vars, -!!status)
    vars.name <- names(vars)
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
               time = time, status = status, strata - strata,
               digits = digits, digits.p = digits.p)

  cox <- Reduce(rbind, temp)

  if (save){
    utils::write.csv(cox, file = paste0(file, "_regression.csv"))
    utils::write.csv(survival, file = paste0(file, "_survival.csv"))
  }

  out <- list(survival = survival, cox = cox)

  return(out)
}

aux_simple_cox <- function(var, var.name, time, status, strata, digits, digits.p){

  var.label <- extract_label(var, var.name)

  if (strata){
    data.model <- bind_cols(time = time, status = status, id = id, var = var)
    temp <- fit_cox(data.model, var.label, strata, digits, digits.p)

    out <- temp %>%
      mutate(Variable = ifelse(.data$`HR (95% CI)` == "Reference" |
                                 .data$Group == "", .data$Variable, ""))

  } else {
    data.model <- bind_cols(time = time, status = status, var = var)
    temp <- fit_cox(data.model, var.label, strata, digits, digits.p)

    out <- temp %>%
      mutate(Variable = ifelse(duplicated(.data$Variable), "", .data$Variable),
             n = ifelse(duplicated(.data$n), "", .data$n))

  }

  return(out)
}

#'@importFrom survival coxph Surv cox.zph
#'@importFrom broom tidy glance
#'@importFrom tidyr separate replace_na
#'@importFrom dplyr select transmute mutate bind_cols
#'@importFrom tibble data_frame
fit_cox <- function(data, var.label, strata, digits, digits.p){

  if (!strata){
    mod <- coxph(Surv(time, status) ~ var, data = data)
    temp <- tidy(mod, exponentiate = TRUE) %>%
      separate(col = .data$term, into = c("var", "Group"), sep = "ar") %>%
      select(-.data$var, -.data$std.error, -.data$statistic) %>%
      transmute(Variable = var.label, .data$Group,
                'HR (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                       round(.data$conf.low, digits), " ; ",
                                       round(.data$conf.high, digits), ")"),
                'p value' = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                   as.character(round(.data$p.value, digits.p))))
  } else {
    mod <- coxph(Surv(time, status) ~ strata(id) + var, data = data)
    temp <- tidy(mod, exponentiate = TRUE) %>%
      separate(col = .data$term, into = c("var", "Group"), sep = "ar") %>%
      select(-.data$var, -.data$std.error, -.data$statistic) %>%
      transmute(Variable = var.label, .data$Group,
                'HR (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                       round(.data$conf.low, digits), " ; ",
                                       round(.data$conf.high, digits), ")"),
                'p value' = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                   as.character(round(.data$p.value, digits.p))))
  }

  zph.table <- cox.zph(mod)$table

  if (!is.numeric(data$var)){
    first_row <- data_frame(Variable = var.label,
                            Group =
                              na.exclude(as.character(unique(data$var)[
                                !(unique(data$var) %in% temp$Group)])),
                            'HR (95% CI)' = "Reference",
                            'p value' = "")

    aux <- glance(mod) %>% select(.data$n, n.event = .data$nevent, .data$concordance,
                                  .data$r.squared, .data$AIC) %>%
      mutate(concordance = round(.data$concordance, 3), r.squared = round(.data$r.squared, 3),
             AIC = round(.data$AIC, 3),
             ph.assumption = round(zph.table[nrow(zph.table), 3], 3))
    first_row <- bind_cols(first_row, aux)

    out <- full_join(first_row, temp, by = c("Variable", "Group", 'HR (95% CI)', 'p value')) %>%
      replace_na(list(n = "", n.event = "", concordance = "", r.squared = "",
                      AIC = "", ph.assumption = ""))

  } else {
    aux <- glance(mod) %>% select(.data$n, n.event = .data$nevent,
                                  .data$concordance, .data$r.squared, .data$AIC) %>%
      mutate(concordance = round(.data$concordance, 3), r.squared = round(.data$r.squared, 3),
             AIC = round(.data$AIC, 3),
             ph.assumption = round(zph.table[nrow(zph.table), 3], 3))

    out <- bind_cols(temp, aux) %>%
      replace_na(list(n = "", n.event = "", concordance = "", r.squared = "",
                      AIC = "", ph.assumption = ""))
  }

  return(out)
}



#'Proportional Hazards Cox regression table
#'
#'@description Tabulating results from fitted Proportional Hazards Cox models.
#'
#'@param fit.list a list of fitted models.
#'@param data a data frame containing the variables used to fit the models listed in fit.list.
#'
#'@examples library(survival)
#'library(magrittr)
#'library(dplyr)
#'
#'data(ovarian)
#'ovarian <- ovarian %>% mutate(resid.ds = ql_var(resid.ds,
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
#'fit.list[[1]] <- coxph(Surv(futime, fustat) ~ ecog.ps, data = ovarian)
#'fit.list[[2]] <- coxph(Surv(futime, fustat) ~ resid.ds + rx, data = ovarian)
#'fit.list[[3]] <- coxph(Surv(futime, fustat) ~ ecog.ps*rx, data = ovarian)
#'
#'nt_table_coxph(fit.list, data = ovarian)
#'
#'@importFrom purrr map
#'@importFrom utils write.csv
#'@export
nt_table_coxph <- function(fit.list, data, save = FALSE, file = "nt_table_cox"){

  temp <- map(fit.list, fit_multiple_cox, data = data)
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
fit_multiple_cox <- function(fit, data){

  vars.name <- names(fit$xlevels)

  group <- data %>% select(vars.name)
  vars.label <- map2(group, vars.name, extract_label)
  nlv <- map(group, function(x) out <- nlevels(x) - 1)

  temp <- tidy(fit, exponentiate = TRUE) %>%
    separate(col = .data$term, into = c("Variable", "Group"),
             sep = rep(vars.name, nlv)) %>%
    mutate(Variable = rep(vars.label, nlv)) %>%
    unite(col = "Group", .data$Variable, .data$Group, sep = ":") %>%
    select(-.data$std.error, -.data$statistic) %>%
    transmute(.data$Group,
              'HR (95% CI)' = paste0(round(.data$estimate, 2), " (",
                                     round(.data$conf.low, 2), " ; ",
                                     round(.data$conf.high, 2), ")"),
              'p value' = ifelse(round(.data$p.value, 3) == 0, "< 0.001",
                                 as.character(round(.data$p.value, 3))))

  if (nrow(temp) > 1){
    ph.assumption <- round(cox.zph(fit)$table[-nrow(temp), 3], 3)
  } else {
    ph.assumption <- round(cox.zph(fit)$table[, 3], 3)
  }

  temp <- bind_cols(temp, ph.assumption = ph.assumption)

  aux <- glance(fit) %>% select(.data$n, n.event = .data$nevent,
                                .data$concordance, .data$r.squared, .data$AIC) %>%
    mutate(concordance = round(.data$concordance, 3), r.squared = round(.data$r.squared, 3),
           AIC = round(.data$AIC, 3))

  first_row <- data_frame(Group = "Reference", 'HR (95% CI)' = "1", 'p value' = "")
  first_row <- bind_cols(first_row, aux)
  last_row <- data_frame(Group = as.character(""), 'HR (95% CI)' = as.character(""),
                         'p value' = as.character(""))

  aux <- full_join(first_row, temp, by = c("Group", "HR (95% CI)", "p value"))
  out <- full_join(aux, last_row, by = c("Group", "HR (95% CI)", "p value")) %>%
    replace_na(list(n = "", n.event = "", concordance = "", r.squared = "",
                    AIC = "", ph.assumption = ""))
  return(out)
}

