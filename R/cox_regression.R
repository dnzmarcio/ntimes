#'Simple Cox Regression
#'
#'@description Performing simple Cox regression.
#'
#'@param data a data frame with the variables.
#'@param time a numeric vector with the follow-up time.
#'@param status a numeric vector indicating status, 0 = censored, 1 = event at time.
#'@param strata a character vector containing the strata.
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
#'@importFrom rlang quo_get_expr
#'
#'@export
nt_simple_cox <- function(data, time, status, ...,
                          cluster = FALSE, strata = NULL,
                          digits = 2, digits.p = 3, save = FALSE,
                          file = "simple_cox", format = TRUE){

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
               strata.var = strata.var, digits = digits,
               digits.p = digits.p)

  cox <- Reduce(rbind, temp)

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
              n, n.event, concordance, r.squared, AIC, ph.assumption) %>%
    rename(`HR (95% CI)` = HR.95CI, `p value` = p.value) %>%
    mutate(Variable = ifelse(duplicated(Variable), "", Variable))
  }

  if (save){
    utils::write.csv(cox, file = paste0(file, "_regression.csv"))
    utils::write.csv(survival, file = paste0(file, "_survival.csv"))
  }

  out <- list(survival = survival, cox = cox)

  return(out)
}

aux_simple_cox <- function(var, var.name, time, status,
                           add, add.name, add.label,
                           strata.var, digits, digits.p){

  var.label <- extract_label(var, var.name)
  aux <- cbind(add, var)
  if(ncol(aux) > 1) {
    var.class <- unlist(map(cbind(add, var), is.numeric))
  } else {
    var.class <- list(var = is.numeric(var))
  }

  for (i in 1:length(var.class)){
    tab.labels <- list()
    tab.levels <- list()

    if (var.class[[i]]){
      tab.labels <- c(add.label, var.label)
      tab.levels <- ""
    } else {
      tab.labels <- paste0(c(add.label, var.label), ": ")
      lv <- levels(var)
      tab.levels <- paste0(lv[2:length(lv)], "/", lv[1])
    }
  }

  if (!is.list(tab.labels))
    tab.labels <- setNames(as.list(tab.labels), "var")

  data.model <- bind_cols(time = time, status = status, var = var, add = add)
  out <- fit_cox(data.model, tab.labels, tab.levels, strata.var)

  if (is.factor(var))
    if (length(levels(var)) > 2)
      out <- out %>%
    mutate(p.value = p.value.lh)

  return(out)
}

#'@importFrom survival coxph Surv cox.zph
#'@importFrom broom tidy glance
#'@importFrom tidyr separate replace_na
#'@importFrom dplyr select transmute mutate bind_cols
#'@importFrom tibble data_frame
fit_cox <- function(data, tab.labels, tab.levels, strata.var){

  if (any(apply(data, 1, is.na)))
    strata.var <- strata.var[-which(apply(data, 1, is.na))]
  data <- na.exclude(data)

  if (is.null(strata.var)){
    fit <- coxph(Surv(time, status) ~ ., data = data)
    temp <- tidy(fit, exponentiate = TRUE) %>%
      mutate(term = str_replace_all(.data$term, unlist(tab.labels))) %>%
      select(-.data$std.error, -.data$statistic) %>%
      separate(term, into = c("term", "group"), sep = ":", fill = "right") %>%
      mutate(group = tab.levels)
  } else {
    fit <- coxph(Surv(time, status) ~ strata(strata.var) + ., data = data)
    temp <- tidy(fit, exponentiate = TRUE) %>%
      mutate(term = str_replace_all(.data$term, unlist(tab.labels))) %>%
      select(-.data$std.error, -.data$statistic) %>%
      separate(term, into = c("term", "group"), sep = ":", fill = "right") %>%
      mutate(group = tab.levels)
  }

  fit0 <- coxph(update.formula(fit$formula, paste0(" ~ . - var")),
                data = data)

  p.value.lh <- anova(fit0, fit)$`P(>|Chi|)`[2]

  zph.table <- cox.zph(fit)$table

  aux <- glance(fit) %>% select(.data$n, n.event = .data$nevent,
                                .data$concordance, .data$r.squared, .data$AIC) %>%
    mutate(p.value.lh = p.value.lh,
           ph.assumption = zph.table[nrow(zph.table), 3])

  out <- merge(data.frame(temp, row.names=NULL), data.frame(aux, row.names=NULL),
               by = 0, all = TRUE)[-1]

  return(out)
}



#'Proportional Hazards Cox regression table
#'
#'@description Tabulating results from fitted Proportional Hazards Cox models.
#'
#'@param fit.list a list of fitted models.
#'@param data a data frame containing the variables used to fit the models listed in fit.list.
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
#'@importFrom purrr map map2
#'@importFrom utils write.csv
#'@export
nt_multiple_cox <- function(fit.list, fit.labels = NULL, format = FALSE, digits = 2, digits.p = 3,
                            save = FALSE, file = "nt_multiple_cox"){

  if (class(fit.list) != "list")
    fit.list <- list(fit.list)
  if (is.null(fit.labels))
    model.labels <- 1:length(fit.list)

  temp <- map2(fit.list, model.labels, aux_multiple_cox,
               format = format)
  tab <- Reduce(rbind, temp)
  adj <- map(fit.list, ~ reference_df(.x)$ref)

  if (format)
    tab <-  tab %>%  transmute(Model = .data$model,
                               Variable = .data$variable, Group = .data$group,
                               'HR (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                                      round(.data$conf.low, digits), " ; ",
                                                      round(.data$conf.high, digits), ")"),
                               'p value' = ifelse(round(.data$p.value, digits.p) == 0, "< 0.001",
                                                as.character(round(.data$p.value, digits.p))))



  if (save)
    write.csv(tab, file = paste0(file, ".csv"))

  out <- list(tab = tab, adj = adj)

  return(out)
}

#'@importFrom stringr str_replace_all
#'@importFrom dplyr select
#'@importFrom tidyr separate
aux_multiple_cox <- function(fit, model.label, format){

  aux <- extract_data(fit)
  temp <- table_fit(fit, exponentiate = TRUE)
  out <- temp %>%
    mutate(model = model.label,
      term = str_replace_all(.data$term, unlist(aux$var.labels))) %>%
    separate(.data$term, into = c("variable", "group"), sep = ":")

  if (format)
    out <- out %>%
    mutate(variable = ifelse(duplicated(.data$variable), "", .data$variable))

  return(out)
}


