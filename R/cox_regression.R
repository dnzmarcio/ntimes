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
#'@param ci.type a character value indicating the procedure to calculate confidence intervals: likelihood ratio (\code{lr}) or wald (\code{wald}).
#'@param tab.type a character value indicating either hazard ratio (\code{hr}) or coefficients (\code{coef}) will be shown in the output.
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
#'
#'fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps*rx, data = ovarian_nt)
#'
#'nt_multiple_cox(fit)
#'
#'@importFrom purrr map2 map
#'@importFrom utils write.csv
#'@importFrom dplyr transmute bind_rows
#'@importFrom tidyr replace_na
#'@export
nt_multiple_cox <- function(fit, ci.type = "lr",
                            format = TRUE, digits = 2, digits.p = 3,
                            save = FALSE, file = "nt_multiple_cox"){

  if (class(fit) != "coxph")
    stop("fit object is not a coxph class")

  out <- aux_multiple_cox(fit, ci.type, format)
  ref <- reference_df(fit)$ref

  if (format)
    out$hr <-  out$hr %>%
    transmute(Variable = .data$variable, HR = .data$hr,
              'Estimate (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                           round(.data$conf.low, digits), " ; ",
                                           round(.data$conf.high, digits), ")"),
              'p value LR' = ifelse(round(.data$p.value.lr, digits.p) == 0, "< 0.001",
                                 as.character(round(.data$p.value.lr, digits.p)))) %>%
    replace_na(list('p value LR' = ""))



  if (save)
    write.csv(out$hr, file = paste0(file, ".csv"))

  out <- list(hr = out$hr, coef = out$coef, ref = ref)

  return(out)
}

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
aux_multiple_cox <- function(fit, ci.type, format){

  aux <- extract_data(fit)

  hr <- effect.coxph(fit, ci.type, exponentiate = TRUE) %>%
    mutate(term = str_replace_all(.data$term, unlist(aux$var.labels))) %>%
    separate(.data$term, into = c("variable", "hr"), sep = ":")

  if (format)
    hr <- hr %>% group_by(.data$variable) %>%
    mutate(aux_variable = ifelse(duplicated(.data$variable), "", .data$variable),
           p.value.lr = ifelse(duplicated(.data$p.value.lr), NA, .data$p.value.lr)) %>%
    ungroup(.data$variable) %>% select(-.data$variable) %>%
    rename(variable = .data$aux_variable)

  temp <- unlist(aux$var.labels)
  labels <- paste0(temp, " ")
  names(labels) <- names(temp)
  coef <- tidy(fit) %>%
    mutate(term = str_replace_all(.data$term, labels),
           term = sub(" $", "", x = term))

  out <- list(hr = hr, coef = coef)

  return(out)
}


#'@importFrom stats model.matrix formula setNames anova vcov update.formula
#'@importFrom survival coxph
#'@importFrom stringr str_split
effect.coxph <- function(fit, type, exponentiate){

  aux <- extract_data(fit)
  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta.var <- as.matrix(vcov(fit))
  term.labels <- attr(fit$terms, "term.labels")

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

  for (i in 1:length(aux$var)){

    if (length(interaction) > 0){
      cond.interaction <- grepl(aux$var[i], x = interaction, fixed = TRUE)
    } else {
      cond.interaction <- FALSE
    }

    if (all(!cond.interaction)){
      temp <- contrast_df(aux$data, aux$var[i], ref)
      design.matrix <- model.matrix(fit, temp$new.data)

      drop <- which(grepl(aux$var[i], x = as.character(term.labels), fixed = TRUE))
      fit0 <- coxph(update.formula(fit$formula, paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
                    data = aux$data)
      p.value.lr <- anova(fit0, fit)$`P(>|Chi|)`[2]


      contrast <- contrast_calc(fit = fit, design.matrix = design.matrix,
                                beta = beta, beta.var = beta.var,
                                p.value = p.value.lr, type = type)

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
        fit0 <- coxph(update.formula(fit$formula,
                                     paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
                      data = aux$data)
        p.value <- anova(fit0, fit)$`P(>|Chi|)`[2]

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

  colnames(out) <- c("term", "estimate", "conf.low", "conf.high",
                     "p.value.lr")

  return(out)
}


