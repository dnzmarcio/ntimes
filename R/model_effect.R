nt_model_effect <- function(fit) {
  UseMethod("model_effect")
}

effect <- function(fit, type, exponentiate) {
  UseMethod("effect")
}

#'@importFrom dplyr mutate group_by ungroup rename
#'@importFrom stringr str_replace_all
#'@importFrom tidyr separate
#'@importFrom broom tidy
nt_model_effect.cox <- function(fit, ci.type = "lr",
                            format = TRUE, digits = 2, digits.p = 3,
                            save = FALSE, file = "nt_model_effect"){

    out <- aux_multiple_cox(fit, ci.type, format)
    ref <- reference_df(fit)$ref

    if (format)
      out$effect <-  out$effect |>
      transmute(Variable = .data$variable, HR = .data$hr,
                'Estimate (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                             round(.data$conf.low, digits), " ; ",
                                             round(.data$conf.high, digits), ")"),
                'p value LR' = ifelse(round(.data$p.value.lr, digits.p) == 0, "< 0.001",
                                      as.character(round(.data$p.value.lr, digits.p)))) |>
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
nt_model_effect.glm <- function(fit, ci.type = "lr",
                                format = TRUE, digits = 2, digits.p = 3,
                                save = FALSE, file = "nt_model_effect"){

  if (fit$family[[1]] == "binomial"){

    out <- aux_multiple_logistic(fit, format, ci.type)
    ref <- reference_df(fit)$ref

    if (format)
      out$effect <- out$effect |>
      transmute(Variable = .data$variable, OR = .data$or,
                'Estimate (95% CI)' = paste0(round(.data$estimate, digits), " (",
                                             round(.data$conf.low, digits), " ; ",
                                             round(.data$conf.high, digits), ")"),
                'p value LR' = ifelse(round(.data$p.value.lr, digits.p) == 0, "< 0.001",
                                      as.character(round(.data$p.value.lr, digits.p)))) |>
      replace_na(list('p value LR' = ""))
  }

  if (save)
    write.csv(out$effect, file = paste0(file, ".csv"))

  out <- list(effect = out$effect, coef = out$coef, ref = ref)

  return(out)
}
