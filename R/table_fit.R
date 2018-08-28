#'@importFrom stats terms get_all_vars
extract_data <- function(fit){

  data <- eval(fit$call$data)
  temp <- get_all_vars(formula(fit), data)
  var.names <- colnames(temp)
  var.labels <- map2(temp, var.names, ~ extract_label(.x, .y))

  suppressWarnings(temp <- temp[!apply(is.na(temp), 1, any), , drop = FALSE])
  var <- sapply(names(temp), grepl,
                x = as.character(formula(fit)[3]),
                fixed = TRUE)
  var <- names(var[var])

  if (!is.null(attr(terms(fit),"specials")$strata))
    var <- var[-(attr(terms(fit),"specials")$strata - 1)]

  out <- list(data = droplevels(temp), var = var, var.labels = var.labels)

  return(out)
}

#'@importFrom stats setNames median
reference_df <- function(fit){

  aux <- extract_data(fit)
  data <- aux$data

  df <- setNames(data.frame(matrix(ncol = ncol(data), nrow = 1)), names(data))

  for (i in 1:length(data)){
    if (is.numeric(data[, i])){
      df[, i] <- median(data[, i], na.rm = TRUE)
    } else if (is.factor(data[, i])){
      df[, i] <- levels(data[, i])[1]
    }
  }

  ref <- df[, which(sapply(aux$var, grepl,
                           x = names(df),
                           fixed = TRUE), arr.ind = TRUE)[, 1]]

  out <- list(df = df, ref = ref)

  return(out)
}

#'@importFrom stats quantile
contrast_df <- function(data, var, ref, interaction = NULL){

  contrast <- ref

  if (class(data[[var]]) == "factor" | class(data[[var]]) == "character"){
    data[[var]] <- as.factor(data[[var]])
    lv <- levels(data[[var]])
    contrast <- contrast[rep(1, each = length(lv)), ]
    contrast[[var]] <- lv
    label <- paste0(var, ": ", contrast[[var]][2:length(lv)], "/", contrast[[var]][1])

  } else if (class(data[[var]]) == "numeric" |
             class(data[[var]]) == "integer") {
    quantiles <- quantile(data[[var]], probs = c(0.25, 0.75))
    contrast <- contrast[c(1, 1), ]
    contrast[[var]] <- round(quantiles, 2)
    label <- paste0(var, ": ", contrast[[var]][2], "/", contrast[[var]][1])
  }

  if (!is.null(interaction)){
    for(k in 1:length(interaction)){
      if (class(data[[interaction[k]]]) == "factor" |
          class(data[[interaction[k]]]) == "character"){
        data[[interaction[k]]] <- as.factor(data[[interaction[k]]])
        lv <- levels(data[[interaction[k]]])
        contrast <- contrast[rep(1:nrow(contrast), each = length(lv)), ]
        contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/length(lv))
        contrast[[interaction[k]]] <-
          factor(contrast[[interaction[k]]], levels = lv)

        if (k == 1){
          label <- as.character(sapply(label, function(x) paste0(x, " at ", interaction[k], " = ", lv)))
        } else {
          label <- as.character(sapply(label, function(x) paste0(x, ", ", interaction[k], " = ", lv)))
        }

      } else if (class(data[[interaction[k]]]) == "numeric" |
                 class(data[[interaction[k]]]) == "integer") {
        quantiles <- quantile(data[[interaction[k]]], probs = c(0.25, 0.75))
        contrast <- contrast[rep(1:nrow(contrast), each = 2), ]
        contrast[[interaction[k]]] <- rep(quantiles, nrow(contrast)/2)

        if (k == 1){
          label <- as.character(sapply(label, function(x) paste0(x, " at ", interaction[k], " = ", round(quantiles, 2))))
        } else {
          label <- as.character(sapply(label, function(x) paste0(x, ", ", interaction[k], " = ", round(quantiles, 2))))
        }

      }
    }

    new.data <- split(contrast,  as.list(contrast)[interaction])
  } else {
    new.data <- contrast
  }





  rownames(new.data) <- NULL
  out <- list(new.data = new.data, label = label)
  return(out)
}


contrast_calc <- function(design.matrix, beta, beta.var,  p.value){

  out <- list()

  pred <- design.matrix%*%beta

  for (j in 2:nrow(design.matrix)){
    estimate <- as.numeric(pred[j] - pred[1])
    diff <- design.matrix[j, ] - design.matrix[1, ]
    pred.se <- sqrt(t(diff)%*%beta.var%*%diff)

    lower <- exp(estimate - 1.96*pred.se)
    upper <- exp(estimate + 1.96*pred.se)
    estimate <- exp(estimate)

    if (j == 2){
      out[[j-1]] <- data.frame(estimate, lower, upper, p.value)
    } else {
      out[[j-1]] <- data.frame(estimate, lower, upper, p.value = NA)
    }
  }

  return(out)
}

#'@importFrom stats model.matrix formula setNames anova vcov update.formula
#'@importFrom survival coxph
effect.coxph <- function(fit){

  aux <- extract_data(fit)
  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta.var <- as.matrix(vcov(fit))
  term.labels <- attr(fit$terms, "term.labels")

  interaction <- colnames(attr(fit$terms, "factors"))[attr(fit$terms, "order") > 1]

  L <- list()

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
      p.value <- anova(fit0, fit)$`P(>|Chi|)`[2]

      contrast <- contrast_calc(design.matrix, beta = beta, beta.var = beta.var,
                                p.value = p.value)

      temp <- setNames(contrast, temp$label)

      if (i > 1)
        temp <- c(L, temp)

      L <- temp

    } else {
      for (k in which(cond.interaction)){
        interaction.vars <- aux$var[sapply(aux$var, grepl, x = as.character(interaction[k]), fixed = TRUE)]
        others <- interaction.vars[interaction.vars != aux$var[i]]
        temp <- contrast_df(aux$data, aux$var[i], ref, others)
        design.matrix <- sapply(temp$new.data, function(x) model.matrix(fit, x), simplify = FALSE)

        drop <- which(grepl(aux$var[i], x = as.character(term.labels), fixed = TRUE))
        fit0 <- coxph(update.formula(fit$formula, paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
                      data = aux$data)
        p.value <- anova(fit0, fit)$`P(>|Chi|)`[2]

        contrast <- sapply(design.matrix, contrast_calc,
                           beta = beta, beta.var = beta.var, p.value = p.value, simplify = FALSE)
        contrast <- flattenlist(contrast)
        temp <- setNames(contrast, temp$label)

        if (i > 1)
          temp <- c(L, temp)

        L <- temp

      }
    }
  }

  term <- names(L)
  L <- Reduce(rbind, L)
  colnames(L) <- c("estimate", "conf.low", "conf.high", "p.value")
  out <- data.frame(term, L)

  return(out)
}

#'@export
table_fit <- function(fit, exponentiate = FALSE){

  out <- effect(fit)
  if (exponentiate)
    out[, apply(out, 2, is.numeric)] <- exp(out[, apply(out, 2, is.numeric)])
  return(out)
}

#'@importFrom stats model.matrix formula setNames anova vcov glm update.formula
effect.glm <- function(fit){














effect.glm <- function(fit){

  aux <- extract_data(fit)
  ref <- reference_df(aux$data)$df
  beta <- as.numeric(fit$coefficients)
  beta.var <- as.matrix(vcov(fit))

  effect <- list()

  for (i in 1:length(aux$var)){
    temp <- contrast_df(aux$data, aux$var[i], ref)
    design.matrix <- model.matrix(terms(fit), temp$newdata)
    pred <- design.matrix%*%beta

    for (j in 2:nrow(design.matrix)){
      estimate <- as.numeric(pred[j] - pred[1])
      diff <- design.matrix[j, ] - design.matrix[1, ]
      pred.se <- sqrt(t(diff)%*%beta.var%*%diff)

      lower <- estimate - 1.96*pred.se
      upper <- estimate + 1.96*pred.se
      effect[[temp$label[j-1]]] <- exp(cbind(estimate, lower, upper))
    }
  }

  Variable <- names(effect)
  effect <- Reduce(rbind, effect)
  colnames(effect) <- c("Effect", "Lower", "Upper")
  out <- data.frame(Variable, effect)

  return(out)
}

