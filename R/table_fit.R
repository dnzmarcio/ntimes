extract_data <- function(fit){

  data <- eval(fit$call$data)
  temp <- get_all_vars(formula(fit), data)
  var.labels <- unlist(lapply(temp, function(x) attr(x, "label")))

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

reference_df <- function(fit){

  aux <- extract_data(fit)
  data <- aux$data

  df <- setNames(data.frame(matrix(ncol = ncol(data), nrow = 1)), names(data))

  for (i in 1:length(data)){
    if (is.numeric(data[, i])){
      df[, i] <- median(data[, i], na.rm = TRUE)
    } else if (is.factor(data[, i])){
      tab <- table(data[, i])
      df[, i] <- names(tab[which.max(tab)])
      df[, i] <- factor(df[, i], levels = levels(data[, i]))
    }
  }

  ref <- df[, which(sapply(aux$var, grepl,
                           x = names(df),
                           fixed = TRUE), arr.ind = TRUE)[, 1]]

  out <- list(df = df, ref = ref)

  return(out)
}

contrast_df <- function(data, var, ref){

  contrast <- ref

  if (class(data[[var]]) == "factor" | class(data[[var]]) == "character"){
    data[[var]] <- as.factor(data[[var]])
    lv <- levels(data[[var]])
    temp <- contrast[rep(1, each = length(lv)), ]
    temp[[var]] <- lv
    label <- paste0(var, ": ", temp[[var]][2:length(lv)], "/", temp[[var]][1])

  } else if (class(data[[var]]) == "numeric" |
             class(data[[var]]) == "integer") {
    quantiles <- quantile(data[[var]], probs = c(0.25, 0.75))
    temp <- contrast[c(1, 1), ]
    temp[[var]] <- round(quantiles, 2)
    label <- paste0(var, ": ", temp[[var]][2], "/", temp[[var]][1])
  }

  rownames(temp) <- NULL
  out <- list(newdata = temp, label = label)
  return(out)
}

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

effect.coxph <- function(fit){

  aux <- extract_data(fit)
  ref <- reference_df(fit)$df
  beta <- as.numeric(fit$coefficients)
  beta.var <- as.matrix(vcov(fit))
  term.labels <- attr(fit$terms, "term.labels")

  L <- list()

  for (i in 1:length(aux$var)){
    temp <- contrast_df(aux$data, aux$var[i], ref)
    design.matrix <- model.matrix(fit, temp$newdata)
    pred <- design.matrix%*%beta

    drop <- which(sapply(aux$var[i], grepl,
           x = as.character(term.labels),
           fixed = TRUE))

    fit0 <- coxph(update.formula(fit$formula, paste0(" ~ . - ", paste(term.labels[drop], collapse = " - "))),
          data = aux$data)

    p.value <- anova(fit0, fit)$`P(>|Chi|)`[2]

    for (j in 2:nrow(design.matrix)){
      estimate <- as.numeric(exp(pred[j] - pred[1]))
      diff <- design.matrix[j, ] - design.matrix[1, ]
      pred.se <- sqrt(t(diff)%*%beta.var%*%diff)

      lower <- exp(estimate - 1.96*pred.se)
      upper <- exp(estimate + 1.96*pred.se)

      if (j == 2){
        L[[temp$label[j-1]]] <- data.frame(estimate, lower, upper, p.value)
      } else {
        L[[temp$label[j-1]]] <- data.frame(estimate, lower, upper, p.value = NA)
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


