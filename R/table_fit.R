table_fit <- function(fit, type, exponentiate = FALSE){

  out <- effect(fit, type)
  if (exponentiate)
    out[, apply(out, 2, is.numeric)] <- exp(out[, apply(out, 2, is.numeric)])

  return(out)
}

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

  if (!is.null(attr(terms(fit),"specials")$random))
    var <- var[-(attr(terms(fit),"specials")$random - 1)]

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
      df[, i] <- factor(levels(data[, i])[1], levels = levels(data[, i]))
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
    label <- paste0(var, ":", contrast[[var]][2:length(lv)], "/", contrast[[var]][1])

  } else if (class(data[[var]]) == "numeric" |
             class(data[[var]]) == "integer") {
    quantiles <- quantile(data[[var]], probs = c(0.25, 0.75))
    contrast <- contrast[c(1, 1), ]
    contrast[[var]] <- round(quantiles, 2)
    label <- paste0(var, ":", contrast[[var]][2], "/", contrast[[var]][1])
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

#'@importFrom multcomp glht
#'@importFrom stats confint
contrast_calc <- function(fit, fit0, design.matrix, beta, beta.var, type){

  contrast_aux <- function(design.matrix){
    out <- sweep(design.matrix, 2, design.matrix[1, ])
    out <- out[-1, , drop = FALSE]
  }

  if (!is.list(design.matrix))
    design.matrix <- list(design.matrix)

  K <- sapply(design.matrix, FUN = contrast_aux, simplify = FALSE)
  K <- Reduce(rbind, K)
  if (!is.matrix(K))
    K <- matrix(K, nrow = 1)
  test <- glht(fit, linfct = K)

  if (type == "wald"){
    estimate <- summary(test)$test$coefficients
    pred.se <- summary(test)$test$sigma

    lower <- exp(estimate - 1.96*pred.se)
    upper <- exp(estimate + 1.96*pred.se)
    estimate <- exp(estimate)
    p.value <- c(p.value, rep(NA, (length(estimate) - 1)))

    out <- data.frame(estimate, lower, upper, p.value)

  } else {
    ci <- confint(test)$confint
    p.value.lr <- 1 - pchisq(-2*(logLik(fit0)[[1]] - logLik(fit)[[1]]),
                             anova(fit0, fit)$Df[2])
    p.value <- c(p.value.lr, rep(NA, (nrow(ci) - 1)))

    out <- data.frame(estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3], p.value)
  }

  return(out)
}


