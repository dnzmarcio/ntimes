table_fit <- function(fit, type, exponentiate = FALSE){

  out <- effect(fit, type)
  if (exponentiate)
    out[, apply(out, 2, is.numeric)] <- exp(out[, apply(out, 2, is.numeric)])

  return(out)
}

#'@importFrom stats terms get_all_vars
extract_data <- function(fit){

  temp <- eval(fit$call$data)
  data <- get_all_vars(formula(fit), temp)
  var.names <- colnames(data)
  var.labels <- mapply(extract_label, data, var.names, SIMPLIFY = FALSE)

  # It is not clear the use of the code below.
  #suppressWarnings(data <- data[!apply(is.na(data), 1, any), , drop = FALSE])
  var <- lapply(setNames(as.list(names(data)), names(data)), grepl,
                x = as.character(formula(fit)[3]),
                fixed = TRUE)
  var <- names(which(unlist(var)))

  if (!is.null(attr(terms(fit),"specials")$strata))
    var <- var[-(attr(terms(fit),"specials")$strata - 1)]

  if (!is.null(attr(terms(fit),"specials")$random))
    var <- var[-(attr(terms(fit),"specials")$random - 1)]

  out <- list(data = droplevels(data), var = var, var.labels = var.labels)

  return(out)
}

#'@importFrom stats setNames median
reference_df <- function(fit){

  aux <- extract_data(fit)
  data <- aux$data

  df <- setNames(data.frame(matrix(NA, ncol = ncol(data), nrow = 1)), names(data))

  for (i in 1:length(data)){
    if (is.numeric(data[, i])){
      df[, i] <- median(data[, i], na.rm = TRUE)
    } else {

      if (!is.factor(data[, i])){
        data[, i] <- as.factor(data[, i])
      }

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
contrast_df <- function(data, var, ref, user.contrast = NULL,
                        interaction = NULL, user.contrast.interaction = NULL){

  contrast <- ref

  if (class(data[[var]]) == "factor" | class(data[[var]]) == "character"){

    if(class(data[[var]]) == "character")
      data[[var]] <- as.factor(data[[var]])

    if (is.null(user.contrast[[var]])) {
      nc <- nlevels(data[[var]])
      contrast <- contrast[rep(1, nc), ]
      lv <- levels(data[[var]])
      contrast[[var]] <- factor(lv, levels = lv)
    } else {
      nc <- length(user.contrast[[var]])
      contrast <- contrast[rep(1, nc), ]
      contrast[[var]] <- factor(user.contrast[[var]],
                                levels = user.contrast[[var]])
    }

    label <- paste0(var, ":",contrast[[var]][1:nc])
    #paste0(var, ":", contrast[[var]][2:nc], "/", contrast[[var]][1])

  } else if (class(data[[var]]) == "numeric" |
             class(data[[var]]) == "integer") {


    if (is.null(user.contrast[[var]])) {
      nc <- 2
      contrast <- contrast[rep(1, nc), ]
      quant <- quantile(data[[var]], probs = c(0.25, 0.75), na.rm = TRUE)
      contrast[[var]] <- round(quant, 2)
    } else {
      nc <- length(user.contrast[[var]])
      contrast <- contrast[rep(1, nc), ]
      contrast[[var]] <- user.contrast[[var]]
    }

    label <- paste0(var, ":", contrast[[var]][1:nc])
      #paste0(var, ":", contrast[[var]][2:nc], "/", contrast[[var]][1])
  }

  if (!is.null(interaction)){

    for(k in 1:length(interaction)){

      if (class(data[[interaction[k]]]) == "factor" |
          class(data[[interaction[k]]]) == "character"){

        if(class(data[[var]]) == "character")
          data[[interaction[k]]] <- as.factor(data[[interaction[k]]])

        if (is.null(user.contrast.interaction[[interaction[k]]])) {
          nc <- nlevels(data[[interaction[k]]])
          lv <- levels(data[[interaction[k]]])
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/nc)
          contrast[[interaction[k]]] <- factor(contrast[[interaction[k]]], levels = lv)
        } else {
          nc <- length(user.contrast.interaction[[interaction[k]]])
          lv <- user.contrast.interaction[[interaction[k]]]
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(user.contrast.interaction[[interaction[k]]], nrow(contrast)/nc)
          contrast[[interaction[k]]] <- factor(contrast[[interaction[k]]],
                                               levels = user.contrast.interaction[[interaction[k]]])
        }

        if (k == 1){
          label <- paste0(rep(label, length(lv)), " at ", interaction[k], " = ", rep(lv, each = length(lv)))
        } else {
          label <- paste0(rep(label, length(lv)), " at ", interaction[k], " = ", rep(lv, each = length(lv)))
        }

      } else if (class(data[[interaction[k]]]) == "numeric" |
                 class(data[[interaction[k]]]) == "integer") {


        if (is.null(user.contrast.interaction[[interaction[k]]])) {
          nc <- 2
          quant <- quantile(data[[interaction[k]]], probs = c(0.25, 0.75), na.rm = TRUE)
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(quant, nrow(contrast)/nc)
        } else {
          nc <- length(user.contrast.interaction[[interaction[k]]])
          quant <- user.contrast.interaction[[interaction[k]]]
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(user.contrast.interaction[[interaction[k]]], nrow(contrast)/nc)
        }

        if (k == 1){
          label <- paste0(rep(label, length(lv)), " at ", interaction[k], " = ", rep(lv, each = length(lv)))
        } else {
          label <- paste0(rep(label, length(lv)), " at ", interaction[k], " = ", rep(lv, each = length(lv)))
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
#'@importFrom stats logLik
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

  if (type == "Wald"){

    sm <- summary(test)
    estimate <- sm$test$coefficients
    pred.se <- sm$test$sigma

    lower <- estimate - 1.96*pred.se
    upper <- estimate + 1.96*pred.se
    estimate <- estimate
    p.value <- sm$test$pvalues

    out <- data.frame(estimate, lower, upper, p.value)

  } else if (type == "profile") {
    ci <- confint(test)$confint
    p.value.lr <- 1 - pchisq(-2*(logLik(fit0)[[1]] - logLik(fit)[[1]]),
                             anova(fit0, fit)$Df[2])
    p.value <- c(p.value.lr, rep(NA, (nrow(ci) - 1)))

    out <- data.frame(estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3], p.value)
  }

  return(out)
}


