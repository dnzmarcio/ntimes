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

  out <- list(data = droplevels(data), var = var)

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
#'@importFrom methods is
contrast_df <- function(data, var, ref, contrast.qt,
                        user.contrast = NULL,
                        interaction = NULL,
                        user.contrast.interaction = NULL,
                        table.reference){

  contrast <- ref

  if (any(is(data[[var]]) == "factor") | any(is(data[[var]]) == "character")){

    if(any(is(data[[var]]) == "character"))
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

    if (table.reference){
      label <- paste0(var, ":",contrast[[var]][1:nc])
    } else {
      label <- paste0(var, ":", contrast[[var]][2:nc], "/", contrast[[var]][1])
    }


  } else if (any(is(data[[var]]) == "numeric") |
             any(is(data[[var]]) == "integer")) {


    if (contrast.qt == "quantiles") {
      nc <- 2
      contrast <- contrast[rep(1, nc), ]
      quant <- quantile(data[[var]], probs = c(0.25, 0.75), na.rm = TRUE)
      contrast[[var]] <- round(quant, 2)

    } else if (contrast.qt == "one-unit"){
      nc <- 2
      contrast <- contrast[rep(1, nc), ]
      quant <- c(data[[var]][1], (data[[var]][1] + 1))
      contrast[[var]] <- quant

    } else if (contrast.qt == "user" & !is.null(user.contrast[[var]])){
      nc <- length(user.contrast[[var]])
      contrast <- contrast[rep(1, nc), ]
      contrast[[var]] <- user.contrast[[var]]

    } else if (contrast.qt == "user" & is.null(user.contrast[[var]])){
      stop('contrast.qt = "user" but user.contrast not found!')
    }

    if (contrast.qt == "one-unit"){
      label <- paste0(var, ":", "one-unit change")

    } else {

      if (table.reference ){
        label <- paste0(var, ":", contrast[[var]][1:nc])
      } else {
        label <- paste0(var, ":", contrast[[var]][2:nc], "/", contrast[[var]][1])
      }
    }



  }

  if (!is.null(interaction)){

    for(k in 1:length(interaction)){

      if (any(is(data[[interaction[k]]]) == "factor") |
          any(is(data[[interaction[k]]]) == "character")){

        if(any(is(data[[var]]) == "character"))
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

        if (!table.reference)
          nc <- 1

        if (k == 1){
          label <- paste0(rep(label, nc), " at ", interaction[k], " = ", rep(lv, each = length(label)))
        } else {
          label <- paste0(rep(label, nc), " at ", interaction[k], " = ", rep(lv, each = length(label)))
        }

      } else if (any(is(data[[interaction[k]]]) == "numeric") |
                 any(is(data[[interaction[k]]]) == "integer")) {


        if (contrast.qt == "quantiles") {

          nc <- 2
          lv <- quantile(data[[interaction[k]]], probs = c(0.25, 0.75), na.rm = TRUE)
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/nc)

        } else if (contrast.qt == "one-unit"){

          nc <- 2
          lv <- c(data[[var]][1], (data[[var]][1] + 1))
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/nc)

        }  else if (contrast.qt == "user" &
                    !is.null(user.contrast.interaction[[interaction[k]]])){

          nc <- length(user.contrast.interaction[[interaction[k]]])
          lv <- user.contrast.interaction[[interaction[k]]]
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(user.contrast.interaction[[interaction[k]]], nrow(contrast)/nc)

        } else if (contrast.qt == "user" &
                   is.null(user.contrast.interaction[[interaction[k]]])){
          stop('contrast.qt = "user" but user.contrast.interaction not found!')
        }

        if (!table.reference)
          nc <- 1

        if (k == 1){
          label <- paste0(rep(label, nc), " at ", interaction[k], " = ", rep(lv, each = nc))
        } else {
          label <- paste0(rep(label, nc), " at ", interaction[k], " = ", rep(lv, each = nc))
        }
    }
  }
    new.data <- split(contrast,  as.list(contrast)[interaction])
  } else {
    new.data <- contrast
  }

  rownames(new.data) <- NULL

  if (table.reference){
    nl <- ifelse(nlevels(data[[var]]) == 0,
                 ifelse(contrast.qt == "one-unit", 1, 2),
                 nlevels(data[[var]]))
    tmp <- ifelse(length(label) > 1, (length(label)-1), 1)
    seq_nl <- seq(1, tmp, by = nl)
    out <- list(new.data = new.data,
                label = label,
                seq_nl = seq_nl)
  } else {
    out <- list(new.data = new.data,
                label = label)
  }
  return(out)
}

#'@importFrom multcomp glht
#'@importFrom stats confint
#'@importFrom stats logLik
contrast_calc <- function(fit, fit0 = NULL, design.matrix, beta, beta.var, type){

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

  if (!is.na(pmatch(type, c("Wald", "wald")))){

    sm <- summary(test)
    estimate <- sm$test$coefficients
    pred.se <- sm$test$sigma

    lower <- estimate - 1.96*pred.se
    upper <- estimate + 1.96*pred.se
    estimate <- estimate
    p.value <- sm$test$pvalues

    out <- data.frame(estimate, lower, upper, p.value)

  } else if (!is.na(pmatch(type, c("profile", "lr")))) {
    ci <- confint(test)$confint
    p.value.lr <- 1 - pchisq(-2*(logLik(fit0)[[1]] - logLik(fit)[[1]]),
                             anova(fit0, fit)$Df[2])
    p.value <- c(p.value.lr, rep(NA, (nrow(ci) - 1)))

    out <- data.frame(estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3], p.value)
  }

  return(out)
}


