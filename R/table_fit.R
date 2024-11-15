table_fit <- function(fit, type, exponentiate = FALSE){

  out <- effect(fit, type)
  if (exponentiate)
    out[, apply(out, 2, is.numeric)] <- exp(out[, apply(out, 2, is.numeric)])

  return(out)
}

#' extract_data <- function(fit){
#'
#'   UseMethod("extract_data", object = fit)
#' }
#'
#' #'@importFrom stats terms get_all_vars
#' extract_data.glm <- function(fit){
#'
#'   temp <- fit$data
#'   data <- get_all_vars(formula(fit), temp)
#'   var.names <- colnames(data)
#'
#'   # It is not clear the use of the code below.
#'   #suppressWarnings(data <- data[!apply(is.na(data), 1, any), , drop = FALSE])
#'   var <- lapply(setNames(as.list(names(data)), names(data)), grepl,
#'                 x = as.character(formula(fit)[3]),
#'                 fixed = TRUE)
#'   var <- names(which(unlist(var)))
#'
#'   if (!is.null(attr(terms(fit),"specials")$strata))
#'     var <- var[-(attr(terms(fit),"specials")$strata - 1)]
#'
#'   if (!is.null(attr(terms(fit),"specials")$random))
#'     var <- var[-(attr(terms(fit),"specials")$random - 1)]
#'
#'   out <- list(data = droplevels(data), var = var)
#'
#'   return(out)
#' }
#'
#'@importFrom stats terms get_all_vars
extract_data <- function(fit, data = NULL){

  if (is.null(data)){
    temp <- eval(fit$call$data)
  } else {
    temp <- data
  }

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
reference_df <- function(fit, data = NULL){

  aux <- extract_data(fit, data)
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
contrast_df <- function(data, var, ref, contrast_qt,
                        user_contrast = NULL,
                        interaction = NULL,
                        user_contrast_interaction = NULL,
                        table_reference){

  contrast <- ref

  if (any(is(data[[var]]) == "factor") | any(is(data[[var]]) == "character")){

    if(any(is(data[[var]]) == "character"))
      data[[var]] <- as.factor(data[[var]])

    if (is.null(user_contrast[[var]])) {
      nc <- nlevels(data[[var]])
      contrast <- contrast[rep(1, nc), ]
      lv <- levels(data[[var]])
      contrast[[var]] <- factor(lv, levels = lv)
    } else {
      nc <- length(user_contrast[[var]])
      contrast <- contrast[rep(1, nc), ]
      contrast[[var]] <- factor(user_contrast[[var]],
                                levels = user_contrast[[var]])
    }

    if (table_reference){
      label <- paste0(var, ":",contrast[[var]][1:nc])
    } else {
      label <- paste0(var, ":", contrast[[var]][2:nc], "/", contrast[[var]][1])
    }


  } else if (any(is(data[[var]]) == "numeric") |
             any(is(data[[var]]) == "integer")) {


    if (contrast_qt == "quartiles") {
      nc <- 2
      contrast <- contrast[rep(1, nc), ]
      quant <- quantile(data[[var]], probs = c(0.25, 0.75), na.rm = TRUE)
      contrast[[var]] <- round(quant, 2)

    } else if (contrast_qt == "one-unit"){
      nc <- 2
      contrast <- contrast[rep(1, nc), ]
      quant <- c(data[[var]][1], (data[[var]][1] + 1))
      contrast[[var]] <- quant

    } else if (contrast_qt == "user" & !is.null(user_contrast[[var]])){
      nc <- length(user_contrast[[var]])
      contrast <- contrast[rep(1, nc), ]
      contrast[[var]] <- user_contrast[[var]]

    } else if (contrast_qt == "user" & is.null(user_contrast[[var]])){
      stop('contrast_qt = "user" but user_contrast not found!')
    }

    if (contrast_qt == "one-unit"){
      label <- paste0(var, ":", "every 1 unit of change")

    } else {

      if (table_reference ){
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

        if (is.null(user_contrast_interaction[[interaction[k]]])) {

          nc <- nlevels(data[[interaction[k]]])
          lv <- levels(data[[interaction[k]]])
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/nc)
          contrast[[interaction[k]]] <- factor(contrast[[interaction[k]]], levels = lv)

          } else {

          nc <- length(user_contrast_interaction[[interaction[k]]])
          lv <- user_contrast_interaction[[interaction[k]]]
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(user_contrast_interaction[[interaction[k]]], nrow(contrast)/nc)
          contrast[[interaction[k]]] <- factor(contrast[[interaction[k]]],
                                               levels = user_contrast_interaction[[interaction[k]]])
        }

        if (!table_reference)
          nc <- 1

        if (k == 1){
          label <- paste0(rep(label, nc), " at ", interaction[k], " = ", rep(lv, each = length(label)))
        } else {
          label <- paste0(rep(label, nc), " at ", interaction[k], " = ", rep(lv, each = length(label)))
        }

      } else if (any(is(data[[interaction[k]]]) == "numeric") |
                 any(is(data[[interaction[k]]]) == "integer")) {


        if (contrast_qt == "quartiles") {

          nc <- 2
          lv <- quantile(data[[interaction[k]]], probs = c(0.25, 0.75), na.rm = TRUE)
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/nc)

        } else if (contrast_qt == "one-unit"){

          nc <- 2
          lv <- c(data[[var]][1], (data[[var]][1] + 1))
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(lv, nrow(contrast)/nc)

        } else if (contrast_qt == "user" &
                    !is.null(user_contrast_interaction[[interaction[k]]])){

          nc <- length(user_contrast_interaction[[interaction[k]]])
          lv <- user_contrast_interaction[[interaction[k]]]
          contrast <- contrast[rep(1:nrow(contrast), each = nc), ]
          contrast[[interaction[k]]] <- rep(user_contrast_interaction[[interaction[k]]], nrow(contrast)/nc)

        } else if (contrast_qt == "user" &
                   is.null(user_contrast_interaction[[interaction[k]]])){
          stop('contrast_qt = "user" but user_contrast_interaction not found!')
        }

        if (!table_reference)
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

  if (table_reference){
    nl <- ifelse(nlevels(data[[var]]) == 0,
                 ifelse(contrast_qt == "one-unit", 1, 2),
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
contrast_calc <- function(fit, fit0 = NULL, design_matrix, beta, beta_var, type){

  contrast_aux <- function(design_matrix){
    out <- sweep(design_matrix, 2, design_matrix[1, ])
    out <- out[-1, , drop = FALSE]
  }

  if (!is.list(design_matrix))
    design_matrix <- list(design_matrix)

  K <- sapply(design_matrix, FUN = contrast_aux, simplify = FALSE)
  K <- Reduce(rbind, K)
  if (!is.matrix(K))
    K <- matrix(K, nrow = 1)
  test <- glht(fit, linfct = K)

  sm <- summary(test)
  p_value_wald <- sm$test$pvalues
  ci <- confint(test)$confint
  p_value_lr <- 1 - pchisq(-2*(logLik(fit0)[[1]] - logLik(fit)[[1]]),
                           anova(fit0, fit)$Df[2])
  p_value_lr <- c(p_value_lr, rep(NA, (nrow(ci) - 1)))

  if (!is.na(pmatch(type, c("Wald", "wald")))){
    estimate <- sm$test$coefficients
    pred.se <- sm$test$sigma

    lower <- estimate - 1.96*pred.se
    upper <- estimate + 1.96*pred.se
    estimate <- estimate

  } else if (!is.na(pmatch(type, c("profile", "lr")))) {

    estimate <- ci[, 1]
    lower <- ci[, 2]
    upper <- ci[, 3]
  }

  out <- data.frame(estimate, lower, upper, p_value_wald, p_value_lr)

  return(out)
}


