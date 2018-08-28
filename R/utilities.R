effect <- function(fit) {
  UseMethod("effect")
}

flattenlist <- function(x){
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){
    Recall(out)
  }else{
    return(out)
  }
}

extract_label <- function(x, y){
  out <- ifelse(is.null(attr(x, "label")), y, attr(x, "label"))
  return(out)
}

extract_unit <- function(x){
  out <- ifelse(is.null(attr(x, "unit")), "", attr(x, "unit"))
  return(out)
}

label_unit <- function(x, y){
  if(x != ""){
    x <- parse(text = x)
    out <- bquote(.(y)~"("*.(x[[1]])*")")
  } else {
    out <- paste(y)
  }
  return(out)
}

keep_levels <- function(x, y){
  x$value <- factor(x$value, levels = y)
  return(x)
}

#https://stackoverflow.com/questions/19410108/cleaning-up-factor-levels-collapsing-multiple-levels-labels
recode_labels <- function(x, code, labels) {

  levels <- levels(x)
  for (i in 1:length(code))
    levels[which(code[i] == levels)] <- labels[i]
  levels(x) <- levels

  return(x)
}

#'@importFrom stats median quantile na.omit
median_iqr <- function (x, limit = "both", type = 7) {
  x <- na.omit(x)
  q25 <- quantile(x, probs = 0.25, type = type)
  q75 <- quantile(x, probs = 0.75, type = type)
  m <- median(x)

  if (limit == "both")
    out <- data.frame(y = m, ymin = q25, ymax = q75)
  if (limit == "lower")
    out <- data.frame(y = m, ymin = q25, ymax = m)
  if (limit == "upper")
    out <- data.frame(y = m, ymin = m, ymax = q75)
  if (limit == "none")
    out <- data.frame(y = m, ymin = m, ymax = m)
  return(out)
}

#'Format qualitative variables
#'
#'@description Recode factors, change the order of levels and
#'add a label to the qualitative variable.
#'
#'@param var vector to be formatted.
#'@param from an optional vector of original levels that \code{var} might have taken. It also requires \code{to}.
#'@param to an optional vector of values to recode the levels of \code{from}. It also requires \code{to}.
#'@param order an optional vector of values to establish the level order using the same values of \code{from}.
#'@param label an optional character to be attributed as a label to \code{var}.
#'
#'@return a vector  of class "factor" recoded if \code{from} and \code{to} are provided,
#'with an attribute label if \code{label} is provided,
#'and reordered if \code{order} is provided.
#'
#'@examples
#'library(dplyr)
#'library(magrittr)
#'data(iris)
#'
#'iris_nt <- iris %>% mutate(Species = ql_var(Species,
#'from = c("setosa", "versicolor", "virginica"),
#'to = c("Setosa", "Versicolor", "Virginica"), label = "Species",
#'order = c("Virginica", "Setosa", "Versicolor")))
#'
#'iris_nt %>% nt_describe(group = Species)
#'
#'@export
ql_var <- function(var, from = NULL, to = NULL, order = NULL, label = NULL){

  if (!is.null(attributes(var)$label) & is.null(label))
    label <- attributes(var)$label

  var <- as.factor(var)

  if (!is.null(from) & !is.null(to))
    var <- recode_labels(var, code = from, labels = to)

  if (!is.null(order))
    var <- factor(var, levels = order)

  if (!is.null(label))
    attributes(var)$label <- label

  out <- var
  return(out)
}


#'Format qualitative variables
#'
#'@description Recode factors, change the order of levels and
#'add a label to the qualitative variable.
#'
#'@param var vector to be formatted.
#'@param label an optional character to be attributed as a label to \code{var}.
#'@param unit an optional character to be attributed as a unit to \code{var}.
#'
#'@return a vector with an attribute label if \code{label} is provided,
#'and an attribute unit if \code{unit} is provided.
#'
#'@examples
#'library(dplyr)
#'library(magrittr)
#'data(iris)
#'
#'iris_nt <- iris %>%
#'mutate(Sepal.Length = qt_var(Sepal.Length, label = "Sepal Length", unit = "cm"),
#'       Sepal.Width = qt_var(Sepal.Width, label = "Sepal Width", unit = "cm"),
#'       Petal.Length = qt_var(Petal.Length, label = "Petal Length", unit = "cm"),
#'       Petal.Width = qt_var(Petal.Width, label = "Petal Length", unit = "cm"))
#'iris_nt %>% nt_describe(group = Species)
#'
#'@export
qt_var <- function(var, label = NULL, unit = NULL){

  if (!is.null(label))
    attributes(var)$label <- label

  if (!is.null(unit))
    attributes(var)$unit <- unit

  out <- var
  return(out)
}

#'Put Together
#'
#'@description Join a descriptive table and a p_value table.
#'
#'@param descriptive.tab a data frame of class "descriptive".
#'@param test.tab a data frame of class "p_value".
#'@param digits a numeric value specifying the number of digits for the p values.
#'@param alpha a numeric value specifying the significance level to indicate statistically significant multiple comparisons.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@examples
#'library(dplyr)
#'library(magrittr)
#'
#'iris_nt <- iris %>% filter(Species != "versicolor")
#'tab01 <- nt_describe(iris_nt, group = Species)
#'tab02 <- nt_compare_tg(iris_nt, group = Species)
#'tab <- put_together(tab01, tab02)
#'
#'@return A data frame similar to the \code{descriptive.tab} with an additional column of
#'p values from the \code{p_value.tab}
#'
#'@importFrom tidyr replace_na
#'@importFrom dplyr mutate select left_join
#'
#'@export
put_together <- function(descriptive.tab, test.tab, digits = 3,
                         alpha = 0.05,
                         save = FALSE,
                         file = "table"){

  if (!any(names(descriptive.tab) == "Variable"))
    stop("'descriptive.tab' does not have any column 'Variable'.")

  if (any(class(test.tab) == "two_groups") |
      any(class(test.tab) == "multiple_groups")){

    if (!any(names(test.tab) == "p value"))
      stop("'test.tab' does not have any column 'p value'.")

    test.tab <- test.tab %>% mutate(Variable = as.character(.data$Variable)) %>%
      select(.data$Variable, .data$`p value`)
    descriptive.tab <- descriptive.tab %>%
      mutate(Variable = as.character(.data$Variable))

    tab <- left_join(descriptive.tab, test.tab, by = "Variable")
    out <- tab %>% mutate(`p value` =
                            ifelse(.data$`p value` < 0.001, "< 0.001",
                                   round(.data$`p value`, digits)))  %>%
      replace_na(list(`p value` = ""))

    if (save)
      write.csv(tab, file = paste0(file, ".csv"))
  }

  if (any(class(test.tab$mc) == "multiple_comparisons")){

    if (!any(names(test.tab$omnibus) == "p value") &
        !any(names(test.tab$mc) == "p value"))
      stop("'test.tab' does not have any column 'p value'.")

    aux <- test.tab$omnibus %>%
      mutate(Variable = as.character(.data$Variable)) %>%
      select(.data$Variable, .data$`p value`)

    temp <- test.tab$mc %>% select(-.data$`95% CI`, -.data$Test, -.data$Group) %>%
      spread(key = .data$Hypothesis, value = "p value", drop = TRUE)

    for(i in 2:ncol(temp)){
      temp[, i] <- ifelse(temp[, i] < alpha, letters[i-1], "")
    }

    mc <- temp %>% unite(col = "Comparisons", -.data$Variable, sep = "")

    test.tab <- left_join(aux, mc, by = "Variable")
    tab <- left_join(descriptive.tab, test.tab, by = "Variable") %>%
      mutate(`p value` =
               ifelse(.data$`p value` < 0.001, "< 0.001",
                      round(.data$`p value`, digits))) %>%
      replace_na(list(`p value` = "", Comparisons = ""))
    comparisons <- setNames(names(temp[-1]), letters[1:(ncol(temp)-1)])

    out <- list(tab = tab, comparisons = comparisons)

    if (save){
      write.csv(tab, file = paste0(file, ".csv"))
      write.csv(comparisons, file = paste0(file, "_legend.csv"))
    }
  }
  return(out)
}

#'Step ribbon statistic
#'
#'Provides stairstep values for ribbon plots
#'
#'@inheritParams ggplot2::geom_ribbon
#'@param geom which geom to use; defaults to "`ribbon`"
#'@param direction \code{hv} for horizontal-veritcal steps, `vh`` for
#'   vertical-horizontal steps
#'@references \url{https://groups.google.com/forum/?fromgroups=#!topic/ggplot2/9cFWHaH1CPs}
#'@importFrom ggplot2 layer
#'@export
stat_stepribbon <- function(mapping=NULL, data=NULL, geom="ribbon",
                            position="identity",
                            na.rm=FALSE, show.legend=NA, inherit.aes=TRUE,
                            direction="hv", ...) {

  layer(
    data = data,
    mapping = mapping,
    stat = StatStepribbon,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      direction = direction,
      ...
    )
  )
}

#'ntimes-ggproto
#'
#'@format NULL
#'@usage NULL
#'@references \url{https://groups.google.com/forum/?fromgroups=#!topic/ggplot2/9cFWHaH1CPs}
#'@importFrom ggplot2 ggproto
#'@export
StatStepribbon <-
  ggproto(
    "StatStepRibbon", Stat,
    required_aes = c("x", "ymin", "ymax"),
    compute_group = function(data, scales, direction="hv",
                             yvars=c("ymin", "ymax"), ...) {
      stairstepn(data=data, direction=direction, yvars=yvars)
    }

  )

stairstepn <- function(data, direction="hv", yvars="y") {

  direction <- match.arg(direction, c("hv", "vh"))

  data <- as.data.frame(data)[order(data$x),]

  n <- nrow(data)

  if (direction == "vh") {
    xs <- rep(1:n, each=2)[-2*n]
    ys <- c(1, rep( 2:n, each=2))
  } else {
    ys <- rep(1:n, each=2)[-2*n]
    xs <- c(1, rep(2:n, each=2))
  }

  data.frame(
    x=data$x[xs],
    data[ys, yvars, drop=FALSE],
    data[xs, setdiff(names(data), c("x", yvars)), drop=FALSE]
  )

}
