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
#'
#'iris <- iris %>% mutate(species = ql_var(Species,
#'from = c("setosa", "Versicolor", "virginica"),
#'to = c("Setosa", "Versicolor", "Virginica"), label = "Species",
#'order = c("Virginica", "Setosa", "Versicolor"))) %>%
#'select(-Species)
#'nt_describe(iris)
#'
#'@export
ql_var <- function(var, from = NULL, to = NULL, order = NULL, label = NULL){

  if (!is.null(attributes(var)$label) & is.null(label))
    label <- attributes(var)$label

  var <- as.factor(var)

  if (!is.null(from) & !is.null(to))
    var <- factor(var, levels = from, labels = to)

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
#'
#'iris <- iris %>% mutate(
#'sepal.length = qt_var(Sepal.Length, label = "Sepal Length", unit = "cm"),
#'sepal.width = qt_var(Sepal.Width, label = "Sepal Width", unit = "cm"),
#'petal.length = qt_var(Petal.Length, label = "Petal Length", unit = "cm")) %>%
#'select(-Sepal.Length, -Sepal.Width, -Petal.Length)
#'nt_describe(iris)
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
#'@importFrom tidyr replace_na
#'@importFrom dplyr mutate select left_join
#'
#'@description Join a descriptive table and a p_value table.
#'
#'@param descriptive.tab a data frame of class "descriptive".
#'@param test.tab a data frame of class "p_value".
#'@param digits a numeric value specifying the number of digits for the p values.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@examples
#'library(dplyr)
#'library(magrittr)
#'library(ntimes)
#'
#'iris_nt <- iris %>% filter(Species != "versicolor")
#'tab01 <- nt_describe(iris_nt, group = Species)
#'tab02 <- nt_compare_tg(iris_nt, group = Species)
#'tab <- put_together(tab01, tab02)
#'
#'@return A data frame similar to the \code{descriptive.tab} with an additional column of
#'p values from the \code{p_value.tab}
#'
#'
#'@export
put_together <- function(descriptive.tab, test.tab, digits = 3,
                         save = FALSE,
                         file = "table"){

  if (!any(names(descriptive.tab) == "Variable"))
    stop("'descriptive.tab' does not have any column 'Variable'.")
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
    write.csv(out, file = paste0(file, ".csv"))

  return(out)
}
