#'Barplot
#'
#'@description Plot barplot for several variables.
#'
#'@param data a data frame with the variables.
#'@param group a character value indicating the group variable.
#'@param labels a list of labels with components given by their variable names.
#'@param ylab a character value specifying the y axis label.
#'@param save a logical value indicating whether the output
#'should be saved as a jpeg file.
#'@param fig.height a numeric value indicating the height (in) of the file.
#'@param fig.width a numeric value indicating the width (in) of the file.
#'@param std_fun a function to plot a barplot when \code{group = NULL}.
#'It must follow the same structure of \code{\link{std_barplot}}.
#'@param std_fun_group a function to plot a dotplot when \code{group}
#'is provided. It must follow the same structure of
#'\code{\link{std_barplot_group}}.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details The functions \code{\link{std_barplot}} and
#'\code{\link{std_barplot_group}} are stardard functions that can be
#'modified by the user in order to customize the barplots a prior.
#'The plots also can be modified a posterior as a regular ggplot object.
#'See \code{\link[ggplot2]{geom_bar}}, \code{\link{std_barplot}} and
#'\code{\link{std_barplot_group}}.
#'
#'@return a list of ggplot objects with each item named by the column names from
#'\code{var}.
#'
#'@examples
#'library(dplyr)
#'data(iris)
#'
#'vars <- iris |>
#'transmute(Species = Species,
#'          Sepal.Length.C = ifelse(Sepal.Length > 5, "> 5", "<= 5"),
#'          Sepal.Width.C = ifelse(Sepal.Width > 3.5, "> 3.5", "<= 3.5"))
#'vars |> nt_barplot(group = Species)
#'
#'@import ggplot2 dplyr
#'@importFrom rlang enquo quo_is_null
#'@importFrom dplyr select
#'@importFrom purrr map2
#'
#'@export
nt_barplot <-  function(data, group = NULL,
                        labels = NULL, ylab = "Percent (%)",
                        save = FALSE, fig.height = 5, fig.width = 5,
                        std_fun = std_barplot,
                        std_fun_group = std_barplot_group, ...) {

  group <- enquo(group)

  if (!quo_is_null(group)){
    vars <- select(.data = data, -!!group)
    group <- select(.data = data, !!group)
    group.name <- names(group)
  } else {
    vars <- data
    group <- NULL
    group.name <- NULL
  }
  vars.name <- names(vars)

  if (!is.null(labels)){
    vars <- data_labeller(vars, labels)
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)

    if (!is.null(group)){
      group <- data_labeller(group, labels)
      group.label <- extract_label(group[[1]], group.name)
    }

  } else {
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)

    if (!is.null(group)){
      group.label <- extract_label(group, group.name)
    }
  }

  out <- map2(.x = vars, .y = vars.name, .f = aux_barplot,
              group = group, group.name = group.name, ylab = ylab,
              fig.height = fig.height, fig.width = fig.width, save = save,
              std_fun = std_fun, std_fun_group = std_fun_group,
              ... = ...)

  return(out)
}

aux_barplot <- function(var, var.name, group, group.name, ylab,
                        fig.height, fig.width, save,
                        std_fun, std_fun_group, ...){

  out <- list()
  var.label <- extract_label(var, var.name)

  if (is.null(group)) {
    gp <- std_fun(var = var,
                  var.label = var.label,
                  ylab = ylab,
                  ... = ...)

    if(save)
      gp <- gp + ggsave(filename = paste0("bar_", var.name, ".jpeg"),
                        height = fig.height, width = fig.width)

    out <- gp

  } else {
    group.label <- extract_label(group[[1]], group.name)

    gp <- std_fun_group(var = var,
                        group = group[[1]],
                        group.label = group.label,
                        var.label = var.label,
                        ylab = ylab,
                        ... = ...)

    if (save)
      gp <- gp + ggsave(filename =
                          paste0("bar_", group.name, "_",
                                 var.name, ".jpeg"),
                        height = fig.height,
                        width = fig.width)

    out <- gp
  }

  return(out)
}

#'Standard barplot
#'
#'@description A function to plot a barplot without groups.
#'
#'@param var a character vector.
#'@param var.label a character value specifying the variable label.
#'@param ylab a character value specifying the y axis label.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function defines the standard barplot without groups to be
#'plotted by the function \code{\link{nt_barplot}}. It can be modified by the
#'user. See more details in \code{\link[ggplot2]{geom_bar}}.
#'
#'@return a ggplot object.
#'
#'@importFrom stats na.omit
#'@importFrom dplyr count mutate
#'
#'@export
std_barplot <- function(var, var.label, ylab, ...){

  ### Data
  data_plot <- data.frame(var = var)
  data_plot <- stats::na.omit(data_plot) |>
    count(var = factor(.data$var)) |>
    mutate(perc = round(prop.table(.data$n) * 100, 2),
           label = paste0(round(.data$perc, 2), '%', " (", .data$n, ")"))

  ### Basic plot
  out <- ggplot(data_plot, aes_string(x = "var", y = "perc")) +
    geom_bar(stat = 'identity', position = position_dodge(width = .9),
             fill = "grey80")

  ### Formatting
  out <- out + labs(y = ylab, x = var.label) +
    scale_y_continuous(limits = c(0, 105)) +
    theme_classic() + theme(legend.position = "none")

  ### Adding labels
  out <- out + geom_text(aes_string(y = "perc + 3", label = "label"),
                         position = position_dodge(width = .9), size = 3.5)

  return(out)
}

#'Standard barplot by group
#'
#'@description A function to plot a barplot with groups.
#'
#'@param var a character vector.
#'@param group a character vector.
#'@param var.label a character value specifying the variable label.
#'@param group.label a character value specifying the group label.
#'@param ylab a character value specifying the y axis label.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function defines the standard barplot with groups to be plotted
#'by the function \code{\link{nt_barplot}}. It can be modified by the user.
#'See more details in \code{\link[ggplot2]{geom_bar}}.
#'
#'@return a ggplot object.
#'
#'@importFrom stats na.omit
#'@importFrom dplyr count group_by mutate
#'
#'@export
std_barplot_group <- function(var, group, var.label, group.label, ylab, ...){

  ### Data
  data_plot <- data.frame(var = var, group = group)
  data_plot <- stats::na.omit(data_plot) |>
    count(group = factor(.data$group), var = factor(.data$var)) |>
    group_by(.data$group) |>
    mutate(perc = round(prop.table(.data$n) * 100, 2),
           label = paste0(round(.data$perc, 2), '%', " (", .data$n, ")"))

  ### Basic plot
  out <- ggplot(data_plot, aes_string(x = "group",
                                      y = "perc",
                                      fill = "var")) +
    geom_bar(stat = 'identity', position = position_dodge(width = .9))

  ### Formatting
  out <- out + labs(y = ylab, x = group.label) +
    scale_y_continuous(limits = c(0, 105)) +
    scale_fill_grey(var.label) +
    theme_classic() + theme(legend.position = "top")

  ### Adding labels
  out <- out + geom_text(aes_string(y = "perc + 3", label = "label"),
                         position = position_dodge(width = .9), size = 3.5)

  return(out)
}
