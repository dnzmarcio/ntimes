#'Dotplot
#'
#'@description Plot dotplot for several variables.
#'
#'@param data a data frame with the variables.
#'@param group an optional data frame with the group variable.
#'@param binwidth a numerical vector specifying the bin width for each variable.
#'@param save a logical value indicating whether the output
#'should be saved as a jpeg file.
#'@param fig.height a numeric value indicating the height (in) of the file.
#'@param fig.width a numeric value indicating the width (in) of the file.
#'@param std_fun a function to plot a dotplot when \code{group = NULL}.
#'It must follow the same structure of the function \code{\link{std_dotplot}}.
#'@param std_fun_group a function to plot a dotplot when \code{group}
#'is provided. It must follow the same structure of the function
#'\code{\link{std_dotplot_group}}.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details The functions std_dotplot and std_dotplot_group can be
#'modified by the user in order to customize the dotplots a prior.
#'The plots also can be modified a posterior as a regular ggplot object.
#'See \code{\link[ggplot2]{geom_dotplot}}, \code{\link{std_dotplot}} and
#'\code{\link{std_dotplot_group}}.
#'
#'@return a list of ggplot objects with each item named by the column names from
#'\code{var}.
#'
#'@examples
#'data(iris)
#'iris |> nt_dotplot(binwidth = 0.05)
#'iris |> nt_dotplot(group = Species, binwidth = 0.05)
#'
#'@import ggplot2 dplyr
#'@importFrom rlang enquo quo_is_null
#'@importFrom dplyr select
#'@importFrom purrr pmap
#'
#'@export
nt_dotplot <-  function(data, group = NULL, binwidth = 1,
                        save = FALSE, fig.height = 5, fig.width = 5,
                        std_fun = std_dotplot,
                        std_fun_group = std_dotplot_group, ...) {

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

  if (length(binwidth) < length(vars.name))
    binwidth[(length(binwidth)+1):length(vars.name)] <- binwidth[length(binwidth)]

  L <- list(vars, vars.name, binwidth)

  out <- pmap(.l = L, .f = aux_dotplot,
              group = group, group.name = group.name,
              fig.height = fig.height, fig.width = fig.width, save = save,
              std_fun = std_fun, std_fun_group = std_fun_group,
              ... = ...)

  return(out)
}

aux_dotplot <- function(var, var.name, binwidth, group, group.name,
                        fig.height, fig.width, save,
                        std_fun, std_fun_group, ...){

  out <- list()
  var.label <- extract_label(var, var.name)
  var.unit <- extract_unit(var)

  if (var.unit != "")
    var.label <- paste(var.label, paste0("(", var.unit, ")"))

  if (is.null(group)) {
    gp <- std_fun(var = var,
                  var.label = var.label,
                  binwidth = binwidth,
                  ... = ...)

    if(save)
      gp <- gp + ggsave(filename = paste0("dot_", var.name, ".jpeg"),
                        height = fig.height, width = fig.width)

    out <- gp

  } else {
    group.label <- extract_label(group[[1]], group.name)

    gp <- std_fun_group(var = var,
                        group = group[[1]],
                        var.label = var.label,
                        group.label = group.label,
                        binwidth = binwidth,
                        ... = ...)

    if (save)
      gp <- gp + ggsave(filename =  paste0("dot_", group.name, "_",
                                 var.name, ".jpeg"),
                        height = fig.height,
                        width = fig.width)

    out <- gp
  }

  return(out)
}

#'Standard dotplot
#'
#'@description A function to plot a dotplot without groups.
#'
#'@param var a numeric vector.
#'@param var.label a character value specifying the variable label.
#'@param binwidth a numerical value specifying the bin width.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function defines the standard dotplot without groups to be
#'plotted by the function \code{\link{nt_dotplot}}. It can be modified by
#'the user. See more details in \code{\link[ggplot2]{geom_dotplot}}.
#'
#'@return a ggplot object.
#'
#'@importFrom stats median
#'@importFrom rlang .data
#'
#'@export
std_dotplot <- function(var, binwidth, var.label, ...){

  ### Data
  data_plot <- data.frame(var = var)

  ### Basic Plot
  out <- ggplot(data_plot, aes(y = .data$var, x = NA)) +
    geom_dotplot(binaxis = "y", stackdir = "center",
                 method = "histodot",
                 fill = "grey80",
                 binwidth = binwidth)

  ### Formatting
  out <- out +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) +
    labs(y = var.label)

  ### Adding summary
  out <- out + stat_summary(fun = stats::median, geom = "crossbar", width = 0.5)

  return(out)
}

#'Standard dotplot by group
#'
#'@description A function to plot a dotplot with groups.
#'
#'@param var a numeric vector.
#'@param group a character vector.
#'@param var.label a character value specifying the variable label.
#'@param group.label a character value specifying the group label.
#'@param binwidth a numerical value specifying the bin width.
#'@param ... additional input arguments that may be used when creating your own function.
#'
#'@details This function defines the standard dotplot with groups to be
#'plotted by the function \code{\link{nt_dotplot}}. It can be modified by
#'the user. See more details in \code{\link[ggplot2]{geom_dotplot}}.
#'
#'@return a ggplot object.
#'
#'@importFrom stats median
#'@importFrom rlang .data
#'
#'@export
std_dotplot_group <- function(var, group, binwidth,
                              var.label, group.label, ...){

  ### Data
  data_plot <- data.frame(var = var, group = group)

  ### Basic Plot
  out <- ggplot(data_plot,
                aes(y = .data$var, x = .data$group, fill = .data$group)) +
    geom_dotplot(binaxis = "y", stackdir = "center", method = "histodot",
                 fill = "grey80",
                 binwidth = binwidth)

  ### Formatting
  out <- out +
    theme_classic() +
    labs(y = var.label, x = group.label) +
    theme(legend.position = "none")

  ### Adding summary
  out <- out + stat_summary(fun = stats::median,
                            geom = "errorbar", width = 0.5, size = 2)

  return(out)
}
