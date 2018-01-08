#'Boxplot
#'
#'@description Plot boxplot for several variables.
#'
#'@importFrom purrr map2
#'@importFrom dplyr select
#'@importFrom magrittr %>%
#'@importFrom rlang .data quo_is_null enquo
#'@importFrom tibble as_data_frame
#'
#'@param data a data frame with the variables.
#'@param group an optional data frame with the group variable.
#'@param save a logical value indicating whether the output
#'should be saved as a jpeg file.
#'@param fig.height a numeric value indicating the height (in) of the file.
#'@param fig.width a numeric value indicating the width (in) of the file.
#'@param std_fun a function to plot a boxplot when \code{group = NULL}.
#'It must follow the same structure of \code{\link{std_boxplot}}.
#'@param std_fun_group a function to plot a boxplot when \code{group}
#'is provided. It must follow the same structure of
#'\code{\link{std_boxplot_group}}.
#'
#'@details The functions \code{\link{std_boxplot}} and
#'\code{\link{std_boxplot_group}} can be modified by the user in order to
#'customize the boxplots a prior.
#'The plots also can be modified a posterior as a regular ggplot object.
#'See \code{\link[ggplot2]{geom_boxplot}}, \code{\link{std_boxplot}} and
#'\code{\link{std_boxplot_group}}.
#'
#'@return a list of ggplot objects with each item named by the column names from
#'\code{var}.
#'
#'@examples
#'library(magrittr)
#'
#'iris %>% nt_boxplot(group = Species)
#'
#'@export
nt_boxplot <-  function(data, group = NULL,
                        save = FALSE, fig.height = 5, fig.width = 5,
                        std_fun = std_boxplot,
                        std_fun_group = std_boxplot_group){

  data <- as_data_frame(data)
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

  out <- map2(.x = vars, .y = vars.name, .f = aux_boxplot,
              group = group, group.name = group.name,
              fig.height = fig.height, fig.width = fig.width, save = save,
              std_fun = std_fun, std_fun_group = std_fun_group)

  return(out)
}

aux_boxplot <- function(var, var.name, group, group.name,
                        fig.height, fig.width, save, std_fun, std_fun_group){

  out <- list()
  var.label <- extract_label(var, var.name)

  if (is.null(group)) {
    gp <- std_fun(var = var,
                  var.label = var.label)

    if(save)
      gp <- gp + ggsave(filename = paste0(var.name, ".jpeg"),
                        height = fig.height, width = fig.width)

    out <- gp

  } else {
    group.label <- extract_label(group[[1]], group.name)

    gp <- std_fun_group(var = var,
                        group = group[[1]],
                        group.label = group.label,
                        var.label = var.label)

    if (save)
      gp <- gp + ggsave(filename =
                          paste0(var.name, "_",
                                 group.name, ".jpeg"),
                        height = fig.height,
                        width = fig.width)

    out <- gp
  }

  return(out)
}

#'Standard boxplot
#'
#'@import ggplot2
#'@description A function to plot a boxplot without groups.
#'
#'@param var a numeric vector.
#'@param var.label a character value specifying the variable label.
#'
#'@details This function defines the standard boxplot without groups to be
#'plotted by the function \code{\link{nt_boxplot}}. It can be modified by the
#'user. See more details in \code{\link[ggplot2]{geom_boxplot}}.
#'
#'@return a ggplot object.
#'
#'@export
std_boxplot <- function(var, var.label){

  ### Data
  data_plot <- data.frame(var = var)

  ### Basic Plot
  out <- ggplot(data_plot, aes_string(y = "var")) +
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot(fill = "grey80", outlier.shape = NA)  +
    geom_jitter(shape = 16, size = 1.5,
                position = position_jitter(width = 0.2))

  ### Formatting
  out <- out + labs(y = var.label) +
    theme_classic()

  return(out)
}

#'Standard boxplot by group
#'
#'@description A function to plot a boxplot with groups.
#'
#'@param var a numeric vector.
#'@param group a character vector.
#'@param var.label a character value specifying the variable label.
#'@param group.label a character value specifying the group label.
#'
#'@details This function defines the standard boxplot with groups to be plotted
#'by the function \code{\link{nt_boxplot}}. It can be modified by the user.
#'See more details in \code{\link[ggplot2]{geom_boxplot}}.
#'
#'@return a ggplot object.
#'
#'@export
std_boxplot_group <- function(var, group, var.label, group.label){

  ### Data
  data_plot <- data_frame(var = var, group = group)

  ### Basic plot
  out <- ggplot(data_plot, aes_string(y = "var", x = "group")) +
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot(outlier.shape = NA, fill = "grey80")  +
    geom_jitter(shape = 16, size = 1.5,
                position = position_jitter(width = 0.2))

  ### Formatting
  out <- out + labs(y = var.label, x = group.label) +
    theme_classic() +
    theme(legend.position = "none")

  return(out)
}