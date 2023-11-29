#'Profile plot
#'
#'@description Plot profile plot for several variables.
#'
#'@param data a data frame with the variables.
#'@param group a character value indicating the group variable.
#'@param time a character value indicating the time variable.
#'@param labels a list of labels with components given by their variable names.
#'@param save a logical value indicating whether the output
#'should be saved as a jpeg file.
#'@param fig.height a numeric value indicating the height (in) of the file.
#'@param fig.width a numeric value indicating the width (in) of the file.
#'@param std_fun a function to plot a profile plot when \code{group = NULL}.
#'It must follow the same structure of \code{\link{std_profileplot}}.
#'@param std_fun_group a function to plot a profile plot when \code{group}
#'is provided. It must follow the same structure of
#'\code{\link{std_profileplot_group}}.
#'
#'@details The functions \code{\link{std_profileplot}} and
#'\code{\link{std_profileplot_group}} can be modified by the user in order to
#'customize the profile plots a prior.
#'The plots also can be modified a posterior as a regular ggplot object.
#'See \code{\link[ggplot2]{geom_errorbar}}, \code{\link{std_profileplot}} and
#'\code{\link{std_profileplot_group}}.
#'
#'@return a list of ggplot objects with each item named by the column names from
#'\code{var}.
#'
#'
#'@import ggplot2 dplyr
#'@importFrom rlang enquo quo_is_null
#'@importFrom dplyr select
#'@importFrom purrr map2
#'
#'@export
nt_profileplot <-  function(data, time = NULL, group = NULL,
                            labels = NULL,
                            save = FALSE, fig.height = 5, fig.width = 5,
                            std_fun = std_profileplot,
                            std_fun_group = std_profileplot_group, ...){

  group <- enquo(group)
  time <- enquo(time)
  #group_time <- list(group, time)

  if (!quo_is_null(group)){
    vars <- select(.data = data, -!!group, -!!time)
    group <- select(.data = data, !!group)
    time <- select(.data = data, !!time)
    group.name <- names(group)
    time.name <- names(time)
  } else {
    vars <- data
    group <- NULL
    group.name <- NULL
    time.name <- NULL
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

    if (!is.null(group)){
      time <- data_labeller(time, labels)
      time.label <- extract_label(time[[1]], group.name)
    }

  } else {
    vars.label <- map2(.x = vars, .y = as.list(vars.name),
                       .f = extract_label)

    if (!is.null(group)){
      group.label <- extract_label(group, group.name)
    }

    if (!is.null(time)){
      time.label <- extract_label(time, time.name)
    }
  }

  out <- pmap(.l = list(vars, vars.name, vars.label),
              .f = aux_profileplot,
              group = group, group.name = group.name, group.label = group.label,
              time = time, time.name = time.name, time.label = time.label,
              fig.height = fig.height, fig.width = fig.width, save = save,
              std_fun = std_fun, std_fun_group = std_fun_group,
              ... = ...)

  return(out)
}

aux_profileplot <- function(var, var.name, var.label,
                            group, group.name, group.label,
                            time, time.name, time.label,
                            fig.height, fig.width, save, std_fun, std_fun_group,
                            ...){

  out <- list()

  if (is.null(group)) {
    gp <- std_fun(var = var,
                  time = time,
                  var.label = var.label,
                  time.label = time.label,
                  ...)

    if(save)
      gp <- gp + ggsave(filename = paste0("profile_", var.name, ".jpeg"),
                        height = fig.height, width = fig.width)

    out <- gp

  } else {
    gp <- std_fun_group(var = var,
                        time = time[[1]],
                        group = group[[1]],
                        var.label = var.label,
                        group.label = group.label,
                        time.label = time.label,
                        ...)

    if (save)
      gp <- gp + ggsave(filename = paste0("profile_", group.name, "_",
                                          var.name, ".jpeg"),
                        height = fig.height,
                        width = fig.width)

    out <- gp
  }

  return(out)
}

#'Standard profile plot
#'
#'@import ggplot2
#'@description A function to plot a profile plot without groups.
#'
#'@param var a numeric vector.
#'@param var.label a character value specifying the variable label.
#'@param time a numeric vector of factor.
#'@param time.label a character value specifying the time label.
#'
#'@details This function defines the standard profile plot without groups to be
#'plotted by the function \code{\link{nt_profileplot}}. It can be modified by the
#'user. See more details in \code{\link[ggplot2]{geom_errorbar}}.
#'
#'@return a ggplot object.
#'
#'@export
std_profileplot <- function(var, time, var.label, time.label, ...){

  ### Data
  dp <- data.frame(var, time) |>
    na.exclude() |>
    group_by(time) |>
    summarize(mean = mean(var, na.rm = TRUE),
              se = sd(var, na.rm = TRUE)/sqrt(n()))


  ### Basic Plot
  out <- ggplot(dp, aes(x = .data$time, y = .data$mean)) +
    geom_point(position = position_dodge(0.05)) +
    geom_errorbar(aes(ymin = .data$mean-.data$se,
                      ymax = .data$mean+.data$se), width=.2,
                  position = position_dodge(0.05))

  ### Formatting
  out <- out +
    labs(y = var.label, x = time.label) +
    theme_classic()

  return(out)
}

#'Standard profile plot by group
#'
#'@description A function to plot a profile plot with groups.
#'
#'@param var a numeric vector.
#'@param group a character vector.
#'@param time a numeric vector.
#'@param var.label a character value specifying the variable label.
#'@param time.label a character value specifying the time label.
#'@param group.label a character value specifying the group label.
#'
#'@details This function defines the standard profile plot with groups to be plotted
#'by the function \code{\link{nt_profileplot}}. It can be modified by the user.
#'See more details in \code{\link[ggplot2]{geom_errorbar}}.
#'
#'@return a ggplot object.
#'
#'@export
std_profileplot_group <- function(var, time, group,
                                  var.label, time.label, group.label,
                                  ...){

  ### Data
  dp <- data.frame(time, group, var) |>
    na.exclude() |>
    group_by(group, time) |>
    summarize(mean = mean(var, na.rm = TRUE),
              se = sd(var, na.rm = TRUE)/sqrt(n()))

  ### Basic plot
  out <- ggplot(dp, aes(x = .data$time, color = .data$group, y = .data$mean)) +
    geom_point(position = position_dodge(0.05)) +
    geom_line(aes(group = group),
              position = position_dodge(0.05)) +
    geom_errorbar(aes(ymin = .data$mean-.data$se,
                      ymax = .data$mean+.data$se), width=.2,
                  position = position_dodge(0.05))

  ### Formatting
  out <- out + labs(y = var.label, x = time.label) +
    theme_classic() +
    theme(legend.position = "top") +
    scale_color_discrete(group.label)

  return(out)
}
