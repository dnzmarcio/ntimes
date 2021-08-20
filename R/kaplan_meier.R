#'Kaplan-Meier plot
#'
#'@description Plot Kaplan-Meier curves for several variables.
#'
#'@param data a data frame with the variables.
#'@param time a data frame with the time to event variable.
#'@param status a data frame with the indicator associated to events.
#'@param labels a list of labels with components given by their variable names.
#'@param xlab a character value specifying the x axis label.
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
#'@param time.points a numeric vector of time points to evaluate the survival curves.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param file a character indicating the name of output file in csv format to be saved.
#'
#'@details The functions \code{\link{std_km}} and
#'\code{\link{std_km_group}} are standard functions that can be
#'modified by the user in order to customize the barplots a prior.
#'The plots also can be modified a posterior as a regular ggplot object.
#'See \code{\link{std_km}} and \code{\link{std_km_group}}.
#'
#'@return a list of ggplot objects with each item named by the column names from
#'\code{var}.
#'
#'@examples
#'library(survival)
#'data(lung)
#'
#'lung_nt <- lung %>% mutate(sex = factor(sex, levels = 1:2,
#'                                     labels = c("Female", "Male")),
#'                        ph.ecog = as.factor(ph.ecog)) %>%
#'                        select(sex, ph.ecog, time, status)
#'lung_nt %>% nt_km(time = time, status = status,
#'                  labels = list(sex = "Sex", ph.ecog = "ECOG"))
#'
#'@import ggplot2
#'@importFrom rlang enquo .data
#'@importFrom dplyr select mutate
#'@importFrom utils write.csv
#'@importFrom purrr map2
#'@importFrom magrittr %>%
#'
#'@export
nt_km <-  function(data, time, status, labels = NULL,
                   xlab = "Time", ylab = "Survival",
                   risktable.title = "n at risk",
                   save = FALSE, fig.height = 5, fig.width = 5,
                   std_fun = std_km,
                   std_fun_group = std_km_group,
                   time.points = NULL, format = TRUE, digits = 2,
                   file = "survival", where = "",
                   ...) {

  time <- enquo(time)
  status <- enquo(status)

  if (ncol(data) > 2){
    vars <- select(.data = data, -!!time)
    vars <- select(.data = vars, -!!status)
    vars.name <- names(vars)
  }

  time <- select(.data = data, !!time)
  time <- time[[1]]
  status <- select(.data = data, !!status)
  status <- as.numeric(as.factor(status[[1]]))

  plot <- list()

  overall <- std_fun(time, status, xlab = xlab, ylab = ylab,
                     time.points = time.points,
                     risktable.title = risktable.title, ...)
  if (!is.null(time.points))
    aux <- tab_km(time, status, time.points, digits = digits)

  if(save)
    ggsave(overall, filename = paste0(where, "km_overall.jpeg"),
                      height = fig.height, width = fig.width)

  if(ncol(data) > 2){
    if (!is.null(labels)){
      vars <- data_labeller(vars, labels)
      vars.label <- map2(.x = vars, .y = as.list(vars.name),
                         .f = extract_label)
    } else {
      vars.label <- map2(.x = vars, .y = as.list(vars.name),
                         .f = extract_label)
    }


    plot <- pmap(.l = list(vars, vars.name, vars.label),
                 .f = aux_km,
                 time = time, status = status,
                 xlab = xlab, ylab = ylab,
                 time.points = time.points,
                 risktable.title = risktable.title,
                 fig.height = fig.height, fig.width = fig.width,
                 save = save, std_fun_group = std_fun_group,
                 ... = ...)
    if (!is.null(time.points)){
      tab <- pmap(.l = list(vars, vars.name, vars.label),
                  .f = tab_km_group,
                  time = time, status = status,
                  time.points = time.points, digits = digits)
     tab <- bind_rows(aux, Reduce(rbind, tab))
    }
  } else {
    if (!is.null(time.points))
      tab <- aux
  }

  plot$overall <- overall

  if (!is.null(time.points)){
    if (format){
      tab <-  tab  %>%
        rename(Variable = .data$variable,
               Group = .data$group,
               Time = .data$time) %>%
        mutate(`Survival (CI 95%)` =
                 paste0(round(.data$survival, digits),
                        " (", round(.data$lower, digits),
                        " - ", round(.data$upper, digits), ")")) %>%
        select(-.data$survival, -.data$lower, -.data$upper) %>%
        select(.data$Time, .data$Variable, .data$Group, .data$`Survival (CI 95%)`)
    }

    if (save)
      write.csv(tab, file = paste0(where, paste0(file, ".csv")))

    out <- list(tab = tab, plot = plot)
  } else {
    plot$overall <- overall
    out <- list(plot = plot)
  }

  return(out)
}

#'@importFrom survival survfit Surv
tab_km <- function(time, status, time.points, digits){

  data.model <- data.frame(time, status)
  fit <- survfit(Surv(time, status) ~ 1, data = data.model)
  temp <- summary(fit, times = time.points)

  out <- data.frame(time = temp$time, variable = "Overall",
                    group = NA, survival = temp$surv,
                    lower = temp$lower, upper = temp$upper)

  return(out)
}

#'@importFrom survival survfit Surv
#'@importFrom tidyr separate
#'@importFrom dplyr mutate
#'@importFrom rlang .data
#'@importFrom magrittr %>%
tab_km_group <- function(var, var.name, var.label, time, status, time.points, digits){

  data.model <- data.frame(time, status)
  fit <- survfit(Surv(time, status) ~ var, data = data.model)
  temp <- summary(fit, times = time.points)

  out <- data.frame(time = temp$time,
                    strata = temp$strata,
                    survival = temp$surv,
                    lower = temp$lower,
                    upper = temp$upper) %>%
    separate(.data$strata, into = c("variable", "group"), sep = "=") %>%
    mutate(variable = var.label)

  return(out)
}


aux_km <- function(var, var.name, var.label, time, status,
                   xlab, ylab, time.points, risktable.title,
                   fig.height, fig.width, save, std_fun_group, ...){

  if (is.character(var))
    var <- as.factor(var)
  if (is.numeric(var))
    stop(paste0(var.label, " is numeric!"))

  if (nlevels(droplevels(var)) >= 2){
    out <- std_fun_group(time = time, status = status,
                         var = var, var.label = var.label,
                         xlab = xlab, ylab = ylab,
                         time.points = time.points,
                         risktable.title = risktable.title,
                         ...)

    if (save)
      ggsave(out, filename = paste0(where, paste0("km_", var.name, ".jpeg")),
             height = fig.height, width = fig.width)

  } else {
    out <- NA
    warning(paste0(var.name, " has only one level."))
  }


  return(out)
}


#'Standard Kaplan-Meier curve
#'
#'@description A function to plot a Kaplan-Meier curve without groups.
#'
#'@param time a numeric vector.
#'@param status a numeric vector of '0' and '1'.
#'@param xlab a character value specifying the x axis label.
#'@param ylab a character value specifying the y axis label.
#'
#'@details This function defines the standard of Kaplan-Meier curves
#'without groups to be plotted by the function \code{\link{nt_km}}.
#'
#'@return a ggplot object.
#'
#'@importFrom survival survfit
#'@importFrom broom tidy
#'@importFrom dplyr bind_rows
#'@importFrom cowplot plot_grid
#'@importFrom magrittr %>%
#'
#'@export
std_km <- function(time, status, xlab, ylab, time.points, risktable.title, ...){

  ### Data
  data.model <- data.frame(time, status)
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = data.model)

  first.row <- data.frame(time = 0, n.risk = nrow(data.model),
                          n.event = 0, n.censor = 0,
                          estimate = 1, std.error = NA,
                          conf.high = 1, conf.low = 1)

  data.plot <- bind_rows(first.row, broom::tidy(fit))

  ### Basic plot
  surv.plot <- ggplot(data.plot, aes_string(x = "time", y = "estimate")) +
    geom_step()

  ### Formatting
  surv.plot <- surv.plot +
    labs(x = xlab, y = ylab) + theme_bw()

  if (!is.null(time.points)){
    surv.plot <- surv.plot + scale_x_continuous(limits = c(0, max(time)),
                                                breaks = time.points)
  } else {
    surv.plot <- surv.plot + scale_x_continuous(limits = c(0, max(time)))
  }

  ### Changing from proportion to percentage
  surv.plot <- surv.plot +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1))

  ### Adding censor marks
  data.censor <- data.plot %>% filter(.data$n.censor > 0)
  surv.plot <- surv.plot +
    geom_point(data = data.censor,
               aes_string(x = "time", y = "estimate"),
               shape = 124)

  ### Adding 95% confidence bands
  surv.plot <- surv.plot +
    geom_ribbon(data = data.plot,
                stat = "stepribbon",
                aes_string(ymin = "conf.low", ymax = "conf.high",
                           fill = "group"),
                alpha = 0.2, size = 0, fill = "grey80")

  ### Adding risk table

  ## Data
  x.ticks <- ggplot_build(surv.plot)$layout$panel_params[[1]]$x$breaks
  table <- summary(fit, times = x.ticks)
  data.table <- data.frame(time = table$time, n.risk = table$n.risk)

  ## Basic plot
  risk.table <- ggplot(data.table, aes_string(x = "time", y = "1")) +
    geom_text(aes_string(label = "n.risk"))

  ## Formatting
  risk.table <- risk.table +
    labs(x = xlab, y = "", title = risktable.title) + theme_bw() +
    theme(title = element_text(size = 9),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


  if (!is.null(time.points)){
    risk.table <- risk.table + scale_x_continuous(limits = c(0, max(time)),
                                                  breaks = time.points)
  } else {
    risk.table <- risk.table + scale_x_continuous(limits = c(0, max(time)))
  }


  ## Combining plots
  out <- cowplot::plot_grid(surv.plot, risk.table, nrow = 2,
                   align = "v", rel_heights = c(3, 1))

  return(out)
}

#'Standard Kaplan-Meier curve by group
#'
#'@description A function to plot a Kaplan-Meier curve with groups.
#'
#'@param time a numeric vector.
#'@param status a numeric vector of '0' and '1'.
#'@param var a character vector.
#'@param xlab a character value specifying the x axis label.
#'@param ylab a character value specifying the y axis label.
#'@param var.label a character value specifying the group label.
#'
#'@details This function defines the standard of Kaplan-Meier curves
#'with groups to be plotted by the function \code{\link{nt_km}}.
#'It can be modified by the user.
#'
#'@return a ggplot object.
#'
#'@importFrom survival survfit survdiff
#'@importFrom broom tidy
#'@importFrom tidyr separate
#'@importFrom dplyr select bind_rows mutate filter
#'@importFrom scales percent
#'@importFrom stats pchisq
#'@importFrom cowplot plot_grid
#'@importFrom magrittr %>%
#'
#'@export
std_km_group <- function(time, status, var, var.label,
                         xlab, ylab, time.points, risktable.title,
                         ...){

  ### Data
  data.model <- data.frame(time, status, var)
  fit <- survival::survfit(survival::Surv(time, status) ~ var, data = data.model)

  first.row <- data.frame(time = 0, n.risk = as.numeric(table(var)),
                          n.event = 0, n.censor = 0,
                          estimate = 1, std.error = NA,
                          conf.high = 1, conf.low = 1,
                          group = levels(var))

  data.plot <- broom::tidy(fit) %>%
    tidyr::separate(.data$strata, into = c("var", "group"), sep = "r=") %>%
    select(-var)
  data.plot <- bind_rows(first.row, data.plot) %>%
    mutate(conf.high =
             ifelse(is.na(.data$conf.high), .data$estimate, .data$conf.high),
           conf.low =
             ifelse(is.na(.data$conf.low), .data$estimate, .data$conf.low)) %>%
    mutate(group = factor(.data$group, levels = levels(var)))

  ### Basic plot
  surv.plot <- ggplot(data.plot, aes_string(x = "time", y = "estimate",
                                            colour = "group")) +
    geom_step()

  ### Formatting
  surv.plot <- surv.plot +
    labs(x = xlab, y = ylab) + theme_bw() +
    theme(legend.position = "top") +
    scale_colour_brewer(var.label, palette = "Set1", drop = FALSE)

  ### Specific time points
  if (!is.null(time.points)){
    surv.plot <- surv.plot + scale_x_continuous(limits = c(0, max(time)),
                                                breaks = time.points)
  } else {
    surv.plot <- surv.plot + scale_x_continuous(limits = c(0, max(time)))
  }

  ### Changing from proportion to percentage
  surv.plot <- surv.plot + scale_y_continuous(labels = scales::percent,
                                              limits = c(0, 1))



  ### Adding censor marks
  data.censor <- data.plot %>% filter(.data$n.censor > 0)
  surv.plot <- surv.plot +
    geom_point(data = data.censor,
               aes_string(x = "time", y = "estimate"), shape = 124)

  ### Adding 95% confidence bands
  surv.plot <- surv.plot +
    geom_ribbon(data = data.plot,
                stat = 'stepribbon',
                aes_string(ymin = "conf.low",
                           ymax = "conf.high",
                           fill = "group"),
                alpha = 0.2, size = 0) +
    scale_fill_brewer(var.label, palette = "Set1", drop = FALSE)


  ### Adding p-values
  test <- survival::survdiff(survival::Surv(time, status) ~ var, data = data.model)
  p <- 1 - stats::pchisq(test$chisq, 1)
  p <- ifelse(round(p, 3) != 0,
              paste0("p = ", round(p, 3)),
              "p < 0.001")

  surv.plot <- surv.plot +
    annotate(geom = "text", label = p,
             x = -Inf, y = -Inf, hjust = -0.2,  vjust = -0.5, size = 3.5)

  ### Adding risk table

  ## Data
  x.ticks <- ggplot_build(surv.plot)$layout$panel_params[[1]]$x$breaks
  table <- summary(fit, times = x.ticks)
  data.table <- data.frame(time = table$time,
                           n.risk = table$n.risk,
                           group = table$strata) %>%
    tidyr::separate(.data$group, into = c("var", "group"), sep = "r=") %>%
    select(-var) %>%
    mutate(group = factor(.data$group, levels = rev(levels(as.factor(.data$group)))))

  ## Basic plot
  risk.table <- ggplot(data.table, aes_string(x = "time", y = "group")) +
    geom_text(aes_string(label = "n.risk"), size = 3.5)

  ## Formatting
  risk.table <- risk.table + theme_bw() +
    labs(x = xlab, y = "", title = risktable.title)

  if (!is.null(time.points)){
    risk.table <- risk.table + scale_x_continuous(limits = c(0, max(time)),
                                                  breaks = time.points)
  } else {
    risk.table <- risk.table + scale_x_continuous(limits = c(0, max(time)))
  }


  ## Changing y axis ticks
  colors <- unique(ggplot_build(surv.plot)$data[[1]]["colour"])

  risk.table <- risk.table +
    scale_y_discrete(labels = rep("-", nlevels(data.table$group))) +
    theme(title = element_text(size = 9),
          axis.text.y = element_text(colour = rev(colors[[1]]),
                                     face = "bold",
                                     size = 48,
                                     vjust = 0.3),
          axis.ticks.y = element_blank())

  ## Combining plots
  out <- cowplot::plot_grid(surv.plot, risk.table, nrow = 2,
                   align = "v", rel_heights = c(3, 1))

  return(out)
}
