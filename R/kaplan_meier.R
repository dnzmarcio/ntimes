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
#'@param time_points a numeric vector of time points to evaluate the survival curves.
#'@param risk_table a logical value indicating whether the risk table should be calculated.
#'@param save a logical value indicating whether the output
#'should be saved as a jpeg file.
#'@param fig_height a numeric value indicating the height (in) of the file.
#'@param fig_width a numeric value indicating the width (in) of the file.
#'@param std_fun a function to plot a barplot when \code{group = NULL}.
#'It must follow the same structure of \code{\link{std_barplot}}.
#'@param std_fun_group a function to plot a dotplot when \code{group}
#'is provided. It must follow the same structure of
#'\code{\link{std_barplot_group}}.
#'@param format a logical value indicating whether the output should be formatted.
#'@param digits a numerical value defining of digits to present the results.
#'@param file a character indicating the name of output file in csv format to be saved.
#'@param ... additional input arguments that may be used when creating your own function.
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
#'library(dplyr)
#'library(survival)
#'data(lung)
#'
#'lung_nt <- lung |> mutate(sex = factor(sex, levels = 1:2,
#'                                     labels = c("Female", "Male")),
#'                        ph.ecog = as.factor(ph.ecog)) |>
#'                        select(sex, ph.ecog, time, status)
#'lung_nt |> nt_km(time = time, status = status,
#'                  labels = list(sex = "Sex", ph.ecog = "ECOG"))
#'
#'@import ggplot2
#'@importFrom rlang enquo .data
#'@importFrom dplyr select mutate
#'@importFrom utils write.csv
#'@importFrom purrr map2
#'
#'@export
nt_km <-  function(data, time, status, labels = NULL,
                   xlab = "Time", ylab = "Survival",
                   base_size = 20,
                   time_points = NULL,
                   risk_table = TRUE, survival_table = NULL,
                   save = FALSE, fig_height = 5, fig_width = 7,
                   std_fun = std_km,
                   std_fun_group = std_km_group,
                   std_survival_table = std_get_survival_table,
                   format = TRUE, digits = 2,
                   file = "survival", where = NULL,
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
  tab <- list()

  overall <- std_fun(time, status,
                     xlab = xlab, ylab = ylab,
                     base_size = base_size,
                     time_points = time_points,
                     risk_table = risk_table,
                     survival_table = survival_table, ...)
  tab$overall <- tab_km(time, status, time_points, digits = digits)

  if(save)
    ggsave(overall$combined_plot, filename = paste0(where, "km_overall.jpeg"),
                      height = fig_height, width = fig_width)

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
                 base_size = base_size,
                 time_points = time_points,
                 risk_table = risk_table,
                 survival_table = survival_table,
                 fig_height = fig_height, fig_width = fig_width,
                 save = save, where = where, std_fun_group = std_fun_group,
                 ... = ...)

      aux <- pmap(.l = list(vars, vars.name, vars.label),
                  .f = tab_km_group,
                  time = time, status = status,
                  time_points = time_points, digits = digits)

      tab <- c(tab, aux)
  }


  plot$overall <- overall

  if (!is.null(time_points)){
    if (format){
      tab <-  map(tab, ~ .x  |>
        rename(Variable = .data$variable,
               Group = .data$group,
               Time = .data$time) |>
        mutate(`Survival (95% CI)` =
                 paste0(round(.data$survival, digits),
                        " (", round(.data$lower, digits),
                        " - ", round(.data$upper, digits), ")")) |>
        select(-.data$survival, -.data$lower, -.data$upper) |>
        select(.data$Variable, .data$Group, .data$Time, .data$`Survival (95% CI)`))

        aux <- bind_rows(tab)

      tab$overall <- tab$overall |> select(-.data$Group)

      aux <- bind_rows(tab) |>
        select(Variable, Group, Time, `Survival (95% CI)`)
    }

  } else {

    if (format){
      tab <-  map(tab, ~ .x  |>
                    rename(Variable = .data$variable,
                           Group = .data$group) |>
                    mutate(`Median (95% CI)` =
                             paste0(round(.data$survival, digits),
                                    " (", round(.data$lower, digits),
                                    " - ", round(.data$upper, digits), ")")) |>
                    select(-.data$survival, -.data$lower, -.data$upper) |>
                    select(.data$Variable, .data$Group, .data$`Median (95% CI)`))

      aux <- bind_rows(tab)

      tab$overall <- tab$overall |> select(-.data$Group)

      aux <- bind_rows(tab) |>
        select(Variable, Group, Time, `Median (95% CI)`)
    }
  }

  if (save){


    write.csv(aux, file = paste0(where, paste0(file, ".csv")))
  }
  rownames(tab) <- NULL

  out <- list(tab = tab, plot = plot)

  return(out)
}

#'@importFrom survival survfit Surv
tab_km <- function(time, status, time_points, digits){

  if (!is.null(time_points)){
    data_model <- data.frame(time, status)
    fit <- survfit(Surv(time, status) ~ 1, data = data_model)
    temp <- summary(fit, times = time_points, extend = TRUE)

    out <- data.frame(variable = "Overall",
                      group = NA, time = temp$time, survival = temp$surv,
                      lower = temp$lower, upper = temp$upper)

  } else {
    data_model <- data.frame(time, status)
    fit <- survfit(Surv(time, status) ~ 1, data = data_model)
    temp <- summary(fit)

    out <- data.frame(variable = "Overall",
                      group = NA,
                      survival = temp$table["median"],
                      lower = temp$table["0.95LCL"],
                      upper = temp$table["0.95UCL"])
  }

  out <- out |>
    mutate(variable = if_else(duplicated(variable), "", variable))
  rownames(out) <- NULL

  return(out)
}

#'@importFrom survival survfit Surv
#'@importFrom tidyr separate
#'@importFrom dplyr mutate
#'@importFrom rlang .data
tab_km_group <- function(var, var.name, var_label,
                         time, status, time_points, digits){

  if (!is.null(time_points)){
    data_model <- data.frame(time, status, var)
    fit <- survfit(Surv(time, status) ~ var, data = data_model)
    temp <- summary(fit, times = time_points, extend = TRUE)

    out <- data.frame(strata = temp$strata,
                      time = temp$time,
                      survival = temp$surv,
                      lower = temp$lower,
                      upper = temp$upper) |>
      separate(.data$strata, into = c("variable", "group"), sep = "=") |>
      mutate(variable = var_label)

  } else {
    data_model <- data.frame(time, status, var)
    fit <- survfit(Surv(time, status) ~ var, data = data_model)
    temp <- summary(fit)

    out <- data.frame(strata = rownames(temp$table),
                      survival = temp$table[, "median"],
                      lower = temp$table[, "0.95LCL"],
                      upper = temp$table[, "0.95UCL"]) |>
      separate(.data$strata, into = c("variable", "group"), sep = "=") |>
      mutate(variable = var_label,
             variable = if_else(duplicated(variable), "", variable))
  }


  rownames(out) <- NULL

  return(out)
}


aux_km <- function(var, var.name, var_label, time, status,
                   xlab, ylab, base_size,
                   time_points, risk_table, survival_table,
                   fig_height, fig_width, save, where, std_fun_group, ...){

  if (is.character(var))
    var <- as.factor(var)
  if (is.numeric(var))
    stop(paste0(var_label, " is numeric!"))

  if (nlevels(droplevels(var)) >= 2){
    out <- std_fun_group(time = time, status = status,
                         var = var, var_label = var_label,
                         xlab = xlab, ylab = ylab,
                         base_size = base_size,
                         time_points = time_points,
                         risk_table = risk_table,
                         survival_table = survival_table,
                         ...)

    if (save)
      ggsave(out$combined_plot,
             filename = paste0(where, paste0("km_", var.name, ".jpeg")),
             height = fig_height, width = fig_width)

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
#'@param time_points a numeric vector of time points to evaluate the survival curves.
#'@param risk_table a logical value indicating whether the risk table should be calculated.
#'@param ... additional input arguments that may be used when creating your own function.
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
#'@importFrom gridExtra tableGrob
#'
#'@export
std_km <- function(time, status, xlab, ylab,
                   time_points, base_size, risk_table, survival_table,
                   ...){

  ### Data
  data_model <- data.frame(time, status)
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = data_model)

  first.row <- data.frame(time = 0, n.risk = nrow(data_model),
                          n.event = 0, n.censor = 0,
                          estimate = 1, std.error = NA,
                          conf.high = 1, conf.low = 1)

  data.plot <- bind_rows(first.row, broom::tidy(fit))

  ### Basic plot
  surv_plot <- ggplot(data.plot, aes_string(x = "time", y = "estimate")) +
    geom_step()

  ### Formatting
  surv_plot <- surv_plot +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = base_size)

  if (!is.null(time_points)){
    surv_plot <- surv_plot +
      scale_x_continuous(limits = c(0, max(time)),
                         breaks = c(0, time_points))
  } else {
    surv_plot <- surv_plot +
      scale_x_continuous(limits = c(0, max(time)))
  }

  ### Changing from proportion to percentage
  surv_plot <- surv_plot +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1))

  ### Adding censor marks
  data.censor <- data.plot |> filter(.data$n.censor > 0)
  surv_plot <- surv_plot +
    geom_point(data = data.censor,
               aes_string(x = "time", y = "estimate"),
               shape = 124)

  ### Adding 95% confidence bands
  surv_plot <- surv_plot +
    geom_ribbon(data = data.plot,
                stat = "stepribbon",
                aes_string(ymin = "conf.low", ymax = "conf.high",
                           fill = "group"),
                alpha = 0.2, size = 0,
                fill = "grey80",
                linetype = "blank")
  tmp <- summary(fit)

  ### Adding a table for median survival
  if (!is.null(survival_table)){

    custom_theme <- ttheme_default(
      core = list(fg_params = list(fontsize = 0.7*base_size,
                                   col = "black"),
                  bg_params = list(fill = "white",
                                   col = "black", lwd = 1)),
      colhead = list(fg_params = list(fontsize = 0.8*base_size,
                                      fontface = "bold",
                                      col = "black"),
                     bg_params = list(fill = "white",
                                      col = "black", lwd = 1))
    )

    surv_table <- survival_table(fit,
                            var_label = NULL)

    if (!is.null(surv_table)){

      table_grob <- tableGrob(surv_table,
                              rows = NULL,
                              theme = custom_theme)

      if (min(data.plot$estimate) < 0.5){
        surv_plot2 <- surv_plot +
          annotation_custom(
            grob = table_grob,
            xmin = max(tmp$time)*0.6,
            xmax = max(tmp$time),
            ymin = 0.8, ymax = 1
          )
      } else {
        surv_plot2 <- surv_plot +
          annotation_custom(
            grob = table_grob,
            xmin = max(tmp$time)*0.6,
            xmax = max(tmp$time),
            ymin = 0.1, ymax = 0.3
          )
      }
    }

  } else {
    surv_table  <- NULL
    surv_plot2 <- surv_plot
  }

  ### Adding risk table
  if (risk_table){
    ## Data
    x_ticks <- ggplot_build(surv_plot)$layout$panel_params[[1]]$x$breaks
    table <- summary(fit, times = x_ticks, extend = TRUE)
    data_table <- data.frame(time = table$time, n.risk = table$n.risk)

    ## Basic plot
    risk_table <- ggplot(data_table, aes(x = .data$time, y = 1)) +
      geom_text(aes_string(label = "n.risk"), size = 0.3*base_size)

    ## Formatting
    risk_table <- risk_table +
      labs(x = xlab, y = "", title = "n at risk") +
      theme_classic(base_size = base_size) +
      theme(title = element_text(size = 0.8*base_size),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank())


    if (!is.null(time_points)){
      risk_table <- risk_table + scale_x_continuous(limits = c(0, max(time)),
                                                    breaks = c(0, time_points))
    } else {
      risk_table <- risk_table + scale_x_continuous(limits = c(0, max(time)))
    }


    ## Combining plots
    combined_plot <- surv_plot2 + risk_table +
      plot_layout(ncol = 1, heights = c(0.8, 0.2))

    out <- list(combined_plot = combined_plot,
                surv_plot = surv_plot,
                risk_table = risk_table,
                surv_table = surv_table)

  } else {
    out <- list(combined_plot = surv_plot2,
                surv_plot = surv_plot,
                surv_table = surv_table)
  }

  return(out)
}

#'Standard Kaplan-Meier curve by group
#'
#'@description A function to plot a Kaplan-Meier curve with groups.
#'
#'@param time a numeric vector.
#'@param status a numeric vector of '0' and '1'.
#'@param var a character vector.
#'@param var_label a character value specifying the group label.
#'@param xlab a character value specifying the x axis label.
#'@param ylab a character value specifying the y axis label.
#'@param time_points a numeric vector of time points to evaluate the survival curves.
#'@param risk_table a logical value indicating whether the risk table should be calculated.
#'@param ... additional input arguments that may be used when creating your own function.
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
#'@importFrom gridExtra tableGrob ttheme_default
#'@importFrom patchwork plot_layout
#'
#'@export
std_km_group <- function(time, status, var, var_label,
                         xlab, ylab, base_size, time_points, risk_table,
                         survival_table,
                         ...){

  ### Data
  data_model <- data.frame(time, status, var)
  fit <- survival::survfit(survival::Surv(time, status) ~ var, data = data_model)

  first.row <- data.frame(time = 0, n.risk = as.numeric(table(var)),
                          n.event = 0, n.censor = 0,
                          estimate = 1, std.error = NA,
                          conf.high = 1, conf.low = 1,
                          group = levels(var))

  data.plot <- broom::tidy(fit) |>
    tidyr::separate(.data$strata, into = c("var", "group"), sep = "r=") |>
    select(-var)
  data.plot <- bind_rows(first.row, data.plot) |>
    mutate(conf.high =
             ifelse(is.na(.data$conf.high), .data$estimate, .data$conf.high),
           conf.low =
             ifelse(is.na(.data$conf.low), .data$estimate, .data$conf.low)) |>
    mutate(group = factor(.data$group, levels = levels(var)))

  ### Basic plot
  surv_plot <- ggplot(data.plot, aes(x = .data$time, y = .data$estimate,
                                     colour = .data$group)) +
    geom_step()

  ### Formatting
  surv_plot <- surv_plot +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "top") +
    scale_colour_brewer(var_label, palette = "Set1", drop = FALSE)

  ### Specific time points
  if (!is.null(time_points)){
    surv_plot <- surv_plot +
      scale_x_continuous(limits = c(0, max(time)),
                         breaks = time_points)
  } else {
    surv_plot <- surv_plot +
      scale_x_continuous(limits = c(0, max(time)))
  }

  ### Changing from proportion to percentage
  surv_plot <- surv_plot +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, 1))

  ### Adding censor marks
  data.censor <- data.plot |> filter(.data$n.censor > 0)
  surv_plot <- surv_plot +
    geom_point(data = data.censor,
               aes_string(x = "time", y = "estimate"), shape = 124)

  ### Adding 95% confidence bands
  surv_plot <- surv_plot +
    geom_ribbon(data = data.plot,
                stat = 'stepribbon',
                aes_string(ymin = "conf.low",
                           ymax = "conf.high",
                           fill = "group"),
                alpha = 0.2, size = 0,
                linetype = "blank") +
    scale_fill_brewer(var_label, palette = "Set1", drop = FALSE)


  ### Adding p-values
  test <- survival::survdiff(survival::Surv(time, status) ~ var, data = data_model)
  p <- 1 - stats::pchisq(test$chisq, 1)
  p <- ifelse(round(p, 3) != 0,
              paste0("Logrank test \n p = ", round(p, 3)),
              "Logrank test \n p < 0.001")

  surv_plot <- surv_plot +
    annotate(geom = "text", label = p,
             x = -Inf, y = -Inf,
             hjust = -0.2,  vjust = -0.5, size = 0.3*base_size)

  tmp <- summary(fit)

  ### Adding survival table
  if (!is.null(survival_table)){

    custom_theme <- ttheme_default(
      core = list(fg_params = list(fontsize = 0.7*base_size,
                                   col = "black"),
                  bg_params = list(fill = "white",
                                   col = "black", lwd = 1)),
      colhead = list(fg_params = list(fontsize = 0.8*base_size,
                                      fontface = "bold",
                                      col = "black"),
                     bg_params = list(fill = "white",
                                      col = "black", lwd = 1))
    )

    surv_table <- survival_table(fit,
                            var_label = {{var_label}})

    table_grob <- tableGrob(surv_table,
                            rows = NULL,
                            theme = custom_theme)

    if (min(data.plot$estimate) < 0.5){
      surv_plot2 <- surv_plot +
        annotation_custom(
          grob = table_grob,
          xmin = max(tmp$time)*0.6,
          xmax = max(tmp$time),
          ymin = 0.8, ymax = 1
        )
    } else {
      surv_plot2 <- surv_plot +
        annotation_custom(
          grob = table_grob,
          xmin = max(tmp$time)*0.6,
          xmax = max(tmp$time),
          ymin = 0.1, ymax = 0.3
        )

    }
  } else {
    surv_table <- NULL
    surv_plot2 <- surv_plot
  }

  ### Adding risk table
  if (risk_table){
  ## Data
    x_ticks <- ggplot_build(surv_plot)$layout$panel_params[[1]]$x$breaks
    table <- summary(fit, times = x_ticks, extend = TRUE)
    data_table <- data.frame(time = table$time,
                             n_risk = table$n.risk,
                             group = table$strata) |>
      tidyr::separate(.data$group, into = c("var", "group"), sep = "r=") |>
      select(-var) |>
      mutate(group = factor(.data$group,
                            levels = rev(levels(var))))

    ## Basic plot
    risk_table <- ggplot(data_table, aes(x = .data$time, y = .data$group)) +
      geom_text(aes(label = .data$n_risk), size = 0.2*base_size)

    ## Formatting
    risk_table <- risk_table +
      theme_classic(base_size = base_size) +
      labs(x = xlab, y = "", title = "n at risk")

    if (!is.null(time_points)){
      risk_table <- risk_table +
        scale_x_continuous(limits = c(0, max(time)),
                           breaks = time_points)
    } else {
      risk_table <- risk_table +
        scale_x_continuous(limits = c(0, max(time)))
    }


    ## Changing y axis ticks
    colors <- unique(ggplot_build(surv_plot)$data[[1]]["colour"])

    risk_table <- risk_table +
      scale_y_discrete(labels = rep("-", nlevels(data_table$group))) +
      theme(title = element_text(size = 1*base_size),
            axis.text.y = element_text(colour = rev(colors[[1]]),
                                       face = "bold",
                                       size = 2.4*base_size,
                                       vjust = 0.3),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank())

    ## Combining plots
    combined_plot <- surv_plot2 + risk_table +
      plot_layout(ncol = 1, heights = c(0.8, 0.2))

    out <- list(combined_plot = combined_plot,
                surv_plot = surv_plot,
                risk_table = risk_table,
                surv_table = surv_table)
  } else {
    out <- list(combined_plot = surv_plot2,
                surv_plot = surv_plot,
                surv_table = surv_table)
  }

  return(out)
}
