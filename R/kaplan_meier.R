#'Kaplan-Meier plot
#'
#'@description Plot Kaplan-Meier curves for several variables.
#'
#'@param data a data frame with the variables.
#'@param time a data frame with the time to event variable.
#'@param status a data frame with the indicator associated to events.
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
#'library(dplyr)
#'data(lung)
#'
#'lung <- lung %>% mutate(sex = ql_var(sex, from = 1:2,
#'                                     to = c("Female", "Male"),
#'                                     label = "Sex"),
#'                        ph.ecog = ql_var(ph.ecog, label = "ECOG"))
#'
#'data_model <- lung %>% dplyr::select(sex, ph.ecog, time, status)
#'data_model %>% nt_km(time = time, status = status)
#'
#'
#'@importFrom broom tidy
#'@importFrom tidyr separate
#'@importFrom cowplot plot_grid
#'@importFrom scales percent
#'@importFrom survival survfit Surv survdiff
#'@importFrom stats pchisq
#'
#'@export
nt_km <-  function(data, time, status,
                   xlab = "Time", ylab = "Survival",
                   save = FALSE, fig.height = 5, fig.width = 5,
                   std_fun = std_km,
                   std_fun_group = std_km_group) {

  data <- as_data_frame(data)
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

  out <- list()

  gp <- std_fun(time, status, xlab = xlab, ylab = ylab)

  if(save)
    gp <- gp + ggsave(filename = "0_overall.jpeg",
                      height = fig.height, width = fig.width)


  if(ncol(data) > 2){
    out <- map2(.x = vars, .y = vars.name, .f = aux_km,
                time = time, status = status,
                xlab = xlab, ylab = ylab,
                fig.height = fig.height, fig.width = fig.width,
                save = save, std_fun_group = std_fun_group)
  }

  out$overall <- gp

  return(out)
}

aux_km <- function(var, var.name, time, status, xlab, ylab,
                   fig.height, fig.width, save, std_fun_group){

  if (is.character(var))
    var <- as.factor(var)
  if (is.numeric(var))
    stop(paste0(var.name, "is numeric!"))

  if (nlevels(droplevels(var)) >= 2){
    var.label <- extract_label(var, var.name)

    out <- std_fun_group(time = time, status = status,
                         var = var, var.label = var.label,
                         xlab = xlab, ylab = ylab)

    if (save)
      out <- out +
      ggsave(filename = paste0("_km_", var.name, ".jpeg"),
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
#'@importFrom survival survfit Surv
#'@importFrom tibble data_frame
#'@importFrom dplyr bind_rows filter
#'@importFrom cowplot plot_grid
#'@import ggplot2
#'@export
std_km <- function(time, status, xlab, ylab){

  ### Data
  data.model <- data.frame(time, status)
  fit <- survfit(Surv(time, status) ~ 1, data = data.model)

  first.row <- data_frame(time = 0, n.risk = nrow(data.model),
                          n.event = 0, n.censor = 0,
                          estimate = 1, std.error = NA,
                          conf.high = 1, conf.low = 1)

  data.plot <- bind_rows(first.row, tidy(fit))

  ### Basic plot
  surv.plot <- ggplot(data.plot, aes_string(x = "time", y = "estimate")) +
    geom_step()

  ### Formatting
  surv.plot <- surv.plot +
    scale_x_continuous(limits = c(0, max(time))) +
    labs(x = xlab, y = ylab) + theme_bw()

  ### Changing from proportion to percentage
  surv.plot <- surv.plot +
    scale_y_continuous(labels = percent, limits = c(0, 1))

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
  x.ticks <- ggplot_build(surv.plot)$layout$panel_ranges[[1]]$x.major_source
  table <- summary(fit, times = x.ticks)
  data.table <- data_frame(time = table$time, n.risk = table$n.risk)
  ## Basic plot
  risk.table <- ggplot(data.table, aes_string(x = "time", y = "1")) +
    geom_text(aes_string(label = "n.risk"))

  ## Formatting
  risk.table <- risk.table +
    scale_x_continuous(limits = c(0, max(time))) +
    labs(x = xlab, y = "", title = "n at risk") + theme_bw() +
    theme(title = element_text(size = 9),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  ## Combining plots
  out <- plot_grid(surv.plot, risk.table, nrow = 2,
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
#'@export
std_km_group <- function(time, status, var, var.label,
                         xlab, ylab){

  ### Data
  data.model <- data.frame(time, status, var)
  fit <- survfit(Surv(time, status) ~ var, data = data.model)

  first.row <- data.frame(time = 0, n.risk = as.numeric(table(var)),
                          n.event = 0, n.censor = 0,
                          estimate = 1, std.error = NA,
                          conf.high = 1, conf.low = 1,
                          group = levels(var))

  data.plot <- tidy(fit) %>%
    separate(.data$strata, into = c("var", "group"), sep = "r=") %>%
    select(-var)
  data.plot <- rbind(first.row, data.plot) %>%
    mutate(conf.high =
             ifelse(is.na(.data$conf.high), .data$estimate, .data$conf.high),
           conf.low =
             ifelse(is.na(.data$conf.low), .data$estimate, .data$conf.low))

  ### Basic plot
  surv.plot <- ggplot(data.plot, aes_string(x = "time", y = "estimate",
                                            colour = "group")) +
    geom_step()

  ### Formatting
  surv.plot <- surv.plot +
    scale_x_continuous(limits = c(0, max(time))) +
    labs(x = xlab, y = ylab) + theme_bw() +
    theme(legend.position = "top") +
    scale_colour_brewer(var.label, palette = "Set1", drop = FALSE)


  ### Changing from proportion to percentage
  surv.plot <- surv.plot + scale_y_continuous(labels = percent,
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
  test <- survdiff(Surv(time, status) ~ var, data = data.model)
  p <- 1 - pchisq(test$chisq, 1)
  p <- ifelse(round(p, 3) != 0,
              paste0("p = ", round(p, 3)),
              "p < 0.001")

  surv.plot <- surv.plot +
    annotate(geom = "text", label = p,
             x = -Inf, y = -Inf, hjust = -0.2,  vjust = -0.5, size = 3.5)

  ### Adding risk table

  ## Data
  x.ticks <- ggplot_build(surv.plot)$layout$panel_ranges[[1]]$x.major_source
  table <- summary(fit, times = x.ticks)
  data.table <- data.frame(time = table$time,
                           n.risk = table$n.risk,
                           group = table$strata) %>%
    separate(.data$group, into = c("var", "group"), sep = "r=") %>% select(-var) %>%
    mutate(group = factor(.data$group, levels = rev(levels(as.factor(.data$group)))))

  ## Basic plot
  risk.table <- ggplot(data.table, aes_string(x = "time", y = "group")) +
    geom_text(aes_string(label = "n.risk"), size = 3.5)

  ## Formatting
  risk.table <- risk.table + theme_bw() +
    scale_x_continuous(limits = c(0, max(time))) +
    labs(x = xlab, y = "", title = "n at risk")


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
  out <- plot_grid(surv.plot, risk.table, nrow = 2,
                   align = "v", rel_heights = c(3, 1))

  return(out)
}
