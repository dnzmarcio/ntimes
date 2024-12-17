

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

data_labeller <- function(data, labels){

  aux_label <- function(x, y){
    attributes(x)$label <- y
    return(x)
  }

  var.names <- colnames(data)
  names(labels) <- paste0("^",  names(labels), "$")
  labels <- str_replace_all(var.names, unlist(labels))

  out <- bind_cols(map2(.x = data, .y = labels, .f = aux_label))
  return(out)
}

#https://stackoverflow.com/questions/19410108/cleaning-up-factor-levels-collapsing-multiple-levels-labels
recode_labels <- function(x, code, labels) {

  x <- as.factor(x)
  levels <- levels(x)
  for (i in 1:length(code))
    levels[which(code[i] == levels)] <- labels[i]
  levels(x) <- levels

  return(x)
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



#'@importFrom dplyr group_by mutate ungroup summarize filter select
#'@importFrom tidyr unnest separate_wider_delim
#'@export
get_survival_table <- function(title = NULL,
                               type = "median", time_points = NULL){

  out  <- function(survfit_obj,
                   var_label) {

    if (is.null(title)){
      if (type == "median"){
        title <- "Median Survival Time (95% CI)"
      } else if (type == "time_points") {
        title <- "Survival (95% CI)"
      }
    }

    s <- summary(survfit_obj)

    if (is.null(s$strata)) {

      if (type == "median" & !is.na(s$table["median"])){
        median <- paste0(round(s$table["median"], 1),
                         ", (",
                         round(s$table["0.95LCL"], 1),
                         " ; ",
                         round(s$table["0.95UCL"], 1),
                         ")")
        events <- paste0(s$table["events"], "/", s$table["n.start"])

        out <-
          tibble({{title}} :=
                   median,
                 `Events/Total` = events)

      } else if (type == "time_points"){
        aux <- summary(survfit_obj, times = time_points)

        survival <- paste0(round(100*aux$surv, 1),
                           ", (",
                           round(100*aux$lower, 1), " ; ",
                           round(100*aux$upper, 1), ")")
        n.events <- sapply(time_points, function(x) sum(s$n.event[s$time < x]))
        events = paste0(n.events, "/", rep(s$table["n.start"], length(n.events)))

        out <-
          tibble(Time = time_points,
                 {{title}} := survival,
                 `Events/Total` = events)
      } else {
        out <- NULL

      }


    } else {

      if (type == "median" & any(!is.na(s$table[, "median"]))) {
        if (all(is.na(s$table[, "median"])))
          stop("There is no median survival time!")

        median <- paste0(round(s$table[, "median"], 1),
                         ", (",
                         round(s$table[, "0.95LCL"], 1),
                         " ; ",
                         round(s$table[, "0.95UCL"], 1),
                         ")")
        events <- paste0(s$table[, "events"], "/", s$table[, "n.start"])

        Groups <- sub(".*=", "", levels(s$strata))
        out <-
          tibble({{var_label}} := Groups,
                 {{title}} := median,
                 `Events/Total` = events)

      } else if (type == "time_points"){

        aux <- function(tp, time) {
          valid_times <- time[time <= tp]
          if (length(valid_times) == 0) NA else valid_times[which.min(abs(valid_times - tp))]
        }

        # Create a data frame with the event data
        event_data <- data.frame(
          time = s$time,
          group = s$strata, # Repeat group labels according to strata sizes
          n_risk = s$n.risk,    # number of subjects at risk
          n_event = s$n.event,   # number of events at each time point
          surv = s$surv,
          lower = s$lower,
          upper = s$upper
        )

        # Compute the cumulative sum of events for each group
        event_data <- event_data |>
          group_by(group) |>
          mutate(cum_event = cumsum(n_event)) |>  # cumulative sum of events by group
          ungroup()

        closest_events <- event_data |>
          group_by(group) |>
          summarize(
            closest_time = list(sapply(time_points, aux, time = time)),
            target_time = list(time_points)
          ) |>
          unnest(c(target_time, closest_time)) |>
          filter(!is.na(closest_time)) |>
          mutate(n.start = unlist(lapply(s$table[, "n.start"], function(x)
            rep(x, each = length(time_points))))) |>
          left_join(event_data, by = c("group", "closest_time" = "time")) |>
          select(group, target_time, cum_event, n.start, surv, lower, upper)

        out <- closest_events |>
          separate_wider_delim(group, "=", names = c("var", "group")) |>
          mutate({{var_label}} := if_else(duplicated(group), "", group),
                 Time = target_time,
                 {{title}} := paste0(round(100*surv, 1),
                                             ", (",
                                             round(100*lower, 1), " ; ",
                                             round(100*upper, 1), ")"),
                 Events = paste0(cum_event, "/", n.start)) |>
          select({{var_label}}, Time, {{title}}, Events)

      } else {
        out <- NULL
      }

    }

    return(out)
  }
}

