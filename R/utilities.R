

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
