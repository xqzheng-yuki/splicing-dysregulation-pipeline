# utilities for plotting
generate_shades <- function(base_color) {
  # Convert hex to RGB
  base_rgb <- col2rgb(base_color) / 255
  shades <- colorRampPalette(c("white", rgb(base_rgb[1], base_rgb[2], base_rgb[3])))(5)
  return(shades[2:length(shades)])
}

get_ylim <- function(first_set,second_set,condition_tag = NULL) {
  ylim <- c(0,max(score(unlist(GRangesList(c(first_set,second_set))))))
  if (!is.null(condition_tag)) {
    ylim <- c(0,max(score(unlist(GRangesList(c(first_set[[condition_tag]],second_set[[condition_tag]]))))))
  }
  return(ylim)
}

enlarge_gr <- function(goi) {
  gr <- geneRanges_GRCm39[goi]
  embiggen_factor <- round(width(gr)*0.1)
  start(gr) <- start(gr) - embiggen_factor
  end(gr) <- end(gr) + embiggen_factor
  return(gr)
}

add_track <- function(new_track,tracks=list()) {
  if (is.list(new_track)) {
    tracks <- append(tracks,new_track)
  } else {
  tracks <- append(tracks,list(new_track))
  }
  return(tracks)
}

#' add_new_plots
#' update the state of the app with the most recent plot
#' @param state
#' @param classic_plot_func 
#' @param dataset_plot_func 
#'
#' @return state
#' @export
#'
#' @examples
add_new_plots <- function(state, classic_plot_func, dataset_plot_func) {
  state$results <- append(state$results, list(
    list(classic = classic_plot_func, dataset = dataset_plot_func)
  ))
  state$current_index <- length(state$results) # Set to the latest added
  # state$results[[state$current_index + 1]] <- recordPlot()
  info(logger,paste0("The state is ",length(state$results)))
  state
}