#' @export
run_hivmappr <- function() {
  appDir <- system.file("shiny", "hivmappr", package = "hivmappr")
  if (appDir == "") {
    stop("Could not find directory. Try re-installing `hivmappr`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}