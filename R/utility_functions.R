#' Refresh Input Data
#'
#' This function conversts the data used in this study from a csv format to Rda format.
#' We do not anticipate this being used but include it for transparancy.
#'
#' @title Refresh Input Data
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords refresh
#'
#' @export
#'
#' @examples
#' refreshInputData()
#'

refreshInputData <- function() {
  pwd <- getwd()
  final.files <- c("kallisto0.05.csv", "kallisto0.2.csv", "RTqPCR.csv", "scSeq.csv")
  output.names <- c("kallisto0.05", "kallisto0.2", "RTqPCR", "scSeq")

  for (i in seq_along(final.files)) {
    x <- read.csv(
      paste0(
        pwd, "/inst/extdata/", final.files[i]
      )
    )
    save(x,
      file = paste0(
        pwd, "/data/", output.names[i], ".Rds"
      )
    )
  }
}
