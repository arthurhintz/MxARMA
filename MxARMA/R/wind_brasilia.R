#' Monthly Wind Speed Data in Brasilia City
#'
#' This dataset contains monthly wind speed data recorded in Brasilia city.
#'
#' @format A data frame with columns:
#'   - \code{Month}: The month of the observation.
#'   - \code{WindSpeed}: The recorded wind speed in meters per second.
#'
#'
#' @return Um dataframe contendo os dados
#' @export
data <- function() {
  data_path <- system.file("data", "wind_speed_brasilia.txt", package = "MxARMA")
  dado <- read.table(data_path, sep = ";", header = TRUE)
  return(dado)
}
