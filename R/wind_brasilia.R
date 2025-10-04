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
#' @keywords datasets
#' @examples
#' data("wind_speed_brasilia")
#' head(wind_speed_brasilia)
#'
#' # Se você distribuiu o .txt junto (inst/extdata), dá para ler assim:
#' # fp <- system.file("extdata", "wind_speed_brasilia.txt", package = "MxARMA")
#' # read.delim(fp, sep = ";")
"wind_speed_brasilia"
