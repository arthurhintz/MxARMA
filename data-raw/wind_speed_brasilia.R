## code to prepare `wind_speed_brasilia` dataset goes here

usethis::use_data(wind_speed_brasilia, overwrite = TRUE)

# data-raw/wind_speed_brasilia.R

# Pacotes que ajudam (opcional)
# install.packages("readr"); install.packages("usethis")
library(readr)
library(usethis)

# Caminho do TXT:
# Opção A: (não será distribuído)
path <- file.path("data-raw", "wind_speed_brasilia.txt")

# Leia o TXT (ajuste separador, decimal e encoding conforme seu arquivo)
wind_speed_brasilia <- read_delim(
  path,
  delim = ";",
  locale = locale(encoding = "UTF-8")
)

# Ajuste de tipos / nomes (exemplo)
# wind_speed_brasilia <- transform(
#   wind_speed_brasilia,
#   Month = as.integer(Month),
#   WindSpeed = as.numeric(WindSpeed)
# )

# Salva como .rda no diretório data/ do pacote
use_data(wind_speed_brasilia, overwrite = TRUE, compress = "xz")
