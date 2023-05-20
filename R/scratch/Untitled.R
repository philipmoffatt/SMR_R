library(gitcreds)
library(usethis)

gitcreds_set()

use_github()


library(whitebox)
dem <- system.file("extdata", "DEM.tif", package="whitebox")
output <- file.path(getwd(), "output.tif")
feature_preserving_denoise(dem, output, filter=9, verbose_mode = FALSE)
