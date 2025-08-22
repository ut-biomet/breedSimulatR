# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of useful functions not related to breeding simulation yet used in
# this package.


# Generate a squence of string with prefix and suffix formated like:
# "01", "02, ... , "10".
.charSeq <- function(prefix = "",
                     seq,
                     suffix = "") {
  sprintf(
    fmt = paste0(
      prefix,
      "%0", floor(log10(max(seq))) + 1, "i",
      suffix
    ),
    seq
  )
}




# Simulate linkage map position from physical position
#
# @param pos vector of physical position
# @param len chromosome length in bp
# @param lenCM chromosome length in centi-morgan
# @param b1 shape parameter 1
# @param b2 shape parameter 2
#
# @return vector of linkage map position
#
# @examples
# len <- 1000
# lenCM <- 100
# n <- 10
# b1 = 10
# b2 = 10
# plotDta <- data.frame(x = seq(0, len, length.out = 200)) %>%
# dplyr::mutate(y1 = simulLinkMapPos(x, len,lenCM,b1 = b1,b2 = b2))
# p <- plot_ly(
#   type = "scatter",
#   mode = "lines",
#   data = plotDta,
#   x = ~ x,
#   y = ~ y1,
#   name = "y1"
# )
.simulLinkMapPos <- function(pos, len, lenCM = 100, b1 = 10, b2 = 10) {
  stopifnot(b1 >= 2)
  stopifnot(b2 >= 2)
  stopifnot(length(len) == 1)
  stopifnot(length(lenCM) == 1)

  posCM <- suppressWarnings((pos / len + b1
    - stats::dbeta((pos / len), 1, b1)
    + stats::dbeta((pos / len), b2, 1))
  * lenCM / (1 + b1 + b2))
  return(posCM)
}
