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
  sprintf(fmt = paste0(prefix,
                       "%0", floor(log10(max(seq))) + 1, "i",
                       suffix),
          seq)
}
