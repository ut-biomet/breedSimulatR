# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# generate package documentation. Stop if any warning are generated.

warn <- 0
x <- withCallingHandlers({
  devtools::document()
}, warning=function(cond) {
  warn <<-  warn+1
})

if (warn > 0) {
  stop("warnings detected in documentation")
}