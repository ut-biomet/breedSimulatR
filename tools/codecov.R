# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# run code coverage and upload report on "codecov.io"


x <- covr::codecov(quiet = FALSE)
print(x)

if (x$meta$status != 200) {
  if (Sys.getenv("CODECOV_TOKEN") == "") {
    warning(paste0("CODECOV_TOKEN not found in the environment variables.\n",
                   "The report will not be upload on \"codecov.io\".\n",
                   "Check that variable in your CI settings.\n"))
  }
  stop(x$error$reason)
}

