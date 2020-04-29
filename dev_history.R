library(devtools)
library(testthat)

#### Only one time ####

# use_git()

# use_r("calculations") # create new functions
# use_r("date") # create new functions

# use_build_ignore("dev_history.R")

# use_mit_license("The University of Tokyo, Laboratory of Biometry and Bioinformatics")

# use_testthat()

# use_spell_check()

# use_readme_rmd()



# use_lifecycle_badge("experimental")
# badger::badge_last_commit()
# badger::badge_codecov()



# use_vignette("example")

# use_gitlab_ci()
# use_github_action_check_full()
#
# use_pkgdown()
#
# use_github_action("pkgdown")

# use_package_doc()




#### Used regularly ####

load_all()

document()

use_tidy_description()
attachment::att_to_description(extra.suggests = c(
  "pkgdown", "covr", "knitr", "rmarkdown", "roxygen2"
))


test()
spell_check()
# spelling::update_wordlist()


check()

goodpractice::goodpractice(checks = goodpractice::all_checks()[-1])
covr::package_coverage()
covr::report()

install()
# build()

# use_github_release()
usethis::use_version()
# use_news_md()



#### pkgdown ####

# pkgdown::build_site()
# pkgdown::template_navbar()
# pkgdown::template_reference()
# pkgdown::clean_site()

