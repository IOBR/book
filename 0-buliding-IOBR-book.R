



library(usethis)
# usethis::use_git()
# use_github(organisation = "DongqiangZeng0808")


# use_build_ignore("5-update-gene-signatures.R")
# use_git_ignore("5-update-gene-signatures.R")


# install.packages("remotes")
# remotes::install_github("rstudio/bookdown")

# remotes::install_github("hadley/emo")

library(future)
options(future.globals.maxSize = 8 * 1024^3)  # 2GB


bookdown::render_book("index.Rmd", output_format="bookdown::gitbook", encoding="UTF-8")


