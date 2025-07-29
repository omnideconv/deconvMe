#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name deconvMestartup
NULL

.onLoad <- function(libname, pkgname){
  cli::cli_alert("checking deconvMe environment and dependencies")
  
  # We ensure to have reticulate
  if (!dir.exists(reticulate::miniconda_path())) {
    message("Setting python version in miniconda to be 3.8")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.8)
    message("Setting up miniconda environment..")
    suppressMessages(reticulate::install_miniconda())
  }
  
  # We ensure to have the r-deconvMe env
  if (!("r-deconvMe" %in% reticulate::conda_list()$name)) {
    message("Create conda evironment 'r-deconvMe' for MethAtlas...")
    reticulate::conda_create("r-deconvMe", python_version = "3.8")
    message("Install all python dependencies...")
    reticulate::py_install(packages = c("numpy", "pandas", "scipy", "matplotlib") , envname = "r-deconvMe",  method = "conda", pip = T)
  }
  
  paths <- reticulate::conda_list()
  path <- paths[paths$name == "r-deconvMe", 2][[1]]
  
  Sys.setenv(RETICULATE_PYTHON = path)
  reticulate::use_miniconda(condaenv = "r-deconvMe", required = TRUE)
  reticulate::py_config()
  reticulate::configure_environment(pkgname, force = TRUE)
}