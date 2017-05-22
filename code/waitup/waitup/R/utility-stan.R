#' Get package's (own) 'src' directory as installed.
src_dir <- function() system.file("src", package="waitup")

#' Get package's (own) directory for Stan files: 
stan_dir <- function() file.path(rappdirs::user_data_dir("waitup", "sakrejda"),'stan')

#' Get package's (own) directory for Stan files: 
cpp_dir <- function() file.path(rappdirs::user_data_dir("waitup", "sakrejda"), 'src')

#' Function for reading in the list of Stan function files.  These are
#' not proper Stan programs as they are missing the 'functions' wrapper
#' as well as the model block.
#' @return vector of strings with paths of Stan function files in the
#' package.
find_stan_functions <- function() {
  function_dir <- system.file('densities', package='waitup')
  function_paths <- list.files(path=function_dir, pattern='\\.stan$',
    full.names=TRUE, recursive=TRUE)
  return(function_paths)
}

#' Function for pulling out declarations from Stan function files.
generate_declarations <- function(s) {
  locs <- gregexpr(
    pattern='\\nreal [a-zA-Z_0-9]+\\([a-zA-Z, ]+\\)',
    text=s)[[1]]
  lens <- attr(locs, 'match.length')
  declarations <- vector(mode='character', length=length(locs))
  for ( i in 1:length(locs)) {
    declarations[i] <- substr(x=s, 
      start=locs[i]+1, stop=locs[i]+lens[i]-1)
  }
  declarations <- paste0('  ', declarations, ';', collapse='\n')
  return(declarations)
}

#' Function for generating a proper dummy Stan functions file with
#' functions block from function files available in the installed
#' package.
merge_stan_functions <- function() {
  function_paths <- find_stan_functions()
  s <- ""
  for ( f in function_paths )
    s <- paste(s, readr::read_file(f), sep = '\n\n')
  header <- 'functions {'
  declarations <- generate_declarations(s)
  footer <- '}\n\nmodel{}\n\n'
  o <- paste(header, declarations, s, footer, sep = '\n\n')
  return(o)
}

#' Merge and transpile Stan functions.
transpile <- function(target_dir = cpp_dir()) { 
  stan_code <- waitup:::merge_stan_functions()
  transpiled <- rstan::stanc(model_code=stan_code,
    model_name='functions', obfuscate_model_name=FALSE)
  if (!dir.exists(target_dir))
    dir.create(target_dir, recursive=TRUE)
  functions_file_path <- file.path(target_dir, 'functions.hpp')
  readr::write_file(transpiled[['cppcode']], path=functions_file_path)
  if (file.exists(functions_file_path)) 
    return(TRUE)
  else
    return(FALSE)
}




