
credentials_file <- function() return('~/credentials/sphhs-bioepi01.json')

connect <- function(credentials) {
  credentials_data <- readr::read_file(credentials) %>% jsonlite::fromJSON()
  conn <- do.call(what=dbConnect, args=c(list('PostgreSQL'), credentials_data))
  return(conn)
}


