get_children <- function(node) {
  children <- EDGES %>%
    filter(from == node) %>%
    dplyr::select(to) %>%
    as.matrix() %>% as.vector()
  return(children)
}

get_parent <- function(node) {
  parent <- EDGES %>%
    filter(to == node) %>%
    dplyr::select(from) %>%
    as.matrix() %>% as.vector()
  return(parent)
}