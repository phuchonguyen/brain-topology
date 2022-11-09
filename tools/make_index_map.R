# Libraries
library(dplyr)

# Load edge list for tree representation
EDGES <- read.csv("data/Desikan_Cortical_Tree_Edge_List.csv")
NAME <- unique(c(as.character(EDGES$from), as.character(EDGES$to)))

# Load ROI order in adjacency matrix
ROI <- read.table("data/Desikan_Cortical_ROI_Order.txt",
                   fill = TRUE,
                   col.names = c("id", "name", "v1", "v2", "v3", "v4")) %>%
  filter(id > 1000) %>%
  separate(name, into = c("layer", "side", "desikan"), sep = "-") %>%
  mutate(roi = paste(side, desikan, sep = "-")) %>%
  select(roi) %>%
  as.matrix() %>% as.vector()

# Create named list to store ROI and corresponding indices in adj matrix
index_map <- vector("list", length(NAME))
names(index_map) <- NAME

# Define function to recursively combine vectors of indices

set_index_map <- function(node, index_map) {
  children <- get_children(node)
  # Base case
  if (length(children) == 0) {
    index_map[[node]] <- which(ROI == node)
    return(index_map)
  }
  # Recursive case
  indices <- c()
  for (ch in children) {
    index_map <- set_index_map(ch, index_map)
    indices <- c(indices, index_map[[ch]])
  }
  index_map[[node]] <- indices
  return(index_map)
}


set_index_map_iter <- function() {
  index_map <- vector("list", length(NAME))
  names(index_map) <- NAME
  
  leaves <- EDGES %>%
    filter(!to %in% from) %>%
    select(to) %>%
    distinct() %>%
    as.matrix() %>% as.vector()
  
  root <- EDGES %>%
    filter(!from %in% to) %>%
    select(from) %>%
    distinct() %>%
    as.matrix() %>% as.vector()
  
  while(length(leaves) > 0) {
    
    for(node in leaves) {
      children <- get_children(node)
      
      # Base case
      if (length(children) == 0) {
        index_map[[node]] <- which(ROI == node)
      }
      
      # Other case
      indices <- c()
      for (ch in children) {
        indices <- c(indices, index_map[[ch]])
      }
      index_map[[node]] <- indices
    }
    
    leaves <- 
      EDGES %>%
      filter(to %in% leaves) %>%
      select(from) %>%
      distinct() %>%
      as.matrix() %>% as.vector()
    print(length(leaves))
  }
}


# Populate index_map starting from root called "ctx"
index_map <- set_index_map("ctx", index_map)

# Save
saveRDS(index_map, file = "data/Desikan_Cortical_Index_Map.Rds")



