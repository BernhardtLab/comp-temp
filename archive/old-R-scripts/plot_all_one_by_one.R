


# Load tidyverse
library(tidyverse) # to pipe (%>%) and map across each file

source(temp_dependences_MacArthur)
# List files and source each
list.files("R-scripts/one-by-one-plots", full.names = TRUE) %>% map(source)
list.files("R-scripts/one-by-one-plots-same-scale", full.names = TRUE) %>% map(source)
