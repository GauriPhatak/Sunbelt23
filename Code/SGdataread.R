library(tidyverse)
library(sqldf)
library(furrr)

file <- "//depot.engr.oregonstate.edu/mime_u1/agor/Safe Graph Data/Weekly Patterns/2020_Weekly_Patterns/2020-01-06-weekly-patterns.csv.gz"
df <- read_csv(file = file, n_max  = 100, skip_empty_rows = FALSE)
#df1 <- read_csv(file = file,skip = 100, n_max  = 200)

#read.csv.sql(file, sep = ",", sql = "select count(*) from file")
df_colNames <- colnames(df)
f <- function(x, pos) {
   return(dplyr::filter(x, .data[["region"]] == "NY"))
  #base::subset(x, city == "NY")
}#subset(x, city =="NY")
df2 <- readr::read_csv2_chunked(file = file,
                         DataFrameCallback$new(f),
                         chunk_size = 100001,
                         col_names = df_colNames,
                         #col_types = "cccccccccTTddccccccddcccc",#str_flatten(rep("c",25)),
                         guess_max = 100001,
                         progress = show_progress(), 
                         skip_empty_rows = FALSE)


readr::read_lines_chunked(file= file, DataFrameCallback$new(f), chunk_size = 10000 )


library(chunked)


df2 <- read_chunkwise(file, chunk_size = 10000) %>%
  dplyr::select(city, region, postal_code) %>%
  dplyr::filter(region == "OR")

######## using furr library

cnt <- length(count.fields(file))

batch_size <- 100001

#plan(multiprocess, workers = availableCores() - 1)

seqCnt  <- seq(from = 0, to = cnt+batch_size+1, by = batch_size)

df2 <- as.data.frame(matrix(nrow = 0, ncol = 25))
colnames(df2) <- df_colNames

for(i in seqCnt){
tmp <- read_csv(file, skip = i, 
                       n_max = batch_size, 
                       col_names = df_colNames, 
                       col_types = cols(.default = "c"),
                       guess_max = batch_size) %>%
                filter(region == "OR")
df2 <- rbind(df2, tmp)

}

  # future_map_dfr(
  #   ~ read_csv(
  #     file,
  #     skip      = .x,
  #     n_max     = batch_size, 
  #     # Can't read in col. names in each chunk since they're only present in the
  #     # 1st chunk
  #     col_names = FALSE,
  #     # This reads in each column as character, which is safest but slowest and
  #     # most memory-intensive. If you're sure that each batch will contain
  #     # enough values in each column so that the type detection in each batch
  #     # will come to the same conclusions, then comment this out and leave just
  #     # the guess_max
  #     col_types = cols(.default = "c"),
  #     guess_max = batch_size
  #   ) %>% 
  #     # This is where you'd insert your filter condition(s)
  #     filter(region == "OR"),
  #   # Progress bar! So you know how many chunks you have left to go
  #   .progress = TRUE
  # ) %>% 
  # # The first row will be the header values, so set the column names to equal
  # # that first row, and then drop it
  # set_names(slice(., 1)) %>% 
  # slice(-1)

















