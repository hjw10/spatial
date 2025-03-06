setwd("C://Users//penguin//Desktop//毕设//SPA//C3/spatial")
x <- read.csv("tissue_positions_list.csv",header = T)
x$pxl_row_in_fullres <- as.integer(x$pxl_row_in_fullres)
x$pxl_col_in_fullres <- as.integer(x$pxl_col_in_fullres)
write.csv(x,"new_tissue_position.csv",row.names = FALSE)


setwd("C://Users//penguin//Desktop//毕设//SPA//C3/spatial")
x <- read.csv("tissue_positions_list.csv",header = T)
x$pxl_row_in_fullres <- as.integer(x$pxl_row_in_fullres)
x$pxl_col_in_fullres <- as.integer(x$pxl_col_in_fullres)
write.csv(x,"new_tissue_position.csv",row.names = FALSE)
