library(VennDiagram)
library(grid)
library(scales)

# Function to read .lst files
read_lst <- function(file) {
  lines <- readLines(file)
  return(lines)
}

draw_venn <- function(set1,set2,titleg=TRUE) {
  grid.newpage()  # Clear previous plot
  venn <- draw.pairwise.venn(
    area1 = length(set1),
    area2 = length(set2),
    cross.area = length(intersect(set1,set2)),
    category = c("Result", "Mashmap"),
    fill = c("blue", "green"),
    lty = "solid",
    cex = 1.2,
    cat.cex = 1.2,
    cat.pos = c(0, 180),
    scaled = TRUE,
    margin = 0.1
  )
  if (titleg) {
  # Add a title with the current CTX identifier
  grid.text(paste("Unique for Results: ", length(setdiff(set1,set2)), " and make up of ",percent(length(setdiff(set1,set2))/length(set1)), sep = ""), 
            x = 0.01, y = 0.9, just="left",
            gp = gpar(fontsize = 10))
  grid.text(paste("Unique for Pure Mashmap: ", length(setdiff(set2,set1))," and make up of ",percent(length(setdiff(set2,set1))/length(set2)), sep = ""), 
            x = 0.01, y = 0.88, just="left",
            gp = gpar(fontsize = 10))
  }
}

file_name <- function(tag) {
  pattern <- paste0(tag,".lst")
  return(pattern)
}
# Directory containing the files
path <- "/Users/xiaoqing.zheng/2025-01-14-list_filter"
file_list <- list.files(path, full.names = TRUE)

# Extract unique CTX identifiers
ctx_ids <- unique(gsub(".*(CTX_\\d+).*", "\\1", file_list))  # Extract CTX_xxx from filenames
tags <- c("d","m1","m2","u")
# ctx_ids <- "CTX_148"
# Loop through each CTX identifier

for (ctx in ctx_ids) {
  # Find files for this CTX identifier
  pattern<-file_name("d")
  mashmap_file <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list) & grepl("mm", file_list)]
  um_file <- file_list[grepl("u.lst", file_list) & grepl(ctx, file_list) & grepl("mm", file_list)]
  result_file <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list)][2]
  
  # Check if both files exist
  if (length(mashmap_file) == 1 && length(result_file) == 1) {
    # Read the sets
    mashmap_set <- read_lst(mashmap_file)
    result_set <- read_lst(result_file)
    unmapmash <- read_lst(um_file)
    
    pdf(file = sprintf("%s.pdf",ctx))
    # Draw the Venn diagram
    draw_venn(result_set,mashmap_set)
    grid.text(paste(pattern," Overlap (", ctx, ")", sep = ""), 
              x = 0.5, y = 0.1, 
              gp = gpar(fontsize = 16, fontface = "italic"))
    
    uniq_re<-setdiff(result_set,mashmap_set)
    
    draw_venn(set1=uniq_re,set2=unmapmash,titleg=FALSE)
    grid.text(paste("There is ", percent(length(intersect(uniq_re,unmapmash))/length(uniq_re)), " of decoy match only with intronic decoy overlap with unmapped read in mashmap.",ctx, sep = ""), 
              x = 0.01, y = 0.9, just="left",
              gp = gpar(fontsize = 10))
    grid.text(paste("Unmapped mashmap with Left Decoy Overlap (", ctx, ")", sep = ""), 
              x = 0.5, y = 0.1, 
              gp = gpar(fontsize = 16, fontface = "italic"))
    dev.off()
  } else {
    warning(paste("Missing files for", ctx, "- Skipping..."))
  }
}


pattern<-file_name("m1")
result_filem1 <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list)][2]
mashmap_filem1 <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list) & grepl("mm", file_list)]
pattern<-file_name("m2")
result_filem2 <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list)][2]
mashmap_filem2 <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list) & grepl("mm", file_list)]
pattern<-file_name("u")
result_fileu <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list)][2]
mashmap_fileu <- file_list[grepl(pattern, file_list) & grepl(ctx, file_list) & grepl("mm", file_list)]

result_m1<-read_lst(result_filem1)
result_m2<-read_lst(result_filem2)
result_u <-read_lst(result_fileu)
mash_m1<-read_lst(mashmap_filem1)
mash_m2<-read_lst(mashmap_filem2)
mash_u<-read_lst(mashmap_fileu)

uniq_m <- setdiff(mashmap_set,result_set)
length(intersect(uniq_m,result_m1))
length(result_m1)
length(intersect(uniq_m,result_m2))
length(result_m2)
length(intersect(uniq_m,result_u))
length(result_u)
length(uniq_m)

uniq_re<-setdiff(result_set,mashmap_set)
length(intersect(mash_m1,uniq_re))
length(mash_m1)
length(intersect(mash_m2,uniq_re))
length(mash_m2)
length(intersect(mash_u,uniq_re))
length(mash_u)
length(uniq_re)
