# Author: Lukas Graz             Date: 2023-11-15
# -----------------------------------------------
# R file: moves all .md files from $QUARTO_PROJECT_OUTPUT_FILES to ../_posts/ and adjusts the paths in the .md files
# -----------------------------------------------

# import environment variable
output.files <- Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES")
output.files <- strsplit(output.files, "\n")[[1]]
# only consider the .md files
md.files  <- output.files[grepl(".md$", output.files)]
stopifnot(length(md.files) > 0)


file <- md.files[1]
for (file in md.files){
  fname_without_extension <- gsub("\\..*$", "", basename(file))
  # rename paths in file
  xfun::gsub_files(file, 
    pattern = paste0("\\]\\(", fname_without_extension),
    replacement = paste0("](../quarto/", dirname(file), "/", fname_without_extension)
    )
  # move file

  newfile <- paste0("../_posts/", basename(file))
  file.rename(file, newfile)
  print(paste("moved", file, "to", newfile, " and adjusted paths"))
}
