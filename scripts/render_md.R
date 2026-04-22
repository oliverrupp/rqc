library(scales)
library(GGally)

template <- "/home/orupp/Projects/rup2/scripts/QCReport.Rmd"
folder <- "/home/orupp/Projects/RanOmics"
plant <- "ND"

for(plant in c("PS", "SP")) {
  rmarkdown::render(
      input = template,
      params = list(results = paste0(folder, "/", plant), min_assigned_read_count = 10000000),
      output_file = paste0(folder, "/", plant, "/", plant,  ".report.html")
  )
}
