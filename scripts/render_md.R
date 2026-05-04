template <- "/home/orupp/Projects/rup2/scripts/QCReport.Rmd"
folder <- "/home/orupp/Projects/trqc"
plant <- "CS"

plants <- list.dirs(folder, recursive = F, full.names = F)

for(plant in plants) {
  rmarkdown::render(
        input = template,
        params = list(results = paste0(folder, "/", plant), min_assigned_read_count = 1000),
        output_file = paste0(folder, "/", plant, "/", plant,  ".report.html")
  )
}
