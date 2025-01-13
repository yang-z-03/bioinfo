
# search history packages
search <- function(pkg) {
  require(dplyr)
  require(pkgsearch)

  cran_package_history(pkg) |> 
  select(Package, Date, Version, dependencies) |> 
  rowwise() |> 
  mutate(Date = Date,
         R_VERSION = filter(dependencies, package == "R")) |> 
  ungroup() |> 
  arrange(desc(Date))
}

search("MASS")