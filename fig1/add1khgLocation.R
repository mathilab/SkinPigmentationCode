## Fill in latitude and longitude for 1000 Genomes samples with population 
#` labels based on a given table with coordinates

add1khgLocation <- function(data.df, loc.df) {
  rownames(loc.df) <- loc.df$Population
  
  for(x in loc.df$Population) {
    data.df$Lat.[which(data.df$Population==x)] <- loc.df[x, "Latitude"]
    data.df$Long.[which(data.df$Population==x)] <- loc.df[x, "Longitude"]
  }
  
  return(data.df)
}