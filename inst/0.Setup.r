#initialize the workspace for outputs

require(bio.lobster)
require(bio.utilities)
require(sf)
require(ggplot2)
require(dplyr)
require(devtools)
require(mgcv)
require(ggeffects)
require(ggforce)

la()

#within a bio_data directory you need to manually create a folder called Framework_LFA35_38

outdir = file.path(project.datadirectory('Framework_LFA35_38'),'outputs')
dir.create(outdir,showWarnings = F)

figdir = file.path(project.datadirectory('Framework_LFA35_38'),'figures')
dir.create(figdir,showWarnings = F)

extra_data = file.path('C:/Users/HowseVJ/Documents/git/bio.lobster.data/')

map_data = file.path(extra_data,'mapping_data')


