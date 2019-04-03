library('data.table')
library('sf')
library("raster")
library('openxlsx')

#load the data
outdir = '/media/dan/variable_selection/v2/'
dir.create(outdir)
pr = read.xlsx('/media/dan/REACT SSA Cities PR data (240717)_Final.xlsx')
setDT(pr)
city_shape = st_read("/home/dan/Documents/react_data/Cities_React/study_areas.shp")

#load rasters to query
sg = readRDS('/media/dan/summary_grid.rds')

#Fix city names
pr[, city_name := Name.of.City]
pr[city_name == "N'Djamena", city_name := 'NDjamena']
pr[city_name == 'Mbuji-Mayi', city_name := 'MbujiMayi']
pr[city_name == 'Dar es salaam', city_name := 'DarEsSalaam']

#Assign rainy/dry
rainy = fread('/media/dan/rainy_seasons.csv')
rainy = melt(rainy, id.vars = 'city', variable.factor = F)
rainy[, month_abb := tolower(variable)]
mth_num = data.table(month_abb = tolower(month.abb), month_num = 1:12 )
rainy = merge(rainy, mth_num, by = 'month_abb')
setnames(rainy, 'city', 'city_name')

pr = merge(pr, rainy, by.x = c('city_name', 'MM'), by.y = c('city_name', 'month_num'))

#Add rasters to the dataset
rasvars = sg[time == 'rainydry' | product == 'srtm',]
uniq_vars = unique(rasvars[, .(funk, time, product, variables)])

assign_ras_var = function(params, pr){
  
  params = params[, .(funk, time, product,variables)]
  
  prcit = copy(pr)
  #for each city in the PR data
  for(cit in unique(prcit$city_name)){
    ras = brick(rasvars[funk == params[, funk] & time == params[,time] & product == params[, product] & variables == params[, variables] & city == cit, raspath])
    rasvals = raster::extract(ras, prcit[city_name == cit, .(Long, Lat)])
    
    if(params[,time] == ""){
      rasvals = cbind(rasvals,rasvals)
    }
    
    prcit = prcit[city_name == cit, (paste0(paste(params, collapse = '-'), c('-dry', '-rainy'))) := list(rasvals[,1], rasvals[,2])]
  }
  
  return(prcit[,  paste0(paste(params, collapse = '-'), c('-dry', '-rainy')), with = F])
}

guides = unique(rasvars[,.(funk,time,product,variables)])
res = mclapply(1:nrow(guides), function(x) assign_ras_var(guides[x, ], pr), mc.cores = 5)

results = do.call(cbind, res)

#generate custom indices

#save adjusted dataset
saveRDS(cbind(pr, results), paste0(outdir, 'extracted_variables.rds'))
saveRDS(names(results), paste0(outdir, 'extracted_var_names.rds'))


