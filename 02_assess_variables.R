#Load libraries
library('data.table')
library('raster')
library('mgcv')

outdir = '/media/dan/variable_selection/v1/'
nfolds = 5
nrounds = 5
#laod files
dat = readRDS(paste0(outdir, 'extracted_variables.rds'))
var_names = readRDS(paste0(outdir, 'extracted_var_names.rds'))

#fill in columns to match by rainy & dry
uniq_vars = gregexpr('-',var_names, fixed = T)
uniq_vars = sapply(uniq_vars, max)
uniq_vars = substr(var_names,1,uniq_vars-1)
uniq_vars = unique(uniq_vars)

for(vvv in uniq_vars){
  dat[value == 'Wet', (vvv) := get(paste0(vvv, '-rainy')) ]
  dat[value == 'Dry', (vvv) := get(paste0(vvv, '-dry')) ]
}
dat[, (var_names) := NULL]

#calculate missingness
miss_assess = dat[, lapply(uniq_vars, function(x) sum(is.na(get(x)))/.N)]
setnames(miss_assess, uniq_vars)
miss_assess = melt(miss_assess, variable.factor = F)
goodies = miss_assess[value<.25, variable]

#create test datasets (SRS & hold out a city or two)
train_vals = sample(1:nrow(dat),size = round(.8 * nrow(dat)))
dat[,id:=.I]
dat[, cases := `PfPR2-10`/100 * Ex]
dat[, not_cases := Ex - cases]
test = dat[!id %in% train_vals, ]
train = dat[train_vals,]

#create kfolds
folds = make_folds(nrow(rhs), nfolds, nrounds)

#create model grid
model_grid = expand.grid(var = goodies, fold = 1:nfolds, round = 1:nrounds, stringsAsFactors = F)
setDT(model_grid)

#for each product and and variable combination, find the best variable
all_mods = sapply(seq(nrow(model_grid)), function(x) fit_model(pr = train,
                                                               iv = model_grid[x,var],
                                                               flagflag = which(folds[,model_grid[x,round]]==model_grid[x,fold]),
                                                               return_rmse = T))
#Compute overall RMSE
model_grid[, rmse := all_mods]

best_prodvar = model_grid[, list(rmse = mean(rmse)), by = 'var']
setnames(best_prodvar, 'var', 'var_name')

best_prodvar[, c('prod', 'var') := tstrsplit(var_name,'-',fixed = T)[3:4]]
unders = sapply(gregexpr("_", best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var], fixed = T), `[[`, 1)
best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var := substr(var, 1, unders-1)]
best_prodvar = best_prodvar[best_prodvar[, .I[which.min(rmse)], by = c('var', 'prod')]$V1]

new_goodies = best_prodvar[, var_name]

#compute interactions
inters = expand.grid(a = new_goodies, b = new_goodies, stringsAsFactors = F)
setDT(inters)

inters[, c('a_prod', 'a_var') := tstrsplit(a,'-',fixed = T)[3:4]]
inters[, c('b_prod', 'b_var') := tstrsplit(b,'-',fixed = T)[3:4]]

#drop rows when the base product and variable are the same
inters = inters[!(a_prod==b_prod & a_var == b_var), ]

#drop two srtm variables
inters = inters[!(b_prod == 'srtm' & a_prod == 'srtm'), ]

#analyze these additional models
imodgrid = expand.grid(id = seq(nrow(inters)), fold = 1:nfolds, round = 1:nrounds)
setDT(imodgrid)
intermods = sapply(seq(nrow(imodgrid)), function(x) fit_model(pr = train,
                                                               iv = unlist(inters[imodgrid[x,id], .(a,b)]),
                                                               flagflag = which(folds[,imodgrid[x,round]]==imodgrid[x,fold]),
                                                               return_rmse = T))



imodgrid[, rmse:= intermods]
inters[, rmse := imodgrid[, list(rmse = mean(rmse)), by = 'id'][,rmse]]

#combine results
res = rbind(best_prodvar[,list(var1 = var_name, rmse)], inters[, .(var1 = a, var2 = b, rmse)], fill = T)

#remove duplicates
res[, id := .I]
res[, v1 := order(unlist(.SD))[1], by = 'id', .SDcols = c('var1','var2')]
res[, v2 := order(unlist(.SD))[2], by = 'id',  .SDcols = c('var1','var2')]

res[v1 == 1, vvv1 := var1]
res[v2 == 1, vvv1 := var2]
res[v1 == 2, vvv2 := var1]
res[v2 == 2, vvv2 := var2]

res = unique(res[,.(vvv1, vvv2, rmse)])
res[, res_id := .I]

#for each of the original prodvars

#take top fraction

#assess OOS fit
