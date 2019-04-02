#Load libraries
library('data.table')
library('raster')
library('mgcv')

outdir = '/media/dan/variable_selection/v1/'
nfolds = 5
nrounds = 5
source('~/Documents/code/umal_varselect/assess_variables_functions.R')

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
train_vals = sample(1:nrow(dat),size = nrow(dat))
dat[,id:=.I]
dat[, cases := `PfPR2-10`/100 * Ex]
dat[, not_cases := Ex - cases]

dat[, YEAR := as.numeric(substr(YEAR, 1, 4))]

#test = dat[!id %in% train_vals, ]
train = dat[train_vals,]

#create kfolds
folds = make_folds(nrow(train), nfolds, nrounds)

#create model grid
model_grid = expand.grid(var = goodies, fold = 1:nfolds, round = 1:nrounds, stringsAsFactors = F)
setDT(model_grid)

#for each product and and variable combination, find the best variable
all_mods = sapply(seq(nrow(model_grid)), function(x) fit_model(pr = train,
                                                               iv = model_grid[x,var],
                                                               flagflag = which(folds[,model_grid[x,round]]==model_grid[x,fold]),
                                                               ret_opt = 'rmse'))
#Compute overall RMSE
model_grid[, rmse := all_mods]

best_prodvar = model_grid[, list(rmse = mean(rmse)), by = 'var']
setnames(best_prodvar, 'var', 'var_name')

best_prodvar[, c('prod', 'var') := tstrsplit(var_name,'-',fixed = T)[3:4]]
unders = sapply(gregexpr("_", best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var], fixed = T), `[[`, 1)
best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var := substr(var, 1, unders-1)]
best_prodvar = best_prodvar[best_prodvar[, .I[which.min(rmse)], by = c('var', 'prod')]$V1]

new_goodies = best_prodvar[, var_name]

#for the "goodies" variables, compute the spline transformation, all insample
trans_vars = lapply(new_goodies, function(x) fit_model(pr = train,
                                                       iv = x,
                                                       ret_opt = 'spline_vals'))

#add to the training dataset
train[, paste0('s_',new_goodies) := trans_vars]

#compute interactions
inters = expand.grid(a = new_goodies, b = new_goodies, stringsAsFactors = F)
setDT(inters)

inters[, c('a_prod', 'a_var') := tstrsplit(a,'-',fixed = T)[3:4]]
inters[, c('b_prod', 'b_var') := tstrsplit(b,'-',fixed = T)[3:4]]

#drop rows when the base product and variable are the same
inters = inters[!(a_prod==b_prod & a_var == b_var), ]

#drop two srtm variables
inters = inters[!(b_prod == 'srtm' & a_prod == 'srtm'), ]

#pass the spline transformed variables
inters[, c('a','b') := list(paste0('s_',a), paste0('s_', b))]

#analyze these additional models
imodgrid = expand.grid(id = seq(nrow(inters)), fold = 1:nfolds, round = 1:nrounds)
setDT(imodgrid)
intermods = sapply(seq(nrow(imodgrid)), function(x) fit_model(pr = train,
                                                               iv = unlist(inters[imodgrid[x,id], .(a,b)]),
                                                               flagflag = which(folds[,imodgrid[x,round]]==imodgrid[x,fold]),
                                                               ret_opt = 'rmse'))



imodgrid[, rmse:= intermods]
inters[, rmse := imodgrid[, list(rmse = mean(rmse)), by = 'id'][,rmse]]

#combine results
res = rbind(best_prodvar[,list(var1 = paste0('s_', var_name), rmse)], inters[, .(var1 = a, var2 = b, rmse)], fill = T)

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
setorder(res, +rmse)
head(res)

#for each unique individual predictor, identify the model that it works best in
res[, v1_best := .I[which.min(rmse)], by = 'vvv1']
res[!is.na(vvv2), v2_best := .I[which.min(rmse)], by = 'vvv2']

#find the best participation by variables
particp = rbind(res[, list(vvv = vvv1, best = v1_best)], res[!is.na(v2_best), list(vvv = vvv2, best = v2_best)])
particp = particp[, list(best = min(best)), by = 'vvv']

#for each of the original prodvars
keepers = res[unique(particp[, best]), ]
keepers[is.na(vvv2), finale := vvv1]
keepers[!is.na(vvv2), finale := paste0(vvv1, '_', vvv2)]

#compute interactions
i_effs = keepers[!is.na(vvv2), ]
for(iii in seq(nrow(i_effs))){
  v1 = i_effs[iii, vvv1]
  v2 = i_effs[iii, vvv2]
  train[, paste0(v1,'_',v2) := get(v1) * get(v2)]
}

#check correlation (they all seem ok)
correlate = cor(na.omit(train[,keepers[,finale],with = F]))
rrr = row.names(correlate)
correlate = data.table(correlate)
correlate[, row_var := rrr]
correlate = melt(correlate, id.vars = 'row_var', variable.factor = F)
correlate = correlate[row_var != variable,]

#compute importance via xgboost
import = expand.grid(fold = 1:nfolds, round = 1:nrounds)

#assess OOS fit
