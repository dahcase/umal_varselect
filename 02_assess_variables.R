#Load libraries
library('data.table')
library('raster')
library('mgcv')
library('ggplot2')

indir = '/media/dan/variable_selection/v2/'
outdir = '/media/dan/variable_selection/v4/'
dir.create(outdir)
nfolds = 5
nrounds = 5
source('~/Documents/code/umal_varselect/assess_variables_functions.R')

#laod files
dat = readRDS(paste0(indir, 'extracted_variables.rds'))
var_names = readRDS(paste0(indir, 'extracted_var_names.rds'))

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
                                                               flagflag = which(folds[,model_grid[x,round]] != model_grid[x,fold]),
                                                               ret_opt = 'rmse'))
#Compute overall RMSE
model_grid[, rmse := all_mods]

best_prodvar = model_grid[, list(rmse = mean(rmse)), by = 'var']
setnames(best_prodvar, 'var', 'var_name')

best_prodvar[, c('prod', 'var') := tstrsplit(var_name,'-',fixed = T)[3:4]]
unders = sapply(gregexpr("_", best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var], fixed = T), `[[`, 1)
best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var := substr(var, 1, unders-1)]
best_prodvar = best_prodvar[best_prodvar[, .I[which.min(rmse)], by = c('var', 'prod')]$V1]

best_prodvar[,id := .I]
best_prodvar[, var := toupper(var)]

#keep best srtm
best_prodvar = best_prodvar[prod != 'srtm' | id == best_prodvar[prod == 'srtm',][which.min(rmse), id],]

#keep best version of indices
best_prodvar = best_prodvar[best_prodvar[, .I[which.min(rmse)], by = c('var')]$V1]

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
intermods = parallel::mclapply(seq(nrow(imodgrid)), function(x) fit_model(pr = train,
                                                               iv = unlist(inters[imodgrid[x,id], .(a,b)]),
                                                               flagflag = which(folds[,imodgrid[x,round]]!=imodgrid[x,fold]),
                                                               ret_opt = 'rmse'), mc.cores = 4)

intermods = unlist(intermods)
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

#get the spline transformed values
trans_inter_vars = lapply(1:nrow(i_effs), function(x) fit_model(pr = train,
                                                       iv = c(i_effs[x,vvv1], i_effs[x,vvv2]),
                                                       ret_opt = 'spline_vals'))

train[, paste0('s_', i_effs[, finale]) := trans_inter_vars]

#adjsut finale
keepers[!is.na(vvv2), finale:= paste0('s_', finale)]

#check correlation
correlate = cor(na.omit(train[,keepers[,finale],with = F]))
rrr = row.names(correlate)
correlate = data.table(correlate)
correlate[, row_var := rrr]
correlate = melt(correlate, id.vars = 'row_var', variable.factor = F)
correlate = correlate[row_var != variable,]

#remove over .8
bad_cor = correlate[abs(value)> .8]
bad_cor[, id:= .I]
vnames = c('row_var','variable')
if(nrow(bad_cor)>0){
  bad_cor[, v1 := get(vnames[order(unlist(.SD))[1]]), by = 'id', .SDcols = vnames]
  bad_cor[, v2 := get(vnames[order(unlist(.SD))[2]]), by = 'id',  .SDcols = vnames]
  bad_cor = unique(bad_cor[, .(v1, v2)])
  
  #do it alphabetically
  rrr = rrr[!rrr %in% bad_cor[,1]]
}

#compute importance via xgboost
ig = expand.grid(fold = 1:nfolds, round = 1:nrounds)

import_res = lapply(1:100, function(x) fit_xgb(pr = train, iv = rrr,
                                                             flagflag = sample(1:nrow(train),size = nrow(train), replace = T), #bootstraps
                                                             ret_opt = 'importance'))

avg_import_res = rbindlist(import_res)
avg_import_res = avg_import_res[,list(Gain = mean(Gain), Cover = mean(Cover), Frequency = mean(Frequency),
                                      lower_Gain = quantile(Gain, .025), upper_Gain = quantile(Gain, .975)), by = 'Feature']
setorder(avg_import_res, -Gain)

# RMSE_res = lapply(1:nrow(ig), function(x) fit_xgb(pr = train, iv = rrr,
#                                                     flagflag = which(folds[, ig[x, 'round']] != ig[x, 'fold']), #id the hold out
#                                                     ret_opt = 'rmse'))
# avg_rmse_xg = mean(unlist(RMSE_res))

full_xg = fit_xgb(pr = train, iv = rrr, flagflag = 1:nrow(train), ret_opt = 'model')

full_xg = xgboost::xgb.load(full_xg)

#graphics for the presentation

#Graph 1: pfpr data over time and city
g1 = ggplot(train, aes(x = YEAR, y = cases /(cases + not_cases), size = (cases + not_cases))) + geom_point() + facet_wrap(~city_name) +
  theme_bw() + scale_size_continuous('N', guide = guide_legend()) + ggtitle('PfPR 2-10 Data by City') +
  xlab('') + ylab('PfPR 2 - 10') + 
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 23, angle = 90),
        axis.text.y = element_text(size = 23),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24))
plot(g1)

ggsave(paste0(outdir, 'pfpr_data_plot.png'), g1, width = 15, height = 12, units = 'in', dpi = 600)


#graph 2 & 3: Scatter plots of variables before and after transform
vartotest = keepers[finale %in% rrr & is.na(vvv2), vvv1 ][1]

train[, v1 := get(substr(vartotest, 3, 100))]
train[, v2 :=  get(vartotest)]
g2 = ggplot(train, aes(x = (v1), y = (v2), color = value)) + geom_point() + theme_bw() + xlab('Original') + ylab('Transformed') +
  ggtitle(substr(vartotest, 3, 100)) + scale_x_continuous(limits = quantile(train[,v1], c(.025,.975), na.rm = T)[1:2]) +
  scale_y_continuous(limits = quantile(train[,v2], c(.025,.975), na.rm = T)) +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        plot.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))
plot(g2)
ggsave(paste0(outdir, 'single_transform.png'), g2, width = 15, height = 12, units = 'in', dpi = 600)


intertotest = keepers[finale %in% rrr & !is.na(vvv2), finale ][1]

train[, v3 := get(substr(intertotest, 3, 100))]
train[, v4 :=  get(intertotest)]
g3 = ggplot(train, aes(x = (v3), y = (v4), color = value)) + geom_point() + theme_bw() + xlab('Pre-transform') + ylab('Post Transform') +
  ggtitle(intertotest) + scale_x_continuous(limits = quantile(train[,v3], c(.025,.975), na.rm = T)[1:2]) +
  scale_y_continuous(limits = quantile(train[,v4], c(.025,.975), na.rm = T)) +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        plot.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))
plot(g3)
ggsave(paste0(outdir, 'inter_transform.png'), g3, width = 15, height = 12, units = 'in', dpi = 600)

#variable importance
vl = rrr #c('EVI', 'NDVI', 'Band7 * EVI', 'NDVI * NDWI', 'Band1', 'NDWI * LST Night',
#        'Band2 * Band7', 'NDVI * Band4', 'NDVI * NTL', 'Band3', 'EVI * Band5', 'NDVI * Aspect', 'EVI * LST Day')
var_cw = data.table(Feature = paste0('ivar', 1:length(rrr)), var_name = rrr,
                                      var_label = vl)
avg_import_res = merge(avg_import_res, var_cw, by = 'Feature', all.x = T)
avg_import_res[is.na(var_label), var_label := Feature]
avg_import_res[, var_label_fact := reorder(var_label, Gain)]

g4 = ggplot(avg_import_res[!var_label %in% unique(dat$city_name), ], aes(x = Gain, y = var_label_fact)) + geom_errorbarh(aes(xmin = lower_Gain, xmax = upper_Gain)) + geom_point() +
  theme_bw() + xlab('Importance/Gain') + ylab('') + ggtitle('Variable Importance via BRT') +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 21),
        plot.title = element_text(size = 30),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
plot(g4)
ggsave(paste0(outdir, 'var_importance.png'), g4, width = 15, height = 12, units = 'in', dpi = 600)

#RMSE ranking
rmse_table = merge(var_cw, keepers, by.x = 'var_name', by.y ='finale', all.x = T)
rmse_table = rmse_table[, .(var_label, rmse)]
write.csv(rmse_table, paste0(outdir, 'var_rmse.csv'), row.names = F)

#transformed variables
#load dakar
ras = brick("/media/dan/processed///summary//MCD43A4/rainydry/max/Bamako_MCD43A4_006_ndwi_nirswi_2_8_2001_2016_rainydry_max.tif")
mod = fit_model(pr = train, iv = 'max-rainydry-MCD43A4-ndwi_nirswi_2_8', ret_opt = 'model' )

pres = as.data.frame(ras, xy = T)
setDT(pres)
pres[,pixel_id := .I]
setnames(pres, c('x', 'y', 'dry', 'rainy', 'id'))

pres = melt(pres, id.vars = c('id', 'x', 'y'), variable.factor = F, value.name = 'ivar')
pres[, wetdry:=as.numeric(variable == 'rainy')]
pres[, city_name_id := 4]

pres[, mod_res := predict(mod, newdata = pres)]
pres[variable== 'dry', variable := 'Dry Season']
pres[variable == 'rainy', variable := 'Rainy Season']
pres[, NDWI := scale(ivar)]
pres[, Season := variable]

minfill = min(pres[,scale(ivar)], pres[, scale(mod_res)])
maxfill = max(pres[,scale(ivar)], pres[, scale(mod_res)])

pre = ggplot(pres, aes(x = x, y = y, fill = NDWI)) + geom_raster() + facet_wrap(~variable ) + theme_bw() + 
  ggtitle('Bamako | NDWI | Pre-transformed') + scale_fill_continuous(limits = c(minfill, maxfill)) +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 22),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +coord_fixed()
plot(pre)
ggsave(paste0(outdir, 'ndwi_pre.png'), pre, width = 15, height = 12, units = 'in', dpi = 600)

pres[,NDWI := scale(mod_res)]
post = ggplot(pres, aes(x = x, y = y, fill = NDWI)) + geom_raster() + facet_wrap(~variable ) + theme_bw() + 
  ggtitle('Bamako | NDWI | Post-transformed') + scale_fill_continuous(limits = c(minfill, maxfill)) + 
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 22),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + coord_fixed()
plot(post)
ggsave(paste0(outdir, 'ndwi_post.png'), post, width = 15, height = 12, units = 'in', dpi = 600)

compare = ggplot(pres, aes(x = ivar, y = mod_res, group = Season, color = Season)) + geom_line(size = 3) + theme_bw() + 
  ggtitle('Bamako: NDWI') + xlab('NDWI') + ylab('Transformed NDWI') + 
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        plot.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
plot(compare)
ggsave(paste0(outdir, 'ndwi_compare.png'), compare, width = 15, height = 12, units = 'in', dpi = 600)
