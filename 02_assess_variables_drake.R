#drake::r_make('~/Documents/code/umal_varselect/02_assess_variables_drake.R')

library('data.table')
library('raster')
library('mgcv')
library('ggplot2')
library('drake')

#Some parameters to get rolling
indir = '/media/dan/variable_selection/v2/'
output = '/media/dan/variable_selection/v4/'
outdir = output
if (!dir.exists(paste0(output, '.drake')))
  c1 = new_cache(path = paste0(output, '.drake'))

nfolds = 5
nrounds = 5
source('~/Documents/code/umal_varselect/assess_variables_functions.R')
source('~/Documents/code/umal_varselect/assess_variables_drake.R')
#do some prework
#load data and var names
dat = readRDS((paste0(indir, 'extracted_variables.rds')))
var_names = readRDS(paste0(indir, 'extracted_var_names.rds'))
var_names_new = gsub('-', '.', var_names, fixed = T)
setnames(dat, var_names, var_names_new)

#Format the dataset and highlight the variables worth processing
uniq_vars = id_uniq_vars(var_names_new)
set.seed(1)
prepped = prepare_data(dat, uniq_vars, nfolds = nfolds, nrounds = nrounds)
goodvars = identify_goodvars(prepped, uniq_vars)
# pvars = unique(prodvars(goodvars)[, paste0(prod, '.', group_var)])
pvars = prodvars(goodvars)

if(file.exists(paste0(outdir, 'prepped.rds'))){
  ttt = readRDS(paste0(outdir, 'prepped.rds'))
  
  if(all.equal(prepped,ttt) != TRUE){
    saveRDS(prepped, paste0(outdir, 'prepped.rds'))    
  }
}else{
  saveRDS(prepped, paste0(outdir, 'prepped.rds'))
}
fcols = paste0('sfold_', 1:nfolds)
  
#todo: reorganize this so it processes by pvars
#that is, for each pvar, run all the permutations, find the best one, continue
#use the tricks found here: https://github.com/ropensci/drake/issues/697

plan = drake_plan(
  
  #load the data
  pr = readRDS(file_in(paste0(outdir, 'prepped.rds'))),
  
  #fit the first round of models
  plain_mods_rmse = target(gam_cross_val(pr = pr,
                                iv,
                                fold_col = !!fcols),
                      transform = map(.data = !!pvars, .id = c(iv))),
  
  #find_best amongst each prod/group_var pair
  best_version = target(select_pvar(plain_mods_rmse),
                          transform = combine(plain_mods_rmse, .by = c(prod, group_var))),
  
  best_version_vars = target(best_version$iv1, transform = map(best_version)),
  
  #get the model objects for these ones
  best_svals = target(fit_model(pr = pr, best_version_vars, ret_opt = 'svals'),
                       transform = map(best_version_vars)),
  
  #add values to the dataset
  pr_s = target(add_cols_to_pr(pr, best_svals),
                transform = combine(best_svals)),
  
  #compute spline interactions
  spline_oos = target(gam_cross_val(pr = pr, iv = c(best_version$iv1, best_version_vars), fold_col = !!fcols),
                       transform = cross(best_version, best_version_vars)),
  
  #find which instance of splines is the best by variable
  by_so = target(select_pvar(spline_oos), transform = combine(spline_oos, .by = c(prod, group_var))),
  
  #for each prod/group_var pair, find which transform is the best
  best_transform = target(best_trans(by_so, best_version), transform = combine(by_so, best_version, .by = c(prod, group_var))),
  
  #And compute the full model
  
  #create new raster bricks for those
  
  trace = T
  
)

a = drake_config(plan, cache_log_file = '~/Documents/code/umal_varselect/cache_log.csv',
                 cache = drake_cache(paste0(output, '.drake')))

a
