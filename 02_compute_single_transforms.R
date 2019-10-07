#drake::r_make('~/Documents/code/umal_varselect/02_assess_variables_drake.R')

library('data.table')
library('raster')
library('mgcv')
library('ggplot2')
library('drake')
library('tidyr')

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
  
  #The first round of models take every spatial and temporal aggregation and compute the rmse
  #the resulting targets are the rmse from the cross validated results
  plain_mods_rmse = target(gam_cross_val(pr = pr,
                                iv,
                                fold_col = !!fcols),
                      transform = map(.data = !!pvars, .id = c(iv))),
  
  #From the first batch of models (plain_mods_rmse), identify by prod and group var, which has the lowest RMSE
  # singlevar_best
  sv_b = target(select_pvar(plain_mods_rmse),
                          transform = combine(plain_mods_rmse, .by = c(prod, group_var))),
   
  # get the splined values
  # The actual splined values (e.g. the predictions) need to be computed so they can be added to the dataset
  # singlevar_splinevalues
  sv_sv = target(fit_model(pr = pr, sv_b$iv1, ret_opt = 'svals'),
                       transform = map(sv_b)),
  # 
  # #add values to the dataset
  pr_s = target(add_cols_to_pr(pr, sv_sv),
                transform = combine(sv_sv)),
  
  #create a dataset that captures which of the individual variables (sv_b) are being passed along
  selected_sv = target(selected_vars(sv_b), transform = combine(sv_b)),
  
  #save the results
  res = saveRDS(list(pr_s, selected_sv), file = file_out(file.path(output, 'sing_var_splines.rds'))),
  
   trace = T
  
)

a = drake_config(plan, cache_log_file = '~/Documents/code/umal_varselect/cache_log.csv',
                 cache = drake_cache(paste0(output, '.drake')))

a
