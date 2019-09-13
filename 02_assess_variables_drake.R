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
pvars = prodvars(goodvars)[, paste0(prod,'.',var)]

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
  
  variable = target(vvv, transform = map(vvv = !!goodvars)),
  
  #fit the first round of models
  plain_mods = target(fit_model(pr = pr,
                                iv = variable,
                                fold_col = fcols),
                      transform = cross(variable, fcols = !!fcols)),
  
  #combine those to calculate rmse
  plain_mods_rmse = target(mean(plain_mods), transform = combine(plain_mods, .by = variable)),
  
  #identify which transformation of variable is best
  pv = target(select_pvar(pvs, plain_mods_rmse), transform = map(pvs = !!pvars)) 
  
  #compute all 
  
)

a = drake_config(plan, cache_log_file = '~/Documents/code/umal_varselect/cache_log.csv',
                 cache = drake_cache(paste0(output, '.drake')))

a
