id_uniq_vars = function(var_names){
  uniq_vars = gregexpr('.',var_names, fixed = T)
  uniq_vars = sapply(uniq_vars, max)
  uniq_vars = substr(var_names,1,uniq_vars-1)
  uniq_vars = unique(uniq_vars)
  
  return(uniq_vars)
  
}
prepare_data = function(dat, uniq_vars, nfolds, nrounds){
  
  for(vvv in uniq_vars){
    dat[value == 'Wet', (vvv) := get(paste0(vvv, '.rainy')) ]
    dat[value == 'Dry', (vvv) := get(paste0(vvv, '.dry')) ]
  }
  dat[, (var_names) := NULL]
  
  
  dat[,id:=.I]
  dat[, cases := `PfPR2-10`/100 * Ex]
  dat[, not_cases := Ex - cases]
  
  dat[, YEAR := as.numeric(substr(YEAR, 1, 4))]
  
  #shuffle
  dat = dat[sample(1:nrow(dat), nrow(dat)),]
  
  #make folds
  folds = make_folds(nrow(dat), nfolds, nrounds)
  
  dat = cbind(dat, folds)
  
  return(dat)
}

identify_goodvars = function(dat, uniq_vars){
  
  miss_assess = dat[, lapply(uniq_vars, function(x) sum(is.na(get(x)))/.N)]
  setnames(miss_assess, uniq_vars)
  miss_assess = melt(miss_assess, variable.factor = F)
  goodies = miss_assess[value<.25, variable]
  
  return(goodies)
  
  
}

prodvars = function(goodies){
  
  best_prodvar = data.table(var_name = goodies)
  best_prodvar[, c('prod', 'var') := tstrsplit(var_name,'.',fixed = T)[3:4]]
  unders = sapply(gregexpr("_", best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var], fixed = T), `[[`, 1)
  best_prodvar[prod == 'MCD43A4' & !grepl('Reflectance', var_name, fixed = T), var := substr(var, 1, unders-1)]
  best_prodvar[,id := .I]
  #best_prodvar[, var := toupper(var)]
  
  best_prodvar[, group_var := var]
  best_prodvar[prod == 'srtm', group_var := '']
  
  #handle srtm variables differently
  return(unique(best_prodvar[,.(prod, group_var, iv = var_name, var)]))
  
}

select_pvar = function(...){

  dot_names = names(pryr::named_dots(...))
  dots = list(...)
  
  best_idx = which.min(unlist(dots))
  
  best = data.table(target = dot_names[best_idx], rmse = dots[[best_idx]])
  
  if(grepl('.', dot_names[best_idx], fixed = T)){
    #find the last _ before the first .
    lastunder = gregexpr('_', dot_names[best_idx], fixed = T)[[1]]
    firstdot = gregexpr('.', dot_names[best_idx], fixed = T)[[1]]
    startpos = max(lastunder[lastunder<firstdot[1]]) + 1
    best[, iv1 := substr(target, startpos, nchar(target))]
    best[, iv2 := NA_character_]
  }else{
    new_name = dot_names[best_idx]
    new_name = gsub('spline_oos_best_version_', "", new_name, fixed = T)
    new_name = gsub('_best_version_vars_best_version_', "|", new_name, fixed = T)
    best[, c('iv1', 'iv2') := tstrsplit(new_name, '|', fixed = T)]
  }
  return(best)
  
}

best_trans = function(...){
  res = rbindlist(list(...))
  
  return(res[which.min(rmse),])
}


add_cols_to_pr = function(pr, ...){
  pr = copy(pr)
  
  dot_names = names(pryr::named_dots(...))
  dots = list(...)
  
  
  pr = pr[, (dot_names) := dots]
  
  return(pr)
  
}

