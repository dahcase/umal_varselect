id_uniq_vars = function(var_names, nfolds, nrounds){
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
