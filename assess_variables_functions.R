#function to create holdouts
make_folds = function(numrows, numfolds = 8, numsets = 1){
  fold_id = lapply(1:numsets, function(x) sample(cut(1:numrows,breaks=as.numeric(numfolds),labels=FALSE)))
  names(fold_id) = paste0('sfold_', 1:numsets)
  
  return(data.frame(fold_id))
}

#For each variable, fit a GAM of the form cbind(cases, not cases) ~ s(variable) + reff(city) & capture AIC
fit_model = function(pr, iv, flagflag = seq(nrow(pr)), return_rmse = F){
  
  lhs = cbind(pr[,cases], pr[,not_cases])
  rhs = pr[, c('city_name', iv, 'value'), with = F]
  
  if(length(iv)==1){
    setnames(rhs, iv, 'ivar')
  }else if(length(iv) == 2){
    rhs[, ivar:= get(iv[1]) * get(iv[2])]
  }else{
    stop('Too many ivs provided')
  }
  
  rhs[, wetdry := as.integer(value == 'Wet')]
  rhs[, city_name_id := as.numeric(as.factor(city_name))]
  setDF(rhs)
  form = as.formula('lhs~1 + wetdry +  s(ivar) + s(city_name_id, bs = "re")')
  
  #subset by provided indices
  lhs_holdout= lhs[-flagflag,]
  lhs = lhs[flagflag,]
  rhs_holdout = rhs[-flagflag,]
  rhs = rhs[flagflag,]
  
  mod = gam(form, family = 'binomial', data =rhs)
  
  if(return_rmse){
    preds = c(predict(mod, rhs_holdout, type = 'response'))
    truuf = lhs_holdout[,1]/rowSums(lhs_holdout)
    rmse = RMSE(preds, truuf)
    return(rmse)
  }
  
  return(mod)
}

#no idea why I always have to look this up (https://stackoverflow.com/questions/26237688/rmse-root-mean-square-deviation-calculation-in-r)
RMSE = function(m, o){
  sqrt(mean((m - o)^2, na.rm = T))
}
