#function to create holdouts
make_folds = function(numrows, numfolds = 8, numsets = 1){
  fold_id = lapply(1:numsets, function(x) sample(cut(1:numrows,breaks=as.numeric(numfolds),labels=FALSE)))
  names(fold_id) = paste0('sfold_', 1:numsets)
  
  return(data.frame(fold_id))
}

#For each variable, fit a GAM of the form cbind(cases, not cases) ~ s(variable) + reff(city) & capture AIC
fit_model = function(pr, iv, fold_col = "", ret_opt = 'rmse'){
  #print(iv)
  
  pr = copy(pr)
  
  if(fold_col == ''){
    pr[, foldy := -1]
  }else{
    pr[, foldy := .SD, .SDcols = fold_col]
  }
  
  lhs = cbind(pr[,cases], pr[,not_cases])
  rhs = pr[, c('city_name', iv, 'value', 'YEAR', 'foldy'), with = F]
  
  if(length(iv)==1){
    setnames(rhs, iv, 'ivar')
  }else if(length(iv) == 2){
    rhs[, ivar:= get(iv[1]) * get(iv[2])]
  }else{
    stop('Too many ivs provided')
  }
  
  rhs[, wetdry := as.integer(value == 'Wet')]
  rhs[, city_name_id := as.numeric(as.factor(city_name))]

  flaggys = lapply(unique(rhs[, foldy]), function(x) which(rhs[,foldy] == x))
  setDF(rhs)
  
  #fit gams
  gams = lapply(flaggys, function(x) gam_fitty(lhs, rhs, x, ret_opt))
  
  if(ret_opt == 'rmse'){
    return(mean(unlist(gams)))
  }else{
    return(gams[[1]])
  }
}

gam_fitty = function(lhs, rhs, flagflag, ret_opt = c('model', 'rmse')){
  
  # #subset by provided indices
  lhs_holdout= lhs[-flagflag,]
  lhs = lhs[flagflag,]
  rhs_holdout = rhs[-flagflag,]
  rhs = rhs[flagflag,]
  form = as.formula('lhs~1 + s(ivar, wetdry, k = 3) + s(city_name_id, bs = "re")')

  mod = suppressWarnings(gam(form, family = 'binomial', data =rhs, method = 'REML'))
  
  if(ret_opt == 'rmse'){
      preds = c(predict(mod, rhs_holdout, type = 'response'))
      truuf = lhs_holdout[,1]/rowSums(lhs_holdout)
      rmse = RMSE(preds, truuf)
      return(rmse)
    } else if(ret_opt == 'model'){
      return(mod)
    } else{
  
      svals = predict(mod, rhs, type = 'terms')
      svals = svals[, which(colnames(svals) == 's(ivar,wetdry)')]
  
      return(svals)
  
    }
  
}


gam_cross_val = function(pr, iv, fold_cols, group_var = "", prod = ""){
  mean(vapply(fold_cols, function(x) fit_model(pr = pr, iv = iv, fold_col = x, ret_opt = 'rmse'), 1))
}

#no idea why I always have to look this up (https://stackoverflow.com/questions/26237688/rmse-root-mean-square-deviation-calculation-in-r)
RMSE = function(m, o){
  sqrt(mean((m - o)^2, na.rm = T))
}

fit_xgb = function(pr, iv, flagflag = seq(nrow(pr)), ret_opt = c('model', 'rmse', 'preds', 'importance')){
  
  #prep dataset
  pr = (pr[, c('city_name', iv, 'value', 'YEAR', 'cases', 'not_cases'), with = F])
  pr[, N := cases + not_cases] 
  pr[, wetdry := as.integer(value == 'Wet')]
  pr[, city_name_id := as.factor(city_name)]
  
  setnames(pr, iv, paste0('ivar', 1:length(iv)))
  ivars = paste0('ivar', 1:length(iv))
  
  cities = unique(pr[, city_name])
  for(ccc in cities[2:length(cities)]){
    pr[, (ccc) := as.numeric(ccc == city_name)]
  }
  
  ptors = c(ivars, 'wetdry',cities[2:length(cities)])
  
  #subset by provided indices
  rhs_holdout = pr[-flagflag, ptors, with = F]
  rhs = pr[flagflag, c(ptors, 'N'), with = F]
  lhs = pr[flagflag, cases/N]
  lhs_holdout = pr[-flagflag, cases/N]
  
  #drop NAs from training
  rhs[,id := .I]
  r_id = rhs[,id]
  rhs = na.omit(rhs)
  drop_r = setdiff(r_id, rhs[,id])
  lhs = lhs[-drop_r]
  
  #drop NAs from test
  rhs_holdout[, id:= .I]
  rh_id = rhs_holdout[, id]
  rhs_holdout = na.omit(rhs_holdout)
  drop_rh = setdiff(rh_id, rhs_holdout[,id])
  lhs_holdout = lhs_holdout[-drop_rh]
  
  dm = xgboost::xgb.DMatrix(data = as.matrix(rhs[, ptors, with = F]),
                            label = lhs,
                            weight = rhs[, N])
  
  mod = xgboost::xgboost(data = as.matrix(rhs[, ptors, with = F]),
                         label = lhs,
                         weight = rhs[, N],
                         params = list(objective = 'reg:logistic', nthread = 4, max_depth = 4, eta = .01,
                                       subsample = .6, colsample_bytree = .75),
                         nrounds = 3000, verbose = F)
  
  #insamp = predict(mod, newdata = xgboost::xgb.DMatrix(data = as.matrix(rhs[,ptors, with = F])))
  #plot(lhs, insamp)
  
  if(ret_opt == 'model'){
    return(xgboost::xgb.save.raw(mod))
  }else if(ret_opt == 'rmse'){
    preds = predict(mod, newdata = xgboost::xgb.DMatrix(data = as.matrix(rhs_holdout[,ptors, with = F])))
    rmse = RMSE(preds, lhs_holdout)
  }else if(ret_opt == 'preds'){
    preds = predict(mod, newdata = xgboost::xgb.DMatrix(data = as.matrix(rhs[,ptors, with = F])))
  }else{
    return(xgboost::xgb.importance(model = mod))
  }
}


