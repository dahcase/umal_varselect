library('data.table')
library('raster')
library('mgcv')
library('ggplot2')
library('drake')
library('tidyr')

#Some parameters to get rolling
indir = '/media/dan/variable_selection/v2/'
output = '/media/dan/variable_selection/v4/'

dc = file.path(output, 'interactions', '.drake')
if (!dir.exists(dc))
  c2 = new_cache(path = dc)

source('~/Documents/code/umal_varselect/assess_variables_functions.R')
source('~/Documents/code/umal_varselect/assess_variables_drake.R')

#Load the results from the single variable stage
res_s1 = readRDS(file.path(output,'sing_var_splines.rds'))
dat = res_s1[[1]]
vars = res_s1[[2]]

#fix some of the column namings
sv_vars = grep('sv_b', names(dat), value = T)
new_names = gsub('sv_sv_', "", sv_vars, fixed = T)
setnames(dat, sv_vars, new_names)

fcols = grep('sfold', names(dat), value = T)
isects = CJ(iv1 = new_names, iv2 = new_names)
isects = isects[iv1 != iv2, ]

setDF(isects)

plan = drake_plan(
  i_s = target(gam_cross_val(pr = pr,
                             iv = c(iv1, iv2), 
                             fold_col = !!fcols),
                             transform = map(.data = !!isects))
  
)


# # Now that all of the best performing single splined variables have been chosen, pass those on and compute all the intersections
# # interaction_splines
# i_s = target(gam_cross_val(pr = pr, iv = c(sv_b$iv1, sv_b_c), fold_col = !!fcols),
#                      transform = cross(sv_b, sv_b_c)),
# 
# i_s_map = target(tidyr::crossing(sv_b$iv1, sv_b_c), transform = map(sv_b, sv_b_c)),
# 
# # Find the best performing spline iteraction by prod and group
# # b_is best_interactionspline
# b_is = target(select_pvar(i_s), transform = combine(i_s, .by = c(prod, group_var))),
# 
# # 
# # for each prod/group_var pair, find which transform is the best
# # bestoveralltransform
# bot = target(best_trans(b_is, sv_b), transform = combine(b_is, sv_b, .by = c(prod, group_var))),
# 
# get the splined values from the interaction spline of the relevant ones
# interactionspline_splinevalues
# is_sv = target(fit_model(pr = pr, b_is, ret_opt = 'svals'),
#                transform = map(b_is)),
# 
# #And compute the full model
# model_select = target(fit_model(pr = pr_s, iv = c(best_transform$iv1, best_transform$iv2), ret_opt = 'model'), transform = map(best_transform)),
# 
# #create new raster bricks for those
# 