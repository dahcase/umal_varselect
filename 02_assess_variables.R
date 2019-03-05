#Load libraries

#fill in columns to match by rainy & dry

#create test datasets (SRS & hold out a city or two)

#For each variable, fit a GAM of the form cbind(cases, not cases) ~ s(variable) + reff(city) & capture AIC

#for the top X % of models, take the interaction between the variables represented

#analyze these additional models

#resort based on aic

#take top fraction

#assess OOS fit
