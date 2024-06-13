library(tidyverse)
library(brms)
library(patchwork)
setwd('D:/OneDrive - UGent/Postdoc Terec/polyploidie/Spirodela experimenten/competition experiment2023')
options(mc.cores = 4, max.print = 10000)

#load data tetraploid invasion experiment (exp1)
comp1_ <- read.csv('FCM/comp2022.csv', header = 1)%>%as_tibble()%>%
  mutate(line = factor(line, labels = c("9242", "9346", "9316", "0013")), rep = factor(rep))%>%unite("pop", line:salt, remove = FALSE, sep = "")%>%
  unite("popt", c(pop,week), remove = FALSE)%>%rename(diploids_ = diploids, tetraploids_ = tetraploids)
str(comp1_)

# here I sum the amount of diploid and polyploid nuclei counted for the 2 instances that each population is measured at a certain week
# I'm pretty sure this is not necessary as each measurement can be interpreted as an independent draw from the same probability determined by week, line and salt 
# and affected by the added error. However, I don't think this is a huge problem as I now  consider it as one bigger sample that is measured instead of two smaller
comp1 <- comp1_%>%group_by(popt)%>%mutate(diploids = sum(diploids_), tetraploids = sum(tetraploids_))%>%distinct(popt, .keep_all=TRUE)%>%
  mutate(source = 'experiment', real_tetraploids = NA, real_total = NA, tetrf = tetraploids/(tetraploids+diploids), realf = real_tetraploids/real_total)%>%ungroup()

#explore
ggplot(comp1, aes(x = week, y = tetrf, color = rep, group = rep))+
  geom_point(size = 2)+
  geom_line()+
  #geom_smooth(method = lm)+
  facet_grid(salt~line)

ggplot(comp1, aes(x = week, y = tetrf))+
  geom_point(size = 2)+
  geom_smooth(method = lm)+
  facet_grid(salt~line)

for (line_ in c("9242", "9346", "9316", "0013")){
  print(ggplot(comp1%>%filter(line == line_), aes(x = week, y = tetrf, color = rep))+
          geom_point(size = 2)+
          geom_line()+
          labs(title = line_)+
          #geom_smooth(method = lm)+
          facet_grid(salt~rep))
}


# load calibration data, add NA-columns for the imputation later. However, I add values to salt (as many 1's as 0's) because BRMS cannot impute categorical variables
# more on this later
cal <- read.csv('FCM/cal.csv', header = 1)
data_cal <- as_tibble(cal)%>%rename(diploids = 'nuclei_res.diploids...Count', tetraploids = 'nuclei_res.tetraploids...Count', real_tetraploids = 'rtet', real_total = 'rtot')%>%
  mutate(source = 'calibration', pop = "Cal", iter = NA, week = NA, salt = rep(c(0,1), length.out = n()), rep = NA, popt = NA, line = factor(line, labels = c('9242', '9346', '9316', '0013')),
         tetrf = tetraploids/(tetraploids+diploids), realf = real_tetraploids/real_total)

#let's combine calibration and experimental data
rdata <- rbind(comp1%>%select(diploids, tetraploids, real_tetraploids, real_total, source, pop, iter, week, salt, line, rep, popt, tetrf, realf), 
               data_cal%>%select(diploids, tetraploids, real_tetraploids, real_total, source, pop, iter, week, salt, line, rep, popt, tetrf, realf))%>%
  #and not the di-, tetraploid counts are important, but their relative frequency. I also immediately add transformed proportions to logit scale.
  mutate(logittetrf = logit_scaled(tetrf), logitrealf = logit_scaled(realf), 
         rep = factor(rep), iter = factor(iter), line = factor(line, labels = c('9242', '9346', '9316', '0013')), salt = factor(salt))
rdata_1_nocal <- rdata%>%filter(source == 'experiment')


#load data diploid invasion experiment data (exp2)
comp2 <- read.csv('FCM/comp2023.csv', header = T)%>%as_tibble()%>%mutate(rep = factor(rep))%>%rename(diploids = colnames(.)[2], tetraploids = colnames(.)[3])%>%
  mutate(totaln = diploids+tetraploids, tetrf = tetraploids/totaln)%>%mutate(line = factor(line, labels = c("9242", "9346", "9316", "0013")), source = 'experiment', real_tetraploids = NA, real_total = NA, iter = NA)%>%
  unite("pop", line:salt, remove = FALSE, sep = "")%>%
  unite("popt", c(pop,week), remove = FALSE)
str(comp2)


ggplot(comp2, aes(x = week, y = tetrf, color = rep, group = rep))+
  geom_point()+
  geom_line()+
  #geom_smooth(aes(fill = rep),alpha = 0.2)+
  labs(title = "diploid invasion: per population measurement", y = "tetraploid proportion")+
  facet_grid(salt~line)

ggplot(comp2, aes(x = week, y = tetrf))+
  geom_point(size = 2)+
  geom_smooth(method = lm)+
  facet_grid(salt~line)

for (line_ in c("9242", "9346", "9316", "0013")){
  print(ggplot(comp2%>%filter(line == line_), aes(x = week, y = tetrf, color = rep))+
          geom_point(size = 2)+
          geom_line()+
          labs(title = line_)+
          #geom_smooth(method = lm)+
          facet_grid(salt~rep))
}

#let's combine calibration and experimental data for exp2
rdata_2 <- rbind(comp2%>%select(diploids, tetraploids, real_tetraploids, real_total, source, pop, iter, week, salt, line, rep, popt), 
                 data_cal%>%select(diploids, tetraploids, real_tetraploids, real_total, source, pop, iter, week, salt, line, rep, popt))%>%
  #and not the di-, tetraploid counts are important, but their relative frequency. I also immediately add transformed proportions to logit scale.
  mutate(tetrf = tetraploids/(tetraploids+diploids), realf = real_tetraploids/real_total,
         logittetrf = logit_scaled(tetrf), logitrealf = logit_scaled(realf), 
         rep = factor(rep), iter = factor(iter), line = factor(line, labels = c('9242', '9346', '9316', '0013')), salt = factor(salt))%>%
  filter(week < 15 | is.na(week))
rdata_2_nocal <- rdata_2%>%filter(source == 'experiment')

####################################################################################################
# Model 1: strain specific calibration of observed tetraploid proportions using imputation in exp1 #
####################################################################################################
#this model imputes 'real proportions' from observed proportions in the experiment and the calibration data that have both real and observed proportions. 
#Here we take the strain into account (in calibration and experimental measurements)

mod1_bfobs <- bf(logittetrf ~ samplef + bias,
                 # I'm using a nonlinear syntax for the observation process so it's easier to map and name priors
                 samplef ~ 0 + line:mi(logitrealf), # this is just some trick since brms doesn't allow mi() variables in the main formula of a non-linear model
                 bias ~ 0+line,
                 nl = TRUE
)

mod1_bfexpe <- bf(logitrealf | mi() ~ 1)

get_prior(mvbf(mod1_bfobs, mod1_bfexpe), data = rdata)
mod1_prior <- c(
  set_prior("normal(0,1)", nlpar = c("bias"), class = "b", resp = "logittetrf"),
  set_prior("normal(0,1)", nlpar = c("samplef"), class = "b", resp = "logittetrf"),
  set_prior("normal(0,0.1)", class = "sigma", resp = c("logittetrf"))+
    set_prior("normal(0,1)", class = "sigma", resp = c("logitrealf"))+
    set_prior("normal(0,2)", class = "Intercept", resp = "logitrealf")
)


mod1 <- brm(mvbf(mod1_bfobs, mod1_bfexpe),
            data = rdata%>%filter(!is.infinite(logitrealf)), backend = "cmdstanr",
            prior = mod1_prior, iter = 4000
)
save(mod1, file = "FCM/fits_/mod1.rda")
load("FCM/fits_/mod1.rda")
pp_check(mod1,resp="logittetrf")

#the pp_check on the (logit)realf is off in this coupled model, 
#because the posterior predictions plot the imputed logitrealfs (all from low frequency), while 'y' represents only the logitrealfs that are not imputed (so only calibration)
pp_check(mod1,resp="logitrealf")
mod1$fit

#Add estimated true proportions to the dataset with a likelihood interval to understand calibration effect (estimated - measured proportion)
mod1_post <- as_draws_df(mod1) %>%tibble()%>%
  select(starts_with("Ymi")) %>%
  gather() %>%
  # this extracts the numerals from the otherwise cumbersome names in `key` and saves them as integers
  mutate(var = str_extract(key, "logit[a-z]{5}"), key = str_extract(key, "\\d+")%>%as.integer())%>%mutate(value = inv_logit_scaled(value))
mod1_post_summ <- mod1_post%>% group_by(key) %>%
  summarise(mean = mean(value), q09 = quantile(value, probs = 0.09), q91 = quantile(value, probs = 0.91))
rdata_est <- rdata%>%filter(source == "experiment")%>%mutate(mod1_post_summ, cal = mean - tetrf)



#plot calibration effect in time. Observed proportions always calibrated to lower proportion. 
#Presumably because of tetraploid nuclei in diploids that results in higher observed tetraploid proportions and the presence of debris that add to the diploid and tetraploid peaks and shifts proportions to 0.5.
ggplot(rdata_est, aes(week, cal))+
  geom_point(color = 'blue')+
  facet_grid(salt~line)
#plot calibration effect in function of proportion, very strong relation with the measured proportion. Larger proportions are reduced more, clearly due to the logit rescaling.
ggplot(rdata_est, aes(tetrf, cal))+
  geom_point(color = 'blue')+
  facet_grid(salt~line)

#new (dummy) dataset for expected values
mod1_nd <- rdata%>%filter(source == 'calibration', is.finite(logitrealf))%>%select(c('logitrealf', line))
#calculate expected values from the posterior
pmod1 <- fitted(mod1, newdata = mod1_nd, re_formula = NA, resp = 'logittetrf', probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(mod1_nd, across(c('Estimate', Q9, Q91, logitrealf), inv_logit_scaled), )%>%rename(realf = logitrealf)

#plot calibration data and estimated relation
ggplot(pmod1, aes(x = realf, y = Estimate, color = line, fill = line))+
  geom_line()+
  geom_point(data = rdata%>%filter(source == "calibration"), aes(x = realf, y = tetrf), size = 1.2)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91), alpha = 0.25, color = NA)+
  geom_abline(slope = 1)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  labs(title = 'calibration of tetraploid proportion of nuclei to individuals', x = 'tetraploid proportion of individuals', y = 'tetraploid proportion of nuclei (fcm)', color = 'strain', fill = 'strain')+
  theme(legend.position = 'top', legend.direction="horizontal")


#construct a list of datasets of imputed realfs that all will be fitted to the model that infers effect of treatments
mod1_imppoints <- as_draws_df(mod1) %>%tibble()%>%
  select(starts_with("Ymi"))%>%slice_head(n = 100)
mod1_impdfs <- list()
for (n in 1:nrow(mod1_imppoints)){
  mod1_impset <- mod1_imppoints%>%slice(n)%>%unlist(use.names = F)
  mod1_impdf <- rdata%>%filter(source == 'experiment')%>%mutate(logitrealf = mod1_impset)
  mod1_impdfs[[n]] <- mod1_impdf
}

#Model that infers effect of line, salt and week on the tetraploid proportions. Posterior of calibrated (imputed) proportions will be used to make inference on: 
#cfr 'imputation before model fitting': https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html)
mod1imp_bf <- bf(logitrealf ~ line*salt*week + (1+week|pop))
mod1impf_bf <- bf(logitrealf ~ 0 + Intercept + line:week*salt:week + (0+week|pop))

get_prior(mod1imp_bf, data = mod1_impdfs)
mod1imp_prior <- c(
  set_prior("normal(-3,1)", class = "Intercept"),
  set_prior("normal(0,0.5)", class = "b", coef = c('line9346', 'line9316', 'line0013', 'salt1', 'line9346:salt1', 'line9316:salt1', 'line0013:salt1')),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("cauchy(0.5,0.5)", class = "sigma"),
  set_prior("cauchy(0,0.2)", class = "sd"),
  set_prior("lkj(4)", class = "cor", group = "pop")
)
mod1impf_prior <- c(
  set_prior("constant(-2.94)", class = "b", coef = "Intercept"),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("cauchy(0.5,0.5)", class = "sigma"),
  set_prior("cauchy(0,0.2)", class = "sd")
)

#brm_multiple models the inference model over a list of data sets and pools results
mod1imp <- brm_multiple(mod1impf_bf,
                         data = mod1_impdfs, backend = "cmdstanr",
                         prior = mod1imp_prior)

save(mod1imp, file = "FCM/fits_/mod1imp.rda")
load('FCM/fits_/mod1imp.rda')

mod1impf <- brm_multiple(mod1impf_bf,
                        data = mod1_impdfs, backend = "cmdstanr",
                        prior = mod1impf_prior)
save(mod1impf, file = "FCM/fits_/mod1impf.rda")
load('FCM/fits_/mod1impf.rda')

pp_check(mod1impf)
## NA warning is normal, remember that it can only check against calibration data here, since it's NA in the input for experiment data

#summary(mod1impf, depth = 2)
mod1impf$fit
#check rhats of separate imputed posterior data sets
round(mod1impf$rhats, 2)

mod1imp_nd <- tibble(unique(rdata_1_nocal[c('week', 'line', 'salt')]))
#calculate expected values from the posterior
mod1imp_pmod <- fitted(mod1impf, newdata = mod1imp_nd, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(mod1imp_nd, across(c(Estimate, Q9, Q91), inv_logit_scaled))

RGR_diff <- read.table("Bafort_growth.txt", header = T)%>%tibble()%>%
  filter(concentration %in% c(0, 2.5))%>%group_by(concentration, strain, ploidy)%>%summarize(r = mean(RGR_fn))%>%pivot_wider(names_from = ploidy, values_from = r, names_prefix = "r_")%>%
  mutate(delta_r = r_tet - r_dip)%>%
  mutate(salt = factor(ifelse(concentration == 0, 0, 1)), line = factor(strain, levels = c("13", "9242", "9316", "9346"), labels = c("0013", "9242", "9316", "9346")))
#From Bafort et al. 2023 data, calculate expected population trajectory, first on odds scale, than on proportion. RGR of Bafort is per day.
expg_sim <- RGR_diff%>%expand_grid(week = 1:12)%>%mutate(tet_inv_odds = 5/95 * exp(week*7*delta_r), dip_inv_odds = 95/5 * exp(week*7*delta_r),
                                                         tet_inv = tet_inv_odds/(1+tet_inv_odds), dip_inv = dip_inv_odds/(1+dip_inv_odds))
ggplot(expg_sim, aes(x = week, y = tet_inv, color = salt))+
  geom_line()+
  facet_grid(.~strain)
ggplot(expg_sim, aes(x = week, y = dip_inv_odds, color = salt))+
  geom_line()+
  facet_grid(.~strain)




#posterior predictive expected value distribution
ggplot(mod1imp_pmod, aes(x = week, y = Estimate, class = salt))+
  geom_line(data = expg_sim, aes(x = week, y = tet_inv, color = salt), size = 1, lty = 2)+
  geom_point(data = rdata_est, aes(x = week, y = mean, color = salt), size = 1.4)+
  geom_line(aes(color = salt), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = salt, fill = salt), alpha = 0.25)+
  labs(title = 'measured proportion of polyploids from flow cytometry', y = 'proportion of polyploids', x = 'week')+
  facet_grid(.~line)+
  ylim(c(0, 0.21))+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  theme(legend.position = c(0.9, 0.9), legend.direction="horizontal")


#posterior differences in slopes (salt-control)
pmod1f_dslopes <- as_draws_df(mod1impf, variable = "week:salt1$", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_week:salt1', '9346'= .$'b_week:salt1'+ .$'b_line9346:week:salt1', '9316'= .$'b_week:salt1' + .$'b_line9316:week:salt1', '0013'= .$'b_week:salt1' + .$'b_line0013:week:salt1')%>%
  select(c('9242', '9346', '9316', '0013'))%>%pivot_longer(everything(), names_to = 'line', values_to = 'dslope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')))

ggplot(pmod1f_dslopes, aes(x = line, y = dslope))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, fill = '#1E64C8')+
  geom_hline(yintercept = 0, lty = 2)+
  #scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'differences in slope (salt - control)', y = 'difference in slope (salt - control)', x = 'strain')+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))


#posterior differences in slopes (salt-control); estimated (unfixed) intercept
pmod1_dslopes <- as_draws_df(mod1imp, variable = "salt1:week$", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_salt1:week', '9346'= .$'b_salt1:week'+ .$'b_line9346:salt1:week', '9316'= .$'b_salt1:week' + .$'b_line9316:salt1:week', '0013'= .$'b_salt1:week' + .$'b_line0013:salt1:week')%>%
  select(c('9242', '9346', '9316', '0013'))%>%pivot_longer(everything(), names_to = 'line', values_to = 'dslope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')))

ggplot(pmod1_dslopes, aes(x = line, y = dslope))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, fill = '#1E64C8')+
  geom_hline(yintercept = 0, lty = 2)+
  #scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'differences in slope (salt - control)', y = 'difference in slope (salt - control)', x = 'strain')+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))


#posterior slopes of all 8 groups (line:salt)
pmod1f_slopes <- as_draws_df(mod1impf, variable = "(_|:)week", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_line9242:week', '9346'= .$'b_line9346:week', '9316'= .$'b_line9316:week', '0013'= .$'b_line0013:week')%>%
  mutate('9242:salt' = .$'9242'+.$'b_week:salt1', '9346:salt'= .$'9346'+.$'b_week:salt1'+ .$'b_line9346:week:salt1', '9316:salt'= .$'9316'+.$'b_week:salt1' + .$'b_line9316:week:salt1', '0013:salt'= .$'0013'+.$'b_week:salt1' + .$'b_line0013:week:salt1')%>%
  select(c('9242', '9346', '9316', '0013', '9242:salt', '9346:salt', '9316:salt', '0013:salt'))%>%pivot_longer(everything(), names_to = c('line', 'salt'), names_sep = ':', values_to = 'slope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')), salt = factor(salt =='salt'&!is.na(salt)))

ggplot(pmod1f_slopes, aes(x = line, y = slope, fill = salt))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3)+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'slopes for lines in control or salt treatment', y = 'slope', x = 'strain')+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))

#posterior slopes of all 8 groups (line:salt); estimated (unfixed) intercept
pmod1_slopes <- as_draws_df(mod1imp, variable = "week$", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_week', '9346'= .$'b_week'+ .$'b_line9346:week', '9316'= .$'b_week' + .$'b_line9316:week', '0013'= .$'b_week' + .$'b_line0013:week')%>%
  mutate('9242:salt' = .$'9242'+.$'b_salt1:week', '9346:salt'= .$'9346'+.$'b_salt1:week'+ .$'b_line9346:salt1:week', '9316:salt'= .$'9316'+.$'b_salt1:week' + .$'b_line9316:salt1:week', '0013:salt'= .$'0013'+.$'b_salt1:week' + .$'b_line0013:salt1:week')%>%
  select(c('9242', '9346', '9316', '0013', '9242:salt', '9346:salt', '9316:salt', '0013:salt'))%>%pivot_longer(everything(), names_to = c('line', 'salt'), names_sep = ':', values_to = 'slope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')), salt = factor(salt =='salt'&!is.na(salt)))

ggplot(pmod1_slopes, aes(x = line, y = slope, fill = salt))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3)+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'slopes for lines in control or salt treatment', y = 'slope', x = 'strain')+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))


##################################
# Model 1.2: same model for exp2 #
##################################

mod1_2_bfobs <- bf(logittetrf ~ samplef + bias,
                   # I'm using a nonlinear syntax for the observation process so it's easier to map and name priors
                   # but the idea should work for the
                   samplef ~ 0 + line:mi(logitrealf), # this is just some trick since brms doesn't allow mi() variables in the main formula of a non-linear model
                   bias ~ 0+line,
                   nl = TRUE
)

mod1_2_bfexpe <- bf(logitrealf | mi() ~ 1)

get_prior(mvbf(mod1_2_bfobs, mod1_2_bfexpe), data = rdata_2)
mod1_2_prior <- c(
  set_prior("normal(0,1)", nlpar = c("bias"), class = "b", resp = "logittetrf"),
  set_prior("normal(0,1)", nlpar = c("samplef"), class = "b", resp = "logittetrf"),
  set_prior("cauchy(0,0.1)", class = "sigma", resp = c("logittetrf"))+
  set_prior("cauchy(0,1)", class = "sigma", resp = c("logitrealf"))+
  set_prior("normal(0,2)", class = "Intercept", resp = "logitrealf")
)


mod1_2 <- brm(mvbf(mod1_2_bfobs, mod1_2_bfexpe),
              data = rdata_2%>%filter(!is.infinite(logitrealf)), backend = "cmdstanr",
              prior = mod1_2_prior, iter = 4000, chains = 4, cores = 4,
)
save(mod1_2, file = "FCM/fits_/mod1_2.rda")
load('FCM/fits_/mod1_2.rda')
pp_check(mod1_2,resp="logittetrf")

pp_check(mod1_2,resp="logitrealf") 
summary(mod1_2, depth = 2)
mod1_2$fit

#calibration effect (estimated - measured proportion)
mod1_2_post <- as_draws_df(mod1_2) %>%tibble()%>%
  select(starts_with("Ymi")) %>%
  gather() %>%
  # this extracts the numerals from the otherwise cumbersome names in `key` and saves them as integers
  mutate(var = str_extract(key, "logit[a-z]{5}"), key = str_extract(key, "\\d+")%>%as.integer())%>%mutate(value = inv_logit_scaled(value))
mod1_2_post_summ <- mod1_2_post%>% group_by(key) %>%
  summarise(mean = mean(value), q09 = quantile(value, probs = 0.09), q91 = quantile(value, probs = 0.91))
rdata_2_est <- rdata_2%>%filter(source == "experiment", week<15)%>%mutate(mod1_2_post_summ, cal = mean - tetrf)


#plot calibration effect in time, always calibrated to lower proportion
ggplot(rdata_2_est, aes(week, cal))+
  geom_point(color = 'blue')+
  facet_grid(salt~line)
#plot calibration effect in function of proportion, very strong relation with the measured proportion. Larger proportions are reduced more, clearly due to the logit rescaling.
ggplot(rdata_2_est, aes(tetrf, cal))+
  geom_point(color = 'blue')+
  facet_grid(salt~line)


mod1_2_nd <- rdata_2%>%filter(source == 'calibration', is.finite(logitrealf))%>%select(c('logitrealf', line))
#calculate expected values from the posterior
pmod1_2 <- fitted(mod1_2, newdata = mod1_2_nd, re_formula = NA, resp = 'logittetrf', probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(mod1_2_nd, across(c('Estimate', Q9, Q91, logitrealf), inv_logit_scaled), )%>%rename(realf = logitrealf)

#plot calibration data and estimated relation
ggplot(pmod1_2, aes(x = realf, y = Estimate, color = line))+
  geom_line()+
  geom_point(data = rdata_2%>%filter(source == "calibration"), aes(x = realf, y = tetrf), size = 1.2)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91), alpha = 0.25)+
  geom_abline(slope = 1)+
  labs(title = 'calibration')+
  theme(legend.position = 'top', legend.direction="horizontal")

#list of datasets with posteriorly drawn realfs
mod1_2_imppoints <- as_draws_df(mod1_2) %>%tibble()%>%
  select(starts_with("Ymi"))%>%slice_tail(n = 50)
mod1_2avimp <- mod1_2_imppoints%>%pivot_longer(cols = everything(), names_to = 'key', values_to = 'impf')%>%
  mutate(key = as.numeric(str_extract(key, "\\d+")), impf = inv_logit_scaled(impf))%>%group_by(key)%>%
  summarise(mean = mean(impf), q09 = quantile(impf, probs = 0.09), q91 = quantile(impf, probs = 0.91))
mod1_2_impdfs <- list()
for (n in 1:nrow(mod1_2_imppoints)){
  mod1_2_impset <- mod1_2_imppoints%>%slice(n)%>%unlist(use.names = F)
  mod1_2_impdf <- rdata_2%>%filter(source == 'experiment', week < 15)%>%mutate(logitrealf = mod1_2_impset)
  mod1_2_impdfs[[n]] <- mod1_2_impdf
}
rdata_2_nocalimp <- rdata_2_nocal%>%mutate(mod1_2avimp)

#calibration effect on tetraploid proportion, illustrated for line 0013 in salt
r_corr <- rdata_2_nocalimp%>%filter(line == "0013", salt == 1)%>%select(popt, week, tetrf, mean, q09, q91)%>%mutate(jweek = runif(n(), week-0.2, week+0.2))%>%rename(real = mean, fcm = tetrf)%>%
  pivot_longer(cols = c(fcm, real), names_to = "m", values_to = "ptetr")
ggplot(r_corr, (aes(x = jweek, y = ptetr)))+
  geom_point(aes(shape = m, color = m, group = popt))+
  geom_line(aes(group = popt), size = 1/4)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  scale_color_brewer(palette = "Dark2", name = "", labels = c("fcm (nuclei)", "real (ind)"))+
  scale_shape(name = "", labels = c("fcm (nuclei)", "real (ind)"))+
  labs(title = "correction of fcm measurement", x = 'time (weeks)', y = "tetraploid proportion")+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14))

#inference model
mod1_2imp_bf <- bf(logitrealf ~ line*salt*week + (1+week|pop))
mod1_2impf_bf <- bf(logitrealf ~ 0 + Intercept + line:week*salt:week + (0+week|pop))

get_prior(mod1_2imp_bf, data = mod1_2_impdfs)
mod1_2imp_prior <- c(
  set_prior("normal(3,1)", class = "Intercept"),
  set_prior("normal(0,0.5)", class = "b", coef = c('line9346', 'line9316', 'line0013', 'salt1', 'line9346:salt1', 'line9316:salt1', 'line0013:salt1')),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("cauchy(0.5,0.5)", class = "sigma"),
  set_prior("cauchy(0,0.2)", class = "sd"),
  set_prior("lkj(4)", class = "cor", group = "pop")
)
mod1_2impf_prior <- c(
  set_prior("constant(2.94, 1)", class = "b", coef = "Intercept"),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("cauchy(0.5,0.5)", class = "sigma"),
  set_prior("cauchy(0,0.2)", class = "sd")
)



mod1_2imp <- brm_multiple(mod1_2imp_bf,
                          data = mod1_2_impdfs, backend = "cmdstanr",
                          prior = mod1_2imp_prior, chains = 4, iter = 2000, cores = 10)
save(mod1_2imp, file = "FCM/fits_/mod1_2imp.rda")
load('FCM/fits_/mod1_2imp.rda')

mod1_2impf <- brm_multiple(mod1_2impf_bf,
                          data = mod1_2_impdfs, backend = "cmdstanr",
                          prior = mod1_2impf_prior, chains = 4, iter = 2000, cores = 10)
save(mod1_2impf, file = "FCM/fits_/mod1_2impf.rda")
load('FCM/fits_/mod1_2impf.rda')

pp_check(mod1_2impf)
## NA warning is normal, remember that it can only check against calibration data here, since it's NA in the input for experiment data

summary(mod1_2impf, depth = 2)
mod1_2impf$fit
round(mod1_2impf$rhats, 2)


mod1_2imp_nd <- tibble(unique(rdata_2_nocal[c('week', 'line', 'salt')]))
#calculate expected values from the posterior
mod1_2imp_pmod <- fitted(mod1_2impf, newdata = mod1_2imp_nd, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(mod1_2imp_nd, across(c(Estimate, Q9, Q91), inv_logit_scaled))

#posterior predictive expected value distribution
ggplot(mod1_2imp_pmod%>%filter(week < 15), aes(x = week, y = Estimate, class = salt))+
  geom_line(data = expg_sim, aes(x = week, y = dip_inv, color = salt), size = 1, lty = 2)+
  geom_point(data = rdata_2_est%>%filter(week < 15), aes(x = week, y = mean, color = salt), size = 1.4)+
  geom_line(aes(color = salt), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = salt, fill = salt), alpha = 0.25)+
  labs(title = 'measured proportion of polyploids from flow cytometry', y = 'proportion of polyploids', x = 'week')+
  ylim(c(0, 1))+
  facet_grid(.~line)+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  theme(legend.position = c(0.13, 0.1), legend.direction="horizontal")




#posterior differences in slopes (salt-control)
pmod1_2f_dslopes <- as_draws_df(mod1_2impf, variable = "week:salt1$", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_week:salt1', '9346'= .$'b_week:salt1'+ .$'b_line9346:week:salt1', '9316'= .$'b_week:salt1' + .$'b_line9316:week:salt1', '0013'= .$'b_week:salt1' + .$'b_line0013:week:salt1')%>%
  select(c('9242', '9346', '9316', '0013'))%>%pivot_longer(everything(), names_to = 'line', values_to = 'dslope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')))

ggplot(pmod1_2f_dslopes, aes(x = line, y = dslope))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, fill = '#1E64C8')+
  geom_hline(yintercept = 0, lty = 2)+
  #scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'differences in slope (salt - control)', y = 'difference in slope (salt - control)', x = 'strain')+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))


#posterior differences in slopes (salt-control); estimated (unfixed) intercept
pmod1_2_dslopes <- as_draws_df(mod1_2imp, variable = "salt1:week$", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_salt1:week', '9346'= .$'b_salt1:week'+ .$'b_line9346:salt1:week', '9316'= .$'b_salt1:week' + .$'b_line9316:salt1:week', '0013'= .$'b_salt1:week' + .$'b_line0013:salt1:week')%>%
  select(c('9242', '9346', '9316', '0013'))%>%pivot_longer(everything(), names_to = 'line', values_to = 'dslope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')))

ggplot(pmod1_2_dslopes, aes(x = line, y = dslope))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, fill = '#1E64C8')+
  geom_hline(yintercept = 0, lty = 2)+
  #scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'differences in slope (salt - control)', y = 'difference in slope (salt - control)', x = 'strain')+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))


#posterior slopes of all 8 groups (line:salt)
pmod1_2f_slopes <- as_draws_df(mod1_2impf, variable = "(_|:)week", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_line9242:week', '9346'= .$'b_line9346:week', '9316'= .$'b_line9316:week', '0013'= .$'b_line0013:week')%>%
  mutate('9242:salt' = .$'9242'+.$'b_week:salt1', '9346:salt'= .$'9346'+.$'b_week:salt1'+ .$'b_line9346:week:salt1', '9316:salt'= .$'9316'+.$'b_week:salt1' + .$'b_line9316:week:salt1', '0013:salt'= .$'0013'+.$'b_week:salt1' + .$'b_line0013:week:salt1')%>%
  select(c('9242', '9346', '9316', '0013', '9242:salt', '9346:salt', '9316:salt', '0013:salt'))%>%pivot_longer(everything(), names_to = c('line', 'salt'), names_sep = ':', values_to = 'slope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')), salt = factor(salt =='salt'&!is.na(salt)))

ggplot(pmod1_2f_slopes, aes(x = line, y = slope, fill = salt))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3)+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'slopes for lines in control or salt treatment', y = 'slope', x = 'strain')+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))

#posterior slopes of all 8 groups (line:salt); estimated (unfixed) intercept
pmod1_2_slopes <- as_draws_df(mod1_2imp, variable = "week$", regex = TRUE)%>%as_tibble()%>%
  mutate('9242' = .$'b_week', '9346'= .$'b_week'+ .$'b_line9346:week', '9316'= .$'b_week' + .$'b_line9316:week', '0013'= .$'b_week' + .$'b_line0013:week')%>%
  mutate('9242:salt' = .$'9242'+.$'b_salt1:week', '9346:salt'= .$'9346'+.$'b_salt1:week'+ .$'b_line9346:salt1:week', '9316:salt'= .$'9316'+.$'b_salt1:week' + .$'b_line9316:salt1:week', '0013:salt'= .$'0013'+.$'b_salt1:week' + .$'b_line0013:salt1:week')%>%
  select(c('9242', '9346', '9316', '0013', '9242:salt', '9346:salt', '9316:salt', '0013:salt'))%>%pivot_longer(everything(), names_to = c('line', 'salt'), names_sep = ':', values_to = 'slope')%>%mutate(line = factor(line, levels = c('9242', '9346', '9316', '0013')), salt = factor(salt =='salt'&!is.na(salt)))

ggplot(pmod1_2_slopes, aes(x = line, y = slope, fill = salt))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3)+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'slopes for lines in control or salt treatment', y = 'slope', x = 'strain')+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))






# plot both experiments together
################################
#combine model posterior expected prediction and the data of tetraploid and diploid invasion experiment
mod1_together <- rbind(mod1_2imp_pmod%>%mutate(minority = 'diploid'),
                       mod1imp_pmod%>%mutate(minority = 'tetraploid'))%>%mutate(minority = factor(minority), salt = factor(salt, labels = c('control', 'salt')))
rdata_together <- rbind(rdata_est%>%mutate(minority = 'tetraploid', salt = factor(salt, labels = c('control', 'salt'))),
                        rdata_2_est%>%mutate(minority = 'diploid', salt = factor(salt, labels = c('control', 'salt'))))#%>%mutate(minority = factor(minority))



expg_together <- expg_sim%>%pivot_longer(c(tet_inv, dip_inv), names_to = "minority", values_to = "exp_prop")%>%mutate(minority = factor(minority, levels = c('dip_inv', 'tet_inv'), labels = c('diploid invasion', 'tetraploid invasion')))
ggplot(mod1_together%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = week, y = Estimate, color = salt))+
  geom_line(data = expg_together%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = week, y = exp_prop, color = salt, class = minority), size = 0.8, lty = 2)+
  geom_point(data = rdata_together, aes(x = week, y = mean, shape = minority), size = 1.3)+
  geom_line(aes(color = salt, class = minority), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = salt, fill = salt, class = minority), alpha = 0.25)+
  labs(title = 'measured proportion of polyploids from flow cytometry', y = 'tetraploid proportion', x = 'week', fill = "", color = "")+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  scale_shape_discrete(name = "", labels = c("diploid invasion", "tetraploid invasion"))+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(.~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14))


#on logit scale
ggplot(mod1_together%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = week, y = logit_scaled(Estimate), color = salt))+
  geom_point(data = rdata_together, aes(x = week, y = logit_scaled(mean), shape = minority), size = 1.4)+
  geom_line(aes(color = salt, class = minority), show.legend=FALSE)+
  geom_ribbon(aes(ymin=logit_scaled(Q9), ymax=logit_scaled(Q91), color = salt, fill = salt, class = minority), alpha = 0.25)+
  labs(title = 'measured proportion of polyploids from flow cytometry (logit scale)', y = 'logit of proportion of polyploids', x = 'week')+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  facet_grid(.~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14))


#posterior of slopes
pmod1_slopestogether <- rbind(pmod1_2f_slopes%>%mutate(minority = 'diploid'),
                              pmod1f_slopes%>%mutate(minority = 'tetraploid'))%>%
  mutate(salttext = ifelse(salt==TRUE, 'salt', 'c'))%>%unite(linesalt, c(line, salttext), remove = F)%>%
  mutate(linesalt = ordered(linesalt, levels = c('9242_c', '9242_salt', '9346_c', '9346_salt', '9316_c', '9316_salt', '0013_c', '0013_salt')))
ggplot(pmod1_slopestogether%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = linesalt, y = slope, fill = salt))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, position=position_dodge(width=0), aes(lty = minority))+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_brewer(palette = 'Dark2')+
  scale_linetype_discrete(name = "", labels = c("diploid invasion", "tetraploid invasion"))+
  labs(title = 'slopes for strains in control or salt treatment', y = 'slope', x = 'strain')+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))

# posterior of differences of slopes
pmod1_dslopestogether <- rbind(pmod1_2f_dslopes%>%mutate(minority = 'diploid'),
                               pmod1f_dslopes%>%mutate(minority = 'tetraploid'))%>%mutate(minority = factor(minority))
ggplot(pmod1_dslopestogether, aes(x = line, y = dslope, fill = minority))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3)+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'differences in slope (salt - control)', y = 'difference in slope (salt - control)', x = 'strain')+
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey98', color = c('#1E64C8')))

# plots with predictive expected per replicate
mod1popimp_nd <- tibble(unique(rdata_1_nocal[c('week', 'line', 'salt', 'pop', 'rep')]))
#calculate expected values from the posterior
mod1popimp_pmod <- fitted(mod1impf, newdata = mod1popimp_nd, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(mod1popimp_nd, across(c(Estimate, Q9, Q91), inv_logit_scaled))

#posterior predictive expected value distribution
(ppp1 <- ggplot(mod1popimp_pmod%>%filter(week < 15)%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = week, y = Estimate, color = salt, class = rep))+
  geom_point(data = rdata_est%>%filter(week < 15)%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = week, y = mean, color = salt), size = 1.4)+
  geom_line(aes(color = salt), show.legend=FALSE)+
  #geom_ribbon(aes(ymin=Q9, ymax=Q91, color = salt, fill = salt), alpha = 0.25)+
  labs(y = 'proportion of polyploids', x = 'week')+
  ylim(c(0, NA))+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(.~line)+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14)))

mod1_2popimp_nd <- tibble(unique(rdata_2_nocal[c('week', 'line', 'salt', 'pop', 'rep')]))
#calculate expected values from the posterior
mod1_2popimp_pmod <- fitted(mod1_2impf, newdata = mod1_2popimp_nd, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(mod1_2popimp_nd, across(c(Estimate, Q9, Q91), inv_logit_scaled))

mod1pop_together <- rbind(mod1_2popimp_pmod%>%mutate(minority = 'diploid'),
                          mod1popimp_pmod%>%mutate(minority = 'tetraploid'))%>%mutate(minority = factor(minority), salt = factor(salt, labels = c('control', 'salt')))

ggplot(mod1pop_together%>%mutate(salt = factor(salt, labels = c('control', 'salt'))), aes(x = week, y = Estimate, color = salt))+
  #geom_line(data = expg_together, aes(x = week, y = exp_prop, color = salt, class = minority), size = 1, lty = 2)+
  geom_point(data = rdata_together, aes(x = week, y = mean, shape = minority), size = 1.3)+
  geom_line(aes(color = salt, class = minority, class = rep), show.legend=FALSE)+
  labs(title = 'measured proportion of polyploids from flow cytometry', y = 'tetraploid proportion', x = 'week', fill = "", color = "")+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  scale_shape_discrete(name = "", labels = c("diploid invasion", "tetraploid invasion"))+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(.~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14))

#####################################################################################################################
# dry weight exp 1 in logistic model                                                                                  #
# model 2.1: calibrate measured dry weight by substracting mean weight and accounting for the added sd by envelopes #
#####################################################################################################################

w_data1_ <- read.csv('FCM/dry weight.csv', header = 1)%>%as_tibble()%>%
  mutate(rep = factor(rep), salt = factor(salt), line = factor(line, labels = c('9242', '9346', '9316', '0013')), size = factor(size), weight = weight*1000)%>%
  unite("pop", line:salt, remove = FALSE, sep = "")

#data on empty envelopes, big and small
envelopes <- read.csv('FCM/envelopes.csv', header = 1)%>%as_tibble()%>%mutate(size = factor(size), weight = weight*1000)
#mean and sd of big and small envelopes
env_sum <- envelopes%>%group_by(size)%>%summarise(mean = mean(weight), sd = sd(weight))

#the proportion of surface area that was sampled each weak
surface_prop = (3.5*0.5)^2*pi/14/9

w_data1 <- w_data1_%>%mutate(envmean = ifelse(size == 'big', env_sum$mean[1], env_sum$mean[2]), envsd = ifelse(size == 'big', env_sum$sd[1], env_sum$sd[2]))%>%
  mutate(dw = weight-envmean, dw_tot = dw/surface_prop, weight_tot = weight/surface_prop)


bf_wK_1 <- bf(dw | mi(sdy = envsd) ~ L/(1+exp(-k*(week-x0))), 
              ##NDmaxime: for other families, such as poisson, it may be best to use an identity link rather than log()-ing the nonlinear formula
              ## but the lognormal family doesn't have a link function independent of the distribution itself
              ##
              nlf(L~exp(logL)),
              nlf(k~exp(logk)),
              logL~ line*salt + (1|pop), 
              logk ~line*salt  + (1|pop), 
              x0 ~ line*salt + (1|pop),
              nl = T)
get_prior(bf_wK_1, data = w_data1)
mod_wK_1 <-brm(data = w_data1, family = gaussian,
              #make use of measurement error capibilities of brms to correct for the variance introduced by 
              formula = bf_wK_1,
              backend = "cmdstanr",
              prior = c(
                set_prior("cauchy(0,2)", class = "sigma"),
                set_prior("normal(0,0.5)", class = "b", nlpar = c("logk"))
                ,set_prior("normal(0,1)", class = "b", nlpar = c("logk"), coef = "Intercept")
                ,set_prior("normal(0,0.5)", class = "b", nlpar = c("logL"))
                ,set_prior("normal(5,1)", class = "b", nlpar = c("logL"), coef = "Intercept")
                ,set_prior("normal(0,0.5)", class = "b", nlpar = c("x0"))
                ,set_prior("normal(0,2)", class = "b", nlpar = c("x0"), coef = "Intercept")
                ,set_prior("cauchy(0,0.1)", class = "sd", nlpar = c("logk"))
                ,set_prior("cauchy(0,0.1)", class = "sd", nlpar = c("x0"))
                ,set_prior("cauchy(0,1)", class = "sd", nlpar = c("logL"))),
              iter = 4000, warmup = 1000,
              save_pars = save_pars(latent = TRUE))

save(mod_wK_1, file = "FCM/fits_/mod_wK_1.rda")
load("FCM/fits_/mod_wK_1.rda")

mod_wK_1$fit
pp_check(mod_wK_1)

ndmod_wK_1 <- tibble(week = rep(1:12, each = 8), salt = factor(rep(0:1, 12, each = 4)), line = rep(unique(rdata$line), 24))%>%
  mutate(envsd = ifelse(week %in% c(1,2,3,4), env_sum$sd[1], env_sum$sd[2]))
#calculate expected values from the posterior
pmod_wK_1 <- fitted(mod_wK_1, newdata = ndmod_wK_1, re_formula = dw ~ 0 + Intercept + line*salt*week, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndmod_wK_1, Estimate = Estimate/surface_prop, Q9 = Q9/surface_prop, Q91 = Q91/surface_prop,salt = factor(salt, levels = c(0, 1), labels = c('control', 'salt')))

(w1 <- ggplot(pmod_wK_1, aes(x = week, y = Estimate, class = salt))+
  geom_point(data = w_data1%>%mutate(salt = factor(salt, levels = c(0, 1), labels = c('control', 'salt'))), aes(x = week, y = dw_tot, color = salt), size = 1.4)+
  geom_line(aes(color = salt), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = salt, fill = salt), alpha = 0.25)+
  labs(title = 'dry weight of the population in tetraploid invasion experiment', y = 'dry weight (mg)', x = 'week', fill = '', color = '')+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(.~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14)))
  
#with group effects (pop)
ndmod_wK_1_pop <- tibble(unique(rdata_1_nocal%>%select(line, rep, salt, week, pop)))%>%mutate(envsd = ifelse(week %in% c(1,2,3,4), env_sum$sd[1], env_sum$sd[2]))
#calculate expected values from the posterior
pmod_wK_1_pop <- fitted(mod_wK_1, newdata = ndmod_wK_1_pop, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndmod_wK_1_pop, Estimate = Estimate/surface_prop, Q9 = Q9/surface_prop, Q91 = Q91/surface_prop)

ggplot(pmod_wK_1_pop, aes(x = week, y = Estimate, class = rep))+
  geom_point(data = w_data1, aes(x = week, y = dw_tot, color = rep), size = 1.4)+
  geom_line(aes(color = rep), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = rep, fill = rep), alpha = 0.25)+
  labs(title = 'dry weight of the population in tetraploid invasion experiment', y = 'dry weight (mg)', x = 'week')+
  facet_grid(salt~line)+
  theme(legend.position = "bottom", legend.direction="horizontal")




#####################################################################################################################
# dry weight exp2 with logistic model                                                                                                #
# model 2.2: calibrate measured dry weight by substracting mean weight and accounting for the added sd by envelopes #
#####################################################################################################################

w_data2 <- read.csv('FCM/dry weight2.csv', header = 1)%>%as_tibble()%>%
  mutate(rep = factor(rep), salt = factor(salt), line = factor(line, labels = c('9242', '9346', '9316', '0013')), size = factor(size), weight = weight*1000)%>%
  unite("pop", line:salt, remove = FALSE, sep = "")
str(w_data2)

w_data2 <- w_data2%>%mutate(envmean = ifelse(size == 'big', env_sum$mean[1], env_sum$mean[2]), envsd = ifelse(size == 'big', env_sum$sd[1], env_sum$sd[2]))%>%
  mutate(dw = weight-envmean, dw_tot = dw/surface_prop, weight_tot = weight/surface_prop)


get_prior(bf_wK_1, data = w_data2)
mod_wK_2 <-brm(data = w_data2, family = gaussian,
               #make use of measurement error capibilities of brms to correct for the variance introduced by 
               formula = bf_wK_1,
               backend = "cmdstanr",
               prior = c(
                 set_prior("cauchy(0,2)", class = "sigma"),
                 set_prior("normal(0,0.5)", class = "b", nlpar = c("logk"))
                 ,set_prior("normal(0,0.5)", class = "b", nlpar = c("logL"))
                 ,set_prior("normal(0,0.5)", class = "b", nlpar = c("x0"))
                 ,set_prior("normal(0,1)", class = "b", nlpar = c("logk"), coef = "Intercept")
                 ,set_prior("normal(4,1)", class = "b", nlpar = c("logL"), coef = "Intercept")
                 ,set_prior("normal(0,2)", class = "b", nlpar = c("x0"), coef = "Intercept")
                 ,set_prior("cauchy(0,0.1)", class = "sd", nlpar = c("logk"))
                 ,set_prior("cauchy(0,0.1)", class = "sd", nlpar = c("x0"))
                 ,set_prior("cauchy(0,1)", class = "sd", nlpar = c("logL"))),
               iter = 4000, warmup = 1000,
               save_pars = save_pars(latent = TRUE))

save(mod_wK_2, file = "FCM/fits_/mod_wK_2.rda")
load("FCM/fits_/mod_wK_2.rda")

mod_wK_2$fit
pp_check(mod_wK_2)

ndmod_wK_2 <- tibble(week = rep(1:12, each = 8), salt = factor(rep(0:1, 12, each = 4)), line = rep(unique(rdata$line), 24))%>%mutate(envsd = ifelse(week %in% c(1,2,3,4), env_sum$sd[1], env_sum$sd[2]))
#calculate expected values from the posterior
pmod_wK_2 <- fitted(mod_wK_2, newdata = ndmod_wK_2, re_formula = dw ~ 0 + Intercept + line*salt*week, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndmod_wK_2, Estimate = Estimate/surface_prop, Q9 = Q9/surface_prop, Q91 = Q91/surface_prop, salt = factor(salt, levels = c(0, 1), labels = c('control', 'salt')))

(w2 <- ggplot(pmod_wK_2, aes(x = week, y = Estimate, class = salt))+
  geom_point(data = w_data2%>%mutate(salt = factor(salt, levels = c(0, 1), labels = c('control', 'salt'))), aes(x = week, y = dw_tot, color = salt), size = 1.4)+
  geom_line(aes(color = salt), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = salt, fill = salt), alpha = 0.25)+
  labs(title = 'dry weight of the population in diploid invasion experiment', y = 'dry weight (mg)', x = 'week', fill = '', color = '')+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(.~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14)))

w1 <- w1+labs(x = "")
(w1 / w2)+ plot_layout(guides = "collect", axis_title = 'collect') & theme(legend.position = "bottom", legend.direction="horizontal")

#with group effects (pop)
ndmod_wK_2_pop <- tibble(unique(rdata_1_nocal%>%select(line, rep, salt, week, pop)))%>%mutate(envsd = ifelse(week %in% c(1,2,3,4), env_sum$sd[1], env_sum$sd[2]))
#calculate expected values from the posterior
pmod_wK_2_pop <- fitted(mod_wK_2, newdata = ndmod_wK_2_pop, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndmod_wK_2_pop, Estimate = Estimate/surface_prop, Q9 = Q9/surface_prop, Q91 = Q91/surface_prop)

ggplot(pmod_wK_2_pop, aes(x = week, y = Estimate, class = rep))+
  geom_point(data = w_data2, aes(x = week, y = dw_tot, color = rep), size = 1.4)+
  geom_line(aes(color = rep), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = rep, fill = rep), alpha = 0.25)+
  labs(title = 'dry weight of the population in diploid invasion experiment', y = 'dry weight (mg)', x = 'week')+
  facet_grid(salt~line)+
  theme(legend.position = "bottom", legend.direction="horizontal")


#######################################
# Model 3: Count to Weight conversion #
#######################################

cw_data <- read.csv('FCM/countweight.csv', header = T)%>%as_tibble()%>%mutate(rep = factor(rep), ploidy = factor(ploidy, labels = c('2n', '4n')))%>%rename(dw = dry.weight)%>%
  mutate(dwmg = dw*1000, line = factor(line))
str(cw_data)

ggplot(cw_data, aes(x = count, y = dwmg))+
  geom_point(aes(color = ploidy))+
  geom_smooth(method = "lm", aes(color = ploidy))+
  facet_grid(.~line)


bf_cw <- bf(dwmg ~ 0+ count*line*ploidy)

get_prior(bf_cw, data = cw_data)
p_cw <- c(
  set_prior("normal(0,1)", class = "b", coef = c('count', 'count:lineB', 'count:lineC', 'count:lineD', 'count:ploidy4n', 'count:lineB:ploidy4n', 'count:lineC:ploidy4n', 'count:lineD:ploidy4n')),
  set_prior("constant(0)", class = "b", coef = c('lineA', 'lineB', 'lineC', 'lineD', 'ploidy4n', 'lineB:ploidy4n', 'lineC:ploidy4n', 'lineD:ploidy4n')),
  set_prior("cauchy(0,1)", class = "sigma"))


mod_cw <- brm(bf_cw,
               data = cw_data, backend = "cmdstanr",
               prior = p_cw,
               iter = 4000, chains = 4
)
save(mod_cw, file = "FCM/fits/mod_cw.rda")
load('FCM/fits/mod_cw.rda')

pp_check(mod_cw)
## NA warning is normal, remember that it can only check against calibration data here, since it's NA in the input for experiment data

summary(mod_cw)
mod_cw$fit

nd_cw <- tibble(unique(cw_data[c('line', 'ploidy', 'count')]))%>%rbind(unique(cw_data[c('line', 'ploidy')])%>%mutate(count = 0))
#calculate expected values from the posterior
pmod_cw <- fitted(mod_cw, newdata = nd_cw, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(nd_cw)

ggplot(pmod_cw%>%mutate(line = factor(line, labels = c("9242", "9346", "9316", "0013"))), aes(x = count, y = Estimate, class = ploidy))+
  geom_point(data = cw_data%>%mutate(line = factor(line, labels = c("9242", "9346", "9316", "0013"))), aes(x = count, y = dwmg, color = ploidy), size = 1.4)+
  geom_line(aes(color = ploidy), show.legend=FALSE)+
  geom_ribbon(aes(ymin=Q9, ymax=Q91, color = ploidy, fill = ploidy), alpha = 0.25)+
  labs(title = 'duckweed weight', y = 'weight (mg)', x = 'number of individuals')+
  facet_grid(.~line)

#posterior draws of weight conversion
cs_long <- as_draws_df(mod_cw)%>%as_tibble()%>%mutate(A2n = b_count, 
                                                       B2n = A2n + .$"b_count:lineB", 
                                                       C2n = A2n + .$"b_count:lineC", 
                                                       D2n = A2n + .$"b_count:lineD", 
                                                       A4n = A2n + .$"b_count:ploidy4n", 
                                                       B4n = A4n + .$"b_count:lineB" + .$"b_count:lineB:ploidy4n", 
                                                       C4n = A4n + .$"b_count:lineC" + .$"b_count:lineC:ploidy4n", 
                                                       D4n = A4n + .$"b_count:lineD" + .$"b_count:lineD:ploidy4n")%>%
  select(A2n, B2n, C2n, D2n, A4n, B4n, C4n, D4n)

#mean weight conversion
cs <- cs_long%>%
  pivot_longer(everything(), names_to = 'linep', values_to = "c")%>%
  group_by(linep)%>%summarise(mean_c = mean(c), c_q09 = quantile(c, prob = 0.09), c_q91 = quantile(c, prob = 0.91))%>%
  mutate(line = c('A', 'A', 'B', 'B', 'C', 'C', 'D', 'D'), ploidy = rep(c('C2n', 'C4n'),4))%>%
  pivot_wider(c(line), names_from = ploidy, values_from = c(mean_c))%>%
  mutate(line = factor(line, labels = c("9242", "9346", "9316", "0013")))
cs

################################################################################
# calculate population size in exp1 from distribution of polyploid proportion, #
# distribution of total dry weight and                                         #
# distribution of count to weight conversion for each line                     #
################################################################################
nsamps <- 50

mod2_1est_prop <- fitted(mod1impf, newdata = ndmod_wK_1, summary = FALSE, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()
mod2_1est_prop <- mod2_1est_prop%>%mutate(across(everything(), inv_logit_scaled))
w4nprops_1 <- slice_sample(mod2_1est_prop, n = nsamps)
w4nprop_1 <- w4nprops_1%>%pivot_longer(everything(), names_to = 'datapoint', values_to = 'p')%>%mutate(ndmod_wK_1%>%slice(rep(1:n(), nsamps)), q = 1-p)


cs_longs <- cs_long%>%slice_sample(n = nsamps)%>%t()%>%as_tibble()%>%
  mutate(line = rep(c('9242', '9346', '9316', '0013'), 2), ploidy = factor(rep(c('2n', '4n'), each = 4)))%>%
  group_by(line)%>%group_split()
names(cs_longs) <- sort(c('9242', '9346', '9316', '0013'))

#weighted mean dry mass per individual (weighted according to cytotype proportion), drawn from 50 posterior cs valuesets over 50 posterior draws of tetraploid proportion (over 96 data points that vary time, ploidy and strain)
csss_1 <- tibble()
for(i in 1:nrow(w4nprop_1)){
  row <- as.matrix(slice(w4nprop_1, i)%>%select(c('q', 'p')), nrow = 1)%*%as.matrix(cs_longs[[as.character(slice(w4nprop_1, i)$line)]][, 1:nsamps])
  csss_1 <- bind_rows(csss_1, row%>%as_tibble())#
}
totalw_1 <- fitted(mod_wK_1, newdata = ndmod_wK_1, re_formula = dw ~ 0 + Intercept + line*salt*week, summary = F)%>%as_tibble()%>%slice_sample(n = nsamps)%>%
  pivot_longer(everything(), names_to = 'datapoint', values_to = 'totalw')%>%mutate(totalw = totalw/surface_prop)
totaln_1 <- bind_cols(w4nprop_1, csss_1)%>%pivot_longer(starts_with('V'), names_to = 'cschain', values_to = 'cstotal')%>%full_join(totalw_1, by = "datapoint")%>%
  mutate(total = totalw/cstotal, tetraploids = total*p, diploids = total*q)%>%pivot_longer(c('tetraploids', 'diploids'), names_to = 'cytotype', values_to = 'pop_size')
tot_1 <- group_by(totaln_1, week, line, salt, cytotype)%>%summarise(mean = mean(pop_size), q09 = quantile(pop_size, probs = 0.09), q91 = quantile(pop_size, probs = 0.91))%>%mutate(salt = factor(salt, levels = c('0', '1'), labels = c('control', 'salt')))

(ptot1 <- ggplot(tot_1, aes(x = week, y = mean, color = cytotype, fill = cytotype))+
  geom_line(show.legend=FALSE)+
  geom_ribbon(aes(ymin=q09, ymax=q91), alpha = 0.25)+
  labs(title = 'estimated population size in tetraploid invasion experiment', y = 'population size (ind)', x = 'week')+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(salt~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14)))
(ptot1_2 <- ggplot(tot_1%>%filter(cytotype == 'tetraploids'), aes(x = week, y = mean, color = salt, fill = salt))+
  geom_line(show.legend=FALSE)+
  geom_ribbon(aes(ymin=q09, ymax=q91), alpha = 0.25)+
  labs(title = 'zoom on estimated tetraploid population size', y = 'population size (ind)', x = 'week')+
  scale_fill_brewer(palette = 'Dark2')+
  scale_color_brewer(palette = 'Dark2')+
  facet_grid(.~line)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14)))

#per population
ndprop_1 <- unique(rdata%>%filter(source == 'experiment')%>%select(c('salt', 'line', 'rep', 'pop')))%>%slice(rep(1:n(), 12))%>%mutate(week = rep(1:12, each = n()/12))
mp2_1 <- fitted(mod1impf, newdata = ndprop_1, summary = FALSE, probs = c(0.09, 0.91), ndraws = nsamps)%>%
  as_tibble()%>%mutate(across(everything(), inv_logit_scaled))

#####################################
# calculate population size in exp2 #
#####################################
#averaged per treatment

mod2_2est_prop <- fitted(mod1_2impf, newdata = ndmod_wK_2, summary = FALSE, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()
mod2_2est_prop <- mod2_2est_prop%>%mutate(across(everything(), inv_logit_scaled))
w4nprops_2 <- slice_sample(mod2_2est_prop, n = nsamps)
w4nprop_2 <- w4nprops_2%>%pivot_longer(everything(), names_to = 'datapoint', values_to = 'p')%>%mutate(ndmod_wK_2%>%slice(rep(1:n(), nsamps)), q = 1-p)


#weighted mean dry mass per individual (weighted according to cytotype proportion), drawn from 50 posterior cs valuesets over 50 posterior draws of tetraploid proportion (over 96 data points that vary time, ploidy and strain)
csss_2 <- tibble()
for(i in 1:nrow(w4nprop_1)){
  row <- as.matrix(slice(w4nprop_1, i)%>%select(c('q', 'p')), nrow = 1)%*%as.matrix(cs_longs[[as.character(slice(w4nprop_1, i)$line)]][, 1:nsamps])
  csss_2 <- bind_rows(csss_2, row%>%as_tibble())#
}
totalw_2 <- fitted(mod_wK_2, newdata = ndmod_wK_2, re_formula = NA, summary = F)%>%as_tibble()%>%slice_sample(n = nsamps)%>%
  pivot_longer(everything(), names_to = 'datapoint', values_to = 'totalw')%>%mutate(totalw = totalw/surface_prop)
totaln_2 <- bind_cols(w4nprop_2, csss_2)%>%pivot_longer(starts_with('V'), names_to = 'cschain', values_to = 'cstotal')%>%full_join(totalw_2, by = "datapoint")%>%
  mutate(total = totalw/cstotal, tetraploids = total*p, diploids = total*q)%>%pivot_longer(c('tetraploids', 'diploids'), names_to = 'cytotype', values_to = 'pop_size')
tot_2 <- group_by(totaln_2, week, line, salt, cytotype)%>%summarise(mean = mean(pop_size), q09 = quantile(pop_size, probs = 0.09), q91 = quantile(pop_size, probs = 0.91))%>%mutate(salt = factor(salt, levels = c('0', '1'), labels = c('control', 'salt')))

(ptot2 <- ggplot(tot_2, aes(x = week, y = mean, color = cytotype, fill = cytotype))+
  geom_line(show.legend=FALSE)+
  geom_ribbon(aes(ymin=q09, ymax=q91), alpha = 0.25)+
  labs(title = 'estimated population size in diploid invasion experiment', y = 'population size (ind)', x = 'week')+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  facet_grid(salt~line)+
  theme(legend.position = "bottom", legend.direction="horizontal",
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14)))
ptot1_ <- ptot1+labs(x ="")
(ptot1_/ptot2)+plot_layout(guides = 'collect', axis_title = 'collect')& theme(legend.position = "bottom", legend.direction="horizontal")

##### adapted figure from grainger
library(wesanderson)
priority_line = tibble(SD = 1:3200*(-0.0005))%>%mutate(FD_e = 1-SD)
coex_line = tibble(SD = 1:1400*(0.0005))%>%mutate(FD_e = 1/(1-SD))%>%rbind(priority_line)%>%mutate(neg = ifelse(SD <0, T, F))
dots = tibble(type = factor(c(1,1,1,3,3,3)), SD = c(-1.3, -1.2, -1.1, -0.8, -0.6, -0.4), FD_e = c(2.1, 2.3, 2.5, 2, 1.8, 1.6), grad = c(0.3, 0.45, 0.7, 0.3, 0.45, 0.7))

ggplot(coex_line, aes(x = SD, y = FD_e))+
  geom_ribbon(aes(ymin = 1, ymax = FD_e, fill = neg), alpha = 0.7)+
  geom_ribbon(aes(ymin = FD_e, ymax = 3.34), fill = wes_palette("AsteroidCity3")[3], alpha = 0.7)+
  geom_line()+
  geom_hline(yintercept = 1, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_text(x = -0.2, y = 2.8, size = 5, label = 'competitive exclusion')+
  geom_text(x = 0.43, y = 1.15, size = 5, label = 'coexistence')+
  geom_text(x = -1, y = 1.15, size = 5, label = 'priority effects')+
  labs(x = "Stabilizing differences", y = "Fitness differences", color = "  environment or genotype")+
  geom_point(data = dots, aes(x = SD, y = FD_e, color = grad, shape = type), size = 3)+
  scale_shape(guide = F)+
  scale_fill_manual(guide = F, values = wes_palette("AsteroidCity3"))+
  guides(color = guide_colourbar(theme = theme(legend.key.width  = unit(10, "lines"),legend.key.height = unit(0.5, "lines"), legend.text = element_blank(), legend.title.position = "bottom")))+
  theme(legend.position = c(0.26, 0.9), legend.direction = "horizontal",legend.background = element_rect(fill=NA), legend.title = element_text(colour="white", face = "bold"),
        panel.background = element_rect(fill = "white", colour = '#1E64C8'),
        strip.background = element_rect(fill = "#1E64C8"), strip.text = element_text(colour = 'white', size = 14))
  