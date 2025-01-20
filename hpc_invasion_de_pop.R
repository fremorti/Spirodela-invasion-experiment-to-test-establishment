library(tidyverse)
library(brms)
load("est_pop_1.rda")
load("est_pop_2.rda")

est_pop_1 <- est_pop_1%>%mutate(pl = ifelse(ploidy == 'tet', 1, 0))
est_pop_2 <- est_pop_2%>%mutate(pl = ifelse(ploidy == 'tet', 1, 0))


#define model
CytocDE <- "
// Sepcify dynamical system (ODEs)
vector ode_cyto(real t,         // time
              vector y,      // the system rate
              real [] alpha) {   // data constant, not used here
  // the outcome
  vector [2] dydt;
  
  // differential equations
  dydt[1] = alpha[1] * (1-alpha[3] * y[1] - alpha[4] * y[2]) * y[1]; // tetraploids
  dydt[2] = alpha[2] * (1-alpha[5] * y[1] - alpha[6] * y[2]) * y[2]; // diploids
  return dydt;  // return a 2-element array
              }

// Integrate ODEs and prepare output
real cyto(real t, real tet0, real dip0, 
        real rtet, real rdip, 
        real atet, real atetdip,
        real adiptet, real adip,
        real pl) {
        
  vector [2] y0;     // Initial values
  real alpha[6];  // Parameters
  vector[2] y[1];   // ODE solution
  // Set initial values
  y0[1] = tet0; 
  y0[2] = dip0;
  // Set parameters
  alpha[1] = rtet; 
  alpha[2] = rdip;
  alpha[3] = atet; 
  alpha[4] = atetdip;
  alpha[5] = adiptet; 
  alpha[6] = adip;
  // Solve ODEs
  y = ode_rk45_tol(ode_cyto,
                y0, 0, rep_array(t, 1),
                0.001, 0.001, 100, // tolerances, steps
                alpha); 
  // Return relevant population values based on our dummy-variable coding method
  return (y[1][1] * pl + 
          y[1][2] * (1-pl));
}
"
de_bf_2 <- bf(mean ~ log(cyto(week, 190, 10, rtet, rdip, atet, atetdip, adiptet, adip, pl)),
            nlf(rtet ~ exp(logrtet)), logrtet ~ line + line:salt + (1|pop),
            nlf(rdip ~ exp(logrdip)), logrdip ~ line + line:salt + (1|pop),
            nlf(atet ~ exp(logatet)), logatet ~ line + line:salt + (1|pop),
            nlf(atetdip ~ exp(logatetdip)), logatetdip ~ line + line:salt + (1|pop), 
            nlf(adiptet ~ exp(logadiptet)), logadiptet ~ line + line:salt + (1|pop), 
            nlf(adip ~ exp(logadip)), logadip ~ line + line:salt + (1|pop),
            nl = T)

de_p_2 <- c(set_prior("exponential(0.5)", class = "sigma"),
          set_prior("normal(0,1)", class = "b", nlpar = c("logrdip", "logrtet")),
          set_prior("normal(0,0.5)", class = "b", nlpar = c("logadip", "logatet", "logadiptet", "logatetdip")),
          set_prior("normal(-7,0.5)", class = "b", coef = "Intercept", nlpar = c("logadip", "logatet", "logadiptet", "logatetdip")),
          set_prior("cauchy(0,0.1)", class = "sd", nlpar = c("logrdip", "logrtet", "logadip", "logatet", "logadiptet", "logatetdip")))

init_list_2 <- list(sigma = 1,
                  b_logrtet = c(0, rep(0, 7)),
                  b_logrdip = c(0, rep(0, 7)),
                  b_logatet = c(-7, rep(0, 7)),
                  b_logadip = c(-7, rep(0, 7)),
                  b_logatetdip = c(-7, rep(0, 7)),
                  b_logadiptet = c(-7, rep(0, 7)),
                  sd_1 = as.array(0.1),
                  sd_2 = as.array(0.1),
                  sd_3 = as.array(0.1),
                  sd_4 = as.array(0.1),
                  sd_5 = as.array(0.1),
                  sd_6 = as.array(0.1),
                  z_1 = list(c(rep(0, 48))),
                  z_2 = list(c(rep(0, 48))),
                  z_3 = list(c(rep(0, 48))),
                  z_4 = list(c(rep(0, 48))),
                  z_5 = list(c(rep(0, 48))),
                  z_6 = list(c(rep(0, 48))))

start_time <- Sys.time()
de_mod_2 <-brm(data = est_pop_2,
               prior = de_p_2,
               formula = de_bf_2,
               family = brmsfamily("lognormal", link_sigma = "identity"),
               init = list(init_list_2, init_list_2, init_list_2, init_list_2),
               iter = 2000, warmup = 1000, chains = 4, cores = 4,
               stanvars = stanvar(scode = CytocDE, block = "functions"), silent = 0)
end_time <- Sys.time()
print(end_time-start_time)
print(de_mod_2$fit)

save(de_mod_2, file = "de_mod_2.rda")
expose_functions(de_mod_2, vectorize = TRUE)

fit_de_2 <- fitted(de_mod_2)%>%as_tibble()
fit_de_2 <- fit_de_2%>%mutate(est_pop_2%>%filter(!is.na(mean)))
save(fit_de_2, file = "fit_de_2.rda")

