# Ariel Boyarsky
# aboyarsky@uchicago.edu
#
# The following code replicates the bounds indicated at the top of Figure 4
# in Mogstad Torgovitsky 2018. We then implement a series of specifications
# producing new bounds.
#
# The procedure allows for partial identification of IV parameters. 
#

# clear enviorment
rm(list = ls())

# make current dir wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load functions
source("MT_FUnctions.R")

# install relevant libraries
library(gurobi)
library(Matrix)

# Save results
res <- data.frame()

#-------------------------------------------------------------------------
# calculate beta_s for tslslike, ivlike, and waldlike estimands
#-------------------------------------------------------------------------

beta_s_tsls <- calc_beta_s(s_tsls, exp_m_true)
beta_s_iv <- calc_beta_s(s_iv, exp_m_true)
beta_s_wald <- calc_beta_s(s_z2_z3, exp_m_true)

#-------------------------------------------------------------------------
# Part 0: Calculate Paper Boounds
#-------------------------------------------------------------------------

gamma_cnstr_iv <- calc_gamma(s_iv,4,bern_4)
gamma_cnstr_tsls <- calc_gamma(s_tsls,4,bern_4)
gamma_star <- calc_gamma_star(4,bern_4)

# Now optimize using Gorubi
A <- matrix(c(c(unlist((gamma_cnstr_iv))),
                    c(unlist((gamma_cnstr_tsls)))),
                  nrow=2,
                  byrow=TRUE
)
rhs <- c(beta_s_iv,beta_s_tsls)
sense <- c('=','=')
obj <- c(unlist((gamma_star)))

r <- lp(A, rhs, sense, obj)

res <- rbind(res, data.frame(r$min, r$max))

#-------------------------------------------------------------------------
# Part 1: Calc Bounds for Polynomial of Degree 2
#-------------------------------------------------------------------------
gamma_cnstr_iv <- calc_gamma(s_iv,2,bern_2)
gamma_cnstr_tsls <- calc_gamma(s_tsls,2,bern_2)
gamma_star <- calc_gamma_star(2,bern_2)

# Now optimize
A <- matrix(c(c(unlist((gamma_cnstr_iv))),
                    c(unlist((gamma_cnstr_tsls)))),
                  nrow=2,
                  byrow=TRUE
)
rhs <- c(beta_s_iv,beta_s_tsls)
sense <- c('=','=')
obj <- c(unlist((gamma_star)))


r <- lp(A, rhs, sense, obj)

res <- rbind(res, data.frame(r$min, r$max))


#-------------------------------------------------------------------------
# Part 2: Calc Bounds for Polynomial of Degree 6
#-------------------------------------------------------------------------

# get constraints
gamma_cnstr_iv <- calc_gamma(s_iv,6,bern_6)
gamma_cnstr_tsls <- calc_gamma(s_tsls,6,bern_6)
gamma_star <- calc_gamma_star(6,bern_6)

# Now optimize
A <- matrix(c(c(unlist((gamma_cnstr_iv))),
                    c(unlist((gamma_cnstr_tsls)))),
                  nrow=2,
                  byrow=TRUE
)
rhs <- c(beta_s_iv,beta_s_tsls)
sense <- c('=','=')
obj <- c(unlist((gamma_star)))


r <- lp(A, rhs, sense, obj)

res <- rbind(res, data.frame(r$min, r$max))


#-------------------------------------------------------------------------
# Part 3: Calc Bounds for MTRs Restricted to Decreasing
#-------------------------------------------------------------------------
gamma_cnstr_iv <- calc_gamma(s_iv,4,bern_4)
gamma_cnstr_tsls <- calc_gamma(s_tsls,4,bern_4)

gamma_star <- calc_gamma_star(4,bern_4)

incr_c1 <- matrix(bandSparse(4,5,0:1,list(rep(1,4), rep(-1,4))), nrow=4)
incr_c2 <- matrix(bandSparse(4,5,0:1,list(rep(1,4), rep(-1,4))), nrow=4)
incr_c <- matrix(bdiag(incr_c1, incr_c2), nrow=8)

A <- matrix(c(c(unlist((gamma_cnstr_iv))),
                    c(unlist((gamma_cnstr_tsls))),
                    c(t(incr_c))),
                  nrow=10,
                  byrow=TRUE
)
rhs <- c(beta_s_iv,beta_s_tsls,rep(0.0,8))
sense <- c('=','=',rep('<=',8))
obj <- c(unlist((gamma_star)))

r <- lp(A, rhs, sense, obj)

res <- rbind(res, data.frame(r$min, r$max))


#-------------------------------------------------------------------------
# Part 4: Calc Bounds for MTRs Restricted to Increasing
#-------------------------------------------------------------------------
gamma_cnstr_iv <- calc_gamma(s_iv,4,bern_4)
gamma_cnstr_tsls <- calc_gamma(s_tsls,4,bern_4)

gamma_star <- calc_gamma_star(4,bern_4)

incr_c1 <- matrix(bandSparse(4,5,0:1,list(rep(1,4), rep(-1,4))), nrow=4)
incr_c2 <- matrix(bandSparse(4,5,0:1,list(rep(1,4), rep(-1,4))), nrow=4)
incr_c <- matrix(bdiag(incr_c1, incr_c2), nrow=8)

A <- matrix(c(c(unlist((gamma_cnstr_iv))),
              c(unlist((gamma_cnstr_tsls))),
              c(t(incr_c))),
            nrow=10,
            byrow=TRUE
)
rhs <- c(beta_s_iv,beta_s_tsls,rep(0.0,8))
sense <- c('=','=',rep('>=',8))
obj <- c(unlist((gamma_star)))

r <- lp(A, rhs, sense, obj)

res <- rbind(res, data.frame(r$min, r$max))

#-------------------------------------------------------------------------
# Part 5: Calculate Bounds with z_0=2 to z_1=3 constraint 
#-------------------------------------------------------------------------

beta_s_tsls <- calc_beta_s(s_tsls, exp_m_true)
beta_s_iv <- calc_beta_s(s_iv, exp_m_true)
beta_s_wald <- calc_beta_s(s_z2_z3, exp_m_true)

gamma_cnstr_iv <- calc_gamma(s_iv,4,bern_4)
gamma_cnstr_tsls <- calc_gamma(s_tsls,4,bern_4)
gamma_cnstr_wald <- calc_gamma(s_z2_z3,4,bern_4)
gamma_star <- calc_gamma_star(4,bern_4)

# Now optimize
A <- matrix(c(c(unlist(gamma_cnstr_iv)),
                    c(unlist(gamma_cnstr_tsls)),
                    c(unlist(gamma_cnstr_wald))),
                  nrow=3,
                  byrow=TRUE
)
rhs <- c(beta_s_iv,beta_s_tsls,beta_s_wald)
sense <- c('=','=','=')
obj <- c(unlist((gamma_star)))

r <- lp(A, rhs, sense, obj)

res <- rbind(res, data.frame(r$min, r$max))

#-------------------------------------------------------------------------
# Print Results
#-------------------------------------------------------------------------

round(res, 3)
