# Ariel Boyarsky
# aboyarsky@uchicago.edu
#
# The following code provides the functions used in Boyarsky_Q2.R
# following the procedure of Mogstad Torgovitsky (2018).

#-------------------------------------------------------------------------
# Specification from Section 2.4 Mogstad Torgovitsky (2018)
#-------------------------------------------------------------------------

# Propensity Score - E[D=1|Z=z]
pz <- function(z){
  return(switch(z,
                0.12,
                0.29,
                0.48,
                0.78))
}

# M_0, D = 0
i_m_0 <- function(lb, ub){
  u <- 0.9 * ub - (1.1 / 2) * ub^2 + (0.3 / 3) * ub^3
  l <- 0.9 * lb - (1.1 / 2) * lb^2 + (0.3 / 3) * lb^3
  return(u-l)
}

# M_1, D = 1
i_m_1 <- function(lb, ub){
  u <- 0.35 * ub - (0.3 / 2) * ub^2 - (0.05 / 3) * ub^3
  l <- 0.35 * lb - (0.3 / 2) * lb^2 - (0.05 / 3) * lb^3
  return(u-l)
}

# E[D]
E_d <- 1/4*(pz(1)+pz(2)+pz(3)+pz(4))

# E[DZ]
E_dz <- 1/4*(pz(1)*1 + pz(2)*2 + pz(3)*3 + pz(4)*4) 

# E[Z]
E_z <- 2.5

#-------------------------------------------------------------------------
# MT2018 Replication Functions - Sections 4 and 5
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# S Functions
#-------------------------------------------------------------------------
s_iv <- function(z){
  num <- z - E_z
  cov_dz <- E_dz - 2.5*E_d
  return(num/cov_dz)
}

s_tsls <- function(z){
  p <- c(1/4 * pz(1), 1/4 * pz(2), 1/4 * pz(3), 1/4 * pz(4))
  XZ <- matrix(rbind(1/4, p), nrow=2)
  ZZ <- diag(0.25, 4)
  pi <- XZ %*%solve(ZZ)
  Z_tilde <- c(1, 2, 3, 4)
  Z <- matrix(Z_tilde, nrow = 4, ncol = 1)
  tslss <- solve(pi%*%t(XZ))%*%pi
  return(tslss[2,z])
}

s_z2_z3 <- function(z){
  denom <- pz(3) - pz(2)
  n1 <- ifelse(z==3,1,0)/0.25
  n2 <- ifelse(z==2,1,0)/0.25
  s <- (n1-n2) / denom
  return(s)
}

#-------------------------------------------------------------------------
# Integral of Bernstein Polynomials
#-------------------------------------------------------------------------

# 6th Degree Bernstein Polynomial
bern_6 <-function(lb,ub,k){
  if(k==0){
    u <- 1/7 * (ub-1)^7
    l <- 1/7 * (lb-1)^7
  }else if(k==1){
    u <- -6*((ub^7 / 7) - (5*ub^6 / 6) + 2*ub^5 - (5/2 * ub^4) + (5/3 * ub^3)
             - (ub^2 / 2))
    l <- -6*((lb^7 / 7) - (5*lb^6 / 6) + 2*lb^5 - (5/2 * lb^4) + (5/3 * lb^3)
             - (lb^2 / 2))
  }else if(k==2){
    u <- 15*((ub^7 / 7) - (2/3 * ub^6) + (6/5 * ub^5) - ub^4 + (ub^3/3))
    l <- 15*((lb^7 / 7) - (2/3 * lb^6) + (6/5 * lb^5) - lb^4 + (lb^3/3))
  }else if(k==3){
    u <- -20*((ub^7 / 7) - (1/2 * ub^6) + (3/5 * ub^5) - (1/4 * ub^4))
    l <- -20*((lb^7 / 7) - (1/2 * lb^6) + (3/5 * lb^5) - (1/4 * lb^4))
  }else if(k==4){
    u <- (15/7 * ub^7) - (5 * ub^6) + (3 * ub^5)
    l <- (15/7 * lb^7) - (5 * lb^6) + (3 * lb^5)
  }else if(k==5){
    u <- -6*((1/7 * ub^7) - (1/6 * ub^6))
    l <- -6*((1/7 * lb^7) - (1/6 * lb^6))
  }else if(k==6){
    u <- (1/7)*ub^7
    l <- (1/7)*lb^7
  }
  return(u-l)
}

# 4th Degree Bernstein Polynomial
bern_4 <- function(lb, ub, k){
  if(k==0){
    l <- (1/5)*(lb-1)^5
    u <- (1/5)*(ub-1)^5
  }else if(k==1){
    l <- -4*((lb^5 / 5) - (3*lb^4 /4) + lb^3 - (lb^2 / 2))
    u <- -4*((ub^5 / 5) - (3*ub^4 /4) + ub^3 - (ub^2 / 2))
  }else if(k==2){
    l <- (6*lb^5 / 5) - 3*lb^4 + 2*lb^3
    u <- (6*ub^5 / 5) - 3*ub^4 + 2*ub^3
  }else if(k==3){
    l <- -4*(lb^5/5 - lb^4/4)
    u <- -4*(ub^5/5 - ub^4/4)
  }else if(k==4){
    l <- 1/5*lb^5
    u <- 1/5*ub^5
  }
  return(u-l)
}

# 3rd Degree Bernstein Polynomial
bern_3 <- function(lb, ub, k){
  if(k == 0){
    l <- (-1/4)*(1-lb)^4
    u <- (-1/4)*(1-ub)^4
  }else if(k==1){
    l <- (3 * lb^4 / 4) - (2 * lb^3) + (3 * lb^2 / 2)
    u <- (3 * ub^4 / 4) - (2 * ub^3) + (3 * ub^2 / 2)
  }else if(k==2){
    l <- lb^3 - (3/4 * lb^4)
    u <- ub^3 - (3/4 * ub^4)
  }else if(k==3){
    l <- (1/4)*(lb)^4
    u <- (1/4)*(ub)^4
  }
  return(u-l)
}

# 2nd Degree Bernstein Polynomial
bern_2 <- function(lb,ub,k){
  if (k==0){
    l <- 1/3*(lb-1)^3
    u <- 1/3*(ub-1)^3
  }else if(k==1){
    l <- -2*((lb^3 / 3) - (lb^2 / 2))
    u <- -2*((ub^3 / 3) - (ub^2 / 2))
  }else if(k==2){
    l <- (1/3)*lb^3
    u <- (1/3)*ub^3
  }
  return(u-l)
}


#-------------------------------------------------------------------------
# Integration Functions
#-------------------------------------------------------------------------

exp_m_bern <- function(k,d,z,bern){
  if(d == 1){
    # do m_1
    lb <- 0
    ub <- pz(z)
    m <- bern(lb, ub, k)
  }else{
    # do m_0
    lb <- pz(z)
    ub <- 1
    m <- bern(lb, ub, k)
  }
  return(m)
}

exp_m_bern_star <- function(k,d,z,bern){
  if(d == 1){
    # do m_1
    lb <- 0
    ub <- pz(z)
    m <- bern(lb, ub, k)
  }else{
    # do m_0
    lb <- 0
    ub <- pz(z)
    m <- bern(lb, ub, k)
  }
  return(m)
}

exp_m_true <- function(z,d){
  if(d == 1){
    # do m_1
    lb <- 0
    ub <- pz(z)
    m <- i_m_1(lb,ub)
  }else{
    # do d == 0
    lb <- pz(z)
    ub <- 1
    m <- i_m_0(lb,ub)
  }
  return(m)
}

#-------------------------------------------------------------------------
# Beta_S Functions
#-------------------------------------------------------------------------

calc_beta_s <- function(s, exp_m){
  # expectation over z
  beta <- exp_z(s, D = 1, exp_m) + exp_z(s, D = 0, exp_m)
  return(beta)
}

exp_z <- function(s, D, exp_m){
  m <- 1/4 *(s(1)*exp_m(1,D) + s(2)*exp_m(2,D) + s(3)*exp_m(3,D) 
             + s(4)*exp_m(4,D))
  return(m)
}

#-------------------------------------------------------------------------
# Gamma Functions
#-------------------------------------------------------------------------

calc_gamma <- function(s,k,bern){
  dat <- matrix(ncol=2, nrow = k+1)
  for(i in 0:k){
    #print(i)
    for(d in 0:1){
      dat[i+1,d+1] <- gamma_dk(d,i,s,bern) 
    }
  }
  return(dat)
}

gamma_dk <- function(d,k,s,bern){
  gamma <- 1/4 * (s(1)*exp_m_bern(k,d,1,bern) + s(2)*exp_m_bern(k,d,2,bern) +
                    s(3)*exp_m_bern(k,d,3,bern) + s(4)*exp_m_bern(k,d,4,bern))
  return(gamma)
}


##
# Gamma Star Functions
##

calc_gamma_star <- function(k,bern){
  dat <- matrix(ncol=2, nrow = k+1)
  for(i in 0:k){
    for(d in 0:1){
      dat[i+1,d+1] <- gamma_dk_star(d,i,bern) 
    }
  }
  return(dat)
}

gamma_dk_star <- function(d,k,bern){
  if(d == 1){
    a <- 1 / E_d
    gamma <- a * 1/4 * (exp_m_bern_star(k,d,1,bern) 
                        + exp_m_bern_star(k,d,2,bern) +
                          exp_m_bern_star(k,d,3,bern) + 
                          exp_m_bern_star(k,d,4,bern))
  }else{
    a <- 1 / E_d
    gamma <- -a * 1/4 * (exp_m_bern_star(k,d,1,bern) + 
                           exp_m_bern_star(k,d,2,bern) +
                           exp_m_bern_star(k,d,3,bern) +
                           exp_m_bern_star(k,d,4,bern))
  }
  return(gamma)
}


#-------------------------------------------------------------------------
# Linear Program Function
#-------------------------------------------------------------------------

lp <- function(A, rhs, sense, obj){
  k <- dim(A)[2]
  
  model <- list()
  model$A <- A
  model$rhs <- rhs
  model$sense <- sense
  model$obj <- obj
  model$modelsense <- 'max'
  model$ub <- c(rep(1,k))
  model$lb <- c(rep(0,k))
  
  # calc max
  result <- gurobi(model)
  max <- result$objval
  
  # calc min
  model_min <- model
  model_min$modelsense <- 'min'
  result2 <- gurobi(model_min)
  min <- result2$objval
  
  return(list('min'= min, 'max'=max))
}