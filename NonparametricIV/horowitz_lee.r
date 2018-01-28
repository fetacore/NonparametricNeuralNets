library(compiler)
enableJIT(3)
# For the Legendre Polynomials ########################################
exclamationMark <- function(x) {
  m <- 1
  if (x==0) {
    return(m)
  } else if (x==1) {
    return(m)
  } else {
    for (i in 1:x) {
      m = m * i
    }
    return(m)
  }
}
combinatorics <- function(up, down) {
  one <- exclamationMark(up)
  two <- exclamationMark(down)
  three <- exclamationMark(up-down)
  c <- one/(two*three)
  return(c)
}
legendrePsi <- function(x, base) {
  s <- rep(0, length(x))
  s0 <- rep(1, length(x))
  s1 <- rep(0, length(x))
  s2 <- rep(0, length(x))
  s3 <- rep(0, length(x))
  s4 <- rep(0, length(x))
  s5 <- rep(0, length(x))
  if (base > 5) {
    for (k in 0:base) {
      #x = x*(2*k+1)/2
      s = s + combinatorics(base, k)*combinatorics(base+k, k)*
        ((2*x*k + x -2)/4)^k
    }
    return(s)
  } else {
    if (base == 1) {
      s1 = combinatorics(base, 1)*combinatorics(base+1, 1)*((2*x + x -2)/4)
      return(s0 + s1)
    } else {
      s1 = combinatorics(base, 1)*combinatorics(base+1, 1)*((2*x + x -2)/4)
      if (base == 2) {
        s2 = combinatorics(base, 2)*combinatorics(base+2, 2)*((4*x + x -2)/4)^2
        return(s0 + s1 + s2)
      } else {
        s2 = combinatorics(base, 2)*combinatorics(base+2, 2)*((4*x + x -2)/4)^2
        if (base == 3) {
          s3 = combinatorics(base, 3)*combinatorics(base+3, 3)*((6*x + x -2)/4)^3
          return(s0 + s1 + s2 + s3)
        } else {
          s3 = combinatorics(base, 3)*combinatorics(base+3, 3)*((6*x + x -2)/4)^3
          if (base == 4) {
            s4 = combinatorics(base, 4)*combinatorics(base+4, 4)*((8*x + x -2)/4)^4
            return(s0 + s1 + s2 + s3 + s4)
          } else {
            s4 = combinatorics(base, 4)*combinatorics(base+4, 4)*((8*x + x -2)/4)^4
            s5 <- combinatorics(base, 5)*combinatorics(base+5, 4)*((10*x + x -2)/4)^5
            return(s0 + s1 + s2 + s3 + s4 + s5)
          }
        }
      }
    }
  }
}
#######################################################################
# In Horowitz(2007) they use j_max = 50 ###############################
# Conditional Densities for Gibbs Sampling ############################
fxGivenW <- function(a, x, w) {
  s_1 <- 0
  s_2 <- 0
  for (j in 1:50) {
    s_1 = s_1 + ((-1)^(j+1))*(j^(-a/2))*sin(j*pi*x)*sin(j*pi*w)
    s_2 = s_2 + ((-1)^(j+2))*(j^(-a/2))*(((-1)^{j}/(j*pi))-(1/(j*pi)))*sin(j*pi*w)
  }
  f_xGivenW <- s_1/s_2
  return(f_xGivenW)
}
fwGivenX <- function(a, w, x) {
  s_1 <- 0
  s_2 <- 0
  for (j in 1:50) {
    s_1 = s_1 + ((-1)^(j+1))*(j^(-a/2))*sin(j*pi*w)*sin(j*pi*x)
    s_2 = s_2 + ((-1)^(j+2))*(j^(-a/2))*(((-1)^{j}/(j*pi))-(1/(j*pi)))*sin(j*pi*x)
  }
  f_wGivenX <- s_1/s_2
  return(f_wGivenX)
}
gibbs <- function(a, x_0, w_0, n) {
  X <- 1:(n+1000)
  W <- 1:(n+1000)
  for (i in 1:(n+1000)) {
    if (i == 1) {
      X[i] = fxGivenW(a, x_0, w_0)
      W[i] = fwGivenX(a, w_0, X[i])
    } else {
      X[i] = fxGivenW(a, X[i-1], W[i-1])
      W[i] = fwGivenX(a, W[i-1], X[i])
    }
  }
  X <- X[300:1000]
  W <- W[300:1000]
  sampleX <- sample(X[which((X <= 1) & (X >= 0))], n, replace = TRUE)
  sampleW <- sample(W[which((W <= 1) & (W >= 0))], n, replace = TRUE)
  return(list(X=sampleX, W=sampleW))
}
conditionalExpXgivenW <- function(a, x, W) {
  E_eachW <- fxGivenW(a, X, w)
  return(sum(E_eachW))
}
dependent <- function(a, X, W, noise, n) {
  E <- 1:n
  for (i in 1:n) {
    E[i] = 1/n * sum(X*fxGivenW(a, X, W[i]))
  }
  return(E + noise)
}
#######################################################################
#######################################################################
# Just want to check if we need to separate scalar x and w or whole ###
# vectors for denominator (for fxw) ###################################
sieve_optimum <- function(psiXj, psiWk, cjk, Y, J_n) {
  n <- length(psiXj[1,])
  ak <- matrix(0, J_n, length(psiXj[1,]))
  ak[1,] = sum(Y*psiWk[1,])
  ak[2,] = sum(Y*psiWk[2,])
  ak[3,] = sum(Y*psiWk[3,])
  if (J_n > 3) {
    ak[4,] = sum(Y*psiWk[4,])
    if (J_n == 5) {
      ak[5,] = sum(Y*psiWk[5,])
    }
  }
  s <- matrix(1/J_n, 1, J_n)%*%((solve(cjk)%*%ak)*psiXj)
  return(s)
}
sieve <- function(psiJx, psiXj, psiWk, cjk, Y, J_n) {
  n <- length(psiXj[1,])
  ak <- matrix(0, J_n, length(psiJx[1,]))
  ak[1,] = sum(Y*psiWk[1,])
  ak[2,] = sum(Y*psiWk[2,])
  ak[3,] = sum(Y*psiWk[3,])
  if (J_n > 3) {
    ak[4,] = sum(Y*psiWk[4,])
    if (J_n == 5) {
      ak[5,] = sum(Y*psiWk[5,])
    }
  }
  s <- matrix(1/J_n, 1, J_n)%*%((solve(cjk)%*%ak)*psiJx)
  return(s)
}
#######################################################################
#### The confidence sets ##############################################
delta_hat_n <- function(psiJx, psiXj, psiWk, cjk, Y, gSieve, J_n, n) {
  sum_k <- matrix(0, n, 100)
  for (k in 1:J_n) {
    sum_k = sum_k + matrix(psiWk[k,], n, 100)*psiJx[k,]
      #A_hat_inverse_psi(psiJx, psiXj, cjk, k, J_n, n)
  }
  # nx100 dimensional
  first_part <- matrix(Y-gSieve, n, 100)
  return(first_part*sum_k)
}
deltas <- function(psiJx1, psiJx2, psiJx3,
                   psiXj,psiWk, cjk,
                   psiXj_bootstrap, psiWk_bootstrap, cjk_bootstrap,
                   Y, bootstrap_Y, gSieve, gSieveBootstrap, J_n, n) {
  # Takes all 3 intervals' grid points so vectorize 3*n x 100
  delta_hat <- matrix(0, 3*n, 100)
  delta_hat_bootstrap <- matrix(0, 3*n, 100)
  # We get 3 sum_deltas
  sum_delta_hat <- matrix(0, 3, 100)
  sum_delta_hat_bootstrap <- matrix(0, 3, 100)
  
  delta_hat[1:n,] = delta_hat_n(psiJx1, psiXj, psiWk, cjk, Y, gSieve, J_n, n)
  delta_hat_bootstrap[1:n,] = delta_hat_n(psiJx1,
                                          psiXj_bootstrap,
                                          psiWk_bootstrap,
                                          cjk_bootstrap,
                                          bootstrap_Y,
                                          gSieveBootstrap,
                                          J_n, n)
  delta_hat[(n+1):(2*n),] = delta_hat_n(psiJx2, psiXj, psiWk, cjk, Y, gSieve, J_n, n)
  delta_hat_bootstrap[(n+1):(2*n),] = delta_hat_n(psiJx2,
                                                  psiXj_bootstrap,
                                                  psiWk_bootstrap,
                                                  cjk_bootstrap,
                                                  bootstrap_Y,
                                                  gSieveBootstrap,
                                                  J_n, n)
  delta_hat[(2*n+1):(3*n),] = delta_hat_n(psiJx3, psiXj, psiWk, cjk, Y, gSieve, J_n, n)
  delta_hat_bootstrap[(2*n+1):(3*n),] = delta_hat_n(psiJx3,
                                                    psiXj_bootstrap,
                                                    psiWk_bootstrap,
                                                    cjk_bootstrap,
                                                    bootstrap_Y,
                                                    gSieveBootstrap,
                                                    J_n, n)
  # We get delta overline n for every point of the grid ################
  ones <- matrix(1, 1, n)
  sum_delta_hat[1,] = ones%*%delta_hat[1:n,]
  sum_delta_hat[2,] = ones%*%delta_hat[(n+1):(2*n),]
  sum_delta_hat[3,] = ones%*%delta_hat[(2*n+1):(3*n),]
  sum_delta_hat_bootstrap[1,] = ones%*%delta_hat_bootstrap[1:n,]
  sum_delta_hat_bootstrap[2,] = ones%*%delta_hat_bootstrap[(n+1):(2*n),]
  sum_delta_hat_bootstrap[3,] = ones%*%delta_hat_bootstrap[(2*n+1):(3*n),]
  # The returned delta hat normal and bootstrap are 3*nx100 ##############
  # The returned delta overline and S are 3x100 ##########################
  return(list(overline = (1/n)*sum_delta_hat,
              overline_bootstrap = (1/n)*sum_delta_hat_bootstrap,
              hat_normal = delta_hat,
              hat_bootstrap = delta_hat_bootstrap)
         )
}
sigmas <- function(delta_overline, delta_hat, A1, A2, A3, n) {
  sum_for_every_grid_point_sq <- matrix(0, 3, 100)
  ones <- matrix(1, 1, n)
  sum_for_every_grid_point_sq[1,] = ones%*%((delta_hat[1:n,] - delta_overline[1,])/A1)^2
  sum_for_every_grid_point_sq[2,] = ones%*%((delta_hat[(n+1):(2*n),] - delta_overline[2,])/A2)^2
  sum_for_every_grid_point_sq[3,] = ones%*%((delta_hat[(2*n+1):(3*n),] - delta_overline[3,])/A3)^2
  # Dimensions 3x100 ###################################################
  return((1/(n^2)*sum_for_every_grid_point_sq))
}
t_n <- function(psiJx1, psiJx2, psiJx3,
                psiXj, psiWk, cjk,
                psiXj_bootstrap, psiWk_bootstrap, cjk_bootstrap,
                Y, bootstrap_Y, gSieve, gSieveBootstrap, J_n, n) {
  d <- list()
  up_onestar <- matrix(0, 3, 100)
  up_twostar <- matrix(0, 3, 100)
  s_one <- matrix(0, 3, 100)
  s_two <- matrix(0, 3, 100)
  ones <- matrix(1, 1, n)
  onesJn <- matrix(1/n, n, J_n)
  onesJnT <- matrix(1, J_n, 1)
  Ainverse1 <- drop(onesJn%*%cjk%*%onesJnT)
  Ainverse2 <- drop(onesJn%*%cjk%*%onesJnT)
  Ainverse3 <- drop(onesJn%*%cjk%*%onesJnT)
  Ainverse1_bootstrap <- drop(onesJn%*%cjk_bootstrap%*%onesJnT)
  Ainverse2_bootstrap <- drop(onesJn%*%cjk_bootstrap%*%onesJnT)
  Ainverse3_bootstrap <- drop(onesJn%*%cjk_bootstrap%*%onesJnT)
  # It takes all intervals' grid points as vector ######################
  d <- deltas(psiJx1, psiJx2, psiJx3,
              psiXj,psiWk, cjk,
              psiXj_bootstrap, psiWk_bootstrap, cjk_bootstrap,
              Y, bootstrap_Y, gSieve, gSieveBootstrap, J_n, n)
  # d$overline_bootstrap 3x100 d$overline 3x100 d$hat_normal 3nx100 d$hat_bootstrap 3nx100
  s_one = sigmas(d$overline, d$hat_normal,
                 Ainverse1, Ainverse2, Ainverse3, n)
  s_two = sigmas(d$overline_bootstrap, d$hat_bootstrap,
                 Ainverse1_bootstrap, Ainverse2_bootstrap, Ainverse3_bootstrap, n)
  # both s_one and s_two are 3x100 so t* and t** are 3x100 #############
  up_onestar[1,] = ones%*%((1/Ainverse1)*(d$hat_bootstrap[1:n,] - d$overline[1,]))
  up_twostar[1,] <- ones%*%((1/Ainverse1_bootstrap)*(d$hat_bootstrap[1:n,] - d$overline_bootstrap[1,]))
  up_onestar[2,] = ones%*%((1/Ainverse2)*(d$hat_bootstrap[(n+1):(2*n),] - d$overline[2,]))
  up_twostar[2,] <- ones%*%((1/Ainverse2_bootstrap)*(d$hat_bootstrap[(n+1):(2*n),] - d$overline_bootstrap[2,]))
  up_onestar[3,] = ones%*%((1/Ainverse3)*(d$hat_bootstrap[(2*n+1):(3*n),] - d$overline[3,]))
  up_twostar[3,] <- ones%*%((1/Ainverse3_bootstrap)*(d$hat_bootstrap[(2*n+1):(3*n),] - d$overline_bootstrap[3,]))
  return(list(one_star = up_onestar/s_one,
              two_star = up_twostar/s_two,
              s_one = s_one,
              s_two = s_two))
}
lipschitz_choice_of_xl <- function(xl_grid_points, interesting_X) {
  # Returns l ##########################################################
  lX <- length(interesting_X)
  grid_matrix <- matrix(xl_grid_points, 100, lX)
  # We get for all 100 grid points the absolute distance between the
  # vector of interesting X and the vector consisting of only xl forall l
  difference_matrix <- abs(grid_matrix - interesting_X)
  difference_vector <- difference_matrix %*% matrix(1, lX, 1)
  return(which.min(difference_vector))
}
############################################################################################
############################################################################################
##### Lipschitz Test #######################################################################
test_lipschitz <- function(intervalLow, intervalHigh, X, xl_grid_points,
                           gSieve, grid_sieves, gReal, s_one,
                           z_a_two_star_10, z_a_two_star_05, z_a_two_star_01) {
  ## INITIALIZATION ########################################################################
  # Suppose C_Lipschitz = 1
  CLip <- 1
  
  lower_bound_twostar_10 <- 0
  lower_bound_twostar_05 <- 0
  lower_bound_twostar_01 <- 0
  
  upper_bound_twostar_10 <- 0
  upper_bound_twostar_05 <- 0
  upper_bound_twostar_01 <- 0
  
  x_choice_low <- 0
  x_choice_high <- 0
  l_index <- 0
  x_low <- 0
  x_high <- 0
  # First 3 rows for one star, last 3 for two star ######################################
  # First 3 columns for pointwise, last 3 columns for uniform ###########################
  results <- matrix(0, 3, 6)
  #######################################################################################
  for (interval in 1:3) {
    x_choice_low <- intervalLow[interval] + (intervalHigh[interval]-intervalLow[interval])/100
    x_choice_high <- intervalLow[interval] + (99*(intervalHigh[interval]-intervalLow[interval]))/100
    # The Xs that we should check because the crappy estimator is meant to be crappy!!
    interesting_X <- X[which((X >= x_choice_low) & (X <= x_choice_high))]
    l_index <- lipschitz_choice_of_xl(xl_grid_points[interval], interesting_X)
    # Within this range because elsewhere the estimator is supercrap ####
    x_low <- xl_grid_points[interval,l_index]-1/100
    x_high <- xl_grid_points[interval,l_index]+1/100
    interesting_i <- which((X >= x_low) & (X <= x_high))
    if (length(interesting_i) == 0) {
      # If the estimator is completelly crap in every point go to next interval #############
      results[interval, 1:6] = matrix(1, 1, 6)
      next()
    } else {
      interesting_sieves <- gSieve[interesting_i]
      lower_bound_twostar_10 = grid_sieves[l_index] -
        z_a_two_star_10[interval] * s_one[interval,l_index] - CLip/100
      upper_bound_twostar_10 = grid_sieves[l_index] +
        z_a_two_star_10[interval] * s_one[interval,l_index] + CLip/100
      lower_bound_twostar_10_vector <- matrix(lower_bound_twostar_10, 1, length(interesting_i))
      upper_bound_twostar_10_vector <- matrix(upper_bound_twostar_10, 1, length(interesting_i))
      #lower_bound_twostar_10_vector = interesting_sieves -
      #  z_a_two_star_10[interval] * s_one[interval,] - CLip/100
      #upper_bound_twostar_10_vector = interesting_sieves +
      #  z_a_two_star_10[interval] * s_one[interval,] + CLip/100
      
      lower_bound_twostar_05 = grid_sieves[l_index] -
        z_a_two_star_05[interval] * s_one[interval,l_index] - CLip/100
      upper_bound_twostar_05 = grid_sieves[l_index] +
        z_a_two_star_05[interval] * s_one[interval,l_index] + CLip/100
      lower_bound_twostar_05_vector <- matrix(lower_bound_twostar_05, 1, length(interesting_i))
      upper_bound_twostar_05_vector <- matrix(upper_bound_twostar_05, 1, length(interesting_i))
      #lower_bound_twostar_05_vector = interesting_sieves -
      #  z_a_two_star_05[interval] * s_one[interval,] - CLip/100
      #upper_bound_twostar_05_vector = interesting_sieves +
      #  z_a_two_star_05[interval] * s_one[interval,] + CLip/100
      
      lower_bound_twostar_01 = grid_sieves[l_index] -
        z_a_two_star_01[interval] * s_one[interval,l_index] - CLip/100
      upper_bound_twostar_01 = grid_sieves[l_index] +
        z_a_two_star_01[interval] * s_one[interval,l_index] + CLip/100
      lower_bound_twostar_01_vector <- matrix(lower_bound_twostar_01, 1, length(interesting_i))
      upper_bound_twostar_01_vector <- matrix(upper_bound_twostar_01, 1, length(interesting_i))
      #lower_bound_twostar_01_vector = interesting_sieves -
      #  z_a_two_star_01[interval] * s_one[interval,] - CLip/100
      #upper_bound_twostar_01_vector = interesting_sieves +
      #  z_a_two_star_01[interval] * s_one[interval,] + CLip/100
      
      ### Pointwise procedure ############################################
      interesting_gReal <- matrix(gReal[interesting_i], 1, length(gReal[interesting_i]))
      if ((all(interesting_gReal >= lower_bound_twostar_10_vector)) & 
          (all(interesting_gReal <= upper_bound_twostar_10_vector))) {
        results[interval, 1] = 1
      } else {
        results[interval, 1] = 0
      }
      if ((all(interesting_gReal >= lower_bound_twostar_05_vector)) & 
          (all(interesting_gReal <= upper_bound_twostar_05_vector))) {
        results[interval, 2] = 1
      } else {
        results[interval, 2] = 0
      }
      if ((all(interesting_gReal >= lower_bound_twostar_01_vector)) & 
          (all(interesting_gReal <= upper_bound_twostar_01_vector))) {
        results[interval, 3] = 1
      } else {
        results[interval, 3] = 0
      }
      ## Uniform Procedure #####################################################################
      distance_bounds_twostar_01 <- uniform_norm(lower_bound_twostar_01_vector,upper_bound_twostar_01_vector)
      distance_bounds_twostar_05 <- uniform_norm(lower_bound_twostar_05_vector,upper_bound_twostar_05_vector)
      distance_bounds_twostar_10 <- uniform_norm(lower_bound_twostar_10_vector,upper_bound_twostar_10_vector)
      if ((uniform_norm(interesting_gReal,lower_bound_twostar_10_vector) <= distance_bounds_twostar_10) & 
          (uniform_norm(interesting_gReal,upper_bound_twostar_10_vector) <= distance_bounds_twostar_10)) {
        results[interval, 4] = 1
      } else {
        results[interval, 4] = 0
      }
      if ((uniform_norm(interesting_gReal,lower_bound_twostar_05_vector) <= distance_bounds_twostar_05) & 
          (uniform_norm(interesting_gReal,upper_bound_twostar_05_vector) <= distance_bounds_twostar_05)) {
        results[interval, 5] = 1
      } else {
        results[interval, 5] = 0
      }
      if ((uniform_norm(interesting_gReal,lower_bound_twostar_01_vector) <= distance_bounds_twostar_01) & 
          (uniform_norm(interesting_gReal,upper_bound_twostar_01_vector) <= distance_bounds_twostar_01)) {
        results[interval, 6] = 1
      } else {
        results[interval, 6] = 0
      }
    }
  }
  return(results)
}
test_monotonicity <- function(X, xl_grid_points,
                              grid_sieves, gReal, s_one, z_a_10,
                              z_a_05, z_a_01) {
  # First 3 rows for one star, last 3 for two star ######################################
  # First 3 columns for pointwise, last 3 columns for uniform ###########################
  results <- matrix(0, 3, 6)
  #######################################################################################
  for (interval in 1:3) {
    # row 1 are the indexes, row 2 are the grid points, row 3 are the upper bounds and row 4 the lower
    # we have to order this from the lower to the higher grid point (row 2) and then we evaluate
    grid_matrix_10 <- matrix(0, 4, 100)
    grid_matrix_05 <- matrix(0, 4, 100)
    grid_matrix_01 <- matrix(0, 4, 100)
    
    grid_matrix_10[1,] = 1:100
    grid_matrix_10[2,] = xl_grid_points[interval,]
    grid_matrix_10[3,] = grid_sieves[interval,] + z_a_10[interval] * s_one[interval,]
    grid_matrix_10[4,] = grid_sieves[interval,] - z_a_10[interval] * s_one[interval,]
    
    grid_matrix_05[1,] = 1:100
    grid_matrix_05[2,] = xl_grid_points[interval,]
    grid_matrix_05[3,] = grid_sieves[interval,] + z_a_05[interval] * s_one[interval,]
    grid_matrix_05[4,] = grid_sieves[interval,] - z_a_05[interval] * s_one[interval,]
    
    grid_matrix_01[1,] = 1:100
    grid_matrix_01[2,] = xl_grid_points[interval,]
    grid_matrix_01[3,] = grid_sieves[interval,] + z_a_01[interval] * s_one[interval,]
    grid_matrix_01[4,] = grid_sieves[interval,] - z_a_01[interval] * s_one[interval,]
    
    # Now everything corresponds to the ascending order of the grid points
    ordered_grid_10 <- grid_matrix_10[1:4,order(grid_matrix_10[2,])]
    ordered_grid_05 <- grid_matrix_05[1:4,order(grid_matrix_05[2,])]
    ordered_grid_01 <- grid_matrix_01[1:4,order(grid_matrix_01[2,])]
    
    interesting_i <- which((X >= min(xl_grid_points[interval,])) & (X <= max(xl_grid_points[interval,])))
    x_with_gReal_matrix <- matrix(0, 3, length(interesting_i))
    x_with_gReal_matrix[1,] = 1:length(interesting_i)
    x_with_gReal_matrix[2,] = X[interesting_i]
    x_with_gReal_matrix[3,] = gReal[interesting_i]
    # This is a matrix with x's in ascending order, matched with their corresponding indexes and gs
    # The ascending xs are in row 2 and indexes in row 1. They are only those inside the interval
    ordered_x_gReal <- x_with_gReal_matrix[1:3, order(x_with_gReal_matrix[2,])]
    
    # For each two grid points (pointwise test)
    # We start assuming that all is good (if something is not in interval these turn to zero)
    r10 <- 1
    r05 <- 1
    r01 <- 1
    lower_unif_bound_10 <- matrix(0, 1, 1)
    upper_unif_bound_10 <- matrix(0, 1, 1)
    lower_unif_bound_05 <- matrix(0, 1, 1)
    upper_unif_bound_05 <- matrix(0, 1, 1)
    lower_unif_bound_01 <- matrix(0, 1, 1)
    upper_unif_bound_01 <- matrix(0, 1, 1)
    test_unif_gReal_10 <- matrix(0, 1, 1)
    test_unif_gReal_05 <- matrix(0, 1, 1)
    test_unif_gReal_01 <- matrix(0, 1, 1)
    for (i in 1:99) {
      if (i == 1) {
        lower_unif_bound_10[1,1] = min(ordered_grid_10[4,i],ordered_grid_10[4,i+1])
        upper_unif_bound_10[1,1] = max(ordered_grid_10[3,i],ordered_grid_10[3,i+1])
        lower_unif_bound_05[1,1] = min(ordered_grid_05[4,i],ordered_grid_05[4,i+1])
        upper_unif_bound_05[1,1] = max(ordered_grid_05[3,i],ordered_grid_05[3,i+1])
        lower_unif_bound_01[1,1] = min(ordered_grid_01[4,i],ordered_grid_01[4,i+1])
        upper_unif_bound_01[1,1] = max(ordered_grid_01[3,i],ordered_grid_01[3,i+1])
      } else if (i < 99) {
        lower_unif_bound_10 <- cbind(lower_unif_bound_10, min(ordered_grid_10[4,i],ordered_grid_10[4,i+1]))
        lower_unif_bound_05 <- cbind(lower_unif_bound_05, min(ordered_grid_05[4,i],ordered_grid_05[4,i+1]))
        lower_unif_bound_01 <- cbind(lower_unif_bound_01, min(ordered_grid_01[4,i],ordered_grid_01[4,i+1]))
        upper_unif_bound_10 <- cbind(upper_unif_bound_10, max(ordered_grid_10[3,i],ordered_grid_10[3,i+1]))
        upper_unif_bound_05 <- cbind(upper_unif_bound_05, max(ordered_grid_05[3,i],ordered_grid_05[3,i+1]))
        upper_unif_bound_01 <- cbind(upper_unif_bound_01, max(ordered_grid_01[3,i],ordered_grid_01[3,i+1]))
      } else {
        lower_unif_bound_10 <- cbind(lower_unif_bound_10, ordered_grid_10[4,100])
        lower_unif_bound_05 <- cbind(lower_unif_bound_05, ordered_grid_05[4,100])
        lower_unif_bound_01 <- cbind(lower_unif_bound_01, ordered_grid_01[4,100])
        upper_unif_bound_10 <- cbind(upper_unif_bound_10, ordered_grid_10[3,100])
        upper_unif_bound_05 <- cbind(upper_unif_bound_05, ordered_grid_05[3,100])
        upper_unif_bound_01 <- cbind(upper_unif_bound_01, ordered_grid_01[3,100])
      }
      # Those xs that are inside xl and xl+1
      interesting_x_10 <- which(
        (ordered_x_gReal[2,] >= ordered_grid_10[2,i]) & 
          (ordered_x_gReal[2,] <= ordered_grid_10[2,i+1])
        )
      interesting_x_05 <- which(
        (ordered_x_gReal[2,] >= ordered_grid_05[2,i]) & 
          (ordered_x_gReal[2,] <= ordered_grid_05[2,i+1])
      )
      interesting_x_01 <- which(
        (ordered_x_gReal[2,] >= ordered_grid_01[2,i]) & 
          (ordered_x_gReal[2,] <= ordered_grid_01[2,i+1])
      )
      # If there is no x inside this interval continue with the other intervals
      if (length(interesting_x_10) == 0) {
        if (i == 1) {
          test_unif_gReal_10[1,1] = min(ordered_grid_10[4,i],ordered_grid_10[4,i+1])
        } else if (i < 99) {
          test_unif_gReal_10 <- cbind(test_unif_gReal_10, min(ordered_grid_10[4,i],ordered_grid_10[4,i+1]))
        } else {
          test_unif_gReal_10 <- cbind(test_unif_gReal_10, min(ordered_grid_10[4,i],ordered_grid_10[4,i+1]))
        }
      } else {
        if (i == 1) {
          test_unif_gReal_10[1,1] = ordered_x_gReal[3, sample(interesting_x_10, 1)]
        } else if (i < 99) {
          test_unif_gReal_10 <- cbind(test_unif_gReal_10, ordered_x_gReal[3, sample(interesting_x_10, 1)])
        } else {
          test_unif_gReal_10 <- cbind(test_unif_gReal_10, ordered_x_gReal[3, sample(interesting_x_10, 1)])
        }
        for (j in interesting_x_10) {
          # Opposite procedure check if something is not there
          if ((ordered_x_gReal[3,j] < min(ordered_grid_10[4,i],ordered_grid_10[4,i+1])) |
              (ordered_x_gReal[3,j] > max(ordered_grid_10[3,i],ordered_grid_10[3,i+1]))) {
            r10 = 0
          }
        }
      }
      if (length(interesting_x_05) == 0) {
        if (i == 1) {
          test_unif_gReal_05[1,1] = min(ordered_grid_05[4,i],ordered_grid_05[4,i+1])
        } else if (i < 99) {
          test_unif_gReal_05 <- cbind(test_unif_gReal_05, min(ordered_grid_05[4,i],ordered_grid_05[4,i+1]))
        } else {
          test_unif_gReal_05 <- cbind(test_unif_gReal_05, min(ordered_grid_05[4,i],ordered_grid_05[4,i+1]))
        }
      } else {
        if (i == 1) {
          test_unif_gReal_05[1,1] = ordered_x_gReal[3, sample(interesting_x_05, 1)]
        } else if (i < 99) {
          test_unif_gReal_05 <- cbind(test_unif_gReal_05, ordered_x_gReal[3, sample(interesting_x_05, 1)])
        } else {
          test_unif_gReal_05 <- cbind(test_unif_gReal_05, ordered_x_gReal[3, sample(interesting_x_05, 1)])
        }
        for (j in interesting_x_05) {
          # Opposite procedure check if something is not there
          if ((ordered_x_gReal[3,j] < min(ordered_grid_05[4,i],ordered_grid_05[4,i+1])) |
              (ordered_x_gReal[3,j] > max(ordered_grid_05[3,i],ordered_grid_05[3,i+1]))) {
            r05 = 0
          }
        }
      }
      if (length(interesting_x_01) == 0) {
        if (i == 1) {
          test_unif_gReal_01[1,1] = min(ordered_grid_01[4,i],ordered_grid_01[4,i+1])
        } else if (i < 99) {
          test_unif_gReal_01 <- cbind(test_unif_gReal_01, min(ordered_grid_01[4,i],ordered_grid_01[4,i+1]))
        } else {
          test_unif_gReal_01 <- cbind(test_unif_gReal_01, min(ordered_grid_01[4,i],ordered_grid_01[4,i+1]))
        }
      } else {
        if (i == 1) {
          test_unif_gReal_01[1,1] = ordered_x_gReal[3, sample(interesting_x_01, 1)]
        } else if (i < 99) {
          test_unif_gReal_01 <- cbind(test_unif_gReal_01, ordered_x_gReal[3, sample(interesting_x_01, 1)])
        } else {
          test_unif_gReal_01 <- cbind(test_unif_gReal_01, ordered_x_gReal[3, sample(interesting_x_01, 1)])
        }
        for (j in interesting_x_01) {
          # Opposite procedure check if something is not there
          if ((ordered_x_gReal[3,j] < min(ordered_grid_01[4,i],ordered_grid_01[4,i+1])) |
              (ordered_x_gReal[3,j] > max(ordered_grid_01[3,i],ordered_grid_01[3,i+1]))) {
            r01 = 0
          }
        }
      }
    }
    
    results[interval, 1] = r10
    results[interval, 2] = r05
    results[interval, 3] = r01
    
    # Uniform confidence bands
    distance_bounds_10 <- uniform_norm(lower_unif_bound_10, upper_unif_bound_10)
    distance_bounds_05 <- uniform_norm(lower_unif_bound_05, upper_unif_bound_05)
    distance_bounds_01 <- uniform_norm(lower_unif_bound_01, upper_unif_bound_01)
    if ((uniform_norm(test_unif_gReal_10, lower_unif_bound_10) <= distance_bounds_10) &
        (uniform_norm(test_unif_gReal_10, upper_unif_bound_10) <= distance_bounds_10)) {
      results[interval, 4] = 1
    }
    if ((uniform_norm(test_unif_gReal_05, lower_unif_bound_05) <= distance_bounds_05) &
        (uniform_norm(test_unif_gReal_05, upper_unif_bound_05) <= distance_bounds_05)) {
      results[interval, 5] = 1
    }
    if ((uniform_norm(test_unif_gReal_01, lower_unif_bound_01) <= distance_bounds_01) &
        (uniform_norm(test_unif_gReal_01, upper_unif_bound_01) <= distance_bounds_01)) {
      results[interval, 6] = 1
    }
  }
  return(results)
}
uniform_norm <- function(vector1, vector2) {
  return(max(abs(vector1-vector2)))
}
########################################################################
########################################################################
# For monte carlo we will use a=1.2 and a=10 ###########################
#################### [a, b] = [0.2, 0.8], [0.1, 0.9], [0.01, c(0.99)] #####
#################### J_n = 3, 4, 5 #####################################
############################################### 1000 simulations #######
monteCarlo <- function(psiJx1, psiJx2, psiJx3, psiWk, psiXj,
                       cjk, gSieve, grid_sieves, g_real,
                       Y, X, W, xl_grid_points, J_n, n,
                       intervalLow, intervalHigh,
                       boot) {
  ### INITIALIZATION ###################################################
  # It took al three intervals as argument #############################
  # So we get a different z_a* for each interval low and high ##########
  z_a_one_star_10 <- 1:3
  z_a_one_star_05 <- 1:3
  z_a_one_star_01 <- 1:3
  z_a_two_star_10 <- 1:3
  z_a_two_star_05 <- 1:3
  z_a_two_star_01 <- 1:3
  # 6x6 -first 3 rows for t* and last 3 rows for t**
  #      first 3 columns for intervals (pointwise) and last 3 for bands (uniform norm)
  res <- matrix(0, 6, 6)
  
  bootstrap_Y <- 1:n
  bootstrap_X <- 1:n
  bootstrap_W <- 1:n
  gSieveBootstrap <- 1:n
  psiWk_bootstrap <- matrix(0, J_n, n)
  psiXj_bootstrap <- matrix(0, J_n, n)
  cjk_bootstrap <- matrix(0, J_n, J_n)
  t <- list()
  zero_vector <- matrix(0, 1, 100)
  # Uniform distance for the two cases of bootstrap for every interval #
  unif_distance <- matrix(0, 3*boot, 2)
  #### END INITIALIZATION ##############################################
  
  for (bootstrap in 1:boot) {
    if (bootstrap == 1) {
      print(paste('J=',J_n,' Bootstrap Sampling Begins'))
    } else if (bootstrap == boot/1.5) {
      print(paste('J=',J_n,' Almost there... Bootstrap Sample no.',bootstrap))
    }
    
    bootstrap_X = sample(X, n, replace = TRUE)
    bootstrap_W = sample(W, n, replace = TRUE)
    bootstrap_Y = sample(Y, n, replace = TRUE)
    psiXj_bootstrap[1,1:n] = legendrePsi(bootstrap_X, 1)
    psiXj_bootstrap[2,1:n] = legendrePsi(bootstrap_X, 2)
    psiXj_bootstrap[3,1:n] = legendrePsi(bootstrap_X, 3)
    psiWk_bootstrap[1,1:n] = legendrePsi(bootstrap_W, 1)
    psiWk_bootstrap[2,1:n] = legendrePsi(bootstrap_W, 2)
    psiWk_bootstrap[3,1:n] = legendrePsi(bootstrap_W, 3)
    cjk_bootstrap[1,1] = sum(psiXj_bootstrap[1,]*psiWk_bootstrap[1,])
    cjk_bootstrap[1,2] = sum(psiXj_bootstrap[2,]*psiWk_bootstrap[1,])
    cjk_bootstrap[1,3] = sum(psiXj_bootstrap[3,]*psiWk_bootstrap[1,])
    cjk_bootstrap[2,1] = sum(psiXj_bootstrap[1,]*psiWk_bootstrap[2,])
    cjk_bootstrap[2,2] = sum(psiXj_bootstrap[2,]*psiWk_bootstrap[2,])
    cjk_bootstrap[2,3] = sum(psiXj_bootstrap[3,]*psiWk_bootstrap[2,])
    cjk_bootstrap[3,1] = sum(psiXj_bootstrap[1,]*psiWk_bootstrap[3,])
    cjk_bootstrap[3,2] = sum(psiXj_bootstrap[2,]*psiWk_bootstrap[3,])
    cjk_bootstrap[3,3] = sum(psiXj_bootstrap[3,]*psiWk_bootstrap[3,])
    if (J_n > 3) {
      psiXj_bootstrap[4,1:n] = legendrePsi(bootstrap_X, 4)
      psiWk_bootstrap[4,1:n] = legendrePsi(bootstrap_W, 4)
      cjk_bootstrap[4,1] = sum(psiXj_bootstrap[1,]*psiWk_bootstrap[4,])
      cjk_bootstrap[4,2] = sum(psiXj_bootstrap[2,]*psiWk_bootstrap[4,])
      cjk_bootstrap[4,3] = sum(psiXj_bootstrap[3,]*psiWk_bootstrap[4,])
      cjk_bootstrap[1,4] = sum(psiXj_bootstrap[4,]*psiWk_bootstrap[1,])
      cjk_bootstrap[2,4] = sum(psiXj_bootstrap[4,]*psiWk_bootstrap[2,])
      cjk_bootstrap[3,4] = sum(psiXj_bootstrap[4,]*psiWk_bootstrap[3,])
      cjk_bootstrap[4,4] = sum(psiXj_bootstrap[4,]*psiWk_bootstrap[4,])
      if (J_n == 5) {
        psiXj_bootstrap[5,1:n] = legendrePsi(bootstrap_X, 5)
        psiWk_bootstrap[5,1:n] = legendrePsi(bootstrap_W, 5)
        cjk_bootstrap[5,1] = sum(psiXj_bootstrap[1,]*psiWk_bootstrap[5,])
        cjk_bootstrap[5,2] = sum(psiXj_bootstrap[2,]*psiWk_bootstrap[5,])
        cjk_bootstrap[5,3] = sum(psiXj_bootstrap[3,]*psiWk_bootstrap[5,])
        cjk_bootstrap[5,4] = sum(psiXj_bootstrap[4,]*psiWk_bootstrap[5,])
        cjk_bootstrap[1,5] = sum(psiXj_bootstrap[5,]*psiWk_bootstrap[1,])
        cjk_bootstrap[2,5] = sum(psiXj_bootstrap[5,]*psiWk_bootstrap[2,])
        cjk_bootstrap[3,5] = sum(psiXj_bootstrap[5,]*psiWk_bootstrap[3,])
        cjk_bootstrap[4,5] = sum(psiXj_bootstrap[5,]*psiWk_bootstrap[4,])
        cjk_bootstrap[5,5] = sum(psiXj_bootstrap[5,]*psiWk_bootstrap[5,])
      }
    }
    gSieveBootstrap = sieve_optimum(psiXj_bootstrap, psiWk_bootstrap, cjk_bootstrap, bootstrap_Y, J_n)
    t = t_n(psiJx1, psiJx2, psiJx3,
            psiXj,psiWk, cjk,
            psiXj_bootstrap, psiWk_bootstrap, cjk_bootstrap,
            Y, bootstrap_Y, gSieve, gSieveBootstrap, J_n, n)
    # t* and t** are 3x100 for every grid point(100) for every interval(3) ###
    # We get the uniform distance of t* and t** from the 0 1x100 vector for
    # every bootstrap sample and this is our dataset for the quantile func
    unif_distance[1:bootstrap,1:2] = c(max(abs(t$one_star[1,])),
                                     max(abs(t$two_star[1,])))
    unif_distance[(bootstrap+1):(2*bootstrap),1:2] = c(max(abs(t$one_star[2,])),
                                       max(abs(t$two_star[2,])))
    unif_distance[(2*bootstrap+1):(3*bootstrap),1:2] = c(max(abs(t$one_star[3,])),
                                                       max(abs(t$two_star[3,])))
  }
  # Prerequisites ########################################################
  z_a_one_star_10[1] = quantile(unif_distance[1:boot,1], .90, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_05[1] = quantile(unif_distance[1:boot,1], .95, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_01[1] = quantile(unif_distance[1:boot,1], .99, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_10[2] = quantile(unif_distance[(boot+1):(2*boot),1], .90, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_05[2] = quantile(unif_distance[(boot+1):(2*boot),1], .95, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_01[2] = quantile(unif_distance[(boot+1):(2*boot),1], .99, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_10[3] = quantile(unif_distance[(2*boot+1):(3*boot),1], .90, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_05[3] = quantile(unif_distance[(2*boot+1):(3*boot),1], .95, name=FALSE, na.rm=TRUE, type=6)
  z_a_one_star_01[3] = quantile(unif_distance[(2*boot+1):(3*boot),1], .99, name=FALSE, na.rm=TRUE, type=6)

  z_a_two_star_10[1] = quantile(unif_distance[1:boot,2], .90, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_05[1] = quantile(unif_distance[1:boot,2], .95, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_01[1] = quantile(unif_distance[1:boot,2], .99, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_10[2] = quantile(unif_distance[(boot+1):(2*boot),2], .90, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_05[2] = quantile(unif_distance[(boot+1):(2*boot),2], .95, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_01[2] = quantile(unif_distance[(boot+1):(2*boot),2], .99, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_10[3] = quantile(unif_distance[(2*boot+1):(3*boot),2], .90, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_05[3] = quantile(unif_distance[(2*boot+1):(3*boot),2], .95, name=FALSE, na.rm=TRUE, type=6)
  z_a_two_star_01[3] = quantile(unif_distance[(2*boot+1):(3*boot),2], .99, name=FALSE, na.rm=TRUE, type=6)
  # We already have s from the t* calculation, it is t$s_one and t$s_two #######
  
  res[1:3, 1:6] <- test_monotonicity(X, xl_grid_points,
                                     grid_sieves, g_real, t$s_one, z_a_one_star_10, 
                                     z_a_one_star_05, z_a_one_star_01)
  res[4:6, 1:6] <- test_monotonicity(X, xl_grid_points,
                                     grid_sieves, g_real, t$s_two, z_a_two_star_10, 
                                     z_a_two_star_05, z_a_two_star_01)
  
  # We need to return an 6x12 dimensional matrix where:
  # -lines 1:2 are for J=3 and columns 1:3 are for intervals and cols 4:6 for bands
  #                            columns 1:6 are for the t* and 7:12 for t**
  # -lines 3:4     for J=4
  # -lines 5:6     for J=5
  return(res)
}
# I propose we do not change the initial values because bad things happen and
# the estimator is even more crap ###########################################
modelVars <- list(experiments = 1000, a = 1.2, x_0=100, w_0=200,
                  J_n=5, n=200, intervalLow = c(0.2, 0.1, 0.01),
                  intervalHigh = c(0.8, 0.9, c(0.99)), boot = 10
)
# Results have 18 rows - every 3 for every 3 different J_n
# every pair of 3 rows is for every grid interval so 3*3 = 9 rows
# every 9 rows are for t* and t** so 2*9 = 18 rows total
# and 6 columns - first 3 for joint intervals and last 3 for uniform bands
# every one of each pair of 3 columns is for a=0.10 0.05 and 0.01 respectively
n <- modelVars$n
allMonte <- matrix(0, 18, 6)
psiWk <- matrix(0, modelVars$J_n, modelVars$n)
psiXj <- matrix(0, modelVars$J_n, modelVars$n)
psiJx1 <- matrix(0, modelVars$J_n, 100)
psiJx2 <- matrix(0, modelVars$J_n, 100)
psiJx3 <- matrix(0, modelVars$J_n, 100)
cjk <- matrix(0, modelVars$J_n, modelVars$J_n)
xl_grid_points <- matrix(0, 3, 100)
grid_sieves <- matrix(0, 3, 100)
indeps <- list()
xl_grid_points[1,] = runif(100, modelVars$intervalLow[1], modelVars$intervalHigh[1])
xl_grid_points[2,] = runif(100, modelVars$intervalLow[2], modelVars$intervalHigh[2])
xl_grid_points[3,] = runif(100, modelVars$intervalLow[3], modelVars$intervalHigh[3])
start_time <- Sys.time()
for (e in 1:modelVars$experiments) {
  print(paste('Experiment no.',e))
  if (e > 995) {
    indeps = list(X=runif(modelVars$n), W=runif(modelVars$n))
  } else {
    indeps = gibbs(modelVars$a, modelVars$x_0, modelVars$w_0, modelVars$n)
  }
  noise <- rnorm(modelVars$n, 0, 0.01)
  Y <- dependent(modelVars$a, 2.2*indeps$X, indeps$W, noise, modelVars$n)
  
  psiWk[1,] = legendrePsi(indeps$W, 1)
  psiWk[2,] = legendrePsi(indeps$W, 2)
  psiWk[3,] = legendrePsi(indeps$W, 3)
  psiXj[1,] = legendrePsi(indeps$X, 1)
  psiXj[2,] = legendrePsi(indeps$X, 2)
  psiXj[3,] = legendrePsi(indeps$X, 3)
  psiJx1[1,1:100] = legendrePsi(xl_grid_points[1,],1)
  psiJx1[2,1:100] = legendrePsi(xl_grid_points[1,],2)
  psiJx1[3,1:100] = legendrePsi(xl_grid_points[1,],3)
  psiJx2[1,1:100] = legendrePsi(xl_grid_points[2,],1)
  psiJx2[2,1:100] = legendrePsi(xl_grid_points[2,],2)
  psiJx2[3,1:100] = legendrePsi(xl_grid_points[2,],3)
  psiJx3[1,1:100] = legendrePsi(xl_grid_points[3,],1)
  psiJx3[2,1:100] = legendrePsi(xl_grid_points[3,],2)
  psiJx3[3,1:100] = legendrePsi(xl_grid_points[3,],3)
  cjk[1,1] = sum(psiXj[1,]*psiWk[1,])
  cjk[1,2] = sum(psiXj[2,]*psiWk[1,])
  cjk[1,3] = sum(psiXj[3,]*psiWk[1,])
  cjk[2,1] = sum(psiXj[1,]*psiWk[2,])
  cjk[2,2] = sum(psiXj[2,]*psiWk[2,])
  cjk[2,3] = sum(psiXj[3,]*psiWk[2,])
  cjk[3,1] = sum(psiXj[1,]*psiWk[3,])
  cjk[3,2] = sum(psiXj[2,]*psiWk[3,])
  cjk[3,3] = sum(psiXj[3,]*psiWk[3,])
  if (modelVars$J_n > 3) {
    psiWk[4,] = legendrePsi(indeps$W, 4)
    psiXj[4,] = legendrePsi(indeps$X, 4)
    psiJx1[4,1:100] = legendrePsi(xl_grid_points[1,],4)
    psiJx2[4,1:100] = legendrePsi(xl_grid_points[2,],4)
    psiJx3[4,1:100] = legendrePsi(xl_grid_points[3,],4)
    cjk[4,1] = sum(psiXj[1,]*psiWk[4,])
    cjk[4,2] = sum(psiXj[2,]*psiWk[4,])
    cjk[4,3] = sum(psiXj[3,]*psiWk[4,])
    cjk[4,4] = sum(psiXj[4,]*psiWk[4,])
    cjk[1,4] = sum(psiXj[4,]*psiWk[1,])
    cjk[2,4] = sum(psiXj[4,]*psiWk[2,])
    cjk[3,4] = sum(psiXj[4,]*psiWk[3,])
    if (modelVars$J_n == 5) {
      psiWk[5,] = legendrePsi(indeps$W, 5)
      psiXj[5,] = legendrePsi(indeps$X, 5)
      psiJx1[5,1:100] = legendrePsi(xl_grid_points[1,],5)
      psiJx2[5,1:100] = legendrePsi(xl_grid_points[2,],5)
      psiJx3[5,1:100] = legendrePsi(xl_grid_points[3,],5)
      cjk[5,1] = sum(psiXj[1,]*psiWk[5,])
      cjk[5,2] = sum(psiXj[2,]*psiWk[5,])
      cjk[5,3] = sum(psiXj[3,]*psiWk[5,])
      cjk[5,4] = sum(psiXj[4,]*psiWk[5,])
      cjk[5,5] = sum(psiXj[5,]*psiWk[5,])
      cjk[1,5] = sum(psiXj[5,]*psiWk[1,])
      cjk[2,5] = sum(psiXj[5,]*psiWk[2,])
      cjk[3,5] = sum(psiXj[5,]*psiWk[3,])
      cjk[4,5] = sum(psiXj[5,]*psiWk[4,])
    }
  }
  g_real <- Y - noise
  for (J_n in 3:modelVars$J_n) {
    gSieve <- sieve_optimum(psiXj[1:J_n,], psiWk[1:J_n,], cjk[1:J_n,1:J_n], Y, J_n)
    grid_sieves[1,] = sieve(psiJx1[1:J_n,], psiXj[1:J_n,], psiWk[1:J_n,], cjk[1:J_n,1:J_n], Y, J_n)
    grid_sieves[2,] = sieve(psiJx2[1:J_n,], psiXj[1:J_n,], psiWk[1:J_n,], cjk[1:J_n,1:J_n], Y, J_n)
    grid_sieves[3,] = sieve(psiJx3[1:J_n,], psiXj[1:J_n,], psiWk[1:J_n,], cjk[1:J_n,1:J_n], Y, J_n)
    #### Awesome graph ######################################################
    plot(sort(indeps$X), sort(Y), xlab="X", ylab="Y",
         main=paste("Nonparametric IV regression with J_n=",J_n),
         xlim=c(min(indeps$X), max(indeps$X)),
         ylim=c(min(Y), max(c(max(Y),max(gSieve))))
         )
    legend(x="bottomright",c("True Curve", "Sieve Estimator"),
           lty=c(1,1),
           lwd=c(1,1),
           col=c("red", "green"),
           bty="n")
    lines(sort(indeps$X), sort(g_real), col="red")
    lines(sort(indeps$X), sort(gSieve), col="green")
    if (J_n == 3) {
      #For confidence intervals we see for each point the distance/for bands the uniform
      # ie ||vecLow - gSieve|| < ||vecLow - vecHigh|| and ||vecHigh - gSieve|| < ||vecLow - vecHigh||
      allMonte[1:6,1:6] = allMonte[1:6,1:6] +
        monteCarlo(psiJx1[1:3,], psiJx2[1:3,], psiJx3[1:3,], psiWk[1:3,], psiXj[1:3,],
                   cjk[1:3,1:3], gSieve, grid_sieves, g_real,
                   Y, indeps$X, indeps$W, xl_grid_points, J_n, modelVars$n,
                   modelVars$intervalLow, modelVars$intervalHigh,
                   modelVars$boot)
    } else if (J_n == 4) {
      allMonte[7:12,1:6] = allMonte[7:12,1:6] +
        monteCarlo(psiJx1[1:4,], psiJx2[1:4,], psiJx3[1:4,], psiWk[1:4,], psiXj[1:4,],
                   cjk[1:4,1:4], gSieve, grid_sieves, g_real,
                   Y, indeps$X, indeps$W, xl_grid_points, J_n, modelVars$n,
                   modelVars$intervalLow, modelVars$intervalHigh,
                   modelVars$boot)
    } else {
      allMonte[13:18,1:6] = allMonte[13:18,1:6] +
        monteCarlo(psiJx1[1:5,], psiJx2[1:5,], psiJx3[1:5,], psiWk[1:5,], psiXj[1:5,],
                   cjk[1:5,1:5], gSieve, grid_sieves, g_real,
                   Y, indeps$X, indeps$W, xl_grid_points, J_n, modelVars$n,
                   modelVars$intervalLow, modelVars$intervalHigh,
                   modelVars$boot)
    }
  }
}
end_time <- Sys.time()
print(allMonte)
print(end_time-start_time)
