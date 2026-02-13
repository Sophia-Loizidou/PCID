cusum_circle <- function(data){
  if (!(is.numeric(data))){
    stop("The input in `data' should be a numeric vector containing the data
         for which the CUSUM function will be calculated.")
  }
  n <- length(data)
  
  lr_test <- rep(NA, n-1)
  
  for(b in 1:(n - 1)){
    c1 <- sum(cos(data[1:b])); c2 <- sum(cos(data[(b+1):n]))
    s1 <- sum(sin(data[1:b])); s2 <- sum(sin(data[(b+1):n]))
    
    lr_test[b] <- sqrt(c1^2 + s1^2) + sqrt(c2^2 + s2^2) - sqrt((c1+c2)^2 + (s1+s2)^2)
  }
  return(lr_test)
}

cusum_circle_alternative <- function(data){
  if (!(is.numeric(data))){
    stop("The input in `data' should be a numeric vector containing the data
         for which the CUSUM function will be calculated.")
  }
  n <- length(data)
  
  # mean under H0
  mu0 <- as.numeric(atan2(sum(sin(data)), sum(cos(data))))
  
  lr_test <- rep(NA, n-1)
  
  for(b in 1:(n - 1)){
    # mean under H1
    mu1 <- as.numeric(atan2(sum(sin(data[1:b])), sum(cos(data[1:b]))))
    mu2 <- as.numeric(atan2(sum(sin(data[(b+1):n])), sum(cos(data[(b+1):n]))))
    
    lr_test[b] <- 2*(sum(cos(data[1:b] - mu1)) + sum(cos(data[(b+1):n] - mu2)) - sum(cos(data - mu0)))
  }
  return(lr_test)
}

s_e_points <- function(r, l, s, e) {
  r <- sort(r)
  l <- sort(l, decreasing = TRUE)
  if (s > e){
    stop("s should be less than or equal to e")
  }
  if (!(is.numeric(c(r, l, s, e))) | (r[1] <= 0) | (l[length(l)] <= 0) | s <= 0 | e <= 0){
    stop("The input arguments must be positive integers")
  }
  if (any(abs(r - round(r)) > .Machine$double.eps ^ 0.5)){
    warning("The input for r should be a vector of positive integers. If there is at least a positive real
            number then the integer part of that number is used.")
  }
  if (any(abs(l - round(l)) > .Machine$double.eps ^ 0.5)){
    warning("The input for l should be a vector of positive integers. If there is at least a positive real
            number then the integer part of that number is used.")
  }
  if (abs(s - round(s)) > .Machine$double.eps ^ 0.5){
    warning("The input for s should be a positive integer. If it is a positive real
            number then the integer part of that number is used.")
  }
  if (abs(e - round(e)) > .Machine$double.eps ^ 0.5){
    warning("The input for e should be a positive integer. If it is a positive real
            number then the integer part of that number is used.")
  }
  r <- as.integer(r)
  l <- as.integer(l)
  e <- as.integer(e)
  s <- as.integer(s)
  e_points <- unique(c(r[which( (r > s) & (r < e))], e))
  s_points <- unique(c(l[which( (l > s) & (l < e))], s))
  return(list(e_points = e_points, s_points = s_points))
}

perm_func <- function(x, B = 1000, test_stat_value, alpha){
  if(factorial(length(x)) < B){
    cpt_yes = F
  } else {
    cut_off <- ceiling(B*alpha)
    ind <- 0; i <- 1
    while((ind < cut_off) && (i <= B)){
      x_perm <- sample(x)
      max_stat <- max(abs(cusum_circle(x_perm)))
      if(test_stat_value <= max_stat){ind <- ind + 1}
      i <- i + 1
    }
    cpt_yes = ifelse(ind == cut_off, F, T)
  }
  return(cpt_yes)
}

PCID_no_win <- function(x, s = 1, e = length(x), points = 5, k_l = 1, k_r = 1, 
                          B = 1000, alpha = 0.005) {
  points <- as.integer(points)
  
  ## draw again the intervals every time
  r_e_points <- tryCatch(seq(s-1+points, e, points), error = function(e) c(s))
  l_e_points <- c(tryCatch(seq(e - points + 1, s, -points), error = function(e) c()), e)

  chp <- 0
  detected_left <- 0
  if (e - s <= 1) {
    cpt <- 0
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- s_e_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ( (chp == 0) & (k_r < min(k_l, rur))) {
        x_temp_r <- x[s:right_points[k_r]]
        ipcr <- cusum_circle(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
        if (perm_func(x_temp_r, B = B, test_stat_value = CUSUM_r[k_r], alpha = alpha)) {
          chp <- pos_r[k_r]
          detected_left <- 0 ## flag if the cpt is detected from left- or right- expanding intervals
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ( (chp == 0) & (k_l < min(k_r, lur))) {
        x_temp_l <- x[left_points[k_l]:e]
        ipcl <- cusum_circle(x_temp_l)
        pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
        CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
        if (perm_func(x_temp_l, B = B, test_stat_value = CUSUM_l[k_l], alpha = alpha)) {
          chp <- pos_l[k_l]
          detected_left <- 1 ## flag if the cpt is detected from left- or right- expanding intervals
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        x_temp_r <- x[s:right_points[k_r]]
        ipcr <- cusum_circle(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
        if (perm_func(x_temp_r, B = B, test_stat_value = CUSUM_r[k_r], alpha = alpha)) {
          chp <- pos_r[k_r]
          detected_left <- 0 ## flag if the cpt is detected from ledt- or right- expanding intervals
        } else {
          x_temp_l <- x[left_points[k_l]:e]
          ipcl <- cusum_circle(x_temp_l)
          pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
          CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
          if (perm_func(x_temp_l, B = B, test_stat_value = CUSUM_l[k_l], alpha = alpha)) {
            chp <- pos_l[k_l]
            detected_left <- 1 ## flag if the cpt is detected from ledt- or right- expanding intervals
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (detected_left == 1) {
        r <- PCID_no_win(x, s = s, e = chp, points = points, k_r = k_r, 
                           k_l = 1, B = B, alpha = alpha)
      } else {
        r <- PCID_no_win(x, s = chp + 1, e = e, points = points, k_r = 1, 
                           k_l = max(1, k_l - 1), B = B, alpha = alpha)
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  return(sort(cpt))
}

PCID_win <- function(x, s = 1, e = length(x), points = 5, win_length = 500, k_l = 1, k_r = 1, 
                       B = 1000, alpha = 0.005, FDR = 0.01, verbose = F, override_default = F){
  l <- length(x)
  
  if(l > win_length){
    end_points <- unique(c(seq(0,l,win_length)[-1], l))
    start_points <- c(1, end_points[-length(end_points)]+1)
    num_intervals <- length(end_points)
    r <- list(); cpt <- c()
    
    if(!is.null(FDR)){
      FDR_used <- round(1 - (1-FDR)^(1/num_intervals), digits = 3)
      if(verbose){
        warning("For the overall type I error to be ", FDR, ", we use ", 
                FDR_used, " to account for multiple testing.")
      }
      FDR <- FDR_used
    }
    
    ## check all disjoint intervals
    for(i in 1:num_intervals){
      length_temp <- end_points[i] - start_points[i] + 1
      
      alpha_B <- choose_alpha(l = length_temp, FDR = FDR, B = B, alpha = alpha, verbose = verbose, override_default = override_default)
      alpha_used <- alpha_B[1]; B_used <- alpha_B[2]
      
      r[[i]] <- PCID_no_win(x = x[start_points[i]:end_points[i]], s = 1, e = length_temp, 
                              points = points, k_l = k_l, k_r = k_r, B = B_used, alpha = alpha_used)
      r[[i]] <- r[[i]] + start_points[i] - 1
    }
    
    ## check the cut-off point
    for(i in 1:(num_intervals-1)){
      if(length(r[[i]]) == 0){
        left_point <- start_points[i] + floor(win_length/2) - 1
      } else {
        left_point <- max(start_points[i] + floor(win_length/2) - 1, r[[i]][length(r[[i]])]+1)
      }
      if(length(r[[i+1]]) == 0){
        right_point <- min(start_points[i+1] + floor(win_length/2) - 1, end_points[i+1])
      } else {
        right_point <- min(min(start_points[i+1] + floor(win_length/2) - 1, end_points[i+1]), r[[i+1]][1])
      }
      
      alpha_B <- choose_alpha(l = right_point - left_point + 1, FDR = FDR, B = B, alpha = alpha, verbose = verbose, override_default = override_default)
      alpha_used <- alpha_B[1]; B_used <- alpha_B[2]
      
      x_temp <- x[left_point:right_point]
      ipc <- cusum_circle(x_temp)
      pos <- which.max(abs(ipc)) + left_point - 1
      CUSUM_temp <- abs(ipc[pos - left_point + 1])
      if (perm_func(x_temp, B = B_used, test_stat_value = CUSUM_temp, alpha = alpha_used)) {
        r[[i]] <- c(r[[i]], pos)
      }
    }
    for(i in 1:num_intervals){
      cpt <- c(cpt, r[[i]])
    }
    
  } else {
    alpha_B <- choose_alpha(l = e-s+1, FDR = FDR, B = B, alpha = alpha, verbose = verbose, override_default = override_default)
    alpha_used <- alpha_B[1]; B_used <- alpha_B[2]
    cpt <- PCID_no_win(x = x, s = s, e = e, points = points, k_l = k_l, k_r = k_r, 
                      B = B_used, alpha = alpha_used)
  }
  return(cpt)
}

alpha_T_table <- data.frame(
  alpha = c(
    # T = 50
    0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0005, 0.0001,
    # T = 100
    0.01, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0005, 0.0001,
    # T = 150
    0.005, 0.003, 0.002, 0.001, 0.0005, 0.0001,
    # T = 200
    0.005, 0.002, 0.001, 0.0005, 0.0003, 0.0002, 0.0001,
    # T = 250
    0.002, 0.001, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001,
    # T = 300
    0.002, 0.001, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001,
    # T = 350
    0.002, 0.001, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001,
    # T = 400
    0.002, 0.001, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001,
    # T = 450
    0.002, 0.001, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001,
    # T = 500
    0.002, 0.001, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001
  ),
  length = c(
    rep(50, 12),
    rep(100, 8),
    rep(150, 6),
    rep(200, 7),
    rep(250, 7),
    rep(300, 7),
    rep(350, 7),
    rep(400, 7),
    rep(450, 7),
    rep(500, 7)
  ),
  FDR = c(
    # T = 50
    0.083, 0.078, 0.066, 0.058, 0.046, 0.041, 0.035, 0.029, 0.008, 0.006, 0.002, 0.000,
    # T = 100
    0.149, 0.083, 0.069, 0.051, 0.037, 0.011, 0.005, 0.001,
    # T = 150
    0.097, 0.055, 0.032, 0.017, 0.010, 0.003,
    # T = 200
    0.131, 0.057, 0.037, 0.017, 0.013, 0.004, 0.003,
    # T = 250
    0.056, 0.034, 0.019, 0.014, 0.012, 0.010, 0.002,
    # T = 300
    0.070, 0.041, 0.021, 0.017, 0.013, 0.009, 0.003,
    # T = 350
    0.068, 0.044, 0.019, 0.018, 0.013, 0.008, 0.007,
    # T = 400
    0.076, 0.045, 0.025, 0.021, 0.013, 0.006, 0.003,
    # T = 450
    0.081, 0.048, 0.020, 0.025, 0.013, 0.009, 0.005,
    # T = 500
    0.096, 0.057, 0.031, 0.028, 0.020, 0.009, 0.002
  )
)


choose_alpha <- function(l, FDR = NULL, B = NULL, alpha = NULL, B_min = 1000, 
                         verbose = F, override_default = F){
  l_round <- round(l/50)* 50 ## round to the closest 50
  if(l_round == 0){l_round <- 50}
  
  if(!override_default){
    if(is.null(alpha)){
      if(l_round > 500){
        stop("The appropriate values for such large intervals have not been calculated.
           To proceed, please either set the value alpha or set the value of 'win_length' to at most 500.")
      }
      ## find the closest value from the table
      alpha_T_table_temp = alpha_T_table[which(alpha_T_table$length == l_round),]
      alpha = alpha_T_table_temp$alpha[which.min(abs(alpha_T_table_temp$FDR - FDR))]
      
      if(is.null(B)){
        B <- max(B_min, 10^nchar(strsplit(format(alpha, scientific = FALSE), "\\.")[[1]][2]))
        if(B > 1000){
          if(verbose){
            warning("For type I error ", FDR, " and length ", l, ", alpha = ", alpha, 
                    " should be used with B = ", B, 
                    ". In order to make the algorithm more efficient, alpha = ", 0.001,
                    ", B = ", 1000, " will be used. If you want to use alpha = ", alpha, 
                    ", set override_default = T.")
          }
          B <- 1000
          alpha <- 0.001
        } else {
          if(verbose){
            warning("For type I error ", FDR, " and length ", l, ", alpha = ", alpha, 
                    " will be used with B = ", B)
          }
        }
      } else { ## B not NULL
        ## force them to have B*alpha to be an integer
        cut_off = B*alpha
        if(abs(cut_off - round(cut_off)) > .Machine$double.eps ^ 0.5){
          B <- max(B_min, 10^nchar(strsplit(format(alpha, scientific = FALSE), "\\.")[[1]][2]))
          if(B > 1000){
            if(verbose){
              warning("For type I error ", FDR, " and length ", l, ", alpha = ", alpha, 
                      " should be used with B = ", B, 
                      ". In order to make the algorithm more efficient, alpha = ", 0.001,
                      ", B = ", 1000, " will be used. If you want to use alpha = ", alpha, 
                      ", set override_default = T.")
            }
            B <- 1000
            alpha <- 0.001
          } else {
            if(verbose){
              warning("For type I error ", FDR, " and length ", l, ", alpha = ", alpha, 
                      " will be used with B = ", B)
            }
          }
        }
      }
    } else { ## alpha not NULL
      stop("If parameter alpha is specified, override_default must be set to TRUE.")
    }
  } else { ## if override_default == T
    if(is.null(alpha)){
      if(l_round > 500){
        stop("The appropriate values for such large intervals have not been calculated.
           To proceed, please either set the value alpha or set the value of 'win_length' to at most 500.")
      }
      ## find the closest value from the table
      alpha_T_table_temp = alpha_T_table[which(alpha_T_table$length == l_round),]
      alpha = alpha_T_table_temp$alpha[which.min(abs(alpha_T_table_temp$FDR - FDR))]
      
      if(is.null(B)){
        B <- max(B_min, 10^nchar(strsplit(format(alpha, scientific = FALSE), "\\.")[[1]][2]))
        if(verbose){
          warning("For type I error ", FDR, " and length ", l, ", alpha = ", alpha, 
                  " will be used with B = ", B)
        }
      } else { ## B not NULL
        ## force them to have B*alpha to be an integer
        cut_off = B*alpha
        if(abs(cut_off - round(cut_off)) > .Machine$double.eps ^ 0.5){
          B <- max(B_min, 10^nchar(strsplit(format(alpha, scientific = FALSE), "\\.")[[1]][2]))
          if(verbose){
            warning("B*alpha must be an integer. For type I error ", FDR, 
                    " and length ", l, ", alpha must be ", alpha," so B will be set to ", B)
          }
        }
      }
    } else { ## alpha not NULL
      warning("When you choose a value for alpha, the overall type I error will not be controlled.")
      if(is.null(B)){
        B <- max(B_min, 10^nchar(strsplit(format(alpha, scientific = FALSE), "\\.")[[1]][2]))
        if(verbose){
          warning("alpha = ", alpha, " will be used with B = ", B)
        }
      } else { ## B not NULL
        ## force them to have B*alpha to be an integer
        cut_off = B*alpha
        if(abs(cut_off - round(cut_off)) > .Machine$double.eps ^ 0.5){
          B <- max(B_min, 10^nchar(strsplit(format(alpha, scientific = FALSE), "\\.")[[1]][2]))
          if(verbose){
            warning("B*alpha must be an integer. Based on your alpha, B will be set to ", B)
          }
        }
      }
    }
  }
  if(B > 5000){
    if(verbose){
      print("B is very large. The computation will be slow.")
    }
  }
  
  return(c(alpha, B))
}

# The user should specify the FDR
# This will determine B and alpha
# If user specifies alpha then the B will be chosen appropriately
PCID <- function(x, FDR = 0.01, s = 1, e = length(x), points = 5, win_length = 500, k_l = 1, k_r = 1, 
                   B = NULL, alpha = NULL, win = T, verbose = F, override_default = F){
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
           which you would like to find change-points.")
  }
  if ((points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    stop("The input for `points' should be a positive integer.")
  }
  
  if(win){
    cpt <- PCID_win(x = x, s = s, e = e, points = points, win_length = win_length, 
                      k_l = k_l, k_r = k_r, B = B, alpha = alpha, FDR = FDR, 
                      verbose = verbose, override_default = override_default)
  } else {
    alpha_B <- choose_alpha(l = e-s+1, FDR = FDR, B = B, alpha = alpha, verbose = verbose, override_default = override_default)
    alpha <- alpha_B[1]; B <- alpha_B[2]
    
    cpt <- PCID_no_win(x = x, s = s, e = e, points = points,
                         k_l = k_l, k_r = k_r, B = B, alpha = alpha)
  }
  return(cpt)
}

## Examples
# library(circular)
# set.seed(1)
# x <- rvonmises(800, circular(0), 2) + c(rep(0, 300), rep(2, 300), rep(0, 200))
# cpt <- PCID(x, FDR = 0.01, verbose = F)
