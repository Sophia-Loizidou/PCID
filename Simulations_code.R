source('~/Dropbox/Sophia-Andreas/ID_circle/PCID.R')

options(expressions = 500000)

# circular 'mse'
# x = noisy signal
# cpt = detected cpt
# signal = true signal
circular_mse <- function(x, cpt, signal){
  require(circular)
  n <- length(x); x <- circular(x)
  len.cpt <- length(cpt)
  if (len.cpt) cpt <- sort(cpt)
  beg <- endd <- rep(0, len.cpt+1)
  beg[1] <- 1
  endd[len.cpt+1] <- n
  if (len.cpt) {
    beg[2:(len.cpt+1)] <- cpt+1
    endd[1:len.cpt] <- cpt
  }
  means <- rep(0, len.cpt+1)
  for (i in 1:(len.cpt+1)) means[i] <- mean.circular(x[beg[i]:endd[i]])
  
  circ_means <- rep(means, endd-beg+1)
  dist <- sqrt(2 - 2*cos(circ_means - signal))
  
  return(dist)
}

# Adjusted rank index
# n = length of signal
# changepoints = true changepoints 
label_timeserie<-function(n,changepoints){ # function to be used for ARI
  df = data.table::data.table(n=c(1:n))
  df[,class := findInterval(n,changepoints,left.open = TRUE)]
  return(df)
}

# Computes ARI
# n = length of timeseries
# true = true changepoints
# est = estimated changepoints

ARI_compute<- function(n,true,est){
  require(MixGHD) ## for ARI
  
  df_true=label_timeserie(n,true)
  df_est=label_timeserie(n,est)
  value=ARI(df_true$class,df_est$class)
  return(value)
}
# example
#ARI_compute(100,c(10,20),c(11,21))

sim <- function(x, B, points, alpha, FDR, verbose, signal, win = T, win_length = 500, override_default){
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric",det_time="numeric", time="numeric", mean="numeric"), 
           prototype(cpt=numeric(0), nocpt=0, det_time=numeric(0), time=numeric(0), mean=numeric(0)))
  
  # if(verbose){print("PCID")}
  PCID <- new("cpt.est")
  start_time <- Sys.time()
  z <- PCID(x, B = B, alpha = alpha, points = points, win = win, win_length = win_length, 
              FDR = FDR, verbose = verbose, override_default = override_default)
  end_time <- Sys.time()
  if(length(z) == 0){PCID@cpt = 0
  PCID@nocpt = 0}
  else{PCID@cpt <- z
  PCID@nocpt <- length(z)}
  PCID@time <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  PCID@mean <- circular_mse(x, PCID@cpt, signal)
  
  list(PCID = PCID)
}

sim.study <- function(signal, true.cpt=NULL, kappa, m = 100, seed = NULL, points = 5, 
                      B = NULL, alpha = NULL, verbose = FALSE, noise = 'vonMises', 
                      win = T, win_length = 500, FDR = 0.05, override_default = F) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",
                                      cptall="numeric",dnc="numeric",ari="numeric",
                                      mse="numeric", time= "numeric", det_delay="numeric"), 
           prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  require(circular)
  
  # if(is.null(points)){
  #   points <- min(ceiling(length(signal)/10), 50)
  # }
  # cat("Parameter 'points' is taken to be", points, "'B' is", B, "and 'alpha' is", alpha, "\n")
  
  ts <- list()
  
  no.of.cpt <- sum(abs(diff(signal)) > 0)
  n <- length(signal)
  ns <- max(c(diff(true.cpt), length(signal)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    if(noise == 'vonMises'){
      set.seed(i)
      epsilon <- rvonmises(n, circular(0), kappa)
    } else if (noise == 'wrpCauchy'){
      if((kappa < 0) | (kappa>1)){stop("'kappa' should be between 0 and 1 if noise = 'wrpCauchy'")}
      epsilon <- rwrappedcauchy(n, circular(0), kappa)
    } else if (noise == 'wrpNormal'){
      if((kappa < 0) | (kappa>1)){stop("'kappa' should be between 0 and 1 if noise = 'wrpNormal'")}
      epsilon <- rwrappednormal(n, circular(0), kappa)
    } else {
      stop("'noise' should be 'vonMises', 'wrpCauchy' or 'wrpNormal'")
    }
    x <- (signal + epsilon) %% (2*pi)
    ts[[i]] <- x
    
    if(i == 1){
      ts.plot(x)
    }
    
    if(i == 1){
      verbose = T
    } else {
      verbose = F
    }
    
    est <- sim(x = x, alpha = alpha, B = B, points=points, verbose = verbose, 
               signal = signal, win = win, win_length = win_length, FDR = FDR,
               override_default = override_default)
    
    if(i == 1){
      for(j in names(est)){
        eval(parse(text=paste(j, " <- new('est.eval')",sep="")))
      }
    }
    
    for(j in names(est)){
      eval(parse(text=paste(j, "@dnc[", i, "] <- est$", j, "@nocpt - no.of.cpt",sep="")))
      eval(parse(text=paste(j, "@mse[", i, "] <- mean((est$", j, "@mean - signal)^2)",sep="")))
      # eval(parse(text=paste(j, "@mse[", i, "] <- mean((est$", j, "@mean - signal)^2)",sep="")))
      eval(parse(text=paste(j, "@cpt[[", i, "]] <- est$", j, "@cpt",sep="")))
      eval(parse(text=paste(j, "@diff <- abs(matrix(est$", j, "@cpt,nrow=no.of.cpt,ncol=length(est$", j, "@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$", j, "@cpt),byr=F))",sep="")))
      eval(parse(text=paste(j, "@dh[i] <- max(apply(", j, "@diff,1,min),apply(", j, "@diff,2,min))/ns",sep="")))
      eval(parse(text=paste(j, "@ari[i] <- ARI_compute(n,true.cpt,est$", j, "@cpt)",sep="")))
      eval(parse(text=paste(j, "@time[i] <- est$", j, "@time",sep="")))
    }
    
    gc()
  }
  
  result <- list(ts = ts)
  for(i in names(est)){
    eval(parse(text=paste("result <- c(result, ", i, " = ", i, ")",sep="")))
  }
  return(result)
}

make_df <- function(x, decimals = 3){
  results <- data.frame('Method' = names(x)[-1],
                        # 'MSE' = rep(NA),
                        'd_h' = rep(NA),
                        'ari' = rep(NA),
                        'time' = rep(NA))
  
  for(i in 1:dim(results)[1]){
    # results$MSE[i] <- eval(parse(text = paste("signif(mean(x$", results$Method[i],"@mse), decimals)", sep='')))
    
    if(!is.null(eval(parse(text = paste("attributes(x$", results$Method[i],")$dh", sep=''))))){
      results$d_h[i] <- eval(parse(text = paste("signif(mean(x$", results$Method[i],"@dh), decimals)", sep='')))
    }
    
    if(!is.null(eval(parse(text = paste("attributes(x$", results$Method[i],")$ari", sep=''))))){
      results$ari[i] <- eval(parse(text = paste("signif(mean(x$", results$Method[i],"@ari), decimals)", sep='')))
    }
    
    if(!is.null(eval(parse(text = paste("attributes(x$", results$Method[i],")$time", sep=''))))){
      results$time[i] <- eval(parse(text = paste("signif(mean(x$", results$Method[i],"@time), decimals)", sep='')))
    }
  }
  
  return(results)
}

cpts_df <- function(x, breaks){
  methods <- names(x)[-1]
  for(i in 1:length(methods)){
    eval(parse(text = paste("x", i,"<- cut(x$", methods[i],"@dnc, breaks = breaks)", sep='')))
  }
  cname <- c('Method', levels(x1))
  df <- data.frame(matrix(ncol = length(cname), nrow=length(methods)))
  colnames(df) <- cname
  df$Method <- methods
  
  for(i in 2:length(cname)){
    for(j in 1:length(methods)){
      eval(parse(text = paste("df[", j,",", i, "] <- sum(x", j,"== cname[", i, "])", sep='')))
    }
  }
  return(df)
}

make_table <- function(x, breaks, decimals = 3){
  df1 <- make_df(x, decimals = 3)
  df2 <- cpts_df(x, breaks)
  df <- cbind(df2, df1[,-1])
  return(df)
}

nice_timeseries_plot <- function(data){
  
  require(ggplot2)
  
  df <- data.frame(
    x = 1:length(data),
    value = data
  )
  
  ggplot(df, aes(x = x, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.7) +
    geom_point(color = "steelblue", size = 1.2) +
    scale_x_continuous(                 # Numeric x-axis
      breaks = waiver(),                # Auto breaks or specify e.g. seq(0, 100, 10)
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "Time",                 # Or "Index", "Iteration", etc.
      y = "Observed signal",
      # title = "Time series plot"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
      axis.title.y = element_text(margin = margin(r = 8)),
      # panel.border = element_rect(
      #   color = "black",     # Border color
      #   fill = NA,           # No fill (transparent)
      #   size = 0.8           # Thickness (adjust as needed)
      # )
    )
  
}
