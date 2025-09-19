library(logger)

f <- function(theta,mu){
  dnorm(theta,mu)
}

f.grad_log <- function(theta,mu){
  -(theta-mu)
}


pdmp.step <- function(x,theta,v,r,delta=0.01){
  t <- 0 
  lambda <- function(v,theta,delta,x,r){
    max(v*(theta+v*delta-x),0)+r
  }
  while(T){
    upper.bound <- max(lambda(v,theta+t*v,0,x,r),lambda(v,theta+t*v,delta,x,r))
    n.p <- rpois(1,upper.bound*delta)
    if (n.p==0){
      t <- t+delta
    } else {
      ts <- runif(n.p,0,delta)
      a <- runif(n.p)<sapply(ts,function(u){lambda(v,theta+t*v,u,x,r)/upper.bound})
      if (sum(a)>0){
        ts <- ts[a]
        return(t+min(ts))
      } else {
        t <- t+delta
      }
    }
  }
  return(t)
}

pdmp <- function(tt,x,theta,v,r){
  t <- 0
  row_guess <- 5000
  skeleton <- matrix(nrow=row_guess,ncol=3)
  skeleton[1,]  <- c(0,theta,v)
  turns <- 1
  while(t<tt){
    turns <- turns+1
    t.step <- min(pdmp.step(x,theta,v,r),tt-t)
    theta.step <- theta+v*t.step
    t <- t+t.step
    theta <- theta.step
    v <- -v
    if (turns<=row_guess) {
      skeleton[turns,] <- c(t,theta.step,v)
    }else {
      skeleton <- rbind(skeleton,c(t,theta.step,v))
    }
    
  }
  return(skeleton[1:turns,])
}

discretize_pdmp <- function(path,tt,n){
  points <- tt*c(1:n)/n
  cur_id <- 1
  path_id <- 1
  out_sample <- matrix(nrow=n,ncol=1)
  while(cur_id<=n){
    if ((points[cur_id]>path[path_id,1])&(points[cur_id]<=path[path_id+1,1])){
      out_sample[cur_id] <- (points[cur_id]-path[path_id,1])*path[path_id,3]+path[path_id,2]
      cur_id <- cur_id+1
      
    } else {
      path_id <- path_id+1
    } 
  }
  return(out_sample)
}


set.seed(1)
n <- 50

rs <- c(0,1,5,10,20,30)
m <- array(dim=c(length(rs),n,2))
run.times <- matrix(nrow=n,ncol=length(rs))
x=1
true_m1 <- x
true_m2 <- x^2+1

pdmp.time <- 1e2
n.sample <- 1e5

for (r.id in 1:length(rs)){
  log_info("r: {rs[r.id]}")
  for (it in 1:n){
    start.time <- Sys.time()
    path <- pdmp(tt=pdmp.time,x=x,theta=0,v=1,r=rs[r.id])
    end.time <- Sys.time()
    run.times[it,r.id] <- end.time-start.time

    sample <- discretize_pdmp(path,pdmp.time,n.sample)

    m[r.id,it,] <- c(mean(sample),mean(sample^2))
  }
}

log_info("sampling complete.")

log_info("x: {x}")
log_info("number of repeats: {n}")
log_info("pdmp sample time: {pdmp.time}")
log_info("sample discretization size: {n.sample}")
log_info("vector of rs used:")
print(rs)
log_info("MSE for first moment:")
mse.m1 <- rowSums((m[,,1]-true_m1)^2)/n
print(mse.m1)
log_info("MSE for second moment:")
mse.m2 <- rowSums((m[,,2]-true_m2)^2)/n
print(mse.m2)
log_info("mean run times:")
mean.run.time <- colMeans(run.times)
print(mean.run.time)

out <- data.frame(
  cbind(
    rs,
    mse.m1,
    mse.m2,
    mean.run.time
  )
)
write.csv(out,"n50t100.csv",row.names = F)
