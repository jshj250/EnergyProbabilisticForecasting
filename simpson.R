simpson <- function(q, fhat)
{
  density  <- FALSE
  gridsize <- length(fhat$eval.points)  
  ## Use Simpson's rule to compute numerical integration
  simp.rule <- rep(0, gridsize-1)
  for (i in 1:(gridsize-1))
  {
    del <- fhat$eval.points[i+1] - fhat$eval.points[i]
    simp.rule[i] <- min(fhat$estimate[i], fhat$estimate[i+1])*del + 1/2*abs(fhat$estimate[i+1] - fhat$estimate[i])*del 
    
  }
  ## add last incomplete trapezoid 
  q.ind <- findInterval(x=q, vec=fhat$eval.points)
  q.prob <- rep(0, length(q))
  i <- 0
  for (qi in q.ind)
  {
    i <- i+1
    
    if (qi==0)
      q.prob[i] <- 0
    else if (qi < gridsize)
    {
      ## linearly interpolate kde 
      fhat.estqi <- (fhat$est[qi+1] - fhat$est[qi])/(fhat$eval[qi+1] - fhat$eval[qi]) * (q[i] - fhat$eval[qi]) + fhat$est[qi]
      delqi <- q[i] - fhat$eval[qi] 
      
      simp.ruleqi <- min(fhat.estqi, fhat$est[qi])*delqi + 1/2*abs(fhat.estqi - fhat$est[qi])*delqi
      q.prob[i] <- sum(simp.rule[1:qi]) + simp.ruleqi
    }
    else
    {
      if (density) q.prob[i] <- 1
      else q.prob[i] <- sum(simp.rule) 
    }
  }
  if (density) q.prob[q.prob>=1] <- 1
  return(q.prob)
}