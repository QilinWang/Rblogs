accel <-
  function (bb, tt,trun=.001, pct = .5,zc=5,sw=0)
  {
    ##  bb is B by p matrix of sufficient vectors, ith row "beta[i]"
    ##  tt is the B-vector of bootstrap values thetahat*[i]
    ## trun does .truncation on d, tt
    ## pct is for the t. calcs
    ## zc for the smoothed z0-type estimate az of 'a'
    
    bb=as.matrix(bb); B=length(tt)
    if (pct > 0) {
      W = scale(bb)
      S = as.vector((W^2)%*%rep(1,ncol(W)))
      s = quantile(S, pct)
      ii = 1:length(tt)
      ii = ii[S < s]
      t.=lm(tt[ii]~bb[ii,])$coef[-1]
      if(sw==1)return(t.)
    }
    tt = as.vector(tt)
    bb=scale(bb,T,F)
    d = as.vector(bb %*% t.)
    
    if(trun>0){
      dlo=quantile(d,trun);dup=quantile(d,1-trun);pmax(pmin(d,dup),dlo)
      tup=quantile(tt,trun);tlo=quantile(tt,1-trun);pmax(pmin(tt,tup),tlo)
    }
    d=d-mean(d); sdd=sd(d)
    tt=tt-mean(tt);sdt=sd(tt)
    dz pmax(zc*d/sdd,-100);dz=pmin(dz,100)#truncate to [-100,100]
    dz=1/(1+exp(dz))
    az=qnorm(mean(dz))
    ad=mean(d^3)/(6*sdd^3)
    Atd=mean((tt^2)*d)/(2*sdt*sdd^2)
    a=c(ad,az,Atd); names(a)=c("a","az","A")
    a
  }