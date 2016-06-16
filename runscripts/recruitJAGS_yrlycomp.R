model{
  for(i in 1:N){
    for(j in 1:Nspp){
      y[i,j]~dnegbin(q[i,j],theta[j])
      q[i,j]<-theta[j]/(theta[j]+lambda[i,j])
      lambda[i,j]<-trueP1[i,j]*fecundity[i,j]
      log(fecundity[i,j])<-intcpt.yr[year[i],j]+intcpt.gr[Group[i],j]+inprod(dd[1:Nspp,j,year[i]],trueP2[i,1:Nspp])
      trueP1[i,j]<-parents1[i,j]*u[j]+parents2[i,j]*(1-u[j])
      trueP2[i,j]<-sqrt(trueP1[i,j])
    }
  }
  
  #priors
  # intercepts
  for(s in 1:Nspp){
    u[s]~dunif(0,1)
    theta[s]~dgamma(0.001,0.001)
    g.tau[s]~dgamma(2,0.5)
    for(g in 1:Ngroups){
      intcpt.gr[g,s]~dnorm(0,g.tau[s])
    }
    intcpt.mu[s]~dnorm(0,0.001)
    intcpt.tau[s]~dgamma(0.001,0.001)
    for(k in 1:Nspp){
      dd.mu[k,s]~dnorm(0,0.001)
      dd.tau[k,s]~dgamma(0.001,0.001)
      for(mm in 1:Nyrs){
        dd[k,s,mm]~dnorm(dd.mu[k,s],dd.tau[k,s])
      }
    }
    for(m in 1:Nyrs){
      intcpt.yr[m,s]~dnorm(intcpt.mu[s],intcpt.tau[s])  
    }
  }
  
}