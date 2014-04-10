estim_herit=function(Y,W)
{
  n=length(Y)
  N=ncol(W)
  a=n/N
  nb_iter=20
  eta_init=0.5
  eta=rep(0,nb_iter+1)
  eta[1]=eta_init
  Z=scale(W,center=TRUE,scale=TRUE)
  M=Z%*%t(Z)/N
  O=eigen(M)$vectors
  lambda=eigen(M)$values
  Y_tilde=t(O)%*%Y
  for (nb in 1:nb_iter)
  {
    A=sum(Y_tilde^2*(lambda-1)/(eta[nb]*(lambda-1)+1)^2)/sum(Y_tilde^2/(eta[nb]*(lambda-1)+1))-1/n*sum((lambda-1)/(eta[nb]*(lambda-1)+1))
    B=((-2*sum(Y_tilde^2*(lambda-1)^2/(eta[nb]*(lambda-1)+1)^3)*sum(Y_tilde^2/(eta[nb]*(lambda-1)+1))+(sum(Y_tilde^2*(lambda-1)/(eta[nb]*(lambda-1)+1)^2))^2)/(sum(Y_tilde^2/(eta[nb]*(lambda-1)+1)))^2 +1/n*sum((lambda-1)^2/(eta[nb]*(lambda-1)+1)^2))
    eta[(nb+1)]=eta[nb]-A/B
  }
  eta_chap=min(max(eta[nb_iter],0.01),0.99)
  w=(lambda-1)/(eta_chap*(lambda-1)+1)
  s_eta=sqrt(2/(sum(w^2)-1/n*sum(w)^2))
  eta_inf=eta_chap -1.96*s_eta
  eta_sup=eta_chap+1.96*s_eta
  list(heritability=eta_chap,IC_inf=eta_inf,IC_sup=eta_sup,standard_deviation=s_eta)
}