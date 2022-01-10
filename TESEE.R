#' Title: The function to estimate the parameters
#' @param T: observed time 
#' @param Delta: censoring indicator 
#' @param X: the threshold covariate
#' @param Z1: covariates for full cohort 
#' @param Z2: covariates for subcohort 
#'
#' @return parameters estimation and their standard errors 
#' @export 
#'
#' @examples
#' # n=800 #sample size
# library(nleqslv)
# set.seed(1005*a)
# Z<-cbind(rbinom(n,1,0.5)-0.5,runif(n,0,1)-0.5)#covariate
# 
# 
# X<-runif(n,-1,1)
# zeta<-0
# ZO<-cbind(Z,(X>zeta)*1,(X>zeta)*Z)
# #dimention of covariate
# beta<-c(0.5,0.5) #coefficient
# eta<-c(0.5,0.5)
# alfa<--0.5
# xi<-c(beta,alfa,eta)
# theta<-c(xi,zeta)
# 
# m0<-function(t) rep(1,length(t)) 
# u <- runif(n, 0, 1)
# T1<-numeric()
# for(i in 1:n){
#   g<-function(t) 1/(m0(t)*exp(as.vector(ZO[i,]%*%xi)))
#   func<-function(t){
#     g(t)/g(0)*exp(-integrate(g,0,t)$value)-u[i]}
#   if(func(100)>=0){T[i]<-100}else{
#     T1[i]<-uniroot(func,interval = c(0,100))$root}
# }
# 
# C<-rexp(n,0.25) #censoring time
# Delta<-(T1<=C)*1 #censoring indicator
# T<-T1*Delta+C*(1-Delta) #observed time
#TESEE(T,Delta,X,Z1,Z2)
TESEE<-function(T,Delta,X,Z1,Z2){
  n<-length(T)
  p1<-ifelse(is.null(Z1),0,ifelse(is.null(dim(Z1)[2]),1,dim(Z1)[2]))
  p2<-ifelse(is.null(Z2),0,ifelse(is.null(dim(Z2)[2]),1,dim(Z2)[2]))
  h<-log(n)*n^(-1/2)*sd(X)
  alpha<-c(0,sort(unique(T[Delta==1])))
  p<-p1+p2+1
  Delta[which(T==max(T))]=1
  alpha<-c(0,sort(unique(T[Delta==1])))
  dif_alpha<-c(diff(alpha),0)
  Y<-numeric()
  for(i in 1:n){
    Y[i]<-1/sum(T>=T[i])
  }
  
  S<-function(t){
    l<-length(t)
    res=numeric()
    for(i in 1:l){
      res[i]<-as.vector(exp(-((T<=t[i])*Delta)%*%Y))
    }
    return(res)
  }
  xi_ini<-rep(0.1,p)
  zeta_ini<-mean(X)
  theta_ini<-c(xi_ini,zeta_ini)
  m_ini<-rep(1,length(alpha))
  iter=1
  convergence=100
  diff<-100
  len<-length(alpha)
  while(iter<=50&diff>1E-8 & convergence>1E-9){
    Zo<-cbind(Z1,pnorm((X-zeta_ini)/h),pnorm((X-zeta_ini)/h)*Z2)
    Sp<-function(t){
      l<-length(t)
      res=numeric()
      for(i in 1:l){
        res[i]<-as.vector(t(exp(-Zo%*%xi_ini))%*%(T>=t[i])/sum(T>=t[i]))
      }
      return(res)
    }
    
    S_int<-S(alpha)*Sp(alpha)*dif_alpha
    int_S<-function(t){
      l<-sum(alpha<=t)-1
      if(l==0){return(t) }else{
        return(sum(S_int[1:l])+(t-alpha[l+1])*S(alpha[l+1])*Sp(alpha[l+1]))
      }
    }
    
    
    
    m<-function(t){
      l=length(t)
      res=numeric()
      for(i in 1:l)
        res[i]<-(int_S(max(alpha))-int_S(t[i]))/S(t[i])
      #res<-res*(res>0)
      return(res)
    }
    plot(alpha,m(alpha),type = 'l')
    
    mx<-m(T)
    esf_xi<-function(u){
      res<-as.vector(t(Zo)%*%(mx*Delta)-t(Zo)%*%(exp(-Zo%*%u)*T)-t(Zo)%*%(mx-m(0)))
      res
    }
    xi_est<-nleqslv(xi_ini,esf_xi)$x
    
    W<-as.vector(xi_ini[p1+1]+xi_ini[(p1+2):(p1+p2+1)]%*%t(Z2))
    up<-quantile(X,0.8)
    low<-quantile(X,0.2)
    esf_zeta<-function(u){
      l<-length(u)
      res<-numeric()
      for(i in 1:l){
        Zo<-cbind(Z1,pnorm((X-u[i])/h)*1,pnorm((X-u[i])/h)*Z2)
        if(u[i]<=up& u[i]>=low ){
          res[i]<-as.vector((W*dnorm((X-u[i])/h))%*%(mx*Delta)-(W*dnorm((X-u[i])/h))%*%(exp(-Zo%*%xi_est)*T)-(W*dnorm((X-u[i])/h))%*%(mx-m(0)))
        }else{ res[i]<-1000}
      }
      res
    }
    
    #plot(seq(low,up,0.01),esf_zeta(seq(low,up,0.01)),type = 'l')
    zeta_est<-nleqslv(zeta_ini,esf_zeta)$x
    
    iter=iter+1
    m_est<-m(alpha)
    convergence = mean((m_est-m_ini)^2)+mean((xi_est-xi_ini)^2+(zeta_est-zeta_ini)^2)
    diff<-mean(abs((m_est-m_ini)[seq(1,len-1,2)]-(m_est-m_ini)[seq(2,len,2)]))
    xi_ini<-xi_est
    zeta_ini<-zeta_est
    m_ini<-m(alpha)
  }
  theta_est<-c(xi_est,zeta_est)
  if(is.null(Z2)){
  W<-xi_est[p1+1]*rep(1,n)
  }else{W<-as.vector(xi_est[p1+1]+xi_est[(p1+2):(p1+p2+1)]%*%t(Z2))}
  Zo<-cbind(Z1,pnorm((X-zeta_est)/h),pnorm((X-zeta_est)/h)*Z2)
  Zw<-cbind(Z1,pnorm((X-zeta_est)/h),pnorm((X-zeta_est)/h)*Z2,W/h)
  
  dif_m<-c(diff(m_ini),0)
  Zw_bar<-matrix(nrow=len,ncol=p+1)
  for(i in 1:len){
    Zw_bar[i,]<- as.vector((T>=alpha[i])%*%Zw/sum(T>=alpha[i]))
  }
  A<-matrix(0,nrow=p+1,ncol = p+1)
  for(i in 1:n){
    l<-sum(T[i]>=alpha)
    if(l==len){
      for(j in 1:len){
        A<-A+(Zw[i,]-Zw_bar[j])%*%t(Zw[i,]-Zw_bar[j])*dif_alpha[j]*(alpha[j]<=T[i])*exp(as.vector(-xi_est%*%Zo[i,]))/n
      } 
    }else{
      for(j in 1:(l+1)){
        A<-A+(Zw[i,]-Zw_bar[j])%*%t(Zw[i,]-Zw_bar[j])*dif_alpha[j]*(alpha[j]<=T[i])*exp(as.vector(-xi_est%*%Zo[i,]))/n
        
      }
    }
  }
  D<-diag(c(rep(1,p),sqrt(h)))
  A=D%*%A%*%D
  V<-matrix(0,nrow=p+1,ncol = p+1)
  for(i in 1:n){
    l<-sum(T[i]>=alpha)
    if(l==len){
      for(j in 1:l){
        V<-V+(Zw[i,]-Zw_bar[j])%*%t(Zw[i,]-Zw_bar[j])*m_ini[j]*(exp(as.vector(-xi_est%*%Zo[i,]))*dif_alpha[j]+dif_m[j])*(alpha[j]<=T[i])/n
      }
    }else{  for(j in 1:(l+1)){
      V<-V+(Zw[i,]-Zw_bar[j])%*%t(Zw[i,]-Zw_bar[j])*m_ini[j]*(exp(as.vector(-xi_est%*%Zo[i,]))*dif_alpha[j]+dif_m[j])*(alpha[j]<=T[i])/n
    }
    }
  }
  V=D%*%V%*%D
  se<-sqrt(abs(diag(solve(A)%*%V%*%solve(A)/n)))%*%D
  return(c(theta_est,se))
}
