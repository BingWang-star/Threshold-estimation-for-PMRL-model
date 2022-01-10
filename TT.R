#' Title The function to calculate the sup-Wald statistics
#' @param T: observed time 
#' @param Delta: censoring indicator 
#' @param X: the threshold covariate
#' @param Z1: covariates for full cohort 
#' @param Z2: covariates for subcohort 
#' @param zetac: threshold candidates 
 
#'
#' @return the sup-Wald statistics
#' @export
#'
#' @examples
TT<-function(T,Delta,X,Z1,Z2,zetac){
  p1<-ifelse(is.null(Z1),0,ifelse(is.null(dim(Z1)[2]),1,dim(Z1)[2]))
  p2<-ifelse(is.null(Z2),0,ifelse(is.null(dim(Z2)[2]),1,dim(Z2)[2]))
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
  len_z<-length(zetac)
  sta<-numeric()
  for(k in 1: len_z){
    xi_ini<-rep(0.2,p)
    m_ini<-rep(1,length(alpha))
    iter=1
    convergence=100
    diff<-100
    zetaf<-zetac[k]
    Zo<-cbind(Z1,(X>zetaf)*1,(X>zetaf)*Z2)
    
    while(iter<=50 &diff>1E-7& convergence>1E-9){
      
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
        res<-res*(res>0)
        return(res)
      }
      mx<-m(T)
      esf_xi<-function(u){
        res<-as.vector(t(Zo)%*%(mx*Delta)-t(Zo)%*%(exp(-Zo%*%u)*T)-t(Zo)%*%(mx-m(0)))
        res
      }
      xi_est<-nleqslv(xi_ini,esf_xi)$x
      iter=iter+1
      m_est<-m(alpha)
      convergence = mean((m_est-m_ini)^2)+mean((xi_est-xi_ini)^2)
      diff<-mean(abs((m_est-m_ini)[seq(1,len-1,2)]-(m_est-m_ini)[seq(2,len,2)]))
      xi_ini<-xi_est
      m_ini<-m_est
    }
    
    dif_m<-c(diff(m_ini),0)
    Z_bar<-matrix(nrow=len,ncol=p-1)
    for(i in 1:len){
      Z_bar[i,]<- as.vector((T>=alpha[i])%*%Zo/sum(T>=alpha[i]))
    }
    A<-matrix(0,nrow=p-1,ncol = p)
    for(i in 1:n){
      for(j in 1:len){
        A<-A+(Zo[i,]-Z_bar[j])%*%t(Zo[i,]-Z_bar[j])*dif_alpha[j]*(alpha[j]<=T[i])*exp(as.vector(-xi_est%*%Zo[i,]))/n
      }
    }
    
    V<-matrix(0,nrow=p-1,ncol = p)
    for(i in 1:n){
      for(j in 1:len)
        V<-V+(Zo[i,]-Z_bar[j])%*%t(Zo[i,]-Z_bar[j])*m_ini[j]*(exp(as.vector(-xi_est%*%Zo[i,]))*dif_alpha[j]+dif_m[j])*(alpha[j]<=T[i])/n
    }
    va<-(solve(A)%*%V%*%solve(A)/n)[(p1+1):(p1+p2+1),(p1+1):(p1+p2+1)]
    sta[k]<-xi_ini[(p1+1):(p1+p2+1)]%*%solve(va)%*%xi_ini[(p1+1):(p1+p2+1)]
  }
  sup_Wald<-max(sta)
  return(sup_Wald)
}