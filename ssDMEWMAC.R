X = matrix(c(-1.347, 25.600, 12.400, 0.256, 58.600,
          -1.386, 24.800, 13.300, 0.285, 58.100,
          -0.478, 24.900, 12.700, 0.248, 57.800,
          -1.772, 25.300, 12.200, 0.265, 58.200,
          -1.309, 25.100, 12.300, 0.236, 58.800,
          -0.580, 24.500, 12.200, 0.252, 56.700,
          -0.892, 26.000, 12.200, 0.270, 58.600,
          -1.470, 24.400, 12.700, 0.248, 57.900,
          -0.673, 25.600, 12.100, 0.218, 58.100,
          -0.315, 25.400, 12.400, 0.212, 57.900), nrow = 10, ncol = 5, byrow = T)

x <- cbind(1, X)

n <-dim(x)[1]
p <- dim(x)[2]

r <- matrix(rep(0,p*p),p,p)
rec <- matrix(rep(0,n*p),n,p)
g <- rep(0,p)
u <- matrix(rep(0,n*p),n,p)

for(i in 1:n){
  eps <- 1e-9
  
  uj<-rep(0,p)
  r1<-r
  xf<-x[i,c(1:p)]
  
  for(j in 1:p){
    if (abs(xf[j]) > eps) {
      d<-sqrt(r[j,j]^2+xf[j]^2)
      cos<-r[j,j]/d
      sin<-xf[j]/d
      r[j,j]<-d
      
    j1<-j+1
    
    if(j1 <= p){
      for(kk in j1:p){
        c1<-r[j,kk]*cos+xf[kk]*sin
        c2<--r[j,kk]*sin+xf[kk]*cos
        r[j,kk]<-c1
        g[kk]<-xf[kk]<-c2
        }
      }
    }
  }
  rec[i,]<-g
  if(i > p){
    
    for(jj in 2:p){
      df<-i-jj
      tj<-g[jj]/(r1[jj,jj]/sqrt(df))
      uj[jj]<-qnorm(pt(tj,df))
    }
  u[i,]<-uj
  
  }
}

rr<-rec[,2:p]
rr