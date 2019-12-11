rm(list=ls())
# Function for generating a Bernoulli network given a random graph model P
rg.sample <- function(P){
  n <-  nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] <- runif(n*(n-1)/2)
  A <- (A + t(A))
  A <- (A < P) + 0 ;
  diag(A) <- 0
  return(A)
}
frob <- function(x){return(sqrt(sum(x^2)))}

# Function for ASE
stfp <- function(A, dim){
  require("irlba")
  
  A.svd <- irlba(A, nu = dim, nv = dim)
  A.svd.values <- A.svd$d[1:dim]
  A.svd.vectors <- A.svd$v[,1:dim]
  if(dim == 1)
    A.coords <- sqrt(A.svd.values) * A.svd.vectors
  else
    A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
  
  return(A.coords)
}

# Function for computing the minimum Frobenius norm after minimizing via Procrustes
procrustes <- function(X,Y, type = "I"){
  if(type == "C"){
    X <- X/norm(X, type = "F")
    Y <- Y/norm(Y, type = "F")
  }
  if(type == "D"){
    tX <- rowSums(X^2)
    tX[tX <= 1e-15] <- 1
    tY <- rowSums(Y^2)
    tY[tY <= 1e-15] <- 1
    X <- X/sqrt(tX)
    Y <- Y/sqrt(tY)
  }
  
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  return(list(error = norm(X%*%W - Y, type = "F"), W = W))
}

# Function for computing the minimum Frobenius norm after minimizing via Procrustes (scaling case)

procrustesC <- function(X,Y){
  
X <- X/sqrt(sum(X^2))
  
Y <- Y/sqrt(sum(Y^2))
  
tmp <- t(X)%*%Y
  
tmp.svd <- svd(tmp)
  
W <- tmp.svd$u %*% t(tmp.svd$v)
  
foo = X%*%W - Y
  
return(list(error = sqrt(sum(foo^2)), W = W))}

# Function for parametric bootstrap as per Tang et al (2017)
T.ase.equal <- function(A1,A2,d,B){
  X = stfp(A1,d)
  Y = stfp(A2,d)
  TX.vec <- numeric(B)
  TY.vec <- numeric(B)
  stat <- procrustes(X,Y, type="I")$error
  
  for(b in 1:B){
    A <- rg.sample(X %*% t(X))
    Xhat0 <- stfp(A,d)
    A <- rg.sample(X %*% t(X))
    Xhat1 <- stfp(A,d)
    
    A <- rg.sample(Y %*% t(Y))
    Yhat0 <- stfp(A,d)
    A <- rg.sample(Y %*% t(Y))
    Yhat1 <- stfp(A,d)
    
    TX.vec[b] <- procrustes(Xhat0, Xhat1, type="I")$error
    TY.vec[b] <- procrustes(Yhat0, Yhat1, type="I")$error
  }
  p1 = (sum(TX.vec > stat) + 0.5)/length(TX.vec)
  p2 = (sum(TY.vec > stat) + 0.5)/length(TY.vec)
  pval = max(p1,p2)
  return(pval)
}


# Function for parametric bootstrap as per Tang et al (2017)
T.ase.scale <- function(A1,A2,d,B){
  X = stfp(A1,d)
  Y = stfp(A2,d)
  TX.vec <- numeric(B)
  TY.vec <- numeric(B)
  stat <- procrustesC(X,Y)$error
  
  for(b in 1:B){
    A <- rg.sample(X %*% t(X))
    Xhat0 <- stfp(A,d)
    A <- rg.sample(X %*% t(X))
    Xhat1 <- stfp(A,d)
    
    A <- rg.sample(Y %*% t(Y))
    Yhat0 <- stfp(A,d)
    A <- rg.sample(Y %*% t(Y))
    Yhat1 <- stfp(A,d)
    
    TX.vec[b] <- procrustesC(Xhat0, Xhat1)$error
    TY.vec[b] <- procrustesC(Yhat0, Yhat1)$error
  }
  p1 = (sum(TX.vec > stat) + 0.5)/length(TX.vec)
  p2 = (sum(TY.vec > stat) + 0.5)/length(TY.vec)
  pval = max(p1,p2)
  return(pval)
}


TWcomm = function(A1,A2,r){
  n = ncol(A1)
  #spectral clustering to find common block structure
  C = (A1+A2)/2;
  d = colSums(C); d[d==0] = 1; d = 1/sqrt(d);
  d = diag(d)
  L = d%*%C%*%d 
  eigen = eigen(L, symmetric = TRUE)	# eigen-decomposition of Laplacian
  val = abs(eigen$values)			# pick the top K 
  id = order(-val)				# eigenvalues of L and put
  id_vec = id[1:r]				# their eigenvectors into a
  X = eigen$vectors[,id_vec]		# N*K matrix
  for (i in 1:n){
    foo = sqrt(sum(X[i,]^2))
    if (foo < 1e-4) {X[i,] = rep(0,r)} else {
      X[i,] = X[i,]/foo}}	# row-normalize
  K = kmeans(X, centers = r, nstart = 10);			# carry out K-meansi
  idx = K$cluster
  return(idx)
}

TWtest = function(A1, A2, r){
  n = ncol(A1)
  # community detection
  if (length(r)==1) {
    idx = TWcomm(A1, A2,r)}  else {
      # if we precompute communities and pass it as r. In this case, we
      # assume that the community id-s start from 1 to a maximum value
      idx = r                
      r = max(idx)}
  
  #computing the scaled difference matrix
  C = A1- A2 
  for (i in 1:r){
    for (j in 1:r){
      if (i==j) {
        foo = A1[idx==i,idx==j]
        Pij = sum(foo)/(nrow(foo)*(nrow(foo)-1))
        foo = A2[idx==i,idx==j]
        Qij = sum(foo)/(nrow(foo)*(nrow(foo)-1))} else {
          foo = A1[idx==i,idx==j]
          Pij = mean(foo);
          foo = A2[idx==i,idx==j]
          Qij = mean(foo)}
      denom = max(sqrt((n-1)*(Pij*(1-Pij)+Qij*(1-Qij))), 1e-5)
      C[idx==i,idx==j] = C[idx==i,idx==j]/denom;
    }
  }  
  # compute test statistic and p-value
  foo = max(abs(eigen(C, symmetric = TRUE)$values))
  testStat = n^(2/3)*(foo-2);
  pval = min(1,2*ptw(q=testStat, beta=1, lower.tail = FALSE))
  # we double p-value for Bonferroni correction. 
  # We know lambda_1 and -lambda_n are TW.
  # Spectral norm = max(lambda_1, -lambda_n), and hence, we need Bonferroni
  return(list(pval=pval, stat=testStat))
}
# Function for test of equality for Chung-Lu
T.frob.chunglu=function(A1,A2,B){  
  n=nrow(A1)
  d_1=rowSums(A1)
  d_2=rowSums(A2)
  P_1hat=(d_1%*%t(d_1))/sum(d_1)
  P_2hat=(d_2%*%t(d_2))/sum(d_2)
  diag(P_1hat)=0
  diag(P_2hat)=0
  foo = P_1hat-P_2hat
  t=sqrt(sum(foo^2))
  P_hat=(P_1hat+P_2hat)/2
  
  ###Bootstrap###
  R=rep(0,B)
  for (k in 1:B) {
    A1hat=rg.sample(P_hat)
    A2hat=rg.sample(P_hat) 
    d_1hat=rowSums(A1hat)
    d_2hat=rowSums(A2hat)
    Phat_1hat=(d_1hat%*%t(d_1hat))/sum(d_1hat)
    Phat_2hat=(d_2hat%*%t(d_2hat))/sum(d_2hat)
    diag(Phat_1hat)=0
    diag(Phat_2hat)=0
    foo = Phat_1hat-Phat_2hat
    R[k]=sqrt(sum(foo^2))}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}


# Function for test of scaling for Chung-Lu
T.scale.chunglu=function(A1,A2,B){  
  n=nrow(A1)
  d_1=rowSums(A1)
  d_2=rowSums(A2)
  P_1hat=(d_1%*%t(d_1))/sum(d_1)
  P_2hat=(d_2%*%t(d_2))/sum(d_2)
  diag(P_1hat)=0
  diag(P_2hat)=0
  rho_1hat = frob(P_1hat)
  rho_2hat = frob(P_2hat)
  foo1 = P_1hat/rho_1hat; foo2 = P_2hat/rho_2hat; foo = foo1-foo2
  t=frob(foo)
  ###Bootstrap###
  P_hat=(foo1+foo2)/2
  P_hat1 = P_hat*rho_1hat;P_hat2 = P_hat*rho_2hat; 
  R=rep(0,B)
  for (k in 1:B) {
    A1hat=rg.sample(P_hat1)
    A2hat=rg.sample(P_hat2) 
    d_1hat=rowSums(A1hat)
    d_2hat=rowSums(A2hat)
    Phat_1hat=(d_1hat%*%t(d_1hat))/sum(d_1hat)
    Phat_2hat=(d_2hat%*%t(d_2hat))/sum(d_2hat)
    diag(Phat_1hat)=0
    diag(Phat_2hat)=0
    rho_1hat = frob(Phat_1hat)
    rho_2hat = frob(Phat_2hat)
    foo1 = Phat_1hat/rho_1hat; foo2 = Phat_2hat/rho_2hat; foo = foo1-foo2
    R[k]=frob(foo)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}

# Function for parametric bootstrap for our proposed method
T.frob.rdpg = function(A1,A2,d,B){
  X1 = stfp(A1,d)
  X2 = stfp(A2,d)
  P1 = X1%*%t(X1)
  P2 = X2%*%t(X2)
  t = norm(P1-P2, type = "F")
  P = (P1+P2)/2
  S = rep(NA,B)
  for(b in 1:B){
    A1star = rg.sample(P)
    A2star = rg.sample(P)
    Xb = stfp(A1star,d)
    Yb = stfp(A2star,d)
    Q1 = Xb%*%t(Xb)
    Q2 = Yb%*%t(Yb)
    f = norm(Q1-Q2, type = "F")
    S[b] = f}
  pval = (sum(S >= t) +0.5)/B
  return(pval)
}


# Function for parametric bootstrap for our proposed method
T.scale.rdpg = function(A1,A2,d,B){
  X1 = stfp(A1,d)
  X2 = stfp(A2,d)
  P1 = X1%*%t(X1)
  P2 = X2%*%t(X2)
  foo1 = P1/sqrt(sum(P1^2)); foo2 = P2/sqrt(sum(P2^2)); foo = foo1-foo2
  t = sqrt(sum(foo^2))
  P = (foo1+foo2)/2
  P1star = P*sqrt(sum(P1^2))
  P2star = P*sqrt(sum(P2^2))
  S = rep(NA,B)
  for(b in 1:B){
    A1star = rg.sample(P1star)
    A2star = rg.sample(P2star)
    Xb = stfp(A1star,d)
    Yb = stfp(A2star,d)
    Q1 = Xb%*%t(Xb)
    Q2 = Yb%*%t(Yb)
    foo1 = Q1/sqrt(sum(Q1^2)); foo2 = Q2/sqrt(sum(Q2^2)); foo = foo1-foo2
    f = sqrt(sum(foo^2))
    S[b] = f}
  pval = (sum(S >= t) +0.5)/B
  return(pval)
}
# Estimation function for DCBM
DCBMfit<-function(A,groups){ 
  if (min(A)<0) stop("network has negative edges")
  if (isSymmetric(A)==FALSE) stop("network is not undirected")
  if (sum(diag(A)) > 0) stop("network has self-loops")
  
  n = nrow(A); tau = mean(rowSums(A))/2; k=groups
  
  #Step 1: estimate P
  
  # Step 1a: community detection via spectral clustering
  community = rep(NA, n)
  d <- colSums(A)			# degree
  null = which(d == 0)		# disconnected nodes
  if (length(null) > 0 ) { 	# remove disconnected nodes and give them random communities
    community[null] = sample(1:k,size=length(null),replace=T)
    A0 = A[-null,-null]
    N = nrow(A0)
    D =  diag(1/sqrt(colSums(A0)+rep(tau,N)))	# regularized degree natrix
    L = D%*%A0%*%D			# regularized Laplacian
    eigen = eigen(L, symmetric = TRUE)	# eigen-decomposition of Laplacian
    val = abs(eigen$values)			# pick the top K 
    id = order(-val)				# eigenvalues of L and put
    id_vec = id[1:k]				# their eigenvectors into a
    X = eigen$vectors[,id_vec]		# N*K matrix
    for (i in 1:N){
      foo = sqrt(sum(X[i,]^2))
      if (foo > 0) {X[i,] = X[i,]/sqrt(sum(X[i,]^2))} else X[i,] = rep(1/sqrt(k),k)}	# row-normalize
    K = kmeans(X, centers = k, nstart = 10);			# carry out K-means
    community[-null] = K$cluster} else {A0=A 
    
    N = nrow(A0)
    D =  diag(1/sqrt(colSums(A0)+rep(tau,N)))	# regularized degree natrix
    L = D%*%A0%*%D			# regularized Laplacian
    eigen = eigen(L, symmetric = TRUE)	# eigen-decomposition of Laplacian
    val = abs(eigen$values)			# pick the top K 
    id = order(-val)				# eigenvalues of L and put
    id_vec = id[1:k]				# their eigenvectors into a
    X = eigen$vectors[,id_vec]		# N*K matrix
    for (i in 1:N){
      foo = sqrt(sum(X[i,]^2))
      if (foo > 0) {X[i,] = X[i,]/sqrt(sum(X[i,]^2))} else X[i,] = rep(1/sqrt(k),k)}	# row-normalize
    K = kmeans(X, centers = k, nstart = 10);			# carry out K-means
    community = K$cluster}
  # End of Step 1a: community detection
  
  omega = matrix(NA, nrow=k, ncol=k)
  theta = colSums(A)
  for (r in 1:k){
    foo1 = which(community==r)
    foo = sum(theta[foo1])
    if (foo>0) {theta[foo1] = theta[foo1]/foo}
    for (s in 1:k){
      foo2 = which(community==s)
      omega[r,s] = sum(A[foo1,foo2])
    }}
  
  M = matrix(0,nrow=n,ncol=k) # membership-degree matrix
  for (i in 1:n){M[i,community[i]] = theta[i]}
  P.hat = M%*%omega%*%t(M)
  diag(P.hat) = 0 
  # End of Step 1: estimate P
  return(P.hat)
}

# T_{frob} function for DCBM
T.frob.dcbm=function(A1,A2,k,B){ 
  n=nrow(A1)
  #2. Fit DCBM model
  P_1.hat = DCBMfit(A=A1,groups = k)
  P_2.hat = DCBMfit(A=A2,groups = k)
  
  #3. Frobenius norm
  t = frob(P_1.hat-P_2.hat)
  
  ##Bootstrapping
  P.hat=(P_1.hat+P_2.hat)/2
  R=rep(NA,B)
  for (i in 1:B) 
  {
    A1hat = rg.sample(P.hat)
    A2hat = rg.sample(P.hat)
    P_1hat.hat = DCBMfit(A1hat,groups = k)
    P_2hat.hat = DCBMfit(A2hat,groups = k)
    
    #3. Frobenius norm
    R[i] = frob(P_1hat.hat-P_2hat.hat)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}


# T_{scale} function for DCBM
T.scale.dcbm=function(A1,A2,k,B){ 
  n=nrow(A1)
  #2. Fit DCBM model
  P_1hat = DCBMfit(A1,groups = k)
  P_2hat = DCBMfit(A2,groups = k)
  rho_1hat = frob(P_1hat)
  rho_2hat = frob(P_2hat)
  foo1 = P_1hat/rho_1hat; foo2 = P_2hat/rho_2hat; foo = foo1-foo2
  t=frob(foo)
  ###Bootstrap###
  P_hat=(foo1+foo2)/2
  P_hat1 = P_hat*rho_1hat;P_hat2 = P_hat*rho_2hat; 
  R=rep(0,B)
  for (i in 1:B) {
    A1hat=rg.sample(P_hat1)
    A2hat=rg.sample(P_hat2) 
    Phat_1hat=DCBMfit(A1hat,groups = k)
    Phat_2hat=DCBMfit(A2hat,groups = k)
    rho_1hat = frob(Phat_1hat)
    rho_2hat = frob(Phat_2hat)
    foo1 = Phat_1hat/rho_1hat; foo2 = Phat_2hat/rho_2hat; foo = foo1-foo2
    R[i]=frob(foo)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}


# Estimation function for PABM
PABMfit<-function(A,groups){ #estimation function for DCBM
  if (min(A)<0) stop("network has negative edges")
  if (isSymmetric(A)==FALSE) stop("network is not undirected")
  if (sum(diag(A)) > 0) stop("network has self-loops")
  
  n = nrow(A); tau = mean(rowSums(A))/2; k=groups
  
  #Step 1: estimate P
  
  # Step 1a: community detection via extreme points
  community = rep(NA, n)
  b.can = EPalgo(A,eps=0) # EP algorithm (no perturbation)
  Q.PA.can = rep(NA, ncol(b.can))	# array to store Q values
  for (i in 1:ncol(b.can)){
    #check if any cluster is empty
    foo = rep(NA, k)
    for (clus in 1:k) {foo[clus]=sum(b.can[,i]==clus)}
    if (min(foo)==0) {stop('Empty groups are not allowed')} 
    Q.PA.can[i] = Q.PA(A, b=b.can[,i])   # fit PABM
  } # end of i for loop
  foo1 = order(-Q.PA.can)[1] 
  b.PA = b.can[,foo1]   # community assignment that maximises Q.PA
  community = b.PA
  # End of Step 1a: community detection
  P.hat = P_PA(A,community)     # prob matrix fitted by PABM by comm detection
  diag(P.hat) = 0 # End of Step 1: estimate P
  return(P.hat)
}

######## Extreme points algo #####
EPalgo<-function(A,eps=0){
  
  ##### perturbed adj matrix ##### (Le pg 15, Amini 2013)
  tau = eps*(mean(colSums(A))/nrow(A))
  A = A + tau
  
  foo<-eigen(A, symmetric = TRUE)
  val = abs(foo$values)			# pick the top 2
  id = order(-val)				# eigenvalues of A and put
  id_vec = id[1:2]				# their eigenvectors into a 2*N matrix
  # columns of foo$vectors = eigenvectors of A
  # we want a 2xn matrix whose rows are the leading eigenvectors
  X = t(foo$vectors[,id_vec])		
  y = X[,1:2]
  
  comms = 1:2
  u <- list(comms)
  v = expand.grid(rep(u, 2))
  v = as.matrix(v)
  # initialize with the parallelogram
  epts = y%*%t(v)	# extreme pts are the columns of this matrix
  b.can = t(v)	# candidate configurations.
  row.names(b.can)=NULL
  
  ptm<-proc.time()
  for (i in 3:ncol(X)){
    b.can1 = rbind(b.can,rep(1,ncol(b.can)))
    b.can2 = rbind(b.can,rep(2,ncol(b.can)))
    b.can = cbind(b.can1,b.can2)
    foo = X[,1:i]%*%b.can
    hull = chull(t(foo))
    epts = foo[,hull]
    b.can = b.can[,hull]}	# next i = next row of X
  proc.time()-ptm
  
  ##### remove invalid candidates
  k = max(b.can)
  foo = b.can
  foo1 = NA
  for (i in 1:ncol(b.can)){
    foo2 = rep(NA,k)
    for (clus in 1:k){foo2[clus]=sum(b.can[,i]==clus)}
    if (min(foo2)==0){foo1 = c(foo1,i)}
  }
  if (length(foo1)>1) {foo1 = foo1[-1]
  b.can = b.can[,-foo1]}
  
  ###### remove eqv candidates
  foo1 = NA
  for (i in 2:ncol(b.can)){
    for (j in 1:i){
      foo4 = abs(b.can[,i] - b.can[,j])
      if (mean(foo4) == 1){ # this means b.can[,i] and b.can[,j] are exactly 1 apart
        foo1 = c(foo1,i)
        break}
    }}
  if (length(foo1)>1){
    foo1 = foo1[-1]
    b.can = b.can[,-foo1]
  }
  return(b.can)}	# end of function

############################################################
############## obtain MLE of probability matrix ########## 
############################################################
P_PA <- function(A, b){
  N<-ncol(A)
  foo<-f.PA(A,b)
  M<-foo$M; O<-foo$O
  lambda <- matrix(NA, nrow=N, ncol = max(b))
  P = matrix(NA,nrow=N,ncol=N)
  for (i in 1:N){
    s <- b[i]
    for (r in 1:max(b)){
      lambda[i,r] <- M[i,r]/sqrt(O[s,r])} }
  for (i in 1:N){  	# populate diagonal entries
    r = b[i]
    P[i,i] = (1/2)*(lambda[i,r]^2)}	# divide by 2 --- see notes section 3.3
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p = lambda[i,b[j]]*lambda[j,b[i]]
      P[i,j] = p
      P[j,i] = P[i,j]
    }}
  return(P)}

P_sim <- function(lambda, b){
  if (nrow(lambda) != length(b)) stop('dimension mismatch between lambda and b')
  N = nrow(lambda)
  k = ncol(lambda)
  b = as.factor(b)
  P = matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p = lambda[i,b[j]]*lambda[j,b[i]]
      P[i,j] = p
      P[j,i] = P[i,j]
    }}
  return(P)}

#####################################
########## PABM Likelihood ##########
#####################################
Q.PA <- function(A, b){
  foo<-f.PA(A,b)
  O=foo$O; M = foo$M
  s1 = sum(M*log(M),na.rm=TRUE) # na.rm = TRUE ignores M=0 cases as log(0) = NA
  s2 = sum(O*log(O),na.rm=TRUE) # na.rm = TRUE ignores O=0 cases as log(0) = NA
  return(2*s1-s2)}

#####################################################################
##### fn to calculate popularity M and block interaction O ##########
#####################################################################
f.PA<-function(A,b){	
  K<-max(b)       # no. of communities
  N<-nrow(A)      # no. of nodes
  M<-matrix(NA,nrow=N,ncol=K)  # popularity matrix
  O<-matrix(NA,nrow=K,ncol=K)  # community interaction matrix
  for (i in 1:N){		# calculate M
    for (r in 1:K){
      nodes = which(b == r)
      M[i,r] = sum(A[i,nodes])
    }}
  for (r in 1:K){		# calculate O
    for (s in r:K){
      nodes1 = which(b == r)
      nodes2 = which(b == s)
      O[r,s] = sum(A[nodes1,nodes2])
      O[s,r] = O[r,s]
    }}
  list(M=M, O=O)}

# T_{frob} function for PABM
T.frob.pabm=function(A1,A2,k,B){ 
  n=nrow(A1)
  #2. Fit DCBM model
  P_1.hat = PABMfit(A=A1,groups = k)
  P_2.hat = PABMfit(A=A2,groups = k)
  
  #3. Frobenius norm
  t = frob(P_1.hat-P_2.hat)
  
  ##Bootstrapping
  P.hat=(P_1.hat+P_2.hat)/2
  R=rep(NA,B)
  for (i in 1:B) 
  {
    A1hat = rg.sample(P.hat)
    A2hat = rg.sample(P.hat)
    P_1hat.hat = PABMfit(A1hat,groups = k)
    P_2hat.hat = PABMfit(A2hat,groups = k)
    
    #3. Frobenius norm
    R[i] = frob(P_1hat.hat-P_2hat.hat)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}


# T_{scale} function for PABM
T.scale.pabm=function(A1,A2,k,B){ 
  n=nrow(A1)
  #2. Fit DCBM model
  P_1hat = PABMfit(A1,groups = k)
  P_2hat = PABMfit(A2,groups = k)
  rho_1hat = frob(P_1hat)
  rho_2hat = frob(P_2hat)
  foo1 = P_1hat/rho_1hat; foo2 = P_2hat/rho_2hat; foo = foo1-foo2
  t=frob(foo)
  ###Bootstrap###
  P_hat=(foo1+foo2)/2
  P_hat1 = P_hat*rho_1hat;P_hat2 = P_hat*rho_2hat; 
  R=rep(0,B)
  for (i in 1:B) {
    A1hat=rg.sample(P_hat1)
    A2hat=rg.sample(P_hat2) 
    Phat_1hat=PABMfit(A1hat,groups = k)
    Phat_2hat=PABMfit(A2hat,groups = k)
    rho_1hat = frob(Phat_1hat)
    rho_2hat = frob(Phat_2hat)
    foo1 = Phat_1hat/rho_1hat; foo2 = Phat_2hat/rho_2hat; foo = foo1-foo2
    R[i]=frob(foo)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}

##### Parallel versions of the test functions #####
# T_{frob} function for PABM in parallel
T.frob.pabm.par=function(A1,A2,k,B,cores){ 
  library("foreach")
  library("doParallel")
  n=nrow(A1)
  #2. Fit DCBM model
  P_1.hat = PABMfit(A=A1,groups = k)
  P_2.hat = PABMfit(A=A2,groups = k)
  
  #3. Frobenius norm
  t = frob(P_1.hat-P_2.hat)
  
  ##Bootstrapping
  P.hat=(P_1.hat+P_2.hat)/2
  registerDoParallel(min(cores,detectCores(),B))
  R = foreach(i = 1:B) %dopar% {
    A1hat = rg.sample(P.hat)
    A2hat = rg.sample(P.hat)
    P_1hat.hat = PABMfit(A1hat,groups = k)
    P_2hat.hat = PABMfit(A2hat,groups = k)
    
    #3. Frobenius norm
    frob(P_1hat.hat-P_2hat.hat)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}

# T_{scale} function for PABM in parallel
T.scale.pabm.par=function(A1,A2,k,B, cores){ 
  library("foreach")
  library("doParallel")
  n=nrow(A1)
  #2. Fit DCBM model
  P_1hat = PABMfit(A1,groups = k)
  P_2hat = PABMfit(A2,groups = k)
  rho_1hat = frob(P_1hat)
  rho_2hat = frob(P_2hat)
  foo1 = P_1hat/rho_1hat; foo2 = P_2hat/rho_2hat; foo = foo1-foo2
  t=frob(foo)
  ###Bootstrap###
  P_hat=(foo1+foo2)/2
  P_hat1 = P_hat*rho_1hat;P_hat2 = P_hat*rho_2hat; 
  registerDoParallel(min(cores,detectCores(),B))
  R = foreach(i = 1:B) %dopar% {
    A1hat=rg.sample(P_hat1)
    A2hat=rg.sample(P_hat2) 
    Phat_1hat=PABMfit(A1hat,groups = k)
    Phat_2hat=PABMfit(A2hat,groups = k)
    rho_1hat = frob(Phat_1hat)
    rho_2hat = frob(Phat_2hat)
    foo1 = Phat_1hat/rho_1hat; foo2 = Phat_2hat/rho_2hat; foo = foo1-foo2
    frob(foo)}
  pval = (sum(R >= t) +0.5)/B
  return(pval)
}
