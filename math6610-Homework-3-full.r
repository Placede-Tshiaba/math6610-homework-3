library(sna)
library(matrixcalc)
library(pracma)
library(Matrix)
library(DescTools)

#### Random permutation matrix######
g<-rgraph(5)
g  
rmperm(g)
rmperm(eye(4,4))


#### LU decomposition ####
A = matrix(c(1,2,2,3), ncol =2)
A
lu.decomposition(A)

R = matrix(c(10^-20, 1, 1,2), ncol =2)
R
lu.decomposition(R)


#### Problem 22.3 ###
### Figure 1 ###
dims=zeros((9+2*15)*12,1);
rhos = dims;

ii = 0
d = 0
for (f in 1:9){
    m=f*10^d;
    for (j in 1:12){
        ii = ii+1;
        A = rand(m)*2-1
        L = expand(lu(A))$L
		U = expand(lu(A))$U
		P = expand(lu(A))$P
        rho = max(max(abs(U)))/max(max(abs(A)))
        rhos[ii] = rho
        dims[ii] = m
    }
}

for (d in 1:2){
    for (f in c(1, 1.1, 1.2, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9)){
        m=floor(f*10^d);
        for (j in 1:12){
            ii = ii+1;
            A = rand(m)*2-1;
            L = expand(lu(A))$L
		    U = expand(lu(A))$U
		    P = expand(lu(A))$P
            rho = max(max(abs(U)))/max(max(A));
            rhos[ii]=rho
            dims[ii] = m
        }
    }
}
dims
rhos

plot(dims, rhos, pch = 20, ylim = c(0,100), xlab = "Matrix dimension", las = 1, ylab = expression(paste("Growth factor ", rho)),)
points(c(1:1000), c(1:1000)^(1/2), type = 'l', lty =2, col = 3)
points(c(1:1000), c(1:1000)^(2/3), type = 'l', lty =2, col = 4)
text(600, 20, lab = expression(m^{1/2}), col = 3)
text(600, 75, lab = expression(m^{3/2}), col = 4)


#### Figure 2 ####
rhos=zeros(10^6,3)
for (m in 1:3){
    for (ii in 1:10^6){
        A = rand(2^(m+2))*2-1;
        L = expand(lu(A))$L
		U = expand(lu(A))$U
		P = expand(lu(A))$P
        rho = max(max(abs(U)))/max(max(abs(A)))
        rhos[ii,m] = rho;
    }
	}
    x=seq(0,10, 0.01);
    N=hist(rhos,x);
    pdf = N/(10^6*0.01);
    plot(x,pdf) 
    }
head(rhos)
dim(rhos)
summary(rhos)

    x=seq(0.01,10, 0.01);
    N1 = hist(rhos[1,])
	str(N1)
    pdf1 = N1/(10^6*0.01)
	length(pdf1)
    plot(x,pdf1)


plot(density(rhos[1,]), grid(10,10), type = "l", las = 1, xlab = expression(rho), ylab = "Probability density", main = "", ylim = c(0,.3), xlim = c(0,15))
points(density(rhos[2,]), col =2, type = "l")
points(density(rhos[3,]), col =3, type = "l")
text(5, 0.06, "m = 8", col = 1)
text(5, 0.18, "m = 16", col = 2)
text(5, 0.14, "m = 32", col = 3)

### Figure 3 and 4  #####
A = rand(128)*2-1
L = as.matrix(expand(lu(A))$L)
U = as.matrix(expand(lu(A))$U)
P = expand(lu(A))$P

qr = qr(L)$qr
q = qr
r = qr
#q[upper.tri(q)] <- 0
#r[lower.tri(r)] <- 0

#### 3.a ####
max(abs(inv(L)))
image(abs(inv(L)>1))

#### 4.a ####
image(abs(q)>1/sqrt(128))


#### Figure 3.b ####
Lp = L*sign(rand(128))
qpr = qr(Lp)$qr
qp = qpr
rp = qpr

### 3b  #####
max(abs(inv(Lp)))
image(abs(inv(Lp))>1)

#### 4b ####
image(abs(qp)>1/sqrt(128))



