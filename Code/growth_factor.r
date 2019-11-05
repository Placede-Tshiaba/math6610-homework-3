library(sna)
library(matrixcalc)
library(pracma)
library(Matrix)
library(DescTools)

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
