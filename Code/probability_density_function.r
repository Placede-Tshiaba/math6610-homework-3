library(sna)
library(matrixcalc)
library(pracma)
library(Matrix)
library(DescTools)

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

