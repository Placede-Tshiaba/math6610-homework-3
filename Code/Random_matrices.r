library(sna)
library(matrixcalc)
library(pracma)
library(Matrix)
library(DescTools)

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


