############################################################################################
## Here are the contraceptive use data from page 46 of the lecture notes (and from the Stata handout), 
## showing the distribution of 1607 currently married and fecund women interviewed in the Fiji 
## Fertility Survey, according to age, education, desire for more children and current use of contraception.
############################################################################################

############################################################################################
## The dataset is also available in the format used in the Stata handout. This version has 32 
## rows corresponding to all possible covariate and response patterns, and includes a weight 
## indicating the frequency of each combination. The file has 5 columns with numeric codes:
##
##    age (four groups, 1=<25, 2=25-29, 3=30-39 and 4=40-49),
##    education (0=none, 1=some),
##    desire for more children (0=more, 1=no more),
##    contraceptive use (0=no, 1=yes), and
##    frequency (number of cases in this category). 
############################################################################################
cuse <- read.table("cuse.dat", header=TRUE)

lrfit <- glm( cbind(using, notUsing) ~ age + education + wantsMore , family = binomial, data = cuse)
lrfit

cuse$noMore <- cuse$wantsMore == "no"
cuse$hiEduc <- cuse$education == "high"

glm( cbind(using,notUsing) ~ age + hiEduc + noMore, family=binomial, data = cuse)

1-pchisq(deviance(lrfit),lrfit$df.residual)

lrfit <- glm( cbind(using,notUsing) ~ age * noMore + hiEduc , family=binomial, data = cuse)
lrfit

1-pchisq(deviance(lrfit),lrfit$df.residual)

residuals(lrfit,type="pearson")
1-pchisq(sum(residuals(lrfit,type="pearson")^2),lrfit$df.residu)

############################################################################################
##  This dataset has information on lung cancer deaths by age and smoking status.
##
##  The file in "raw" format, smoking.raw, has four columns:
##
##    age: in five-year age groups coded 1 to 9 for 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80+.
##    smoking status: coded 1 = doesn't smoke, 2 = smokes cigars or pipe only, 3 = smokes cigarrettes and cigar or pipe, and 4 = smokes cigarrettes only,
##    population: in hundreds of thousands, and
##    deaths: number of lung cancer deaths in a year. 
############################################################################################
smoke<- read.table("smoking.dat", header=TRUE)

lrfit1 <- glm( dead ~ age + smoke, family = poisson, data = smoke)
lrfit1

1-pchisq(deviance(lrfit1),lrfit1$df.residual)

lrfit2 <- glm( dead ~ age + smoke, offset = log(pop), family = poisson, data = smoke)
lrfit2

1-pchisq(deviance(lrfit2),lrfit2$df.residual)

residuals(lrfit2,type="pearson")
1-pchisq(sum(residuals(lrfit2,type="pearson")^2),lrfit2$df.residu)

#################################################################
## Envelope de Probabilidade - Modelo Poisson
fit.model <- lrfit2	## MUITO IMPORTANTE !!!!
X <- model.matrix(fit.model)
n <- nrow(X)
p <- ncol(X)
w <- fit.model$weights
W <- diag(w)
H <- solve(t(X)%*%W%*%X)
H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
h <- diag(H)

# Desvio Residual - "Ajustado"
td <- resid(fit.model,type="deviance")/sqrt(1-h)

# Simulacao do envelope
m <- 1000
e <- matrix(0,n,m)	# Armazena as simulacoes
for(i in 1:m){
  # Gera novas saidas a partir dos valores ajustados (originais)
  nresp <- rpois(n, fitted(fit.model))
  fit   <- glm(nresp ~ X, poisson(link = "log"),offset = log(smoke$pop))
  w     <- fit$weights
  W     <- diag(w)
  H     <- solve(t(X)%*%W%*%X)
  H     <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h     <- diag(H)
  e[,i] <- sort( resid(fit,type="deviance")/sqrt(1-h) )
}

# Construcao do envelope com 95% de confianca
e1 <- numeric(n)		# limites inferiores
e2 <- numeric(n)		# limites superiores
for(i in 1:n){
  e0    <- sort(e[i,])	# ordena cada percentil
  e1[i] <- e0[floor(m*.05)]		# Percentil 5
  e2[i] <- e0[floor(m*.95)]		# Percentil 95
}

med   <- apply(e,1,mean)	# Define as medias para cada percentil
faixa <- range(td,e1,e2)	# Define os minimos e maximos
par(pty="s")			# square plotting region
qqnorm(td,xlab="Percentis da N(0,1)",ylab="Componente do Desvio Padronizado",ylim=faixa)
par(new=T)									# Permite inserir novas linhas
qqnorm(e1,xlab="",ylab="",type="l",ylim=faixa,lty=1)		# Limite Inferior
par(new=T)
qqnorm(e2,xlab="",ylab="",type="l",ylim=faixa,lty=1)		# Limite Superior
par(new=T)
qqnorm(med,xlab="",ylab="",type="l",ylim=faixa,lty=2)	# Linha Central
