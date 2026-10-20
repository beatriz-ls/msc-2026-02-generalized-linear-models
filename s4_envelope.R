## ----setup, include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, size="footnotesize")



## ----eval=TRUE-----------------------------------------------------------------
# Dados (Contraceptive Use - Fiji)
cuse <- read.table("cuse.dat", header=TRUE)
cuse$noMore <- cuse$wantsMore == "no"
cuse$hiEduc <- cuse$education == "high"

# Modelo MLG Binomial
lrfit <- glm(cbind(using,notUsing) ~ age * noMore + hiEduc, 
             family=binomial, data = cuse)

# Extração Matemática do Resíduo Padronizado
X <- model.matrix(lrfit); n <- nrow(X)
W <- diag(lrfit$weights)
H <- sqrt(W) %*% X %*% solve(t(X)%*%W%*%X) %*% t(X) %*% sqrt(W)
h <- diag(H) # Leverage
fi <- deviance(lrfit)/(nrow(cuse)-ncol(X))
td <- resid(lrfit, type="deviance") * sqrt(1/(fi*(1-h)))


## ----eval=FALSE----------------------------------------------------------------
# m <- 100 # Simulações
# e <- matrix(0, n, m)
# size <- (cuse$using + cuse$notUsing)
# 
# for(i in 1:m){
#   resp <- rbinom(n, size=size, prob=fitted(lrfit))
#   fit_sim <- glm(cbind(resp, size-resp) ~ X, family=binomial)
# 
#   w_sim <- fit_sim$weights
#   H_sim <- sqrt(diag(w_sim)) %*% X %*%
#     solve(t(X)%*%diag(w_sim)%*%X) %*% t(X) %*% sqrt(diag(w_sim))
#   h_sim <- diag(H_sim)
# 
#   phi <- deviance(lrfit)/(nrow(cuse)-ncol(X))
#   res_sim <- resid(fit_sim, type="deviance") * sqrt(1/(phi*(1-h_sim)))
#   e[,i] <- sort(res_sim) # Ordenação para o envelope
# }


## ----eval=FALSE----------------------------------------------------------------
# m <- 1000
# e <- matrix(0, n, m)
# mu_orig <- fitted(lrfit); fam <- binomial()
# phi <- deviance(lrfit) / (nrow(cuse) - ncol(X)) # Fora do loop!
# 
# for(i in 1:m){
#   resp <- rbinom(n, size=size, prob=mu_orig)
#   prop <- resp / size
#   # Ajuste sem overhead de fórmula
#   fit_sim <- glm.fit(x = X, y = prop, weights = size, family = fam)
#   # Leverage rápido via matriz ortogonal Q
#   h_sim <- rowSums(qr.Q(fit_sim$qr)^2)
#   # Cálculo vetorizado do resíduo deviance
#   dev_sq <- fam$dev.resids(prop, fit_sim$fitted.values, size)
#   sgn <- sign(prop - fit_sim$fitted.values)
#   res_sim <- sgn * sqrt(dev_sq) * sqrt(1 / (phi * (1 - h_sim)))
#   e[, i] <- sort(res_sim)
# }


## ----eval=FALSE----------------------------------------------------------------
# library(hnp)
# 
# # O pacote calcula internamente |r_i| e extrai os quantis
# # resid.type = "deviance" define o tipo de resíduo basal
# hnp_obj <- hnp(lrfit,
#                resid.type = "deviance",
#                how.many.out = TRUE,
#                scale = TRUE,
#                sim = 99)
# 
# # A saída principal é puramente gráfica (Gráfico de envelopes)
# # Valores fora da banda sugerem sobredispersão ou falta de ajuste


## ----eval=FALSE----------------------------------------------------------------
# library(DHARMa)
# 
# # simulateResiduals executa a simulação e calcula o PIT via eCDF
# # NOTA: O modelo original NÃO é reajustado M vezes.
# # Isso torna o método computacionalmente muito mais leve.
# res_dharma <- simulateResiduals(fittedModel = lrfit,
#                                 n = 1000)
# 
# # Plota QQ-plot Uniforme e Resíduos vs Previstos
# plot(res_dharma)
# 
# # Testes formais baseados na Uniformidade (u_i)
# testDispersion(res_dharma)
# testUniformity(res_dharma)


## ----eval=TRUE, echo=FALSE-----------------------------------------------------
m <- 1000 
e <- matrix(0, n, m)
mu_orig <- fitted(lrfit); fam <- binomial()
phi <- deviance(lrfit) / (nrow(cuse) - ncol(X)) # Fora do loop!
size <- (cuse$using+cuse$notUsing)

for(i in 1:m){
  resp <- rbinom(n, size=size, prob=mu_orig)
  prop <- resp / size
  # Ajuste sem overhead de fórmula
  fit_sim <- glm.fit(x = X, y = prop, weights = size, family = fam)
  # Leverage rápido via matriz ortogonal Q
  h_sim <- rowSums(qr.Q(fit_sim$qr)^2)
  # Cálculo vetorizado do resíduo deviance
  dev_sq <- fam$dev.resids(prop, fit_sim$fitted.values, size)
  sgn <- sign(prop - fit_sim$fitted.values)
  res_sim <- sgn * sqrt(dev_sq) * sqrt(1 / (phi * (1 - h_sim)))
  e[, i] <- sort(res_sim)
}

library(hnp)
hnp_obj <- hnp(lrfit, 
               resid.type = "deviance", 
               how.many.out = TRUE, 
               scale = TRUE,
               sim = 999, plot.sim=FALSE, print.on=FALSE)

library(DHARMa)
res_dharma <- simulateResiduals(fittedModel = lrfit,n = 1000)


## ----eval=TRUE, echo=FALSE, results='hide', fig.width=10, fig.height=3.5, out.width="100%"----
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# 1. Manual (Desvio Padronizado)
faixa <- range(td, e)
qqnorm(td, main = "Manual (Desvio)", ylab = "Resíduos", ylim = faixa)
lines(qnorm(ppoints(n)), apply(e, 1, quantile, 0.025), lty = 1)
lines(qnorm(ppoints(n)), apply(e, 1, quantile, 0.975), lty = 1)

# 2. hnp (Half-Normal Plot)
plot(hnp_obj, main = "hnp (Meio-Normal)", ylab = "Resíduos Absolutos")

# 3. DHARMa (Uniform QQ-Plot)
plotQQunif(res_dharma, main = "DHARMa (PIT Uniforme)", 
           testUniformity = FALSE, testOutliers = FALSE, 
           testDispersion = FALSE)

