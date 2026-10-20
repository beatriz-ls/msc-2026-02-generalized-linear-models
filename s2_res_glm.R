################################################################################
## ANÁLISE DE DIAGNÓSTICO E INFLUÊNCIA LOCAL PARA MODELOS LINEARES GENERALIZADOS
## Dados: Uso de Contraceptivos (Fiji Fertility Survey)
################################################################################

# 1. Leitura e Preparação dos Dados --------------------------------------------
cuse <- read.table("cuse.dat", header=TRUE)

# Criando variáveis indicadoras booleanas
cuse$noMore <- cuse$wantsMore == "no"
cuse$hiEduc <- cuse$education == "high"

# Ajuste do Modelo Binomial (Interação idade * não quer mais filhos + educação)
lrfit <- glm(cbind(using, notUsing) ~ age * noMore + hiEduc, 
             family = binomial, data = cuse)

# Resumo do modelo
summary(lrfit)


# 2. Matriz de Projeção (H) e Medidas Iniciais ---------------------------------
X <- model.matrix(lrfit)     # Matriz do modelo (Covariáveis)
w <- lrfit$weights           # Pesos do MLG na convergência
W <- diag(w)                 # Matriz diagonal de pesos W
n <- nrow(cuse)

# Cálculo da Matriz de Projeção (H) aproximada para MLG
# H = W^{1/2} * X * (X'WX)^{-1} * X' * W^{1/2}
H <- sqrt(W) %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% sqrt(W)
h <- diag(H)                 # Elementos da diagonal principal (Alavancagem - hii)

# Parâmetro de dispersão (phi estimado)
fi <- deviance(lrfit)/(nrow(cuse) - ncol(X))


# 3. Análise de Resíduos (Deviance e Pearson) ----------------------------------
# Resíduos originais
rd <- residuals(lrfit, type = "deviance")
rp <- residuals(lrfit, type = "pearson")

# Resíduos Padronizados (Studentizados)
td <- rd / sqrt(fi * (1 - h))        # Deviance Padronizado
ts <- rp / sqrt(fi * w * (1 - h))    # Pearson Padronizado


# 4. Medida de Influência Global: Distância de Cook (LD) -----------------------
# Fórmula aproximada baseada na expansão de Taylor de segunda ordem
LD <- h * ts^2 / (1 - h)

# Gráfico da Distância de Cook
par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(1:n, LD, type = "h", lwd = 2, col = "blue",
     main = "Distância de Cook (LD)", xlab = "Índice da Observação", ylab = "LD")
points(1:n, LD, pch = 16, col = "blue")
# Adicionando um limiar prático de atenção e textos
limiar_cook <- 4/n
abline(h = limiar_cook, col = "red", lty = 2)
outliers_cook <- which(LD > limiar_cook)
if(length(outliers_cook) > 0) text(outliers_cook, LD[outliers_cook], labels=outliers_cook, pos=3, cex=0.8, col="red")


# 5. Influência Local: Avaliação Individual (Análogo ao DFFITS) ----------------
# Mede o impacto ao perturbar a i-ésima observação isoladamente
# Ci = 2 * hii * rPi
Ci <- 2 * h * rp

# Limiares recomendados (Média +/- 2 desvios padrões)
limiar_sup <- mean(Ci) + 2*sd(Ci)
limiar_inf <- mean(Ci) - 2*sd(Ci)

plot(1:n, Ci, type = "h", lwd = 2, col = "darkgreen", ylim = c(min(Ci,limiar_inf),max(Ci,limiar_sup))*1.1,
     main = "Curvatura Ci (Análogo DFFITS)", xlab = "Índice da Observação", ylab = "Ci")
points(1:n, Ci, pch = 16, col = "darkgreen")
abline(h = c(limiar_sup, limiar_inf), col = "red", lty = 2)

outliers_ci <- which(Ci > limiar_sup | Ci < limiar_inf)
if(length(outliers_ci) > 0) text(outliers_ci, Ci[outliers_ci], labels=outliers_ci, pos=3, cex=0.8, col="red")


# 6. Influência Local: Direção de Maior Curvatura Global (amax) ----------------
# A direção a_max detecta onde uma pequena perturbação causa mais danos globais
# B = D(rp) * H * D(rp)
B <- diag(rp) %*% H %*% diag(rp)
Cmax <- eigen(B)$val[1]              # Maior autovalor (Maior Curvatura)
amax <- abs(eigen(B)$vec[, 1])       # Autovetor correspondente

plot(1:n, amax, type = "h", lwd = 2, col = "purple", ylim=c(0,max(amax)*1.1),
     main = "Direção de Maior Influência Global (|amax|)", xlab = "Índice da Observação", ylab = "|amax|")
points(1:n, amax, pch = 16, col = "purple")
# Destaque empírico para observações muito proeminentes
limiar_amax <- mean(amax) + 2*sd(amax)
abline(h = limiar_amax, col = "red", lty = 2)
outliers_amax <- which(amax > limiar_amax)
if(length(outliers_amax) > 0) text(outliers_amax, amax[outliers_amax], labels=outliers_amax, pos=3, cex=0.8, col="red")


# 7. Influência Local Específica: Avaliando Coeficientes (Análogo ao DFBETAS) --
# Perturbação focada em um parâmetro específico. 
# Exemplo: Analisar a influência na variável "hiEducTRUE"
# Acessando os nomes das colunas para achar o índice automaticamente
nome_var <- "hiEducTRUE"
i <- which(colnames(X) == nome_var) 

if(length(i) > 0) {
  Xi <- as.matrix(X[, i])     # Coluna da variável de interesse
  X_i <- X[, -i]              # Restante das matriz X
  
  # Regressão linear da variável X_i nas demais (ortogonalização)
  # v = W^{1/2}Xi - W^{1/2}X_{-i}(X'_{-i}WX_{-i})^{-1}X'_{-i}W Xi
  proj_Xi <- sqrt(W) %*% X_i %*% solve(t(X_i) %*% W %*% X_i) %*% t(X_i) %*% W %*% Xi
  v <- sqrt(W) %*% Xi - proj_Xi
  
  # Recalculando o amax específico para a direção ortogonalizada
  amax_var <- abs(v * rp / sqrt(Cmax))

  limiar_var <- mean(amax_var) + 2*sd(amax_var)
  
  plot(1:n, amax_var, type = "h", lwd = 2, col = "darkorange", ylim=c(0,max(amax_var,limiar_var)*1.1),
       main = paste("Influência em", nome_var, "(Análogo DFBETAS)"), 
       xlab = "Índice da Observação", ylab = "|amax| Específico")
  points(1:n, amax_var, pch = 16, col = "darkorange")
  abline(h = limiar_var, col = "red", lty = 2)
  outliers_var <- which(amax_var > limiar_var)
  if(length(outliers_var) > 0) text(outliers_var, amax_var[outliers_var], labels=outliers_var, pos=3, cex=0.8, col="red")
} else {
  plot(0, 0, type="n", main="Variável não encontrada", xlab="", ylab="")
  text(0, 0, paste("Variável", nome_var, "não encontrada em X"))
}

# Restaurar layout original dos gráficos
par(mfrow=c(1,1))