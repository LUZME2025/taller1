install.packages("tidyverse")
install.packages("rmdformats")
install.packages("ggExtra")
install.packages("pracma")
install.packages("ggpairs")
install.packages("plotly")
install.packages("plot3D")
install.packages("car")
install.packages("glmnet")

library(tidyverse)
library(ggpairs)
library(GGally) # para ggpairs
library(pracma) # para la función meshgrid
library(plotly) # para gráfico 3D
library(plot3D) # para gráfico 3D
library(car)    # para calcular VIF
library(glmnet) # sirve para ajustar modelos lineales generalizados (GLM)

--'1.'

taller1 <- read.csv("C:/Users/User/Downloads/taller1.txt/taller1.txt", header=TRUE)
View(taller1)
eigen_vals <- eigen(cor(taller1[, -1]))$values
min(eigen_vals)
max(eigen_vals)
max(eigen_vals) / min(eigen_vals)

--'2.'

set.seed(123)  # guardamos la semilla

# Índices de entrenamiento
idx_train <- sample(1:nrow(taller1), 1000)

# Separar conjuntos
train <- taller1[idx_train, ]
test  <- taller1[-idx_train, ]

# Verificar
dim(train)  # 1000 5001
dim(test)   # 200 5001


--'3.'

# Preparación
X_train <- as.matrix(train[, -1])  # genes
y_train <- train[, 1]              # variable

# Ridge (alpha = 0)
set.seed(123)
cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0, nfolds = 10)
lambda_r <- cv_ridge$lambda.min
lambda_r

# Lasso (alpha = 1)
set.seed(123)
cv_lasso <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 10)
lambda_l <- cv_lasso$lambda.min
lambda_l

--'4.'

# Ajuste Ridge
modelo_ridge <- glmnet(X_train, y_train, alpha = 0, lambda = lambda_r)

# Ajuste Lasso
modelo_lasso <- glmnet(X_train, y_train, alpha = 1, lambda = lambda_l)

# Verificar
coef(modelo_ridge)
coef(modelo_lasso)

coef_lasso <- coef(modelo_lasso)
sum(coef_lasso != 0) - 1  # el -1 es por el intercepto

--'5.'

# Preparar datos de prueba
X_test <- as.matrix(test[, -1])
y_test <- test[, 1]

# Predicciones
pred_ridge <- predict(modelo_ridge, newx = X_test)
pred_lasso <- predict(modelo_lasso, newx = X_test)

# ECM en datos de prueba
ecm_ridge <- mean((y_test - pred_ridge)^2)
ecm_lasso <- mean((y_test - pred_lasso)^2)

ecm_ridge
ecm_lasso

--'6.'

# Preparar todos los datos
X_full <- as.matrix(taller1[, -1])
y_full <- taller1[, 1]

# Ajuste Lasso con todos los 1200 datos
modelo_final <- glmnet(X_full, y_full, alpha = 1, lambda = lambda_l)

# Verificar cuántos genes seleccionó
coef_final <- coef(modelo_final)
sum(coef_final != 0) - 1

--'7'

# Traza de coeficientes Lasso para todos los 1200 datos
modelo_trazas <- glmnet(X_full, y_full, alpha = 1)
plot(modelo_trazas, xvar = "lambda", label = FALSE)
abline(v = -log(lambda_l), col = "red", lty = 2)


--'8'

