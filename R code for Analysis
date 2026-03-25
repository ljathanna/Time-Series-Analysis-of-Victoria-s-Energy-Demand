rm(list = ls())

library(TSA)

setwd("C:/Users/DELL/OneDrive/Desktop/Sem 3/Time Series Analysis")

vic_energy <- read.csv("Victoria_Energy_Data.csv", header= TRUE)

#Checking whether our data is in the right form
head(vic_energy)
str(vic_energy)

colnames(vic_energy)[2] <- "Demand"
vic_energyTS <- ts(vic_energy$Demand, start=1,frequency=1)

plot(vic_energyTS, ylab='Demand (GW/100)',xlab='Minute', type='o',
     main="Time Series Plot for the Energy Demand in Victoria")

#Fitting the model
t <- time(vic_energyTS)

#Building the function for linear/quadratic models
model_diagnostics <- function(model, model_name) {
  std_res <- rstudent(model)
  
  plot(ts(std_res),
       ylab = "Standardized Residuals",
       xlab = "Time",
       type = "o",
       main = paste("Standardized Residuals Time Series Plot -", model_name))
  
  hist(std_res,
       main = paste("Histogram -", model_name),
       xlab = "Standardized Residuals")
  
  qqnorm(std_res,
         main = paste("QQ Plot -", model_name))
  qqline(std_res, col = 2, lwd = 2, lty = 2)
  
  print(shapiro.test(std_res))
  
  acf(std_res,
      main = paste("ACF of Standardized Residuals -", model_name))
}

#Building the functions for cosine models
fit_cosine_model <- function(P_value, series, t_var)
  {
  cos.t <- cos(2 * pi * t_var / P_value)
  sin.t <- sin(2 * pi * t_var / P_value)
  model <- lm(series ~ cos.t + sin.t)
  plot(ts(fitted(model)),
       ylab = "Demand (GW/100)",
       xlab = "Minute",
       main = paste("Cosine Trend Fit (P =", P_value, ")"),
       ylim = c(min(c(fitted(model), as.vector(series))),
                max(c(fitted(model), as.vector(series)))))
  lines(as.vector(series), type = "o")
  return(model)
}

#Building the function for combined models
fit_combined_model <- function(P_value, series, t_var, t2_var)
  {
  cos.t <- cos(2 * pi * t_var / P_value)
  sin.t <- sin(2 * pi * t_var / P_value)
  model <- lm(series ~ t_var + t2_var + cos.t + sin.t)
  plot(ts(fitted(model)),
       ylab = "Demand (GW/100)",
       xlab = "Minute",
       main = paste("Combined Quadratic & Cosine Trend (P =", P_value, ")"),
       ylim = c(min(c(fitted(model), as.vector(series))),
                max(c(fitted(model), as.vector(series)))))
  lines(as.vector(series), type = "o")
  return(model)
}

#Linear model
model_linear <- lm(vic_energyTS ~t)
summary(model_linear)
plot(vic_energyTS, ylab='Demand (GW/100)', xlab='Minute',type='o',
     main="Linear Trend Fit")
abline(model_linear, lwd=2)

model_diagnostics(model_linear, "Linear Model")
#not best fit because its oscillating

#Quadratic Model
t2 = t^2
model_quadratic= lm(vic_energyTS ~ t+t2)
summary(model_quadratic)
plot(ts(fitted(model_quadratic)), ylab = "Demand (GW/100)", xlab = "Minute",
     main = "Quadratic Trend Fit",
     ylim = c(min(c(fitted(model_quadratic), as.vector(vic_energyTS))),
              max(c(fitted(model_quadratic), as.vector(vic_energyTS)))))
lines(as.vector(vic_energyTS), type = "o")

model_diagnostics(model_quadratic, "Quadratic Model")

#Seasonal Model
#We create a dummy model for this data set and assume a 15 min cycle
vic_energyTS_season <- ts(vic_energy$Demand, start = 1, frequency = 15)
season_factor <- factor(cycle(vic_energyTS_season))
model_seasonal_15 <- lm(vic_energyTS_season ~ season_factor - 1)
summary(model_seasonal_15)

plot(ts(fitted(model_seasonal_15)),
     ylab = "Demand (GW/100)",
     xlab = "Minute",
     main = "Seasonal Model ",
     ylim = c(min(c(fitted(model_seasonal_15), as.vector(vic_energyTS_season))),
              max(c(fitted(model_seasonal_15), as.vector(vic_energyTS_season)))),
     col = "red")
lines(as.vector(vic_energyTS_season), type = "o")

model_diagnostics(model_seasonal_15, "Seasonal Model")


#Cosine / Harmonic Model
model_cos_10 <- fit_cosine_model(10, vic_energyTS, t)
summary(model_cos_10)

model_cos_12 <- fit_cosine_model(12, vic_energyTS, t)
summary(model_cos_12)

model_cos_15 <- fit_cosine_model(15, vic_energyTS, t)
summary(model_cos_15)

model_cos_20 <- fit_cosine_model(20, vic_energyTS, t)
summary(model_cos_20)

#From the summary data we can see that P=15 gives us the best data so we continue the rest of the modelling with P=15
model_diagnostics(model_cos_15, "Cosine Model")

#Combined Harmonic Trend (Quadratic + Cosine)
model_combined_15 <- fit_combined_model(15, vic_energyTS, t, t2)
summary(model_combined_15)

model_combined_13 <- fit_combined_model(13, vic_energyTS, t, t2)
summary(model_combined_13)

model_combined_17 <- fit_combined_model(17, vic_energyTS, t, t2)
summary(model_combined_17)

# Best combined model chosen for diagnostics
model_diagnostics(model_combined_15, "Combined Model")

#Comparison Table
comparison_models <- data.frame(
  Model = c("Linear", "Quadratic","Seasonal", "Cosine", "Combined"),
  Adjusted_R2 = c(summary(model_linear)$adj.r.squared,
                  summary(model_quadratic)$adj.r.squared,
                  summary(model_seasonal_15)$adj.r.squared,
                  summary(model_cos_15)$adj.r.squared,
                  summary(model_combined_15)$adj.r.squared),
  Residual_SE = c(summary(model_linear)$sigma,
                  summary(model_quadratic)$sigma,
                  summary(model_seasonal_15)$sigma,
                  summary(model_cos_15)$sigma,
                  summary(model_combined_15)$sigma),
  AIC = c(AIC(model_linear),
          AIC(model_quadratic),
          AIC(model_seasonal_15),
          AIC(model_cos_15),
          AIC(model_combined_15))
)

comparison_models

selected_model <- model_combined_15

t_forecast <- 153:167
t2_forecast <- t_forecast^2

fitted_values <- fitted(selected_model)

P<- 15
cos_forecast <-  cos(2*pi*t_forecast/P)
sin_forecast <- sin(2*pi*t_forecast/P)

forecast_15 <- predict(selected_model,
                       newdata = data.frame(
                         t_var = t_forecast,
                         t2_var = t2_forecast,
                         cos.t = cos_forecast,
                         sin.t = sin_forecast
                       ),
                       interval = "prediction")

forecast_table <- data.frame(Minute = t_forecast,
                  Forecast_demand = forecast_15)

forecast_table

par(mar = c(5, 4, 4, 8))
plot(c(1:length(vic_energyTS), t_forecast),
     c(as.vector(vic_energyTS), forecast_15[,1]),
     type = "n",
     xlab = "Minute",
     ylab = "Demand (GW/100)",
     main = "15-Minute Forecast")
lines(1:length(vic_energyTS), vic_energyTS,
      type = "o", col = "black")
lines(t_forecast, forecast_15[,1],
      col = "red", lwd = 2)
lines(t_forecast, forecast_15[,2],
      col = "blue", lty = 2)
lines(t_forecast, forecast_15[,3],
      col = "blue", lty = 2)
legend("bottomleft",
       legend = c("Actual", "Forecast", "Prediction Interval"),
       col = c("black", "red", "blue"),
       lty = c(1,1,2),
       lwd = c(1,2,1),
       pch=1,
       text.width=18,
       bty = "n")
       
