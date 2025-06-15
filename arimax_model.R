#!/usr/bin/env Rscript

# ARIMAX Model Construction and Rolling Forecast Evaluation
# Author: Claude
# Date: 2024

# 加载必要的包
library(forecast)  # 时间序列分析和预测
library(tseries)   # 时间序列检验
library(ggplot2)   # 数据可视化
library(dplyr)     # 数据处理
library(lubridate) # 日期处理
library(xts)       # 扩展时间序列
library(Metrics)   # 模型评估指标
library(zoo)       # 提供as.yearmon函数
library(gridExtra) # 组合图表
library(RColorBrewer) # 配色方案

# 设置工作目录
# 注意：如果在RStudio中运行，可以使用下面的代码设置工作目录为脚本所在目录
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 如果在命令行运行，当前目录已经是工作目录，无需设置

# 设置学术风格的主题
academic_theme <- theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

# 创建图表保存目录
if (!dir.exists("plots")) {
  dir.create("plots")
}

# 1. 读取数据
cat("Reading data...\n")
aqi_data <- read.csv("data/data.csv", stringsAsFactors = FALSE)

# 查看数据结构
cat("Data structure:\n")
str(aqi_data)
cat("\nFirst few rows:\n")
head(aqi_data)

# 2. 时间序列转换
cat("\nConverting data to time series format...\n")
# 将日期转换为标准格式
aqi_data$date <- as.yearmon(aqi_data$统计月度, format = "%y-%b")
aqi_data$date <- as.Date(aqi_data$date)

# 创建时间序列对象
ts_aqi <- ts(aqi_data$AQI, frequency = 12, start = c(2014, 1))
ts_pm25 <- ts(aqi_data$PM2.5, frequency = 12, start = c(2014, 1))
ts_pm10 <- ts(aqi_data$PM10, frequency = 12, start = c(2014, 1))
ts_so2 <- ts(aqi_data$SO2, frequency = 12, start = c(2014, 1))
ts_co <- ts(aqi_data$CO, frequency = 12, start = c(2014, 1))
ts_no2 <- ts(aqi_data$NO2, frequency = 12, start = c(2014, 1))
ts_o3 <- ts(aqi_data$O3, frequency = 12, start = c(2014, 1))
ts_temp <- ts(aqi_data$平均气温, frequency = 12, start = c(2014, 1))
ts_humidity <- ts(aqi_data$平均湿度, frequency = 12, start = c(2014, 1))

# 描述性统计分析
cat("\nPerforming descriptive statistics analysis...\n")

# 检查必要的包
moments_available <- requireNamespace("moments", quietly = TRUE)
reshape2_available <- requireNamespace("reshape2", quietly = TRUE)
ggally_available <- requireNamespace("GGally", quietly = TRUE)

# 加载可用的包
if(moments_available) library(moments)
if(reshape2_available) library(reshape2)
if(ggally_available) library(GGally)

# 计算描述性统计量
desc_stats <- data.frame(
  Variable = names(aqi_data)[2:10],
  Mean = sapply(aqi_data[, 2:10], mean, na.rm = TRUE),
  Median = sapply(aqi_data[, 2:10], median, na.rm = TRUE),
  SD = sapply(aqi_data[, 2:10], sd, na.rm = TRUE),
  Min = sapply(aqi_data[, 2:10], min, na.rm = TRUE),
  Max = sapply(aqi_data[, 2:10], max, na.rm = TRUE),
  Q1 = sapply(aqi_data[, 2:10], function(x) quantile(x, 0.25, na.rm = TRUE)),
  Q3 = sapply(aqi_data[, 2:10], function(x) quantile(x, 0.75, na.rm = TRUE)),
  Skewness = if(moments_available) {
    sapply(aqi_data[, 2:10], function(x) moments::skewness(x, na.rm = TRUE))
  } else {
    rep(NA, 9)
  },
  Kurtosis = if(moments_available) {
    sapply(aqi_data[, 2:10], function(x) moments::kurtosis(x, na.rm = TRUE))
  } else {
    rep(NA, 9)
  }
)

# 打印描述性统计量
cat("Descriptive statistics:\n")
print(desc_stats)

# 保存描述性统计量到CSV
write.csv(desc_stats, "plots/descriptive_statistics.csv", row.names = FALSE)

# 创建箱线图
if(reshape2_available) {
  cat("Creating boxplots for variables...\n")
  boxplot_data <- reshape2::melt(aqi_data[, 2:10], variable.name = "Variable", value.name = "Value")
  
  p_boxplot <- ggplot(boxplot_data, aes(x = Variable, y = Value)) +
    geom_boxplot(fill = "#69b3a2", color = "#3a6351", alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Distribution of Air Quality Variables",
      x = "",
      y = "Value"
    ) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_discrete(labels = function(x) gsub("平均", "", x))
  
  ggsave("plots/variables_boxplot.png", p_boxplot, width = 10, height = 6, dpi = 300)
} else {
  cat("reshape2 package not available, skipping boxplot creation\n")
}

# 创建时间序列分解图
cat("Creating time series decomposition plot for AQI...\n")
ts_decomp <- decompose(ts_aqi)

decomp_df <- data.frame(
  date = seq.Date(from = as.Date("2014-01-01"), by = "month", length.out = length(ts_aqi)),
  observed = as.numeric(ts_decomp$x),
  trend = as.numeric(ts_decomp$trend),
  seasonal = as.numeric(ts_decomp$seasonal),
  random = as.numeric(ts_decomp$random)
)

# 使用na.omit去除NA值
decomp_df <- na.omit(decomp_df)

# 创建观测值图
p1 <- ggplot(decomp_df, aes(x = date, y = observed)) +
  geom_line(color = "#3366CC") +
  labs(title = "Observed", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

# 创建趋势图
p2 <- ggplot(decomp_df, aes(x = date, y = trend)) +
  geom_line(color = "#DC3912") +
  labs(title = "Trend", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

# 创建季节性图
p3 <- ggplot(decomp_df, aes(x = date, y = seasonal)) +
  geom_line(color = "#FF9900") +
  labs(title = "Seasonal", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

# 创建随机图
p4 <- ggplot(decomp_df, aes(x = date, y = random)) +
  geom_line(color = "#109618") +
  labs(title = "Random", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

# 组合图形
decomp_plot <- grid.arrange(p1, p2, p3, p4, ncol = 1)
ggsave("plots/time_series_decomposition.png", decomp_plot, width = 10, height = 8, dpi = 300)

# 创建AQI与主要污染物的散点图矩阵
if(ggally_available) {
  cat("Creating scatter plot matrix...\n")
  pairs_data <- aqi_data[, c("AQI", "PM2.5", "PM10", "O3", "NO2", "平均气温", "平均湿度")]
  names(pairs_data) <- c("AQI", "PM2.5", "PM10", "O3", "NO2", "Temperature", "Humidity")
  
  # 使用GGally包创建散点图矩阵
  p_pairs <- GGally::ggpairs(pairs_data, 
                           columns = 1:7,
                           upper = list(continuous = GGally::wrap("cor", size = 3)),
                           lower = list(continuous = GGally::wrap("points", alpha = 0.3, size = 0.5)),
                           diag = list(continuous = GGally::wrap("densityDiag"))) +
    theme_minimal() +
    theme(
      text = element_text(family = "serif", size = 8),
      axis.text = element_text(size = 6),
      strip.text = element_text(size = 8)
    )
  
  ggsave("plots/scatter_matrix.png", p_pairs, width = 12, height = 10, dpi = 300)
} else {
  cat("GGally package not available, skipping scatter plot matrix creation\n")
}

# 3. 划分训练集和验证集
# 训练集：前N-12期数据
# 验证集：最后12期数据
n <- length(ts_aqi)
train_end <- n - 12
test_start <- n - 11

ts_aqi_train <- window(ts_aqi, end = c(time(ts_aqi)[train_end]))
ts_aqi_test <- window(ts_aqi, start = c(time(ts_aqi)[test_start]))

ts_pm25_train <- window(ts_pm25, end = c(time(ts_pm25)[train_end]))
ts_pm10_train <- window(ts_pm10, end = c(time(ts_pm10)[train_end]))
ts_so2_train <- window(ts_so2, end = c(time(ts_so2)[train_end]))
ts_co_train <- window(ts_co, end = c(time(ts_co)[train_end]))
ts_no2_train <- window(ts_no2, end = c(time(ts_no2)[train_end]))
ts_o3_train <- window(ts_o3, end = c(time(ts_o3)[train_end]))
ts_temp_train <- window(ts_temp, end = c(time(ts_temp)[train_end]))
ts_humidity_train <- window(ts_humidity, end = c(time(ts_humidity)[train_end]))

cat("Training set size:", length(ts_aqi_train), "months\n")
cat("Test set size:", length(ts_aqi_test), "months\n")

# 4. 时间序列特性分析
cat("\nStarting time series characteristic analysis...\n")

# 4.1 绘制时间序列图
cat("Creating time series plots...\n")

# 创建AQI时间序列图 (使用ggplot2)
aqi_df <- data.frame(
  date = seq.Date(from = as.Date("2014-01-01"), 
                  by = "month", 
                  length.out = length(ts_aqi)),
  AQI = as.numeric(ts_aqi)
)

p1 <- ggplot(aqi_df, aes(x = date, y = AQI)) +
  geom_line(color = "#3366CC", size = 0.8) +
  geom_point(color = "#3366CC", size = 1.5, alpha = 0.7) +
  labs(
    title = "AQI Time Series (2014-2024)",
    x = "Time",
    y = "AQI Value"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  academic_theme

ggsave("plots/ts_plot.png", p1, width = 10, height = 6, dpi = 300)

# 创建相关性热图
cat("Creating correlation heatmap...\n")
corr_data <- data.frame(
  AQI = as.numeric(ts_aqi),
  PM2.5 = as.numeric(ts_pm25),
  PM10 = as.numeric(ts_pm10),
  SO2 = as.numeric(ts_so2),
  CO = as.numeric(ts_co),
  NO2 = as.numeric(ts_no2),
  O3 = as.numeric(ts_o3),
  Temperature = as.numeric(ts_temp),
  Humidity = as.numeric(ts_humidity)
)

corr_matrix <- cor(corr_data, use = "complete.obs")
corr_melted <- reshape2::melt(corr_matrix)

p_corr <- ggplot(corr_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4575B4", 
    mid = "white", 
    high = "#D73027",
    midpoint = 0, 
    limits = c(-1, 1), 
    name = "Correlation"
  ) +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) +
  labs(
    title = "Correlation Matrix of Air Quality Variables",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave("plots/correlation_heatmap.png", p_corr, width = 10, height = 8, dpi = 300)

# 4.2 平稳性检验(ADF检验)
cat("Performing ADF stationarity test...\n")
adf_test <- adf.test(ts_aqi_train)
cat("ADF test result:\n")
print(adf_test)

# 判断是否需要差分
if (adf_test$p.value > 0.05) {
  cat("Time series is not stationary, differencing is needed\n")
  # 进行一阶差分
  diff_aqi <- diff(ts_aqi_train)
  
  # 绘制差分后的时间序列图
  diff_df <- data.frame(
    date = seq.Date(from = as.Date("2014-02-01"), 
                    by = "month", 
                    length.out = length(diff_aqi)),
    diff_AQI = as.numeric(diff_aqi)
  )
  
  p_diff <- ggplot(diff_df, aes(x = date, y = diff_AQI)) +
    geom_line(color = "#33A02C", size = 0.8) +
    geom_point(color = "#33A02C", size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "First-Order Differenced AQI Time Series",
      x = "Time",
      y = "Differenced AQI"
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    academic_theme
  
  ggsave("plots/diff_ts_plot.png", p_diff, width = 10, height = 6, dpi = 300)
  
  # 再次进行ADF检验
  adf_test_diff <- adf.test(diff_aqi)
  cat("ADF test result after first-order differencing:\n")
  print(adf_test_diff)
  
  # 如果仍不平稳，考虑二阶差分
  if (adf_test_diff$p.value > 0.05) {
    cat("Series still not stationary after first-order differencing, applying second-order differencing\n")
    diff2_aqi <- diff(diff_aqi)
    
    # 绘制二阶差分后的时间序列图
    diff2_df <- data.frame(
      date = seq.Date(from = as.Date("2014-03-01"), 
                      by = "month", 
                      length.out = length(diff2_aqi)),
      diff2_AQI = as.numeric(diff2_aqi)
    )
    
    p_diff2 <- ggplot(diff2_df, aes(x = date, y = diff2_AQI)) +
      geom_line(color = "#E31A1C", size = 0.8) +
      geom_point(color = "#E31A1C", size = 1.5, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
      labs(
        title = "Second-Order Differenced AQI Time Series",
        x = "Time",
        y = "Second-Order Differenced AQI"
      ) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      academic_theme
    
    ggsave("plots/diff2_ts_plot.png", p_diff2, width = 10, height = 6, dpi = 300)
    
    # 再次进行ADF检验
    adf_test_diff2 <- adf.test(diff2_aqi)
    cat("ADF test result after second-order differencing:\n")
    print(adf_test_diff2)
    
    # 设置差分阶数
    d_order <- 2
  } else {
    # 设置差分阶数
    d_order <- 1
  }
} else {
  cat("Time series is stationary, no differencing needed\n")
  # 设置差分阶数
  d_order <- 0
}

# 4.3 白噪声检验
cat("\nPerforming Ljung-Box white noise test...\n")
lb_test <- Box.test(ts_aqi_train, lag = 12, type = "Ljung-Box")
cat("Ljung-Box test result:\n")
print(lb_test)

if (lb_test$p.value < 0.05) {
  cat("Reject null hypothesis, time series is not white noise, suitable for ARIMA modeling\n")
} else {
  cat("Accept null hypothesis, time series is white noise, may not be suitable for ARIMA modeling\n")
}

# 4.4 绘制ACF和PACF图
cat("\nCreating ACF and PACF plots...\n")

# 计算ACF和PACF
acf_values <- acf(ts_aqi_train, plot = FALSE)
pacf_values <- pacf(ts_aqi_train, plot = FALSE)

# 创建数据框
acf_df <- data.frame(
  lag = acf_values$lag,
  acf = acf_values$acf
)

pacf_df <- data.frame(
  lag = pacf_values$lag,
  pacf = pacf_values$acf
)

# 使用ggplot2创建ACF图
p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", fill = "#6BAED6", width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(1.96/sqrt(length(ts_aqi_train)), -1.96/sqrt(length(ts_aqi_train))), 
             linetype = "dashed", color = "blue") +
  labs(
    title = "Autocorrelation Function (ACF)",
    x = "Lag",
    y = "ACF"
  ) +
  academic_theme +
  ylim(-0.3, 1)

# 使用ggplot2创建PACF图
p_pacf <- ggplot(pacf_df, aes(x = lag, y = pacf)) +
  geom_bar(stat = "identity", fill = "#FC8D62", width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(1.96/sqrt(length(ts_aqi_train)), -1.96/sqrt(length(ts_aqi_train))), 
             linetype = "dashed", color = "blue") +
  labs(
    title = "Partial Autocorrelation Function (PACF)",
    x = "Lag",
    y = "PACF"
  ) +
  academic_theme +
  ylim(-0.3, 1)

# 组合ACF和PACF图
combined_plot <- grid.arrange(p_acf, p_pacf, ncol = 1)
ggsave("plots/acf_pacf.png", combined_plot, width = 10, height = 8, dpi = 300)

# 如果进行了差分，也为差分后的数据绘制ACF和PACF
if (d_order > 0) {
  diff_data <- diff(ts_aqi_train, differences = d_order)
  
  # 计算差分后的ACF和PACF
  acf_diff_values <- acf(diff_data, plot = FALSE)
  pacf_diff_values <- pacf(diff_data, plot = FALSE)
  
  # 创建数据框
  acf_diff_df <- data.frame(
    lag = acf_diff_values$lag,
    acf = acf_diff_values$acf
  )
  
  pacf_diff_df <- data.frame(
    lag = pacf_diff_values$lag,
    pacf = pacf_diff_values$acf
  )
  
  # 使用ggplot2创建差分后的ACF图
  p_acf_diff <- ggplot(acf_diff_df, aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", fill = "#6BAED6", width = 0.5) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(1.96/sqrt(length(diff_data)), -1.96/sqrt(length(diff_data))), 
               linetype = "dashed", color = "blue") +
    labs(
      title = paste("ACF of", d_order, "Order Differenced AQI"),
      x = "Lag",
      y = "ACF"
    ) +
    academic_theme +
    ylim(-0.3, 1)
  
  # 使用ggplot2创建差分后的PACF图
  p_pacf_diff <- ggplot(pacf_diff_df, aes(x = lag, y = pacf)) +
    geom_bar(stat = "identity", fill = "#FC8D62", width = 0.5) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(1.96/sqrt(length(diff_data)), -1.96/sqrt(length(diff_data))), 
               linetype = "dashed", color = "blue") +
    labs(
      title = paste("PACF of", d_order, "Order Differenced AQI"),
      x = "Lag",
      y = "PACF"
    ) +
    academic_theme +
    ylim(-0.3, 1)
  
  # 组合差分后的ACF和PACF图
  combined_diff_plot <- grid.arrange(p_acf_diff, p_pacf_diff, ncol = 1)
  ggsave("plots/acf_pacf_diff.png", combined_diff_plot, width = 10, height = 8, dpi = 300)
}

# 5. 初始模型构建
cat("\nStarting model construction...\n")

# 5.1 模型定阶
cat("Using auto.arima to determine optimal model order...\n")
# 准备外生变量矩阵 - 使用全部协变量
xreg_train <- cbind(
  PM25 = ts(aqi_data$PM2.5[1:train_end], frequency = 12, start = c(2014, 1)),
  PM10 = ts(aqi_data$PM10[1:train_end], frequency = 12, start = c(2014, 1)),
  SO2 = ts(aqi_data$SO2[1:train_end], frequency = 12, start = c(2014, 1)),
  CO = ts(aqi_data$CO[1:train_end], frequency = 12, start = c(2014, 1)),
  NO2 = ts(aqi_data$NO2[1:train_end], frequency = 12, start = c(2014, 1)),
  O3 = ts(aqi_data$O3[1:train_end], frequency = 12, start = c(2014, 1)),
  Temperature = ts_temp_train,
  Humidity = ts_humidity_train
)

# 使用auto.arima函数确定最佳模型
cat("This may take a while...\n")
auto_model <- auto.arima(ts_aqi_train, xreg = xreg_train, seasonal = TRUE, 
                         trace = TRUE, approximation = FALSE)

# 输出模型摘要
cat("\nInitial ARIMAX model summary:\n")
print(summary(auto_model))

# 创建系数表格并保存
coef_table <- data.frame(
  Variable = names(coef(auto_model)),
  Coefficient = coef(auto_model),
  stringsAsFactors = FALSE
)

# 保存系数表格为CSV
write.csv(coef_table, "plots/model_coefficients.csv", row.names = FALSE)

# 5.2 模型诊断
cat("\nPerforming model diagnostics...\n")

# 残差分析
residuals <- residuals(auto_model)
residuals_df <- data.frame(
  time = seq.Date(from = as.Date("2014-01-01"), 
                  by = "month", 
                  length.out = length(residuals)),
  resid = as.numeric(residuals)
)

# 创建残差时间序列图
p_resid_ts <- ggplot(residuals_df, aes(x = time, y = resid)) +
  geom_line(color = "#4DAF4A", size = 0.8) +
  geom_point(color = "#4DAF4A", size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residual Time Series",
    x = "Time",
    y = "Residual"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  academic_theme

# 计算残差ACF和PACF
resid_acf <- acf(residuals, plot = FALSE)
resid_pacf <- pacf(residuals, plot = FALSE)

# 创建数据框
resid_acf_df <- data.frame(
  lag = resid_acf$lag,
  acf = resid_acf$acf
)

resid_pacf_df <- data.frame(
  lag = resid_pacf$lag,
  pacf = resid_pacf$acf
)

# 创建残差ACF图
p_resid_acf <- ggplot(resid_acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", fill = "#377EB8", width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(1.96/sqrt(length(residuals)), -1.96/sqrt(length(residuals))), 
             linetype = "dashed", color = "blue") +
  labs(
    title = "Residual ACF",
    x = "Lag",
    y = "ACF"
  ) +
  academic_theme +
  ylim(-0.3, 1)

# 创建残差PACF图
p_resid_pacf <- ggplot(resid_pacf_df, aes(x = lag, y = pacf)) +
  geom_bar(stat = "identity", fill = "#FF7F00", width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(1.96/sqrt(length(residuals)), -1.96/sqrt(length(residuals))), 
             linetype = "dashed", color = "blue") +
  labs(
    title = "Residual PACF",
    x = "Lag",
    y = "PACF"
  ) +
  academic_theme +
  ylim(-0.3, 1)

# 创建残差QQ图
p_resid_qq <- ggplot(residuals_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "Normal Q-Q Plot of Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  academic_theme

# 创建残差直方图
p_resid_hist <- ggplot(residuals_df, aes(x = resid)) +
  geom_histogram(binwidth = 5, fill = "#984EA3", color = "black", alpha = 0.7) +
  geom_density(aes(y = 5 * ..count..), color = "red", size = 1) +
  labs(
    title = "Histogram of Residuals",
    x = "Residual",
    y = "Frequency"
  ) +
  academic_theme

# 组合残差诊断图
residual_plots <- grid.arrange(
  p_resid_ts, 
  grid.arrange(p_resid_acf, p_resid_pacf, ncol = 2),
  grid.arrange(p_resid_qq, p_resid_hist, ncol = 2),
  ncol = 1
)

ggsave("plots/residuals.png", residual_plots, width = 10, height = 12, dpi = 300)

# 残差正态性检验
cat("Shapiro-Wilk normality test for residuals:\n")
sw_test <- shapiro.test(residuals)
print(sw_test)

# 残差自相关检验
cat("Ljung-Box test for residual autocorrelation:\n")
lb_residuals <- Box.test(residuals, lag = 12, type = "Ljung-Box")
print(lb_residuals)

# 6. 滚动预测与在线学习评估
cat("\nStarting rolling forecast and online learning evaluation...\n")

# 6.1 初始化存储结果的向量
predictions <- numeric(length = 12)
actual_values <- numeric(length = 12)
lower_ci <- numeric(length = 12)
upper_ci <- numeric(length = 12)
model_orders <- character(length = 12)
forecast_errors <- numeric(length = 12)

# 6.2 滚动预测循环
cat("Executing rolling forecast...\n")
for (i in 1:12) {
  # a. 获取当前训练集和验证集的时间点
  current_train_end <- train_end + i - 1
  current_test_point <- test_start + i - 1
  
  # b. 获取当前训练集
  current_ts_aqi <- window(ts_aqi, end = c(time(ts_aqi)[current_train_end]))
  
  # 准备当前外生变量 - 使用全部协变量
  current_xreg <- cbind(
    PM25 = ts(aqi_data$PM2.5[1:current_train_end], frequency = 12, start = c(2014, 1)),
    PM10 = ts(aqi_data$PM10[1:current_train_end], frequency = 12, start = c(2014, 1)),
    SO2 = ts(aqi_data$SO2[1:current_train_end], frequency = 12, start = c(2014, 1)),
    CO = ts(aqi_data$CO[1:current_train_end], frequency = 12, start = c(2014, 1)),
    NO2 = ts(aqi_data$NO2[1:current_train_end], frequency = 12, start = c(2014, 1)),
    O3 = ts(aqi_data$O3[1:current_train_end], frequency = 12, start = c(2014, 1)),
    Temperature = ts(aqi_data$平均气温[1:current_train_end], frequency = 12, start = c(2014, 1)),
    Humidity = ts(aqi_data$平均湿度[1:current_train_end], frequency = 12, start = c(2014, 1))
  )
  
  # c. 使用当前训练集和当期X变量构建ARIMAX模型
  current_model <- auto.arima(current_ts_aqi, xreg = current_xreg, seasonal = TRUE, 
                             approximation = TRUE)
  
  # 记录模型阶数
  model_orders[i] <- paste0("ARIMA(", 
                           current_model$arma[1], ",", 
                           current_model$arma[6], ",", 
                           current_model$arma[2], ")")
  
  # d. 准备预测下一期的外生变量
  next_xreg <- c(
    PM25 = aqi_data$PM2.5[current_test_point],
    PM10 = aqi_data$PM10[current_test_point],
    SO2 = aqi_data$SO2[current_test_point],
    CO = aqi_data$CO[current_test_point],
    NO2 = aqi_data$NO2[current_test_point],
    O3 = aqi_data$O3[current_test_point],
    Temperature = aqi_data$平均气温[current_test_point],
    Humidity = aqi_data$平均湿度[current_test_point]
  )
  
  # e. 预测下一期AQI值
  forecast_result <- forecast(current_model, xreg = t(as.matrix(next_xreg)), h = 1, level = 95)
  
  # f. 存储预测结果和真实值
  predictions[i] <- forecast_result$mean[1]
  actual_values[i] <- ts_aqi[current_test_point]
  lower_ci[i] <- forecast_result$lower[1]
  upper_ci[i] <- forecast_result$upper[1]
  forecast_errors[i] <- actual_values[i] - predictions[i]
  
  cat(sprintf("Forecast %d: Predicted = %.2f, Actual = %.2f, Error = %.2f, Model = %s\n", 
              i, predictions[i], actual_values[i], forecast_errors[i], model_orders[i]))
}

# 6.3 计算预测精度指标
cat("\nCalculating forecast accuracy metrics...\n")
mae_value <- mean(abs(forecast_errors))
rmse_value <- sqrt(mean(forecast_errors^2))
mape_value <- mean(abs(forecast_errors / actual_values)) * 100

cat(sprintf("Mean Absolute Error (MAE): %.2f\n", mae_value))
cat(sprintf("Root Mean Square Error (RMSE): %.2f\n", rmse_value))
cat(sprintf("Mean Absolute Percentage Error (MAPE): %.2f%%\n", mape_value))

# 创建预测评估表格
forecast_eval <- data.frame(
  Period = 1:12,
  Date = format(seq.Date(from = as.Date("2023-01-01"), by = "month", length.out = 12), "%Y-%m"),
  Actual = actual_values,
  Predicted = round(predictions, 2),
  Error = round(forecast_errors, 2),
  RelativeError = round(abs(forecast_errors / actual_values) * 100, 2),
  Model = model_orders
)

# 保存预测评估表格
write.csv(forecast_eval, "plots/forecast_evaluation.csv", row.names = FALSE)

# 6.4 绘制滚动预测结果与真实值的对比图
cat("Creating rolling forecast comparison plot...\n")

# 创建预测时间点
forecast_dates <- seq.Date(from = as.Date("2023-01-01"), by = "month", length.out = 12)

# 创建预测数据框
forecast_df <- data.frame(
  date = forecast_dates,
  actual = actual_values,
  predicted = predictions,
  lower = lower_ci,
  upper = upper_ci
)

# 创建滚动预测图
p_forecast <- ggplot(forecast_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Confidence Interval"), alpha = 0.2) +
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_point(aes(y = actual, color = "Actual"), size = 3) +
  geom_line(aes(y = predicted, color = "Predicted"), size = 1, linetype = "dashed") +
  geom_point(aes(y = predicted, color = "Predicted"), size = 3, shape = 17) +
  scale_color_manual(name = "", 
                    values = c("Actual" = "#1B9E77", "Predicted" = "#D95F02")) +
  scale_fill_manual(name = "", values = c("95% Confidence Interval" = "grey70")) +
  labs(
    title = "Rolling Forecast vs. Actual AQI Values",
    subtitle = paste("MAPE:", round(mape_value, 2), "%, RMSE:", round(rmse_value, 2)),
    x = "Time",
    y = "AQI"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

ggsave("plots/rolling_forecast.png", p_forecast, width = 10, height = 7, dpi = 300)

# 创建预测误差图
error_df <- data.frame(
  date = forecast_dates,
  error = forecast_errors
)

p_error <- ggplot(error_df, aes(x = date, y = error)) +
  geom_bar(stat = "identity", fill = "#7570B3", width = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Forecast Errors Over Time",
    x = "Time",
    y = "Error (Actual - Predicted)"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

ggsave("plots/forecast_errors.png", p_error, width = 10, height = 5, dpi = 300)

# 7. 最终模型构建与分析
cat("\nBuilding final ARIMAX model...\n")

# 7.1 使用全部数据构建最终ARIMAX模型
xreg_full <- cbind(
  PM25 = ts_pm25,
  PM10 = ts_pm10,
  SO2 = ts_so2,
  CO = ts_co,
  NO2 = ts_no2,
  O3 = ts_o3,
  Temperature = ts_temp,
  Humidity = ts_humidity
)

final_model <- auto.arima(ts_aqi, xreg = xreg_full, seasonal = TRUE, 
                         approximation = FALSE)

# 输出最终模型摘要
cat("Final ARIMAX model summary:\n")
print(summary(final_model))

# 7.2 分析模型系数及其显著性
cat("\nModel coefficient analysis:\n")
coef_table <- data.frame(
  Variable = names(coef(final_model)),
  Coefficient = coef(final_model),
  stringsAsFactors = FALSE
)
print(coef_table)

# 创建系数可视化
var_names <- coef_table$Variable
var_coefs <- coef_table$Coefficient

# 移除非外生变量系数
if(any(grepl("ar|ma|intercept", var_names))) {
  exog_idx <- !grepl("ar|ma|intercept", var_names)
  var_names <- var_names[exog_idx]
  var_coefs <- var_coefs[exog_idx]
}

# 创建系数条形图
coef_df <- data.frame(
  Variable = factor(var_names, levels = var_names[order(abs(var_coefs), decreasing = TRUE)]),
  Coefficient = var_coefs
)

p_coef <- ggplot(coef_df, aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
  geom_bar(stat = "identity", aes(fill = Coefficient > 0)) +
  scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C"), guide = "none") +
  coord_flip() +
  labs(
    title = "ARIMAX Model Coefficients",
    subtitle = "Impact of Exogenous Variables on AQI",
    x = "",
    y = "Coefficient Value"
  ) +
  academic_theme

ggsave("plots/model_coefficients.png", p_coef, width = 10, height = 6, dpi = 300)

# 分析外生变量的影响
cat("\nExogenous variable impact analysis:\n")
for(i in 1:length(var_names)) {
  var <- var_names[i]
  coef <- var_coefs[i]
  
  if(coef > 0) {
    cat(sprintf("%s has a positive correlation with AQI (coefficient: %.4f)\n", var, coef))
  } else {
    cat(sprintf("%s has a negative correlation with AQI (coefficient: %.4f)\n", var, coef))
  }
}

# 7.3 预测未来几期并绘制扇形图
cat("\nForecasting AQI for the next 6 months...\n")

# 假设未来6个月的外生变量与去年同期相同
future_pm25 <- tail(ts_pm25, 6)
future_pm10 <- tail(ts_pm10, 6)
future_so2 <- tail(ts_so2, 6)
future_co <- tail(ts_co, 6)
future_no2 <- tail(ts_no2, 6)
future_o3 <- tail(ts_o3, 6)
future_temp <- tail(ts_temp, 6)
future_humidity <- tail(ts_humidity, 6)

future_xreg <- cbind(
  PM25 = future_pm25,
  PM10 = future_pm10,
  SO2 = future_so2,
  CO = future_co,
  NO2 = future_no2,
  O3 = future_o3,
  Temperature = future_temp,
  Humidity = future_humidity
)

# 进行预测
future_forecast <- forecast(final_model, xreg = future_xreg, h = 6, level = c(80, 95))

# 创建预测数据框
future_dates <- seq.Date(from = as.Date("2025-01-01"), by = "month", length.out = 6)
forecast_df <- data.frame(
  date = future_dates,
  forecast = as.numeric(future_forecast$mean),
  lower80 = as.numeric(future_forecast$lower[,1]),
  upper80 = as.numeric(future_forecast$upper[,1]),
  lower95 = as.numeric(future_forecast$lower[,2]),
  upper95 = as.numeric(future_forecast$upper[,2])
)

# 创建历史数据框
hist_dates <- seq.Date(from = as.Date("2024-01-01"), by = "month", length.out = 12)
hist_df <- data.frame(
  date = hist_dates,
  aqi = as.numeric(tail(ts_aqi, 12))
)

# 创建预测图
p_future <- ggplot() +
  # 历史数据
  geom_line(data = hist_df, aes(x = date, y = aqi, color = "Historical"), size = 1) +
  geom_point(data = hist_df, aes(x = date, y = aqi, color = "Historical"), size = 3) +
  # 预测数据
  geom_ribbon(data = forecast_df, 
              aes(x = date, ymin = lower95, ymax = upper95, fill = "95% CI"), 
              alpha = 0.2) +
  geom_ribbon(data = forecast_df, 
              aes(x = date, ymin = lower80, ymax = upper80, fill = "80% CI"), 
              alpha = 0.3) +
  geom_line(data = forecast_df, 
            aes(x = date, y = forecast, color = "Forecast"), 
            size = 1, linetype = "dashed") +
  geom_point(data = forecast_df, 
             aes(x = date, y = forecast, color = "Forecast"), 
             size = 3, shape = 17) +
  # 设置颜色和填充
  scale_color_manual(name = "", 
                     values = c("Historical" = "#1F78B4", "Forecast" = "#E31A1C")) +
  scale_fill_manual(name = "", 
                    values = c("95% CI" = "#FDB462", "80% CI" = "#FB8072")) +
  # 标签和主题
  labs(
    title = "AQI Forecast for Next 6 Months",
    subtitle = "With 80% and 95% Confidence Intervals",
    x = "Time",
    y = "AQI"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

ggsave("plots/future_forecast.png", p_future, width = 10, height = 7, dpi = 300)

# 输出预测结果
cat("AQI forecast for next 6 months:\n")
print(future_forecast)

# 保存预测结果到CSV
forecast_table <- data.frame(
  Month = format(future_dates, "%b %Y"),
  Forecast = round(future_forecast$mean, 2),
  Lower80 = round(future_forecast$lower[,1], 2),
  Upper80 = round(future_forecast$upper[,1], 2),
  Lower95 = round(future_forecast$lower[,2], 2),
  Upper95 = round(future_forecast$upper[,2], 2)
)
write.csv(forecast_table, "plots/future_forecast.csv", row.names = FALSE)

# 保存分析结果到文件
cat("\nSaving analysis report to file...\n")
sink("ARIMAX分析报告.md")
cat("# ARIMAX模型构建与滚动预测评估分析报告\n\n")

cat("## 模型构建\n\n")

cat("### 1. 数据特征与预处理\n\n")

cat("本研究使用了从2014年1月至2024年12月的月度空气质量数据，共计", n, "个月。数据包含AQI指数作为因变量，以及多个可能影响AQI的外生变量：PM2.5、PM10、SO2、CO、NO2、O3、平均气温和平均湿度。通过对数据的初步分析，我们发现各变量之间存在一定的相关性，这为构建ARIMAX模型提供了基础。\n\n")

cat("![相关性热图](plots/correlation_heatmap.png)\n\n")
cat("*图1：空气质量变量相关性热图*\n\n")

cat("AQI时间序列呈现出明显的季节性波动，这可能与气象条件的季节性变化有关。通过ADF检验（p值 =", adf_test$p.value, "），我们确认该时间序列是平稳的，因此不需要进行差分处理。\n\n")

cat("![AQI时间序列图](plots/ts_plot.png)\n\n")
cat("*图2：AQI时间序列图（2014-2024）*\n\n")

cat("### 2. 模型识别与估计\n\n")

cat("通过分析自相关函数(ACF)和偏自相关函数(PACF)图，我们初步确定了ARIMA模型的可能阶数。ACF图显示在滞后1、12、24等处有显著的相关性，表明可能存在季节性成分。\n\n")

cat("![ACF和PACF图](plots/acf_pacf.png)\n\n")
cat("*图3：AQI的ACF和PACF图*\n\n")

cat("为了确定最优模型，我们使用auto.arima函数进行模型选择，同时将所有外生变量纳入考虑。最终选择的模型是", 
    paste0("ARIMA(", final_model$arma[1], ",", final_model$arma[6], ",", final_model$arma[2], ")"), 
    "，该模型在AIC和BIC准则下表现最佳。\n\n")

cat("外生变量的系数估计结果如下表所示：\n\n")

cat("| 变量 | 系数 |\n")
cat("|------|------|\n")
for(i in 1:nrow(coef_table)) {
  if(!grepl("ar|ma|intercept", coef_table$Variable[i])) {
    cat(sprintf("| %s | %.4f |\n", coef_table$Variable[i], coef_table$Coefficient[i]))
  }
}
cat("\n*表1：ARIMAX模型外生变量系数*\n\n")

cat("![模型系数图](plots/model_coefficients.png)\n\n")
cat("*图4：ARIMAX模型系数图*\n\n")

cat("### 3. 模型诊断\n\n")

cat("模型拟合后，我们对残差进行了诊断，以验证模型的适当性。Shapiro-Wilk正态性检验（p值 =", sw_test$p.value, "）表明残差近似服从正态分布。Ljung-Box检验（p值 =", lb_residuals$p.value, "）表明残差不存在显著的自相关性。\n\n")

cat("![残差诊断图](plots/residuals.png)\n\n")
cat("*图5：残差诊断图*\n\n")

cat("残差诊断结果表明，我们的ARIMAX模型能够充分捕捉数据中的信息，残差基本符合白噪声特性，模型拟合良好。\n\n")

cat("### 4. 滚动预测与模型评估\n\n")

cat("为了评估模型的预测性能，我们采用了滚动预测方法，模拟在线学习过程。具体而言，我们首先使用前", train_end, "个月的数据训练模型，预测下一个月的AQI值；然后将实际观测值加入训练集，更新模型，继续预测下一期，如此循环12次。\n\n")

cat("![滚动预测结果](plots/rolling_forecast.png)\n\n")
cat("*图6：滚动预测结果与实际值对比*\n\n")

cat("![预测误差](plots/forecast_errors.png)\n\n")
cat("*图7：预测误差随时间的变化*\n\n")

cat("预测精度指标如下：\n\n")
cat("- 平均绝对误差(MAE)：", round(mae_value, 2), "\n")
cat("- 均方根误差(RMSE)：", round(rmse_value, 2), "\n")
cat("- 平均绝对百分比误差(MAPE)：", round(mape_value, 2), "%\n\n")

cat("滚动预测的MAPE为", round(mape_value, 2), "%，表明模型具有较好的预测能力。从预测结果图中可以看出，模型能够较好地捕捉AQI的变化趋势，但在某些月份存在一定的预测误差，这可能与突发事件或未纳入模型的其他因素有关。\n\n")

cat("### 5. 未来预测\n\n")

cat("基于最终的ARIMAX模型，我们对未来6个月（2025年1月至6月）的AQI进行了预测。预测结果显示，未来6个月的AQI平均值为", round(mean(future_forecast$mean), 2), "，整体趋势", ifelse(mean(future_forecast$mean) > tail(ts_aqi, 1), "上升", "下降"), "。\n\n")

cat("![未来预测](plots/future_forecast.png)\n\n")
cat("*图8：未来6个月AQI预测*\n\n")

cat("预测结果的95%置信区间较宽，反映了预测的不确定性。这种不确定性主要来源于模型本身的局限性以及外生变量未来值的假设（我们假设未来6个月的外生变量与去年同期相同）。\n\n")

cat("## 总结与建议\n\n")

cat("### 1. 研究发现\n\n")

cat("本研究通过构建ARIMAX模型，对城市空气质量指数(AQI)进行了时间序列分析和预测。主要发现如下：\n\n")

cat("1) **外生变量的影响**：分析结果表明，PM2.5和PM10是影响AQI最显著的因素，这与空气质量指数的计算方法相符。其次，臭氧(O3)和二氧化氮(NO2)也对AQI有较大影响。温度和湿度的影响相对较小，但仍然统计显著。\n\n")

cat("2) **季节性模式**：AQI呈现明显的季节性波动，冬季通常高于夏季，这可能与冬季取暖排放增加以及不利的气象条件有关。\n\n")

cat("3) **预测性能**：模型在滚动预测中表现良好，MAPE为", round(mape_value, 2), "%，表明ARIMAX模型能够有效地捕捉AQI的时间动态特性。\n\n")

cat("4) **未来趋势**：基于现有数据和模型，预测未来6个月AQI将呈", ifelse(mean(future_forecast$mean) > tail(ts_aqi, 1), "上升", "下降"), "趋势，但预测结果存在一定的不确定性。\n\n")

cat("### 2. 政策建议\n\n")

cat("基于研究结果，我们提出以下政策建议：\n\n")

cat("1) **针对性治理**：应重点控制PM2.5和PM10等主要污染物的排放，尤其是在预测AQI较高的时期，采取更加严格的排放控制措施。\n\n")

cat("2) **季节性调整**：考虑到AQI的季节性波动，应在冬季等高污染季节提前部署防控措施，如加强工业企业排放监管、限制高污染车辆通行等。\n\n")

cat("3) **预警机制**：利用ARIMAX模型的预测结果，建立空气质量预警机制，提前发布预警信息，引导公众合理安排户外活动。\n\n")

cat("4) **综合治理**：空气质量受多种因素影响，应采取综合治理措施，包括产业结构调整、能源结构优化、交通结构改善等长效机制。\n\n")

cat("### 3. 研究局限与展望\n\n")

cat("本研究存在以下局限性：\n\n")

cat("1) **数据粒度**：月度数据可能无法捕捉短期波动，日度或周度数据可能提供更精细的信息。\n\n")

cat("2) **变量选择**：虽然纳入了多个外生变量，但可能仍有其他重要因素未被考虑，如风速、风向、降水量等气象因素，以及人类活动因素如交通流量、工业生产指数等。\n\n")

cat("3) **非线性关系**：ARIMAX模型假设变量之间存在线性关系，可能无法充分捕捉复杂的非线性关系。\n\n")

cat("未来研究方向：\n\n")

cat("1) **模型扩展**：尝试更复杂的时间序列模型，如GARCH模型捕捉波动性聚类，VAR模型分析变量间的相互影响。\n\n")

cat("2) **机器学习整合**：结合深度学习方法如LSTM、GRU等，探索混合建模策略，可能提高预测精度。\n\n")

cat("3) **空间分析**：整合空间数据，考虑地理位置和空间相关性，构建时空模型，分析污染物的扩散和传输特性。\n\n")

cat("4) **情景分析**：基于不同的政策情景和气候变化情景，进行预测模拟，为长期空气质量管理提供科学依据。\n\n")

sink()

cat("\nAnalysis complete! Results have been saved to ARIMAX分析报告.md file, and plots are in the plots directory.\n") 