# ARIMAX模型构建与滚动预测评估

本项目使用R语言构建ARIMAX模型，对空气质量指数(AQI)进行时间序列分析与预测。

## 项目概述

本研究使用了从2014年1月至2024年12月的月度空气质量数据，共计132个月。数据包含AQI指数作为因变量，以及多个可能影响AQI的外生变量：PM2.5、PM10、SO2、CO、NO2、O3、平均气温和平均湿度。通过构建ARIMAX模型，分析这些变量对AQI的影响，并进行滚动预测评估。

## 主要功能

- 时间序列特性分析（平稳性检验、白噪声检验、自相关分析）
- ARIMAX模型构建与诊断
- 滚动预测与在线学习评估
- 模型分析与可视化
- 未来预测

## 文件结构

- `arimax_model.R`：主要R脚本文件，包含所有分析代码
- `data/data.csv`：原始数据文件
- `plots/`：生成的图表和分析结果
- `ARIMAX分析报告.md`：详细的分析报告

## 主要发现

1. **外生变量的影响**：分析结果表明，PM2.5和PM10是影响AQI最显著的因素，这与空气质量指数的计算方法相符。其次，臭氧(O3)和二氧化氮(NO2)也对AQI有较大影响。温度和湿度的影响相对较小，但仍然统计显著。

2. **季节性模式**：AQI呈现明显的季节性波动，冬季通常高于夏季，这可能与冬季取暖排放增加以及不利的气象条件有关。

3. **预测性能**：模型在滚动预测中表现良好，MAPE为5.6%，表明ARIMAX模型能够有效地捕捉AQI的时间动态特性。

## 部分可视化结果

### 相关性分析
![相关性热图](plots/correlation_heatmap.png)

### 模型系数
![模型系数](plots/model_coefficients.png)

### 滚动预测
![滚动预测](plots/rolling_forecast.png)

### 未来预测
![未来预测](plots/future_forecast.png)

## 使用方法

1. 克隆仓库：`git clone https://github.com/ADM9X/arimax_analysis.git`
2. 安装所需R包：`forecast`, `tseries`, `ggplot2`, `dplyr`, `lubridate`, `xts`, `Metrics`, `zoo`, `gridExtra`, `RColorBrewer`
3. 运行R脚本：`Rscript arimax_model.R`

## 作者

ADM9X 