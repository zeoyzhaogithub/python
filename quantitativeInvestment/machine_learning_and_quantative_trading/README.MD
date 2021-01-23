# python 机器学习与量化投资

[http://www.topquant.vip](http://www.topquant.vip)

python 功能：
大数据
机器学习
自然语言处理
人脸识别
游戏开发
字体设计
电脑设计

1.1 scikit-learn 模块库

1.2 scikit-learn 六大功能

分类

回归

聚类

数据降维

模型选择

数据预处理

test_version 显示相关库及其对应的版本号

> 注：
若手动升级pandas 之后提示缺乏某种属性，无法使用，则需升级Dash为最新版本

经典的机器学习算法：

线性回归，函数名：LinearRegression

朴素贝叶斯算法，函数名：Multinomialnb

KNN近邻算法，函数名：KNeighborsClassifier

逻辑回归算法，函数名：LogisticRegression

随机森林算法，函数名：RandomForestClassifier

决策树算法，函数名：tree.DecisionTreeClassifier

迭代决策树算法,又叫MATA(Mulitple Additive Regression Tree)，函数名：GradientBostingClassifier

SVM算法，函数名：SVC

# 金融数据预处理

数据归一化：
把数变为（0，1）之间的小数
把有量纲到达式变为无量纲表达式

pandas 在其独立出来的pandas_datareader模块库中集成了多个著名的金融数据网站，一站式提供了股票期货、外汇黄金等多种金融数据。

pandas_datareader模块库，内置的金融数据网站包括：Yahoo、Enigma、FRED、Google、MOEX、Morningstar、Quandl、Stooq、NASDAQ、IEX等。

pandas_datareader模块库[!源码](https://pydata.github.io/pandas-datareader/stable/)

量化分析三大经典指标

夏普比率(Sharpe Ratio):简称SR，投资回报与风险的比率
詹森指数(Jensen):综合考虑基金收益和风险因素
特雷诺指数(Treynor Ratio):简称TR，每单位风险获得收益的溢价。
