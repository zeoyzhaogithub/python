#! /usr/bin/python
# -*- coding:utf-8 -*-

# 使用pandas_datareader模块库下载最新的外汇交易数据

from pandas_datareader import data as web

jpy = web.DataReader(
        #'USDEUR',   # 外汇交易对应代码
        'BTCUSD',
        'stooq'     # 数据源名称
    )
jpy.to_csv('./data/usdjpy2019.csv')
print(jpy.tail())