#! /usr/bin/python
# -*-coding: utf-8-*-


import sys

import datetime

import pandas
import pandas_datareader as web

reload(sys)
sys.setdefaultencoding("utf-8")

class StockAna:
    def __init__(self): 
        print("stock analysis")
        print(sys.version)

    def getStockInfo(self):
        '''
        获取股票数据
        :return:
        '''
        start = datetime.datetime(2017, 6, 5)
        end   = datetime.date.today()

        #goog = web.data.DataReader('GOOG', "yahoo", start, end)
        #print(goog)
        print(end)


if __name__ == '__main__':
    stockAna = StockAna()

    stockAna.getStockInfo()
    