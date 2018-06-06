#! /usr/bin/python
# -*-coding: utf-8-*-


import sys

#import pandas_datareader

reload(sys)
sys.setdefaultencoding("utf-8")

class StockAna:
    def __init__(self): 
        print("stock analysis")
        print(sys.version)

    def getStockInfo(self):
        print()


if __name__ == '__main__':
    stockAna = StockAna()
    