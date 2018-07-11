#! /usr/bin/python
# -*-coding: utf-8-*-


import sys

import six
from abc import ABCMeta,abstractmethod




reload(sys)
sys.setdefaultencoding("utf-8")

class TradeStrategyBase(six.with_metaclass(ABCMeta, object)):
    '''
    交易策略抽象基类
    '''

    @abstractmethod
    def buy_strategy(self, *args, **kwargs):
        '''
        卖出策略基类
        :param args: 任意数量参数 list
        :param kwargs: 任意数量的关键字参数 dict
        :return:
        '''
        pass

    @abstractmethod
    def sell_strategy(self, *args, **kwargs):
        '''
        卖出策略基类
        :param args:
        :param kwargs:
        :return:
        '''
        pass


class TradeStrategy1(TradeStrategyBase):
    '''
    交易策略1：追涨策略，当股价上涨一个阈值默认为7%时
    买入股票并持有 s_keep_stock_threshold(20)天
    '''
    s_keep_stock_threshold = 20

    def __init__(self):
        self.keep_stock_day = 0
        
