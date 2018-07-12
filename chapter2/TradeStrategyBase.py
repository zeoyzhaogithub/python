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
        # 持有时间
        self.keep_stock_day = 0
        # 7% 上涨策略作为买入阈值
        self.__buy_change_threshold = 0.07

    def buy_strategy(self, trade_ind, trade_day, trade_days):
        '''
        买入
        :param trade_ind:
        :param trade_day:
        :param trade_days:
        :return:
        '''
        if self.keep_stock_day == 0 and trade_day.change > self.__buy_change_threshold:
            # 买入
            self.keep_stock_day += 1
        elif self.keep_stock_day > 0:
            # 已经持有，追加天数
            self.keep_stock_day += 1


    def sell_strategy(self, trade_ind, trade_day, trade_days):
        '''
        卖出
        :param trade_ind:
        :param trade_day:
        :param trade_days:
        :return:
        '''
        if self.keep_stock_day >= TradeStrategy1.s_keep_stock_threshold:
            # 当持有股票天数超过阈值s_keep_stock_threshold,卖出股票
            self.keep_stock_day = 0


    '''
    property 属性
    '''
    @property
    def buy_change_threshold(self):
        '''
        私有变量
        getter 函数
        :return:
        '''
        return self.__buy_change_threshold

    @buy_change_threshold.setter  # 添加 setter函数
    def buy_change_threshold(self, buy_change_threshold):
        '''
        给示例增加除访问与修改之外的其他处理逻辑
        :param buy_change_threshold:
        :return:
        '''
        if not isinstance(buy_change_threshold, float):
            '''
            上涨阈值需要为float类型
            '''
            raise TypeError('buy_change_threshold must be float!')

        # 上涨阈值只取小数点后两位
        self.__buy_change_threshold = round(buy_change_threshold, 2)
        
