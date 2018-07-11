#! /usr/bin/python
# -*-coding: utf-8-*-


import sys

from collections import namedtuple
from collections import OrderedDict
from collections import Iterable
from abupy import ABuSymbolPd
from functools import partial




reload(sys)
sys.setdefaultencoding("utf-8")

class StockTradDays(object):
    def __init__(self, price_array, start_date, date_array=None):
        '''
        初始化变量
        :param price_array:
        :param start_date: 开始时间
        :param date_array:
        '''
        # 私有价格序列
        self.__price_array = price_array
        # 私有日期序列
        self.__date_array = self._init_days(start_date, date_array)
        # 私有涨跌幅序列
        self.__change_array = self.__init_change()
        # 进行OrderedDict封装
        self.stock_dict = self._init_stock_dict()



    def __init_change(self):
        '''
        private
        从price_array生成change_array
        :return:
        '''
        price_float_array = [
            float(price_str) for price_str in self.__price_array
        ]
        # 开盘价和收盘价
        pp_array = [
            (price1, price2) for price1, price2 in zip(
                price_float_array[:-1], price_float_array[1:]
            )
        ]
        change_array = map(
            lambda pp: reduce(
                lambda a, b: round((b - a)/a, 3),
                pp
            ),
            pp_array
        )
        change_array.insert(0, 0)
        return change_array


    def _init_days(self, start_date, date_array):
        '''
        protect方法
        :param start_date: 初始日期
        :param date_array: 给定日期序列
        :return:
        '''

        if date_array is None:
            # 由 start_date 和 self.__price_array确定日期序列
            date_array = [
                str(start_date + ind) for ind, _ in enumerate(self.__price_array)
            ]
        else:
            date_array = [
                str(date) for date in date_array
            ]

        return date_array

    def _init_stock_dict(self):
        '''
        使用 namedtuple, OrderedDict将结果合并
        :return:
        '''
        stock_namedtuple = namedtuple(
            'stock', ('date', 'price', 'change')
        )

        # 使用以被赋值的 __date_array 进行组装
        stock_dict = OrderedDict(
            (date, stock_namedtuple(
                date, price, change
            ))
            for date, price, change in zip(
                self.__date_array, self.__price_array, self.__change_array
            )
        )

        return stock_dict


    def filter_stock(self, want_up=True, want_calc_sum=False):
        '''
        筛选结果子集
        :param want_up: 是否筛选上涨
        :param want_calc_sum: 是否计算涨跌幅和
        :return:
        '''

        filter_func = (lambda day: day.change > 0) if want_up else (lambda day: day.change < 0)
        # 使用filter_func 作为筛选函数
        want_days = filter(
            filter_func, self.stock_dict.values()
        )

        if not want_calc_sum:
            return want_days
        change_sum = 0.0
        for day in want_days:
            change_sum += day.change

        return change_sum



    def __str__(self):
        '''
        自定义__repr__ 和 __str__ 是为了简化调试和实例化输出
        :return:
        '''
        return str(self.stock_dict)
    __repr__ = __str__



    def __iter__(self):
        '''
        通过代理stock_dict的迭代，yield元素
        :return:
        '''
        for key in self.stock_dict:
            yield self.stock_dict[key]


    def __getitem__(self, ind):
        '''

        :param ind:
        :return:
        '''
        date_key = self.__date_array[ind]
        return self.stock_dict[date_key]


    def __len__(self):
        '''
        通过代理self.stock_dict的len()方法，实现对象长度查询
        :return:
        '''
        return len(self.stock_dict)


if __name__ == '__main__':

    # price_str = '30.14,29.58,26.36,32.56,32.82'
    # price_array = price_str.split(',')
    # # print(price_array)
    # date_base = 20170118
    #
    # trad_days = StockTradDays(price_array, date_base)
    # print(trad_days)
    # print('trad_days 对象的长度：{}'.format(len(trad_days)))
    # # 判断 trad_days是否可迭代
    # if isinstance(trad_days, Iterable):
    #     for day in trad_days:
    #         print(day)
    #
    # filter_stock = trad_days.filter_stock()
    # print(filter_stock)

    # 实战数据
    # 两年的TSLA 收盘数据
    price_array = ABuSymbolPd.make_kl_df('TSLA', n_folds=2).close.tolist()
    # print(price_array)
    date_base = 20170118
    date_array = ABuSymbolPd.make_kl_df('TSLA', n_folds=2).date.tolist()
    print(price_array[:5])
    print(date_array[:5])


