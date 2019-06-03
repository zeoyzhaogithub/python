#! /usr/bin/python
# -*-coding: utf-8-*-


import sys

from collections import namedtuple
from collections import OrderedDict
from functools import partial



reload(sys)
sys.setdefaultencoding("utf-8")

date_array1 = []
date_base1 = 20170118
# 这里for只是为了计数，无用的变量PYTHON 建议用'_'声明

price_str = '30.14,29.58,26.36,32.56,32.82'
price_array = price_str.split(',')
for _ in xrange(0, len(price_array)):
    date_array1.append(str(date_base1))
    date_base1 += 1
# print(date_array1)

date_base = 20170118
date_array = [str(date_base + ind) for ind, _ in enumerate(price_array)]
#print(date_array)

# Zip 同时迭代多个序列
stock_tuple_list = [
    (date, price) for date, price in zip(date_array, price_array)
]

# print(stock_tuple_list[1])
# print(stock_tuple_list)

# 命名元祖
stock_namedtuple = namedtuple('stock', ('date', 'price'))
stock_tuple_list1 = [
    stock_namedtuple(date, price) for date, price in zip(date_array, price_array)
]

# print(stock_tuple_list1)
# print(stock_tuple_list1[1])

# 字典推导式
stock_dict = {
    date: price for date, price in zip(date_array, price_array)
}
# print(stock_dict)
# print(stock_dict.keys())
# 2 有序字典
stock_dict1 = OrderedDict(
    (date, price) for date, price in zip(date_array, price_array)
)
# print(stock_dict1)

# 3
price_min = min(zip(stock_dict1.values(), stock_dict1.keys()))
# print(price_min)

# 4第二大数据
find_second_max_lambda = lambda dict_array: sorted(
    zip(dict_array.values(), dict_array.keys())
)[-2]
stock_second = find_second_max_lambda(stock_dict1)
# print(stock_second)

# 5从收盘价格推导出每天的涨跌幅度
# 将字符串的价格通过列表推导式转换为float
price_float_array = [
    float(price_str) for price_str in stock_dict1.values()
]

# 通过将时间平移形成两个错开的收盘价序列
pp_array = [
    (price1, price2) for price1, price2 in zip(
        price_float_array[:-1], price_float_array[1:]
    )
]

# print(pp_array)
change_array = map(
    lambda pp: reduce(
        lambda a, b: round((b - a)/a, 3),
        pp
    ),
    pp_array
)

# 将第一天的涨跌幅设置为0
change_array.insert(0, 0)
# print(change_array)

# 将涨跌幅数据加入，重新构建stocj_dict
stock_namedtuple = namedtuple('stock', (
    'date', 'price', 'change'
))

# 通过Zip分别拿数据以date作为key组成OrderedDict
stock_dict = OrderedDict(
    (date, stock_namedtuple(date, price, change))
    for date, price, change in zip(
        date_array, price_array, change_array
    )
)
# print(stock_dict)
print(stock_dict.values())
# 上涨的交易日
up_days = filter(
    lambda day: day.change > 0,
    stock_dict.values()
)
# print(up_days)

# 筛选出上涨，下跌，相应数值
def filter_stock(stock_array_dict, want_up=True, want_calc_sum=False):
    # 类型不正确，抛出错误
    if not isinstance(stock_array_dict, OrderedDict):
        raise TypeError('stock_array_dict must be OrderedDict!')
    filter_func = (lambda day: day.change > 0) if want_up else (lambda day: day.change < 0)
    # 使用filter_func 作为筛选函数
    want_days = filter(
        filter_func, stock_array_dict.values()
    )

    if not want_calc_sum:
        return want_days
    change_sum = 0.0
    for day in want_days:
        change_sum += day.change

    return change_sum

print(filter_stock(stock_dict, want_up=False))

# 筛选上涨交易
filter_stock_up_days = partial(
    filter_stock, want_up=True, want_calc_sum=False
)
# 筛选下跌交易
filter_stock_down_days = partial(
    filter_stock, want_up=False, want_calc_sum=False
)
# 筛选上涨交易涨幅和
filter_stock_up_sums = partial(
    filter_stock, want_up=True, want_calc_sum=True
)
# 筛选下跌交易跌幅和
filter_stock_down_days = partial(
    filter_stock, want_up=False, want_calc_sum=True
)

print('所有上涨的交易日：{}'.format(filter_stock_up_days(stock_dict)))
