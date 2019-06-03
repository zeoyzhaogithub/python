#! /usr/bin/python
# -*-coding: utf-8-*-


import sys

import logging
import pdb
from collections import namedtuple
from collections import OrderedDict
from functools import partial



reload(sys)
sys.setdefaultencoding("utf-8")
logging.basicConfig(level=logging.DEBUG)

# -------------------------
# logging 级别：
# debug: 调试信息
# info: 输出信息
# warning: 警告信息
# error: 错误信息
# critical: 严重错误
# --------------------------

def gen_buy_cahnge_lust():
    buy_change_list = []
    for buy_change in range(-5, -16, -1):
        buy_change = buy_change / 100.0
        buy_change_list.append(buy_change)
    return buy_change_list

print(gen_buy_cahnge_lust())


# print log
def gen_buy_cahnge_lust1():
    logging.info("gen_buy_cahnge_lust1 begin")
    buy_change_list = []
    for buy_change in range(-5, -16, -1):
        # logging.debug(buy_change)
        buy_change = buy_change / 100.0
        buy_change_list.append(buy_change)
    return buy_change_list

print(gen_buy_cahnge_lust1())


