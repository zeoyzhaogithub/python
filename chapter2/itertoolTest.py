#! /usr/bin/python
# -*- coding:utf-8 -*-

import sys

import itertools



reload(sys)
sys.setdefaultencoding("utf-8")

items = [1, 2, 3]
'''
permutations() 顺序组合元素
'''
for item in itertools.permutations(items):
    print(item)

'''
combinations() 不考虑顺序，不放回元素
'''
print()
for item in itertools.combinations(items, 2):
    print(item)

'''
combinations_with_replacement() 不考虑顺序，放回元素
'''
print()
for item in itertools.combinations_with_replacement(items, 2):
    print(item)


'''
product() 笛卡尔积，针对多个输入序列进行排列组合
'''
print()
a1 = ['a', 'b', 'c']
a2 = [1, 2, 3]
a3 = [3, 2, 1]
for item in itertools.product(a1, a2, a3):
    print(item)


