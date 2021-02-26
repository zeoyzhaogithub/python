#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-17 21:04:48
# @Author: zeoy
# rdkit Chem 和AllChem的区别

from rdkit.Chem import AllChem as Chem

Chem ： 负责基础常用的化学功能 （ 如 ： 读写分子 ， 子结构搜索 ， 分子美化等 ）
AllChem: 负责高级但不常用的化学功能 。

区分它们的目的是为了加速载入速度 。

也可以简化使用
