#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-01 12:58:47
# @Author: zeoy
# rdkit 子结构搜索

# 一、引入所需库
from rdkit import Chem

# 二、子结构搜索
>子结构搜索可以通过SMARTS匹配符完成。
## 2.1 判断是否有子结构
首先创建分子对象，然后定义匹配模式，最后判断是否有子结构。
匹配模式：支持Chem.MolFromSmarts() 和 Chem.MolFromSmiles() 两种形式。 Smarts的表现形式更加丰富。

mol = Chem.MolFromSmiles('c1ccccc1OC')
patt = Chem.MolFromSmarts('OC')
flag = mol.HasSubstructMatch(patt)
if flag:
    print('分子中包含基团 -OCH3')
else:
    print('分子中不包含基团 -OCH3')
    

# 分子中包含基团 -OCH3
## 2.2 获取第一个子结构对应的原子编号
m = Chem.MolFromSmiles('c1ccccc1OC')
patt = Chem.MolFromSmarts('OC')
flag = m.HasSubstructMatch(patt)
if flag:
    atomids = m.GetSubstructMatch(patt)
    print("匹配到的原子id:",atomids)
else:
    print("分子中不包含基团 -OCH3")
# 匹配到的原子id: (6, 7)
# 根据输出结果可知，6号原子对应的是甲氧基原子的O原子，7号原子对应的是甲氧基原子的C原子

## 2.3 获取对应所有子结构
m = Chem.MolFromSmiles('c1ccc(OC)cc1OC')
patt = Chem.MolFromSmarts('OC')
flag = m.HasSubstructMatch(patt)
if flag:
    atomids = m.GetSubstructMatches(patt)
    print("匹配到的原子id:",atomids)
else:
    print('分子中不包含基团 -OCH3')

# 匹配到的原子id: ((4, 5), (8, 9))

## 2.4 子结构搜索考虑手性

# > 不考虑手性的时候，有手性的分子可以匹配无手性的模式和错误手性的模式
# 考虑手性的时候，有手性信息的分子可以匹配无手性的模式
# 无手性信息的分子不能匹配有手性的模式。


### 2.4.1 不考虑手性
m = Chem.MolFromSmiles('CC[C@H](F)Cl')
print(m.HasSubstructMatch(Chem.MolFromSmiles('C[C@H](F)Cl')))
print(m.HasSubstructMatch(Chem.MolFromSmiles('C[C@@H](F)Cl')))
print(m.HasSubstructMatch(Chem.MolFromSmiles('CC(F)Cl')))
# True
# True
# True

### 2.4.2 考虑手性
m = Chem.MolFromSmiles('CC[C@H](F)Cl')
print(m.HasSubstructMatch(Chem.MolFromSmiles('C[C@H](F)Cl'), useChirality=True))
print(m.HasSubstructMatch(Chem.MolFromSmiles('C[C@@H](F)Cl'), useChirality=True))
print(m.HasSubstructMatch(Chem.MolFromSmiles('CC(F)Cl'), useChirality=True))
# True
# False
# True

# 三、smarts
# smarts在子结构匹配、化学反应等方面发挥着重要的作用

1. C c 大写小写C是不一样的，大写代表脂肪碳；小写代表芳香碳
2. 冒号后面的数字为atom map id




