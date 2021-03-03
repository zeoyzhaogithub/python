#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-24 21:18:04
# @Author:
# rdkit 化学反应处理Chemical Reaction Handling

# 一、引入所需库
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# 二、反应SMARTS（reaction smarts）
# rdkit中的反应SMARTS是基于SMARTS表达，但是又不同于SMILES和reaction SMILES，称为 reaction SMARTS。
# reaction SMARTS语法如下：
# reaction   ::= reactants">>"
# reactants  ::= molecules
# products   ::= molecules
# molecules  ::= molecules
#                molecules "." molecule
# molecule   ::= a valid SMARTS string without "." characters
#               "(" a valid SMARTS string without "." characters ")"

# 三、特性
# # 3.1 特性1 当模板中带有map id的虚原子时，生成的产物中的虚原子会被反应中对应的原子替代
# 创建具体反应规则的引擎对象
rxn = AllChem.ReactionFromSmarts('[C:1]=[O,N:2]>>[C:1][*:2]')
img = Draw.ReactionToImage(
    rxn
)
img.save('/drug_development/studyRdkit/st_rdcit/img/mol57.jpg')

# 产生化学反应
reactant = rxn.RunReactants((Chem.MolFromSmiles('CC=O'),))

products_r1 = [Chem.MolToSmiles(x, 1) for x in rxn.RunReactants(
    (Chem.MolFromSmiles('CC=O'),))[0]]
print(products_r1)  # ['CCO']

products_r2 = [Chem.MolToSmiles(x, 1) for x in rxn.RunReactants(
    (Chem.MolFromSmiles('CC=N'),))[0]]
print(products_r2)  # ['CCN']

# acid = Chem.MolFromSmiles('CC(=O)O')
# base = Chem.MolFromSmiles('CC(=O)NCCN')
# mols = [acid, base]
# for index, m in enumerate(mols):
#     img_file = str(index) + '.jpg'
#     Draw.MolToImageFile(
#         m,
#         img_file
#     )

# 四、手性
# 五、其他反应规则和警告
