#! /usr/bin/python
# coding: utf-8
# @Time: 2020-05-27 15:44:21
# @Author: zeoy
# rdkit 读写分子

# rdkit支持从Smiles、mol、sdf文件中读入分子获取分子对象。
# Smiles、mol通常用于保存单个分子；而sdf格式是作为分子库形式设计的。
# 因此读入sdf得到的是分子迭代器，读入Smiles、mol文件得到分子对象。

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# 一、读分子操作
# 1.1、读入smiles
smi = 'CC(C)OC(=O)C(C)NP(=O)(OCC1C(C(C(O1)N2C=CC(=O)NC2=O)(C)F)O)OC3=CC=CC=C3'
mol = Chem.MolFromSmiles(smi)   # 将Smiles转换为mol对象
# 将Mol分子画出结构图，并存储在相应地址
Draw.MolToImageFile(
    mol,    # mol分子对象
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol2.jpg"  # 分子结构图存储地址
)
print('mol的类型=', type(mol))   # mol的类型=<class 'rdkit.Chem.rdchem.Mol'>


# 1.2、读入mol文件

# 将mol文件转换为mol对象
mol3 = Chem.MolFromMolFile(
    '/Users/zeoy/st/drug_development/st_rdcit/952883.mol')
# 将Mol分子画出结构图，并存储在相应地址
Draw.MolToImageFile(
    mol3,    # mol分子对象
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol3.jpg"  # 分子结构图存储地址
)
print('mol3的类型=', type(mol))  # mol3的类型=<class 'rdkit.Chem.rdchem.Mol'>

# 1.3 读入sdf文件
mols_suppl = Chem.SDMolSupplier(
    '/Users/zeoy/st/drug_development/st_rdcit/2d.sdf')
print('类型=', type(mols_suppl))
for _mol in mols_suppl:
    print('类型=', type(_mol))  # mol3的类型=<class 'rdkit.Chem.rdchem.Mol'>

# 二、写分子操作
# RDKit可以把分子对象保存成Smiles、molBlock、mol文件

# 2.1 写将分子对象存储为mol文件

smi4 = 'CC(C)OC(=O)C(C)NP(=O)(OCC1C(C(C(O1)N2C=CC(=O)NC2=O)(C)F)O)OC3=CC=CC=C3'
mol4 = Chem.MolFromSmiles(smi4)
smi4 = Chem.MolToSmiles(mol4)
print('mol分子对象=', smi4)
molblock = Chem.MolToMolBlock(mol4)
print(molblock)
# print(*objects, sep=’ ‘, end=’n’, file=sys.stdout, flush=False)
# python的print函数file参数支持定义输出位置。
print(molblock, file=open('/Users/zeoy/st/drug_development/st_rdcit/stock.mol', 'w+'))
