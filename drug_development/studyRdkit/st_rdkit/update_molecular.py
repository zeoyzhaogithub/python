#! /usr/bin/python
# coding: utf-8
# @Time: 2020-05-29 14:36:04
# @Author: zeoy
# rdkit 修改分子

# 一、引入所需库
from rdkit import Chem
from rdkit.Chem import Draw


# 二、增删H原子
mol = Chem.MolFromSmiles('OC1C2C1CC2')
# 画分子结构
Draw.MolToImageFile(
    mol,
    '/drug_development/studyRdkit/st_rdcit/img/mol5.jpg'
)

# 2.1 增加H原子函数解析

# 将氢添加到分子图上
rdkit.Chem.rdmolops.AddHs(
    (Mol)mol   # 要修饰的分子
    [, (bool) explicitOnly=False  # （可选）如果设置了此切换，则仅将显式Hs添加到分子中。默认值为0（添加隐式和显式Hs）。
     [, (bool) addCoords=False  # (可选) 如果设置了此开关，则Hs将设置3D坐标。默认值为0（无3D坐标）。
      [, (AtomPairsParameters) onlyOnAtoms=None  # （可选）如果提供了此序列，则仅将这些原子视为具有添加的Hs
       [, (bool)addResidueInfo=False  # （可选）如果为true，则将残基信息添加到氢原子（对PDB文件有用）。
        ]]]]
)

# 2.2 增加H原子
mol2 = Chem.AddHs(mol)
# 画分子结构
Draw.MolToImageFile(
    mol2,
    '/drug_development/studyRdkit/st_rdcit/img/mol6.jpg'
)
print('mol smiles:', Chem.MolToSmiles(mol))  # mol smiles: OC1C2CCC12
# mol2 smiles: [H]OC1([H])C2([H])C([H])([H])C([H])([H])C12[H]
print('mol2 smiles:', Chem.MolToSmiles(mol2))
print('num ATOMs in mol:', mol.GetNumAtoms())  # num ATOMs in mol: 6
print('num ATOMs in mol:', mol2.GetNumAtoms())  # num ATOMs in mol: 14
# 注：3D构象优化的时候，需要采用显式H原子


# 2.3 删除H原子函数解析

从分子图中除去所有氢。
rdkit.Chem.rdmolops.RemoveHs(
    (Mol)mol  # 要修饰的分子
    [，（bool）implicitOnly=False  # hiddenOnly ：（可选）如果设置了此切换，则只会从图中删除隐式Hs。默认值为0（删除隐式和显式Hs）。
     [，（bool）updateExplicitCount=False  # （可选）如果设置了此切换，则将更新具有Hs的原子的显式H计数。默认值为0（不更新显式H计数）。
      [，（bool）sanitize=True  # （可选）如果设置了此切换开关，则去除Hs后将对分子进行消毒。缺省值为1（进行消毒）。
       ]]])

笔记：
未与重原子连接的氢将不会被去除。这样可防止分子[H][H]除去所有原子。
标记的氢（例如原子序数 = 1，但同位素 > 1的原子）将不会被去除。
两个坐标Hs，例如C[H-] C中的中心H，将不会被删除
连接到虚拟原子的Hs将不会被移除
属于双键立体化学定义的Hs将不会被移除
未连接到其他任何对象的HS将不会被删除


# 2.4 删除H原子
mol3 = Chem.RemoveHs(mol2)
print('mol3 smiles:', Chem.MolToSmiles(mol3))  # mol3 smiles: OC1C2CCC12
print('num ATOMs in mol3:', mol3.GetNumAtoms())  # num ATOMs in mol3: 6


# 三、芳香共轭键和库里单双键
# RDKit 默认把芳香体系的键的类型识别为芳香键。

# 以苯为例
mol4 = Chem.MolFromSmiles('c1ccccc1')
# 画分子结构
Draw.MolToImageFile(
    mol4,
    '/drug_development/studyRdkit/st_rdcit/img/mol7.jpg'
)

# 键类型
for bond in mol4.GetBonds():
    print(bond.GetBondType())
# AROMATIC
# AROMATIC
# AROMATIC
# AROMATIC
# AROMATIC
# AROMATIC

# 3.1 将芳香键的类型修改为单双建的类型
Chem.Kekulize(mol4)
# 画分子结构
Draw.MolToImageFile(
    mol4,
    '/drug_development/studyRdkit/st_rdcit/img/mol8.jpg'
)
for bond in mol4.GetBonds():
    print(bond.GetBondType())

# DOUBLE
# SINGLE
# DOUBLE
# SINGLE
# DOUBLE
# SINGLE

# bond 1 is aromatic True
print('bond 1 is aromatic', mol4.GetBondWithIdx(1).GetIsAromatic())
# atom 1 is aromatic True
print('atom 1 is aromatic', mol4.GetAtomWithIdx(1).GetIsAromatic())
