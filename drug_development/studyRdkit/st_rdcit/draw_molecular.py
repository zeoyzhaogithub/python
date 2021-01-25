#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-01 11:01:04
# @Author: zeoy
# rdkit 绘制分子【可视化分子】

# 一、引入所需库
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
# from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults


# rdkit 内置了Draw模块，用于绘图，把一些经常用到的方法直接放在Draw下面。

# 一、分子对象转化为图片
## 1.1 分子对象转图片文件函数解析
Draw.MolToFile(
    mol,  # mol对象
    'filename.png',  # 图片存储地址
    size=(300, 300), 
    kekulize=True, 
    wedgeBonds=True, 
    imageType=None, 
    fitImage=False, 
    options=None, 
    **kwargs
)

## 1.2 分子对象转图片函数解析
MolToImage(
    mol, 
    size=(300, 300), 
    kekulize=True, 
    wedgeBonds=True, 
    fitImage=False, 
    options=None, 
    canvas=None, 
    **kwargs
)

## 1.3 分子对象转图片
opts = DrawingOptions()
m = Chem.MolFromSmiles('OC1C2C1CC2')
opts.includeAtomNumbers=True
opts.bondLineWidth=2.8
draw = Draw.MolToImage(m, options=opts)
draw.save('/Users/zeoy/st/drug_development/st_rdcit/img/mol10.jpg')

## 1.4 多个分子按照grid显示
smis=[
    'COC1=C(C=CC(=C1)NS(=O)(=O)C)C2=CN=CN3C2=CC=C3',
    'C1=CC2=C(C(=C1)C3=CN=CN4C3=CC=C4)ON=C2C5=CC=C(C=C5)F',
    'COC(=O)C1=CC2=CC=CN2C=N1',
    'C1=C2C=C(N=CN2C(=C1)Cl)C(=O)O',
]
mols = []
for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    mols.append(mol)

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=4,
    subImgSize=(200,200),
    legends=['' for x in mols]
)

img.save('/Users/zeoy/st/drug_development/st_rdcit/img/mol11.jpg')

## 1.5 多个分子基于公共骨架按照grid显示
template = Chem.MolFromSmiles('c1nccc2n1ccc2')
AllChem.Compute2DCoords(template)
mols = []
for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    # 生成一个分子的描述，其中一部分 分子被约束为具有与参考相同的坐标。
    AllChem.GenerateDepictionMatching2DStructure(mol, template)
    mols.append(mol)

# 基于分子文件输出分子结构
img = Draw.MolsToGridImage(
    mols,   # mol对象
    molsPerRow=4,
    subImgSize=(200,200),
    legends=['' for x in mols]
)
img.save('/Users/zeoy/st/drug_development/st_rdcit/img/mol12.jpg')


