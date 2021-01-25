#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-24 16:40:05
# @Author: zeoy
# rdkit SMARTS支持和扩展

# 一、引入所需库
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# rdkit 支持Daylight定义的SMARTS的绝大部分标准特性以及一些有用的拓展

# 二、SMARTS 不支持的特性
# 1.非四面体手性轴手性
# 2. @ ？ 操作符号
# 3.显式原子质量 （ 支持同位素定义查询 ）
# 4.片段组的查询(C).(C)

# 三、SMARTS 支持的扩展

# # 3.1 杂化方式查询
# 杂化方式在SMARTS 中通过^符号进行定义。 如：

# 1.^0 匹配S 杂化的原子
# 2.^1 匹配SP 杂化的原子
# 3.^2 匹配SP2 杂化的原子
# 4.^3 匹配SP3 杂化的原子
# 5.^4 匹配SP3D 杂化的原子
# 6.^5 匹配SP3D2 杂化的原子

aspirin = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
Draw.MolToImageFile(
    aspirin,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol52.jpg',
    legend='aspirin'
)
# 阿司匹林
# sp2杂化的原子
sp2_atoms = aspirin.GetSubstructMatches(Chem.MolFromSmarts('[^2]'))
# sp3杂化的原子
sp3_atoms = aspirin.GetSubstructMatches(Chem.MolFromSmarts('[^3]'))
print('sp2 atoms', sp2_atoms)
# sp2 atoms ((1,), (2,), (3,), (4,), (5,), (6,), (7,), (8,), (9,), (10,), (11,), (12,))
print('sp3 atoms', sp3_atoms)  # sp3 atoms ((0,),)

# 对于分子阿司匹林，只有0号原子是sp3杂化，其他原子都是sp2杂化
# >注：苯酚中的氧都是sp2杂化，所以羟基氧才具有更强的酸性。COO中的两个氧都是sp2杂化 羧基中也有类似苯环的共轭体系，并且羧酸中羟基氢酸性更强，共轭更明显，更应该是sp2。 醇羟基中的氧是sp3杂化。
# # 3.2 配位键
# rdkit的SMARTS通过 -> 和 < -符号表示配位键 ， 箭头的方向代表电子转译的方向
m = Chem.MolFromSmiles('C1=CC=CC=N1->[Fe]')
Draw.MolToImageFile(
    m,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol53.jpg',
    legend=''
)
match7_ = m.GetSubstructMatches(Chem.MolFromSmarts('[#7]->*'))
match_7 = m.GetSubstructMatches(Chem.MolFromSmarts('*<-[#7]'))
print(match7_)  # ((5, 6),)
print(match_7)  # ((6, 5),)
# # 3.3 邻居杂原子查询
# 根据碳原子相连杂原子的个数进行查询
# 1.zn代表匹配相连n个杂原子 （ 非CH原子 ） 的碳原子 ， 比如z2代表相连2个杂原子的任意碳原子
# 2.Zn代表匹配相连n个脂肪族杂原子（ 非CH原子 ） 的碳原子，比如Z2代表相连2个脂肪杂原子的碳原子

# >注：在有机化学中，一般将有机物分为三类：
# 1.开链化合物 ， 分子中的碳原子连接成链状 ， 又称为 “ 脂肪族化合物 ”
# 2.碳环化合物 ， 分子汇总的碳原子连接成环状 ， 包括脂肪族和芳香族化合物
# 3.杂环化合物，即分子中含有其他原子（如O、N、S、P等）的环状化合物。
mol = Chem.MolFromSmiles('O=C(O)c1nc(O)ccn1')
Draw.MolToImageFile(
    m,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol54.jpg',
    legend=''
)
z2_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts('[z2]'))
Z2_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts('[Z2]'))
Z1_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts('[Z1]'))
print(z2_atoms)  # ((1,), (3,), (5,))
print(Z2_atoms)  # ((1,),)
print(Z1_atoms)  # ((5,),)

# 该示例分子中和2个杂原子相连的碳原子编号是1，3，5； 和2个脂肪杂原子相连的碳原子编号是1； 和1个脂肪族原子相连的碳原子编号是5.

# # 3.4 范围查询
# rdkit SMARTS语法支持指定范围查询 ， 支持的类型有 ： D ， h, r, R, v, x, X, z, Z, +, -
# D代表度
# 1.D{2-4}匹配原子的度为2 - 4的原子 ， 链接2 - 4个原子的原子
# 2.D{-3}匹配原子的度小于等于3的原子
# 3.D{2-}匹配原子的度大于等于2的原子
mol = Chem.MolFromSmiles('CC(=O)OC')
Draw.MolToImageFile(
    mol,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol55.jpg',
    legend='methyl acetate'
)
# 乙酸甲酯
z1234 = mol.GetSubstructMatches(Chem.MolFromSmarts('[z{1-}]'))
D23 = mol.GetSubstructMatches(Chem.MolFromSmarts('[D{2-3}]'))
D012 = mol.GetSubstructMatches(Chem.MolFromSmarts('[D{-2}]'))
print('z1234', z1234)  # z1234 ((1,), (4,))
print('D23', D23)  # D23 ((1,), (3,))
print('D012', D012)  # D012 ((0,), (2,), (3,), (4,))

# # 3.5 SMARTS语法参考
Bond
Primitive	Property	Notes
“”	“single or aromatic”	“unspecified bonds”
-	single
=	double
#	triple
:	aromatic
~	“any bond”
@	“ring bond”
/	“directional”
\	“directional”
-> “dative right” extension
<- “dative left” extension

# 四、子结构匹配
# # 4.1 具有SMARTS查询的功能组
sucrose = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O"
sucrose_mol = Chem.MolFromSmiles(sucrose)
Draw.MolToImageFile(
    sucrose_mol,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol56.jpg',
    legend=''
)
primary_alcohol = Chem.MolFromSmarts("[CH2][OH1]")
match_alcohol = sucrose_mol.GetSubstructMatches(primary_alcohol)
print(match_alcohol)  # ((0, 22), (13, 14), (17, 18))
secondary_alcohol = Chem.MolFromSmarts('[CH1][OH1]')
match_secondary = sucrose_mol.GetSubstructMatches(secondary_alcohol)
print(match_secondary)  # ((2, 21), (3, 20), (4, 19), (9, 16), (10, 15))

# # 4.2 具有SMARTS查询的宏周期
erythromycin = Chem.MolFromSmiles(
    "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O")

Draw.MolToImageFile(
    erythromycin,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol57.jpg',
    legend=''
)
# 定义一个环原子大于12的SMARTS模型
macro = Chem.MolFromSmarts('[r{12-}]')
match_macro = erythromycin.GetSubstructMatches(macro)
print(match_macro)
