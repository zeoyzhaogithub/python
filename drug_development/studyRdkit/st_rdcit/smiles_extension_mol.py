#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-23 13:08:08
# @Author: zeoy
# rdkit smiles支持和扩展

# 一、引入所需库
import os

from rdkit import Chem

from rdkit.Chem import Draw
# rdkit涵盖了Daylight SMILES所有的标准功能以及一些有用的扩展，下面是扩展的部分内容

# 二、芳香性
# 和氧同族的Te(碲 ， 拼音 ： dì ， 原子序数52 ， 是银白色的类金属 ） 元素也可能具有芳香性 ， 当其连接2个芳香原子时 ， 它贡献2个pi电子 。
m = Chem.MolFromSmiles('OC(=O)c1[te]ccc1')
Draw.MolToImageFile(
    m,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol48.jpg',
    legend='tellurophene-2-carboxylic acid'
)
# 碲吩-2甲酸分子
# Te原子的编号是4，下面检查其芳香性
aromatic_atom4 = m.GetAtomWithIdx(4).GetIsAromatic()
print('atom4 is aromatic', aromatic_atom4)  # atom4 is aromatic True
# 三 配位键
rdkit通过 -> 和 < -来支持配位键表示 ， 箭头的方向非常重要 ， 代表了谁提供电子
配位键不会影响起始原子的价态 ， 只会影响指向原子的价态
cu_mol = Chem.MolFromSmiles('[Cu](Cl)Cl')
bipy = Chem.MolFromSmiles('C1(C2=NC=CC=C2)=CC=CC=N1')
bipycu = Chem.MolFromSmiles('c1cccn->2c1-c1n->3cccc1.[Cu]23(Cl)Cl')
mols = [cu_mol, bipy, bipycu]
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(400, 400),
    legends=['CuCl2', 'bipy', 'bipycu']
)
img.save(
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol49.jpg'
)

# 获取Cu原子
cu_atom = cu_mol.GetAtoms()[0]
# Cu原子在分子中的价位
cu_valence = cu_atom.GetExplicitValence()
print('Cu atom valence in the first molecule', cu_valence)
# Cu atom valence in the first molecule 2
n_atom = bipy.GetAtoms()[2]
n_valence = n_atom.GetExplicitValence()
print('N atom valence in the second molecule', n_valence)
# N atom valence in the second molecule 3
n_atom = bipycu.GetAtoms()[4]
n_valence = n_atom.GetExplicitValence()
print('N atom valence in the third molecule', n_valence)
# N atom valence in the third molecule 3
cu_atom = bipycu.GetAtoms()[12]
cu_valence = cu_atom.GetExplicitValence()
print('Cu atom valence in the third molecule', cu_valence)
# Cu atom valence in the third molecule 4

有上述的数据我们看到，在配位后Cu的价态由二价变为四价，N的价态不变 还是三价。

# 四、闭环
rdkit除了支持[atom]n的方式表示闭环外，还支持 %（n）的表示方法。N的方位是0-99999；n可以不从0开始，可以是任意数字。
mol1 = Chem.MolFromSmiles('C%(1000)OC%(1000)')
mol2 = Chem.MolFromSmiles('C2OC2')
mol3 = Chem.MolFromSmiles('C1OC1')
Draw.MolToImageFile(
    mol1,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol50.jpg',
    legend='Propylene Oxide'
)
# 环氧丙烷

# 五、通过atomic number 指定原子
# rdkit 除了直接指定原子symbol也支持通过atomic number来指定原子[atomic number]。 atomic number 默认是在SMARTS中使用的，Smiles 也支持这种形式。
mol1 = Chem.MolFromSmiles('C1OC1')
mol2 = Chem.MolFromSmiles('[#6]1[#8][#6]1')

# 六、ChemAxon SMILES 拓展 CXSMILES extensions
RDKit 支持部分ChemAxon 拓展的SMILES 语法功能。

1.atomic coordinates 原子坐标
2.atomic values  原子值
3.atomic labels  原子标签
4.atomic properties  原子性质
5.coordinate bonds(these are translated into double bonds)
6.radicals
7.enhanced stereo(these are converted into StereoGroups)
# rdkit.Chem.rdmolfiles.MolToCXSmiles>`_ 导出ChemAxon 格式的SMILES.
前6个性质可通过方法` rdkit.Chem.rdmolfiles.MolToCXSmiles() < https: // www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html

m = Chem.MolFromSmiles('OC')
m.GetAtomWithIdx(0).SetProp('p1', '2')
m.GetAtomWithIdx(1).SetProp('p1', '5')
m.GetAtomWithIdx(1).SetProp('p2', 'A1')
m.GetAtomWithIdx(0).SetProp('atomLabel', 'O1')
m.GetAtomWithIdx(1).SetProp('atomLabel', 'C2')
sm = Chem.MolToCXSmiles(m)
print(sm)  # CO |$C2;O1$,atomProp:0.p1.5:0.p2.A1:1.p1.2|

# 七、使用坐标和异构SMILES从molfile中识别手性中心
mol = Chem.MolFromMolBlock("""
     RDKit          2D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  2  4  1  0
  4  5  2  0
  4  6  1  0
M  END""")
Draw.MolToImageFile(
    mol,
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol51.jpg',
    legend='L-alanin'
)

# L-丙氨酸

# 根据结构指定原子手性标记
Chem.AssignAtomChiralTagsFromStructure(mol)
# 找到分子的手性中心
chiral_center = Chem.FindMolChiralCenters(mol)
print(chiral_center)  # [(1, 'S')]

# 在smiles中也有该性质
sm = Chem.MolToSmiles(mol)
print(sm)  # C[C@H](N)C(=O)O
m2 = Chem.MolFromSmiles(sm)
Chem.AssignAtomChiralTagsFromStructure(m2)
chiral_center = Chem.FindMolChiralCenters(m2)
print(chiral_center)  # [(1, 'S')]

# 当以非异构体的形式输出读取时，因为分子不在有构象所以手性信息会丢失
m3 = Chem.MolToSmiles(mol, isomericSmiles=False)
print(m3)  # CC(N)C(=O)O
m4 = Chem.MolFromSmiles(m3)
Chem.AssignAtomChiralTagsFromStructure(m4)
chiral_center = Chem.FindMolChiralCenters(m4)
print(chiral_center)  # []
