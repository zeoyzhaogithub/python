#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-20 17:14:41
# @Author: zeoy
# rdkit 分子片段、R集团分解、片段指纹与指纹重要性分析

# 一、引入所需库
import os

from rdkit import Chem


from rdkit.Chem import Draw
# from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions


# 二、芳香性
芳香性是一个简单又复杂的话题。无论是实验化学家还是理论化学家都无法给出一个明确的定义

rdkit除了通过模式匹配已知的芳香体系进行识别 ， 还定义了一些识别方向体系的规则 。
这些具体规则如下 ：

1.芳香性仅仅和环上的原子和键有关 。
2.芳香性键必须由2两个芳香性原子组成 。
3.由2个芳香性原子组成的键不一定就是芳香键。

m = Chem.MolFromSmiles('C1=CC2=C(C=C1)C1=CC=CC=C21')
opts = DrawingOptions()
opts.atomLabelFontSize = 400

opts.includeAtomNumbers = True
opts.dblBondLengthFrac = 0.8
opts.includeAtomNumbers = True
opts.dotsPerAngstrom = 1000
img = Draw.MolToImage(m, options=opts, size=(500, 500))
img.save(
    '/Users/zeoy/st/drug_development/st_rdcit/img/mol41.jpg'
)

# 基于RDKit方向模型判断其上的原子和键的芳香性
atom3 = m.GetAtomWithIdx(3).GetIsAromatic()
atom6 = m.GetAtomWithIdx(6).GetIsAromatic()
bond_36 = m.GetBondBetweenAtoms(3, 6).GetIsAromatic()

print('atom 3 aromatic is', atom3)
print('atom 6 aromatic is', atom6)
print('bond 36 aromatic is', bond_36)
# atom 3 aromatic is True
# atom 6 aromatic is True
# bond 36 aromatic is False

# 三、RDKit 芳香性模型

# # 3.1 芳香性原子类型
芳香性是针对环体系原子和键而言的 ， 如果该环上的pi电子满足4N + 2规则 ， 则组成该环的原子和键具有芳香性
PI电子的计算根据原子的类型和环境 ， 具体如下表所示 ：
匹配模式	pi electrons 数目	解释
c(a)a	1	芳香碳原子，链接两个芳香原子，则该碳原子贡献一个PI电子
n(a)a	1	芳香氮原子，链接两个芳香原子，则该氮原子贡献一个PI电子
An(a)a	2	芳香氮原子，链接2个芳香原子和1个任意原子的时候，贡献2个PI电子
o(a)a	2	芳香氧原子，链接2个芳香原子，贡献2个PI电子
s(a)a	2	同上
se(a)a	2	同上
te(a)a	2	同上
O = c(a)a	0	芳香碳原子链接2个芳香原子，同时参与形成羰基，贡献0个PI电子
N = c(a)a	0	同上
*(a)a	0, 1, or 2
>a: 芳香原子 A: 任何原子，包含H * : 虚拟原子

# # 3.2 简单芳香模型
这个非常简单：只有五元和六元简单环才被视为具有芳香性。使用上面列出的相同的电子贡献计数。
# # 3.3 MDL芳香模型

这个没有找到很好的文档（至少不是公开的），因此我们尝试重现chem文档（https: // docs.eyesopen.com/toolkits/python/oechemtk/aromaticity.html）中提供的内容。

. 稠环（即ie）可以是芳香族的

.五元环不是芳香族的（尽管它们可以是稠合芳香族系统的一部分）

.只有C和N可以是芳香族

.仅接受一个电子供体

.具有环外双键的原子不是芳族的

>注：出于计算方便的原因，芳香性感知仅在所有成员的大小最多为24个原子的稠环系统中进行。

# 四、分子中的环
# #4.1 计算分子中的环系统


def GetRingSystems(mol, includeSpiro=False):
    '''
    计算分子中的环系统
    '''
    if (mol == ''):
        return False
    # 直接获取环的信息
    ri = mol.GetRingInfo()
    systems = []
    # 各个环上的原子信息
    print(ri.AtomRings())
    for ring in ri.AtomRings():
        # 每一个环的信息
        ringAts = set(ring)
        print(ringAts)
        nSystems = []
        for system in systems:
            print(system)
            print(ringAts)
            # 同时存在于两个环上的原子个数
            nInCommon = len(ringAts.intersection(system))
            print(nInCommon)
            # 如果两个环上有公用的原子
            if nInCommon and (includeSpiro or nInCommon > 1):
                # 公用的原子，说明两个环是并在一起的
                # 将两个环的原子去重合并
                ringAts = ringAts.union(system)
                print(ringAts)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


mol = Chem.MolFromSmiles('CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3')
Draw.MolToImageFile(
    mol,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol42.jpg",
)
ringInfo = GetRingSystems(mol)
print(ringInfo)

# # 4.2 环外原子对芳香环（Aromatic）的影响
# >注：环键上连接的负电性原子会“窃取”环原子的价电子，且这些亚原子，提供了使环芳香性所必须的元素。
# 使用稠环来增加芳香度可能导致单个环不是芳香的情况，但稠环系统是芳香性的。其中一个例子就是azulene（甘菊蓝）
# 下面的例子，展示了两个稠环和环外双键的影响
m = Chem.MolFromSmiles('O=C1C=CC(=O)C2=C1OC=CO2')
Draw.MolToImageFile(
    m,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol43.jpg",
)
isAromatic6 = m.GetAtomWithIdx(6).GetIsAromatic()
print(isAromatic6)  # True
isAromatic7 = m.GetAtomWithIdx(7).GetIsAromatic()
print(isAromatic7)  # True
isAromatic67 = m.GetBondBetweenAtoms(6, 7).GetIsAromatic()
print(isAromatic67)  # False

# # 4.3 判别芳香环的特殊情况
# # #4.3.1 带有自由基的杂原子不被视为芳香环
m = Chem.MolFromSmiles('C1=C[N]C=C1')
Draw.MolToImageFile(
    m,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol44.jpg",
)

isAromatic0 = m.GetAtomWithIdx(0).GetIsAromatic()
print(isAromatic0)  # False
isAromatic2 = m.GetAtomWithIdx(2).GetIsAromatic()
print(isAromatic2)  # False
radical_electrons2 = m.GetAtomWithIdx(2).GetNumRadicalElectrons()
print(radical_electrons2)  # 1

# ##4.3.2 带自由基的带电荷碳不被视为芳香环
m = Chem.MolFromSmiles('C1=CC=CC=C[C+]1')
isAromatic0 = m.GetAtomWithIdx(0).GetIsAromatic()
print(isAromatic0)  # False
isAromatic6 = m.GetAtomWithIdx(6).GetIsAromatic()
print(isAromatic6)  # False
formal_charge6 = m.GetAtomWithIdx(6).GetFormalCharge()
print(formal_charge6)  # 1
radical_electrons6 = m.GetAtomWithIdx(6).GetNumRadicalElectrons()
print(radical_electrons6)  # 1

# ##4.3.3 带自由基的中性碳可视为芳香环
m = Chem.MolFromSmiles('C1=[C]NC=C1')
Draw.MolToImageFile(
    m,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol45.jpg",
)
isAromatic0 = m.GetAtomWithIdx(0).GetIsAromatic()
print(isAromatic0)  # True
isAromatic1 = m.GetAtomWithIdx(1).GetIsAromatic()
print(isAromatic1)  # True
radical_electrons1 = m.GetAtomWithIdx(1).GetNumRadicalElectrons()
print(radical_electrons1)  # 1

# ## 4.3.4 识别芳香环
m = Chem.MolFromSmiles('c1cccc2c1CCCC2')
Draw.MolToImageFile(
    m,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol46.jpg",
)
# 获取分子的环信息
ri = m.GetRingInfo()
# 组成每个环的原子
atom_rings = ri.AtomRings()
print(atom_rings)  # ((0, 5, 4, 3, 2, 1), (6, 7, 8, 9, 4, 5))
# 组成每个环的键
bond_ring = ri.BondRings()
print(bond_ring)  # ((9, 4, 3, 2, 1, 0), (6, 7, 8, 10, 4, 5))


def isRingAromatic(mol, bondRings):
    '''
    将每个环上的键遍历一遍，如果所有键都是芳香性的，则将环标记为芳香环
    '''
    is_ring_aromatic = []
    for bondRing in bondRings:
        for id in bondRing:
            if not mol.GetBondWithIdx(id).GetIsAromatic():
                is_ring_aromatic.append(False)
                break
        else:
            # for循环完整执行过之后才会被执行，如果for循环被break，则else块将不会被执行。
            is_ring_aromatic.append(True)
    return is_ring_aromatic


print(isRingAromatic(m, bond_ring))  # [True, False]

# ## 4.3.5 识别芳香原子
m = Chem.MolFromSmiles('c1ccccc1C=CCC')
Draw.MolToImageFile(
    m,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol47.jpg",
)

aromatic_carbon = Chem.MolFromSmarts('c')
ring_aromatic_carbon = m.GetSubstructMatches(aromatic_carbon)
print(ring_aromatic_carbon)  # ((0,), (1,), (2,), (3,), (4,), (5,))

# rdkit包含一个查询杂化类型的SMARTS扩展，下面我们查询脂肪碳的sp2杂化
olefinic_carbon = Chem.MolFromSmarts('[C^2]')
olefinic_carbon_sp2 = m.GetSubstructMatches(olefinic_carbon)
print(olefinic_carbon_sp2)  # ((6,), (7,))
