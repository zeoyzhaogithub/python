#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-01 16:54:16
# @Author: zeoy
# rdkit 化学转换

# 一、引入所需库
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS

# > RDKit包含许多用于修饰分子的功能。注意，这些变换功能旨在提供一种简单的方法，可以对分子进行简单的修饰。

# 二、基于子结构的转换
# 使用RDKit的子结构匹配机器进行快速分子转化的功能多种多样

# 2.1 删除子结构
m = Chem.MolFromSmiles('CC(=O)O')
patt = Chem.MolFromSmiles('C(=O)[OH]')
rm = AllChem.DeleteSubstructs(m, patt)
smi = Chem.MolToSmiles(rm)
print(smi)  # C
mols = [m, rm]
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=2,
    subImgSize=(200, 200),
    legends=['' for x in mols]
)
img.save('/drug_development/studyRdkit/st_rdcit/img/mol12.jpg')

# 2.2 取代基替换
m = Chem.MolFromSmiles('COc1c(Br)cccc1OC')
patt = Chem.MolFromSmarts('OC')
repsmis = ['F', 'Cl', 'Br', 'O']

mols = []
mols.append(m)
for r in repsmis:
    rep = Chem.MolFromSmarts(r)
    res = AllChem.ReplaceSubstructs(m, patt, rep)
    mols.extend(res)

print(mols)
smis = [Chem.MolToSmiles(mol) for mol in mols]
print(smis)  # ['COC(C)=O']
mols = [Chem.MolFromSmiles(smi) for smi in smis]
print(mols)
img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(
    200, 200), legends=['' for x in mols])
img.save('/drug_development/studyRdkit/st_rdcit/img/mol13.jpg')

# > 注：ReplaceSubstructs()替换操作返回的是分子对象操作列表，如果分子只有一个地方匹配到，则返回一个分子的列表。
# 如果分子中有2个地方匹配到，则返回2个分子的列表。为了标准化smiles，可以将得到的分子mol化->smiles->mol，然后在进行可视化

# 2.3 SAR分析-core可视化
# Chem.ReplaceSidechains(m1,core) : 我们需要定义分子对象，骨架分子； 然后执行ReplaceSidechains函数，删除侧链就能得到骨架可视化。
# 定义嘧啶为核心结构，对其骨架进行可视化
m1 = Chem.MolFromSmiles('BrCCc1cncnc1C(=O)O')
core = Chem.MolFromSmiles('c1cncnc1')
tmp = Chem.ReplaceSidechains(m1, core)
Chem.MolToSmiles(tmp)

Draw.MolToImageFile(
    tmp, '/drug_development/studyRdkit/st_rdcit/img/mol14.jpg')

# 2.4 SAR分析-sidechain可视化
m1 = Chem.MolFromSmiles('BrCCc1cncnc1C(=O)O')
core = Chem.MolFromSmiles('c1cncnc1')
tmp = Chem.ReplaceCore(m1, core)
Draw.MolToImageFile(
    tmp, '/drug_development/studyRdkit/st_rdcit/img/mol15.jpg')

# >注：侧链的编号默认是从1开始的，这取决于算法找到侧链的先后顺序。
# 也可以根据侧链连接到骨架上的原子进行编号tmp=CHem.ReqlaceCore(m1, core)
tmp = Chem.ReplaceCore(m1, core, labelByIndex=True)
Draw.MolToImageFile(
    tmp, '/drug_development/studyRdkit/st_rdcit/img/mol16.jpg')

# 2.5 拆分手段
rs = Chem.GetMolFrags(tmp, asMols=True)
print(len(rs))  # 2
smi0 = Chem.MolToSmiles(rs[0])
print(smi0)  # *CCBr
smi1 = Chem.MolToSmiles(rs[1])
print(smi1)  # [5*]C(=O)O

# 三、Murcho分解
# 把分子中环结构提取出来，然后保留连接环结构的最少的键，如果该结构上的原子是双键连接，则保留双键，得到的结构称为Murcho骨架
# 3.1 获取分子骨架Murcho
m1 = Chem.MolFromSmiles(
    'C=CC(=O)N1CCC(CC1)C2CCNC3=C(C(=NN23)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=O)N')
m2 = Chem.MolFromSmiles(
    'CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4')
m3 = Chem.MolFromSmiles('CNC1(CCCCC1=O)C1=CC=CC=C1Cl')

core_1 = MurckoScaffold.GetScaffoldForMol(m1)
core_2 = MurckoScaffold.GetScaffoldForMol(m2)
core_3 = MurckoScaffold.GetScaffoldForMol(m3)

core_mols = [core_1, core_2, core_3]
img = Draw.MolsToGridImage(
    core_mols,
    molsPerRow=3,
    subImgSize=(300, 300),
    legends=['' for x in core_mols]
)
img.save('/drug_development/studyRdkit/st_rdcit/img/mol17.jpg')

# 四、最大公共分子
# RDKit内置了计算最大公共子结构的函数FindMCS
# 4.1 最大公共子结构FindMCS函数解析
FindMCS(
    (AtomPairsParameters)mols   # mol对象
    [, (bool)maximizeBonds=True
    [, (float)threshold=1.0
    [, (int)timeout=3600
    [, (bool)verbose=False
    [, (bool)matchValences=False
    [, (bool)ringMatchesRingOnly=False
    [, (bool)completeRingsOnly=False
    [, (bool)matchChiralTag=False
    [, (AtomCompare) atomCompare=rdkit.Chem.rdFMCS.AtomCompare.CompareElements
    [, (BondCompare) bondCompare=rdkit.Chem.rd.FMCS.BondCompare.CompareOrder
    [, (RingCompare) ringCompare=rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion
    [, (str) seedSmarts=''
    ]]]]]]]]]]]
)


# 该函数返回一个MCSResult实例，其中包含有关MCS中的原子数和键数，与所标识的MCS匹配的SMARTS字符串以及有关的算法是否超时的标志信息。
# 如果未找到MCS，则将原子和键数设置为0，SMARTS设置为''。
# 默认atomCompare和bondCompare采用(AtomCompare)AtomCompare=rdkit.Chem.rdFMCS.AtomCompare.CompareElements 要求他们元素相同和 (BondCompare)bondCompare=rdkit.Chem.rdFMCS.BondCompare.CompareOrder 有相同的键的类型。

# atomCompare也有其他的内置函数如：
rdkit.Chem.rdFMCS.AtomCompare.CompareElements
CompareAny = rdkit.Chem.rdFMCS.AtomCompare.CompareAny
CompareAnyHeavyAtom = rdkit.Chem.rdFMCS.AtomCompare.CompareAnyHeavyAtom
CompareElements = rdkit.Chem.rdFMCS.AtomCompare.CompareElements
CompareIsotopes = rdkit.Chem.rdFMCS.AtomCompare.CompareIsotopes

bondCompare也有其他的内置函数如：
CompareAny = rdkit.Chem.rdFMCS.BondCompare.CompareAny
CompareOrder = rdkit.Chem.rdFMCS.BondCompare.CompareOrder
CompareOrderExact = rdkit.Chem.rdFMCS.BondCompare.CompareOrderExac

MCS算法搜索公共子结构的通常花费几秒钟的时间，如果遇到复杂结构，通常需要几分钟甚至更长时间。 默认timeout= 3600秒。如果超过默认时间，则这个res.canceled属性会被设置成True。


# 4.2 最大公共分子结构FindMCS
# 以下三个分子为例，计算它们的最大公共分子

mol1 = Chem.MolFromSmiles("O=C(NCc1cc(OC)c(O)cc1)CCCC/C=C/C(C)C")
mol2 = Chem.MolFromSmiles("CC(C)CCCCCC(=O)NCC1=CC(=C(C=C1)O)OC")
mol3 = Chem.MolFromSmiles("c1(C=O)cc(OC)c(O)cc1")
mols = [mol1, mol2, mol3]

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(300, 300),
    legends=['' for x in mols]
)
img.save('/drug_development/studyRdkit/st_rdcit/img/mol18.jpg')

res = rdFMCS.FindMCS(mols)
print('原子个数=', res.numAtoms)  # 原子个数= 10
print('键个数=', res.numBonds)  # 键个数= 10

common_mol = Chem.MolFromSmarts(res.smartsString)

Draw.MolToImageFile(
    common_mol, '/drug_development/studyRdkit/st_rdcit/img/mol19.jpg')


# 默认情况下，如果两个原子是相同的元素，则匹配；如果两个原子具有相同的键类型，则键匹配。指定atomCompare和bondCompare使用不同的比较功能。

mols = (
    Chem.MolFromSmiles('NCC'),
    Chem.MolFromSmiles('OC=C')
)
res1 = rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareAny).smartsString
print(res1)  # [#7,#8]-[#6]
res2 = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareAny).smartsString
print(res2)   # [#6]-,=[#6]

# atomCompare参数的选项为：CompareAny表示任何原子都与任何其他原子匹配，
# CompareElements按元素类型进行比较，Comparelsotopes根据同位素标签进行匹配。同位素标记可用于实现用户定义的原子类型。
# CompareAny的bondCompare表示任何键都与任何其他键匹配，CompareOrderExact表示，当且仅当键具有相同的键类型时，它们才是等价的，并且CompareOrder允许单个键和芳族键相互匹配，但要求精确的顺序匹配：

mols = (
    Chem.MolFromSmiles('c1ccccc1'),
    Chem.MolFromSmiles('C1CCCC=C1')
)
res_bond_any = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareAny).smartsString
print(res_bond_any)
res_bond_order_exact = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString
print(res_bond_order_exact)
res_bond_order = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareOrder).smartsString
print(res_bond_order)
