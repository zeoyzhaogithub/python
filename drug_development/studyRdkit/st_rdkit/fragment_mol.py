#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-14 08:21:48
# @Author: zeoy
# rdkit 分子片段、R集团分解、片段指纹与指纹重要性分析

# 一、引入所需库
import os

from rdkit import Chem
from rdkit import RDConfig

from rdkit.Chem import Draw

from rdkit.Chem import FragmentCatalog

from rdkit.Chem import rdRGroupDecomposition as rdRGD
# from rdkit.Chem.Pharm2D.SigFactory import SigFactory
# from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.ML import InfoTheory


# 二、分子片段(molecular fragments)
# 分子片段(molecular fragments)是一组可能具有相关功能的相连的原子组成的。RDKit中包含了大量把分子分解成片段的方法和处理片段的工具。

# # 2.1 分子片段获取
# 获取官能团库
fName = os.path.join(
    RDConfig.RDDataDir,
    '/drug_development/studyRdkit/st_rdcit/data/FunctionalGroups.txt'
)

# 根据官能团库实例化一个参数器
fparams = FragmentCatalog.FragCatParams(1, 6, fName)
# 查看官能团库中包含的官能团数量
fparams_num = fparams.GetNumFuncGroups()

mols = []
# 查看每个官能团对应的集团
for i in range(fparams_num):
    mols.append(fparams.GetFuncGroup(i))

# 可视化官能团库
img = Draw.MolsToGridImage(mols, molsPerRow=8)
img.save(
    '/drug_development/studyRdkit/st_rdcit/img/mol39.jpg'
)

# 根据这份官能团列表，我们就可以分析分子在有多少个官能团

# 传入参数器，创建一个片段存储器，产生的分子片段都会存储在该对象中
fcat = FragmentCatalog.FragCatalog(fparams)

# 创建一个片段生成器，通过该对象生成片段
fcgen = FragmentCatalog.FragCatGenerator()

m = Chem.MolFromSmiles('OCC=CC(=O)O')
# 计算分子片段
fcgen.AddFragsFromMol(m, fcat)
# 查看分子片段数量
num_entries = fcat.GetNumEntries()
# 通过存储器查看片段
print(fcat.GetEntryDescription(0))
print(fcat.GetEntryDescription(1))
# C<-O>C
# C=C<-C(=O)O>

# 注 ：
# 尖括号中的内容 ： 表示与片段相连的官能团 ， 以上面的结果为例
# 第0号片段中 ， 对应着一个乙基片段 ， 该乙基片段与一个羟基相连
# 第1号片段中 ， 对应着一个乙烯片段 ， 该乙烯片段与一个羧基相连


# 关于官能团的详细信息 ， 可以通过下述方法获取 ：

# 向存储器传入分子片段id ， 获取片段中所包含的官能团标号 ： GetEntryFuncGroupIds
entries_group_ids = fcat.GetEntryFuncGroupIds(num_entries-1)

print('matched the function group ids is', list(entries_group_ids))
# matched the function group ids is [34, 1]

# 向参数器传入官能团编号，获取官能团对应的mol对象
fg1 = fparams.GetFuncGroup(1)
fg34 = fparams.GetFuncGroup(34)

print(fg1)  # <rdkit.Chem.rdchem.Mol object at 0x112b1b8f0>
print('name of group 1', fg1.GetProp('_Name'))  # name of group 1 -C(=O)O
print('name of group 34', fg34.GetProp('_Name'))  # name of group 34 -O

# 可视化官能团
mols = [fg1, fg34]
img = Draw.MolsToGridImage(mols, molsPerRow=2)
img.save(
    '/drug_development/studyRdkit/st_rdcit/img/mol40.jpg'
)

# 用该方法提取到的片段是层级结构 ， 小片段在最底层 ， 逐渐合并形成大片段 。
# 可以查看一个小片段形成了那些大片段

# 根据id获取片段
print(fcat.GetEntryDescription(0))
print(fcat.GetEntryDescription(1))
# C<-O>C
# C=C<-C(=O)O>

# 获取上级片段id
fcat_down_ids = list(fcat.GetEntryDownIds(0))
print(fcat_down_ids)  # [2]

# 根据上级片段获取上级片段信息
fcat_down_desc = fcat.GetEntryDescription(fcat_down_ids[0])
print(fcat_down_desc)
# C<-C(=O)O>=CC<-O>

# # 2.2 片段指纹生成

# 先将多个分子的片段汇总到一个片段存储器中
ms = [
    Chem.MolFromSmiles('OCC(NC1CC1)CCC'),
    Chem.MolFromSmiles('OCC=CC(=O)O')
]
# 片段存储器
fcat = FragmentCatalog.FragCatalog(fparams)

# 片段生成器
for m in ms:
    fcgen.AddFragsFromMol(m, fcat)

# 查看分子片段数量
num_entries = fcat.GetNumEntries()
print(num_entries)  # 17

# 存储器收集完所有片段后 ， 再用它来生成分子指纹

# 创建一个片段指纹生成器：FragFPGenerator()
fpgen = FragmentCatalog.FragFPGenerator()
# 传入分子和存储器用于生成指纹:GetFPForMol(mol,fcat)
fp1 = fpgen.GetFPForMol(ms[1], fcat)
# 以字符串形式查看指纹:ToBitString()
print(fp1.ToBitString())  # 10000000000000011

# 查看指纹中哪些位是有效的:GetOnBits()
print(list(fp1.GetOnBits()))  # [0, 15, 16]

# 可以用处理一般分子指纹的方法来处理片段分子指纹，例如寻找相同的片段

# 先对分子指纹做“&”位运算，两个指纹结果都为1时，结果为1，否则为0
# 获取两个指纹中都出现的片段：GetOnBits()
# 查看片段信息：GetEnteyDescription()

fp0 = fpgen.GetFPForMol(ms[0], fcat)
andfp = fp0 & fp1
onbit = list(andfp.GetOnBits())
print(onbit)  # [0]
onbit_0_desc = fcat.GetEntryDescription(onbit[0])
print(onbit_0_desc)  # C<-O>C

# 也可以查看一下，哪些片段导致了分子的差异

# 对分子指纹做“^”运算，两个指纹相同时，结果为0，否则为1。在做一步“&”运算
# 按以上方法查看相异片段

dis = fp0 ^ fp1
print(list(dis.GetOnBits()))
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
combinedfp = fp0 & dis
onbit = list(combinedfp.GetOnBits())
print(onbit)   # [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
fcat_entry_desc = fcat.GetEntryDescription(onbit[-1])
print(fcat_entry_desc)  # CCCC(C<-O>)N<-cPropyl>

# # 2.3 指纹重要性分析
# 以下内容主要分析指纹对离散标签的重要性。在rdkit.ML.InfoTheory.rdInfoTheory.InfoBitRanker
# 中提供了对指纹分析的功能。这个类可以根据分子指纹和离散标签，对特征进行计算和排序。

# 生成指纹片段
fName = os.path.join(
    RDConfig.RDDataDir,
    '/drug_development/studyRdkit/st_rdcit/data/bzr.sdf'
)
supple = Chem.SDMolSupplier(fName)
sdms = [x for x in supple]
print('指纹片段数', len(sdms))  # 指纹片段数 163

# 获取官能团库
fName = os.path.join(
    RDConfig.RDDataDir,
    '/drug_development/studyRdkit/st_rdcit/data/FunctionalGroups.txt'
)

# 根据官能团库实例化一个参数器
fparams = FragmentCatalog.FragCatParams(1, 6, fName)
# 查看官能团库中包含的官能团数量
fparams_num = fparams.GetNumFuncGroups()

# 片段存储器
fcat = FragmentCatalog.FragCatalog(fparams)
# 片段生成器
fcgen = FragmentCatalog.FragCatGenerator()
# 片段指纹生成器
fpgen = FragmentCatalog.FragFPGenerator()
# 汇总所有片段
for m in sdms:
    fcgen.AddFragsFromMol(m, fcat)

# 生成指纹片段
fps = [fpgen.GetFPForMol(x, fcat) for x in sdms]

# 信息增益（infoGain）分析，实例化一个排序对象：
# InfoBitRanker(
# nBits,    nBits：指纹长度
# nClasses, nClasses：类别数量，需要和标签满足的关系：0 <= 标签 < 类别数量
# infoType  infoType：度量指标。默认使用rdInfoTheory.InfoType.ENTROPY，即信息增益作为比较标准，它反映了使用某个特征进行分类后，系统混乱程度降低的多少，数值越大表明特征越重要。
# )

ranker = InfoTheory.InfoBitRanker(len(fps[0]), 2)

# 获取每个分子的活性信息 ： GetDoubleProp('ACTIVITY')
# 以7为标准对活性离散化 ： 大于7为1 ， 小于7为0
# 根据指纹和类别进行投票 ： AccumulateVotes(fp, act)
# 获取前5个重要特征 ： GetTopN(5)
# 依次输出特征id、信息增益、特征为0类别中的无活性分子数、特征为1类别中的有活性分子数。

acts = [x.GetDoubleProp('ACTIVITY') for x in sdms]
print(acts)
[6.87, 7.7, 7.74, 6.45, 6.89, 8.74, 7.23, 8.74, 6.51, 6.68, 7.47, 8.09, 8.07, 8.51, 8.42, 7.04, 8.46, 8.92, 6.06, 8.32, 8.0, 8.03, 7.74, 7.03, 6.96, 6.77, 5.0, 7.43, 5.0, 7.89, 6.46, 7.4, 6.41, 5.0, 7.06, 5.0, 5.0, 6.49, 8.46, 5.0, 5.0, 6.34, 7.68, 8.82, 7.85, 6.42, 7.12, 7.77, 8.13, 8.29, 5.0, 5.0, 7.31, 6.37, 8.08, 7.61, 8.8, 8.39, 7.72, 5.0, 5.0, 7.55, 7.15, 8.41, 8.03, 6.85, 8.85, 8.46, 8.48, 8.4, 8.15, 8.54, 7.43, 6.34, 7.14, 6.82, 6.52, 7.38, 8.19, 7.61, 7.15, 5.0,
    8.25, 8.28, 7.82, 6.64, 7.72, 8.19, 7.21, 5.0, 6.39, 8.52, 8.6, 8.17, 5.0, 8.38, 8.82, 5.0, 5.0, 8.57, 7.8, 8.77, 7.82, 6.14, 7.6, 8.3, 8.49, 7.85, 8.6, 7.34, 8.49, 8.85, 5.52, 8.11, 8.47, 8.48, 8.51, 8.0, 8.4, 8.55, 8.2, 8.54, 8.21, 8.77, 8.64, 8.28, 6.59, 8.46, 8.62, 6.21, 7.19, 7.44, 7.52, 7.74, 7.37, 7.62, 8.28, 7.02, 8.44, 7.85, 7.72, 8.13, 8.46, 7.59, 7.96, 7.89, 7.96, 7.12, 8.4, 7.92, 8.55, 8.72, 8.28, 6.38, 8.55, 6.52, 7.4, 7.8, 7.47, 8.4, 8.37, 8.35, 8.38]
for i, fp in enumerate(fps):
    act = int(acts[i] > 7)
    ranker.AccumulateVotes(fp, act)

top5 = ranker.GetTopN(5)
for id, gain, n0, n1 in top5:
    print(int(id), '%.3f' % gain, int(n0), int(n1))
# 698 0.081 20 17
# 222 0.073 23 25
# 196 0.073 30 43
# 378 0.073 30 43
# 1207 0.073 0 25

# 加入偏置 ， 以信息增益为例 ， 重新设置infoType
# 设置偏置类别 ： SetBiasList()
# 在这种模式下 ， 一个特征与所设置了偏置类别的相关性要高于所有非偏置类别 ，
# 例如设置偏置类别为4，某位特征为为1对应的标签中，类别为4的数量应该大于其他类别的数量
ranker = InfoTheory.InfoBitRanker(
    len(fps[0]), 2, InfoTheory.InfoType.BIASENTROPY)
ranker.SetBiasList((0,))
acts = [x.GetDoubleProp('ACTIVITY') for x in sdms]

for i, fp in enumerate(fps):
    act = 0 if acts[i] < 7 else 1
    ranker.AccumulateVotes(fp, act)

top5 = ranker.GetTopN(5)
for id, gain, n0, n1 in top5:
    print(int(id), '%.3f' % gain, int(n0), int(n1))
# 698 0.081 20 17
# 222 0.073 23 25
# 196 0.073 30 43
# 378 0.073 30 43
# 2375 0.062 5 0
# 三、R集团分解
# rdkit中由自动提取R集团的算法rdRGroupDecomposition 。

smis = [
    'OC(C1=CC(O)=CC=C1)=O',
    'OC(C2=CC(F)=CC=C2)=O',
    'OC(C3=C(Br)C=CC=C3)=O'
]

mols = []
for smi in smis:
    mols.append(Chem.MolFromSmiles(smi))

print('number of mols', len(mols))  # number of mols 3

core = Chem.MolFromSmarts('[*:1]c1cc([*:2])ccc1')

res, unmatched = rdRGD.RGroupDecompose([core], mols, asSmiles=True)
print(res)
[
    {
        'Core': 'c1cc([*:1])c([*:3])c([*:2])c1',
        'R1': 'O=C(O)[*:1]',
        'R2': 'O[*:2]',
        'R3': '[H][*:3]'
    }, {
        'Core': 'c1cc([*:1])c([*:3])c([*:2])c1',
        'R1': 'O=C(O)[*:1]',
        'R2': 'F[*:2]',
        'R3': '[H][*:3]'
    }, {
        'Core': 'c1cc([*:1])c([*:3])c([*:2])c1',
        'R1': 'O=C(O)[*:1]',
        'R2': '[H][*:2]',
        'R3': 'Br[*:3]'
    }
]
