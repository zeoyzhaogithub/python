#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-08 23:23:04
# @Author: zeoy
# rdkit 分子性质（描述符）

# 一、引入所需库
from rdkit import Chem
from rdkit import DataStructs

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw

from rdkit.Chem.Draw import SimilarityMaps

# 二、性质描述符计算
# 分子性质也被称为描述符。 RDKit中内置了大量的分子描述符的计算方法， 这些方法主要位于`rdkit.Chem.Descriptors <https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html#module-rdkit.Chem.Descriptors>`_ 也有些常用的性质在AllChem模块下面。

# 计算分子的The topological polar surface area (TPSA) descriptor 、logP、电荷等性质

m = Chem.MolFromSmiles('c1ccccc1C(=O)O')

tpsa_m = Descriptors.TPSA(m)
logp_m = Descriptors.MolLogP(m)

AllChem.ComputeGasteigerCharges(m)
charge_atm0 = float(m.GetAtomWithIdx(0).GetProp('_GasteigerCharge'))

print('the TPSA of m is', tpsa_m)  # the TPSA of m is 37.3
print('the logP of m is', logp_m)  # the logP of m is 1.3848
print('the gasteigerCharge of the first atom', charge_atm0)
# the gasteigerCharge of the first atom - 0.04769375004654255

# 三、原子对性质的贡献可视化

# 相似性地图也可用于性质的可视化，只要性质可以分解到原子上就可以进行可视化
mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')

AllChem.ComputeGasteigerCharges(mol)
contribs = [float(mol.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
            for i in range(mol.GetNumAtoms())]
# print(contribs)
# [
#     0.07771844728655561, -0.492820796983233,
#     0.16152918919178583, -0.01602219709262404,
#     - 0.05770529265149885, -0.050816379353719424,
#     0.024978185790941746, -0.00237222532795198,
#     0.19244618980318173, 0.28654036365436475,
#     - 0.2657003945289995, -0.34935691936673163,
#     0.017840495848308394, -0.03446611044002525,
#     - 0.03873971707569541, -0.001658098696669069,
#     - 0.29982407384979287, 0.015855299823588617,
#     0.030662879547007613, -0.36715829261179367,
#     0.06502826940121241, -0.03636566023234588,
#     - 0.05790012527553211, -0.03404808267884825,
#     0.09076826966888336, -0.25292830777068237,
#     0.04589736871527414, 0.04598776260803381,
#     - 0.2508095590227565, 0.11193305580241227,
#     0.030662879547007613, 0.015855299823588617,
#     - 0.4470028418642525, 0.1764058687346651
# ]
d = Draw.MolDraw2DSVG(400, 400)
d.ClearDrawing()
fig = SimilarityMaps.GetSimilarityMapFromWeights(
    mol,
    contribs,
    colorMap='jet',
    contourLines=10,
    draw2d=d
)
d.FinishDrawing()
with open('/Users/zeoy/st/drug_development/st_rdcit/img/mol29.svg', 'w+') as outf:
    outf.write(d.GetDrawingText())

# RDKit中内置Crippen方法计算原子logP。
mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')
contribs = rdMolDescriptors._CalcCrippenContribs(mol)

d = Draw.MolDraw2DSVG(400, 400)
d.ClearDrawing()
fig = SimilarityMaps.GetSimilarityMapFromWeights(
    mol,
    [x for x, y in contribs],
    colorMap='jet',
    contourLines=10,
    draw2d=d
)
d.FinishDrawing()
with open('/Users/zeoy/st/drug_development/st_rdcit/img/mol30.svg', 'w+') as outf:
    outf.write(d.GetDrawingText())
