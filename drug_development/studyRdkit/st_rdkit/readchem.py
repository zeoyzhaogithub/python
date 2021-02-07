#! /usr/bin/python
# coding: utf-8
# @Time: 2020-05-20 14:53
# @Author: zeoy

from rdkit.Chem import AllChem as Chem
from rdkit import rdBase
from rdkit.Chem import Draw

print('rdkit版本号=', rdBase.rdkitVersion)

# 西地那非（Sildenafil)化学式
smi = 'CCCc1nn(C)c2C(=O)NC(=Nc12)c3cc(ccc3OCC)S(=O)(=O)N4CCN(C)CC4'
m = Chem.MolFromSmiles(smi)

# 画结构式
Draw.MolToImageFile(m, "/drug_development/studyRdkit/st_rdcit/img/mol1.jpg")
print('化学结构式Sildenafil')

# vardenafil
smi2 = 'CCCc1nc(C)c2C(=O)N=C(Nn12)c3cc(ccc3OCC)S(=O)(=O)N4CCN(CC)CC4'
m2 = Chem.MolFromSmiles(smi2)
fp2 = Chem.GetMorganFingerprintAsBitVect(m2, 2, nBits=2048)
