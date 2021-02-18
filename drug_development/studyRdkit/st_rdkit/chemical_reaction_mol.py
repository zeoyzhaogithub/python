#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-09 11:58:19
# @Author: zeoy
# rdkit 化学反应

# 一、引入所需库
from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem import Draw


# RDKit提供化学反应引擎，其中化学反应模板是基于smarts构建。反应物+反应引擎就可以生成产物。

# 二、化学反应实例
# 示例反应
rxn = AllChem.ReactionFromSmarts(
    '([Cl;H0;D1;+0:1]-[c;H0;D3;+0:2](:[c:3]):[n;H0;D2;+0:4]:[c:5])>>(C-[n;H0;D3;+0:4](:[c:5]):[c;H0;D3;+0:2](=O):[c:3]).(Cl-P(-Cl)(=O)-[Cl;H0;D1;+0:1])')
img = Draw.ReactionToImage(
    rxn
)
img.save(
    '/drug_development/studyRdkit/st_rdcit/img/mol31.jpg'
)
# 反应模板如下图所示：

# 从反应模板中，我们看到主要的变化是Cl变成羰基氧，N上多了一个甲基
# >注：这是一个逆反应模板
# 反应物如下图所示 ：
mol = Chem.MolFromSmiles(
    'CC(C)(Nc1nc(Cl)c(-c2ccc(F)cc2)c(-c2ccncc2)n1)c1ccccc1')
Draw.MolToImageFile(
    mol,
    "/drug_development/studyRdkit/st_rdcit/img/mol32.jpg",
    size=(350, 300),
    legend='CC(C)(Nc1nc(Cl)c(-c2ccc(F)cc2)c(-c2ccncc2)n1)c1ccccc1'
)
# .创建具体反应规则的引擎对象rxn = AllChem.ReactionFromSmarts(tem)
# .输入反应物，借助引擎产生反应rxn.RunReactants([productmol])


def getrxns(rxn, product_smi):
    """
    获取反应规则的引擎对象
    product_smi 反应物
    """
    product_mol = Chem.MolFromSmiles(product_smi)
    reactions = rxn.RunReactants([product_mol])
    rxns = []
    for reaction in reactions:
        smis = []
        for compound in reaction:
            smi = Chem.MolToSmiles(compound)
            smis.append(smi)

        rxnstr = '.'.join(smis) + '>>' + product_smi
        newr = canon_reaction(rxnstr)
        rxns.append(newr)
    return rxns


tem = '([Cl;H0;D1;+0:1]-[c;H0;D3;+0:2](:[c:3]):[n;H0;D2;+0:4]:[c:5])>>(C-[n;H0;D3;+0:4](:[c:5]):[c;H0;D3;+0:2](=O):[c:3]).(Cl-P(-Cl)(=O)-[Cl;H0;D1;+0:1])'
rxn = AllChem.ReactionFromSmarts(tem)
product_smi = 'CC(C)(Nc1nc(Cl)c(-c2ccc(F)cc2)c(-c2ccncc2)n1)c1ccccc1'
reactions = getrxns(rxn, product_smi)
for reaction in reactions:
    img = Draw.ReactionStringToImage(reaction)
    display(img)

# 三、化学反应模板
# 化学反应模板主要可通过两种方法获取：1. 自动抽取；2. 化学家编码。
# 1. 自动提取 软件：RDKit中没有提供自动提取反应模板的方法。 ASKCOS 开源了自动提取模板的方法。 原理就是化学环境发生变化的原子定义为反应中心，然后基于反应中心拓展一定的半径或者基团。
# 对反应模板进行高亮显示 ：
# 2. 化学家编码 化学家通过对反应进行归纳，整理出的反应规则，典型代表Chemtic 软件。
# 四、化学反应注意事项

# 利用化学反应引擎产生的产物，不一定符合化学规则，因此需要进行对产物进行检查。
tem = '([Cl;H0;D1;+0:1]-[c;H0;D3;+0:2](:[c:3]):[n;H0;D2;+0:4]:[c:5])>>(C-[n;H0;D3;+0:4](:[c:5]):[c;H0;D3;+0:2](=O):[c:3]).(Cl-P(-Cl)(=O)-[Cl;H0;D1;+0:1])'
rxn = AllChem.ReactionFromSmarts(tem)


def getrxns_reactants(rxn, product_smi):
    product_mol = Chem.MolFromSmiles(product_smi)
    reactions = rxn.RunReactants([product_mol])
    rxns = []
    for reaction in reactions:
        smis = []
        for compound in reaction:
            smi = Chem.MolToSmiles(compound)
            smis.append(smi)

        newr = '.'.join(smis)
        rxns.append(newr)
    return rxns


prosmi = "COC(=O)c1cccc(-c2nc(Cl)cc3c(OC)nc(C(C)C)n23)c1"
rs = getrxns_reactants(rxn, prosmi)

print(rs)  # ['COC(=O)c1cccc(-c2n(C)c(=O)cc3c(OC)nc(C(C)C)n23)c1.O=P(Cl)(Cl)Cl']
print(len(rs))  # 1

smi = rs[0]
m = Chem.MolFromSmiles(smi, sanitize=False)

if m is None:
    print('invalid Smiles')
else:
    try:
        mi = Chem.SanitizeMol(m)
        print(mi)
        print('smiles is ok')
    except:
        print('invalid chemistry')  # invalid chemistry
