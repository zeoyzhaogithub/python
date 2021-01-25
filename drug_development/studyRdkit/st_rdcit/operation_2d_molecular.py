#! /usr/bin/python
# coding: utf-8
# @Time: 2020-05-29 16:19:04
# @Author: zeoy
# rdkit 处理2D、3D分子

# 一、引入所需库
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import pickle

# Smiles 可以看成分子的1D形式，分子的平面结构可以看成分子的2D形式。
# 该算法能够减少分子中原子在平面内的碰撞，使得绘制的分子更加清晰。

# 二、处理2D分子
# 可用于作为绘制分子的模板。

# 如果有多个分子共享一个骨架，我们希望绘图的时候能够保证在这些分子中的骨架方向是一致的。 首先把骨架提取成绘图模板，然后在绘图上添加基团。
# 我们可以看到这4个结构都含有吡咯并嘧啶（c1nccc2n1ccc2）的子结构， 但是在图片上他们的子结构取向不一致，不便于比较他们的结构。 我们可以采用模板绘制法，固定公共子结构的取向。

smis = [
    'COC1=C(C=CC(=C1)NS(=O)(=O)C)C2=CN=CN3C2=CC=C3',
    'C1=CC2=C(C(=C1)C3=CN=CN4C3=CC=C4)ON=C2C5=CC=C(C=C5)F',
    'COC(=O)C1=CC2=CC=CN2C=N1',
    'C1=C2C=C(N=CN2C(=C1)Cl)C(=O)O',
]
template = Chem.MolFromSmiles('c1nccc2n1ccc2')
## 2.1 计算分子的2D坐标函数解析
Compute2DCoords（
    （Mol）mol   # 目标分子
    [,（bool）canonOrient = True #以规范的方式定向分子
    [,（bool）clearConfs = True # 如果为真，则该分子上的所有现有构象将被清除
    [,（dict）coordMap = {}  # 映射原子ID的字典-> Point2D对象与应该锁定其位置的原子的起始坐标。
    [,（int）nFlipsPerSample = 0  # 可旋转债券的数量 一次随机翻转。
    [,（int）nSample = 0   # 可旋转键的随机采样数。
    [,（int ）sampleSeed = 0  # 随机抽样过程的种子。
    [,（bool）permuteDeg4Nodes = False  # 允许键的排列为4度
    [,（float）bondLength = -1.0   # 更改描述的默认键长
    [,（bool）forceRDKit = False 
    ] ] ] ] ] ] ] ] ] 
)

## 2.2 生成一个分子的描述函数解析
GenerateDepictionMatching2DStructure(
    （Mol）mol, # mol-要排列的分子，它会回来一个单一的conformer。
    （Mol）reference  # 与参考原子对齐的分子；这应该有一个描述。
    [,（int）confId = -1  # confId-（可选）要使用reference的参考构象的idPattern-（可选）要用于的查询分子生成分子和参考之间的原子映射。
    [,（AtomPairsParameters）refPatt = None 
    [,（bool）acceptFailure = False  # （可选）如果为True，则将生成标准描述  对于没有与参考匹配的子结构的分子；如果为False，则抛出DepictException。
    [,（bool）forceRDKit = False 
    ] ] ] ] 
)
## 2.3 计算分子的2D坐标  结果坐标存储在分子的每个原子上
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
img.save('/Users/zeoy/st/drug_development/st_rdcit/img/mol9.jpg')

# 三、处理3D分子
# RDKit可以使用两种不同的方法为分子生成构象。原始方法使用距离几何。请注意，此过程产生的构象往往很难看。应使用力场对其进行清理。可以在RDKit中使用通用力场（UFF）的实现来完成。还有一种Riniker和Landrum方法的实现，该方法使用剑桥几何结构数据库（CSD）中的扭转角首选项在已经使用距离几何图形生成它们之后校正矫正器。使用这种方法，无需使用最小化步骤来清理结构。
# 从RDKit的2018.09版本开始，ETKDG是默认的构象生成方法。

# 3.1 产生3D构象
# 默认情况下，RDKit分子在图中没有明确存在H原子，但是H原子对于获得逼真的几何形状很重要，因此通常在产生3D构象前，需要为分子添加H原子。
mol = Chem.MolFromSmiles('c1nccc2n1ccc2')
m2 = Chem.AddHs(mol) # 加氢原子
AllChem.EmbedMolecule(m2) # 2D->3D化
AllChem.MMFFOptimizeMolecule(m2)   # 使用MMFF94最小化RDKit生成的构象
m3 = Chem.RemoveHs(m2) # 删除氢原子

# 3.2 产生多个3D构象
AllChem.EmbedMultipleConfs(m2, numConfs=10) # 为m2分子产生10个构象，保存在m2分子中。
m2.GetConformer(1)  # 访问指定构象
m2.GetConformers()  # 获取构象

# 产生多个构象比较耗时，支持多线程加速。
cids = AllChem.EmbedMultipleConfs(m2,numConfs=10, numThreads=4)
print('构象个数=',len(cids))


# 四、 计算构象的RMS值
## 4.1 计算其他构象和第一个构象的RMS值
rmslist = []
AllChem.AlignMolConformers(m2, RMSlist=rmslist)
print('第一个构想者与所有其他构想者之间的RMS值')
print(rmslist)
# [0.031617482029182735, 0.027100784678672763, 0.03211558617566241, 0.03876593076230512, 0.03981170998684297, 0.023608843156566642, 0.023619922568325093, 0.037621957801700204, 0.0532861818792842]

## 4.2 计算指定两个构象的RMS值
rms = AllChem.GetConformerRMS(m2, 1, 9, prealigned=True)
# print('计算指定两个构象的RMS值=', rms) # 计算指定两个构象的RMS值= 0.049340647238401286

## 4.3 MMFF立场对构象进行优化
res = AllChem.MMFFOptimizeMoleculeConfs(m2)
print(res)
# [(0, 44.63060837612333), (0, 44.6306083325595), (0, 44.63060829914013), (0, 44.63060834780363), (0, 44.63060829995848), (0, 44.63060839000017), (0, 44.630608353905146), (0, 44.6306086577808), (0, 44.630608361427775), (0, 44.63060856588376)]
# 构象优化比较耗时，支持多线程加速
res = AllChem.MMFFOptimizeMoleculeConfs(m2, numThreads=10)

# 五、保存分子对象
### RDKit内置了2种保存分子对象的方法。
## 5.1使用pickling机制进行保存。
pkl = pickle.dumps(mol)
print(pkl)
m2 = pickle.loads(pkl)
mol2 = Chem.MolToSmiles(m2)
print(mol2)  # c1nccc2n1ccc2
## 5.2 保存为普通二进制文件
binStr = mol.ToBinary()
m2 = Chem.Mol(binStr)
mol2 = Chem.MolToSmiles(M2)
### rdkit的pickling文件非常紧凑，从pickling文件中加载分子比从mol文件和SMILES字符串开很多。
### 将经常使用的分子，保存为picklingz是一个很好的办法。
## 二进制文件的大小和分子大小有关系，而pickle文件和分子大小关系不明显。推荐用pickle保存文件。
