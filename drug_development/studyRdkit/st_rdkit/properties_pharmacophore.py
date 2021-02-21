#! /usr/bin/python
# coding: utf-8
# @Time: 2020-06-12 18:46:11
# @Author: zeoy
# rdkit 化学性质和药效团

# 一、引入所需库
import os

from rdkit import Chem
from rdkit import RDConfig

from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D


# 二、化学性质
# 建立一个化学性质对象,通过该对象可以得到分子的化学性质

fdefName = os.path.join(
    RDConfig.RDDataDir,
    '/Users/zeoy/st/drug_development/st_rdcit/BaseFeatures.fdef'
)
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

smi = 'C=CC(=O)N1CCC(CC1)C2CCNC3=C(C(=NN23)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=O)N'
m = Chem.MolFromSmiles(smi)
Draw.MolToImageFile(
    m,
    "/Users/zeoy/st/drug_development/st_rdcit/img/mol38.jpg",
)
# 使用特征工厂搜索特征
feats = factory.GetFeaturesForMol(m)
print(len(feats))  # 16

# 搜索到的每个特征都包含了改特征家族（例如受体、供体等）特征类别、该特征对应的原子、特征对应的序号等
for f in feats:
    print(
        f.GetFamily(),  # 特征家族信息
        f.GetType(),    # 特征类型信息
        f.GetAtomIds()  # 特征对应原子
    )

# Donor SingleAtomDonor (4,)
# Donor SingleAtomDonor (13,)
# Donor SingleAtomDonor (34,)
# Acceptor SingleAtomAcceptor (3,)
# Acceptor SingleAtomAcceptor (17,)
# Acceptor SingleAtomAcceptor (25,)
# Acceptor SingleAtomAcceptor (33,)
# Aromatic Arom5 (14, 15, 16, 17, 18)
# Aromatic Arom6 (19, 20, 21, 22, 23, 24)
# Aromatic Arom6 (26, 27, 28, 29, 30, 31)
# Hydrophobe ThreeWayAttach (7,)
# Hydrophobe ThreeWayAttach (15,)
# Hydrophobe ThreeWayAttach (19,)
# Hydrophobe ChainTwoWayAttach (1,)
# LumpedHydrophobe RH6_6 (19, 20, 21, 22, 23, 24)
# LumpedHydrophobe RH6_6 (26, 27, 28, 29, 30, 31)

# 三、化学特征文件介绍
# 以上都是基于现有的特征库进行分析和提取 ， 而特征库就是一个特征定义文件 （ FeatureDefinitionFile, FDef ） 。
# 该文件包含了一系列的化学特征及他们的所有信息 ， 并通过SMARTS来表示 。 除了化学特征 ， FDef文件也有对原子类型的定义及解释 ， 让特征更容易理解 。

# 每个化学特征对应的SMARTS，参考文档
# key为“特征家族.特征类型”，值为SMARTS
fdefName = os.path.join(
    RDConfig.RDDataDir, '/Users/zeoy/st/drug_development/st_rdcit/BaseFeatures.fdef')
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)

print(list(featFactory.GetFeatureDefs().items())[:2])
[
    (
        'Donor.SingleAtomDonor',
        '[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),$([$(n[n;H1]),$(nc[n;H1])])]),$([O,S;H1;+0])]'
    ),
    (
        'Acceptor.SingleAtomAcceptor',
        '[$([O;H1;v2]),$([O;H0;v2;!$(O=N-*),$([O;-;!$(*-N=O)]),$([o;+0])]),$([n;+0;!X3;!$([n;H1](cc)cc),$([$([N;H0]#[C&v4])]),$([N&v3;H0;$(Nc)])]),$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]'
    )
]
# # 3.1 化学特征（chemical features）
# 化学特征由一个特征类型 （ FeatureType ） 和特征家族 （ FeatureFamily ） 共同定义 。
# 特征家族是对特征总体上的分类 ， 例如氢键供体 、 芳香性等 。 药效团匹配就是根据特征家族来实现的 。
# 而特征类型提供了特征的一些更详细的信息 。 每个特征类型包含了以下信息 ：
# .一个smarts表达式 ， 定义了该特征类型所匹配的原子 。
# .特征位置的权重，特征位置由原子位置决定。

# # 3.2 FDef文件语法
# FDef中包含两个概念 ： 原子类型和特征 。 可以大致这么认为 ： 原子类型是最底层的抽象 ， 对相似的原子做归类 。
# 特征是更高级的抽象，对相似的原子和原子类型再进行归类。以下都围绕原子类型定义和特征定义的语法展开。[原文]('http://www.rdkit.org/docs/RDKit_Book.html#the-feature-definition-file-format')

# # # 3.2.1 原子类型（Atom Type）定义
# .原子类型定义相当于给原子的SMARTS起了一个别名 ， 可以让FDef更具有可读性 。
# 例如 ， 可以用以下方式来定义一个非极性碳原子 ， Carbon_NonPolar就是一个原子类型定义名称 ：

# **AtomTypeCarbon_NonPolar[C & ! $ (C=[O, N, P, S])]**
# .因此，可以创建一个有利于理解和使用SMARTS的名称。要引用某个原子类型，可以写下该名称并用大括号括起来。例如，定义另一个原子类型，取名为Hphobe，让Hphobe包含Carbon_NonPolar：
# **
# AtomType Carbon_NonPolar[C &!$(C=[O, N, P, S])]
# AtomType Hphobe[{Carbon_NonPolar}, c, s, S & H0 & v2, F, Cl, Br, I]
# **

# 重复写下一个原子类型时 ， 则意味着把他们合并在一起 ， 相当于用 “ ， ” 连接 （ 在SMARTS中表示或 ） ， 例如下面两个语句 ：
# AtomType d1[N &!HO]
# AtomType d1[N &!HO]
# 等价于 ：
# AtomType d1[N &!H0, O &!H0]

# 更简洁的写法 ：
# AtomType d1[N, O; !H0]

# 要注意“&”和“；”都表示“与”，“，”表示“或”，但“&”的优先级高于“，”，而“；”的优先级低于“，”

# 类似于SMARTS，原子特征定义也可以用“！”来表示“非”，而“！”会与自身的表达式结合，例如：
# AtomType d1[N, O, S]
# AtomType !d1[H0]]

# "!d1"等价于：
# AtomType d1[!H0; N, O, S]

# # # 3.2.2 特征(Feature)定义
# 特征定义比原子类型定义更复杂 ， 并且需要多行实现 ， 例如 ：
# DefineFeatureHDonor1[N, O!H0]
# Family HBondDonor
# Weights 1.0
# EndFeature

# 特征定义的第一行包含了特征类型和所规定的SMARTS，第二行和第三行（这两个行没有先后顺序）定义了特征家族和原子权重（
# 权重值的数量和SMARTS中包含的原子数量相同，有多个权重时，用逗号隔开）。原子权重用来计算特征的位置。最后一行表示结束定义，必须是“EndFeature”。原子类型定义和特征定义可以混在一起使用，只要在引用前定义好就行了。

# # # 3.2.3 其它语法

# 井号“  # ”开头表示注释，程序运行时会忽略该行。
# 反斜杠“\”结尾表示续行，也就是该行没有结束。
# 开头的空格将被忽略。

# 四、2D药效团指纹
# # 4.1 参数设置
# 对上面计算的性质进行组合可以用作分子的2D药效团。 药效团可以进一步转化为分子药效团指纹。
# 参考文件
fdefName = os.path.join(
    RDConfig.RDDataDir, '/Users/zeoy/st/drug_development/st_rdcit/BaseFeatures.fdef')
# 实例化特征工厂
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)

# 使用特征工厂再来构建指纹工厂signature，factory用于设置指纹参数
# 构建指纹工厂 ：
SigFactory(
    featFactory,      # 特征工厂
    useCounts=False,  # 默认False。False不考虑指纹频数，并生成SparseBitVect
    minPointCount=2,  # 默认为2.生成指纹时包括的最少的药效团数量。
    maxPointCount=3,  # 默认为3。生成指纹时包括的最多的药效团数量。
    ...
)
sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3)
# 对拓扑距离进行分段
sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
# 每次修改参数后，都要初始化一下
sigFactory.Init()
# 计算指纹的长度
print('指纹长度=', sigFactory.GetSigSize())   # 指纹长度= 2988

# # 4.2 生成2D药效团指纹
# 指纹工厂中的参数设置完毕，接下来就可以生成2D指纹了。
# 计算2D药效团指纹 ：
Gen2DFingerprint(
    mol,  # 要计算指纹的mol对象
    sigFactory,  # 设置了参数的指纹工厂
    bitinfo   # 获取指纹id及对应的原子
)

mol = Chem.MolFromSmiles('OCC(=O)CCCN')
fp = Generate.Gen2DFingerprint(mol, sigFactory)
print(len(fp))   # 2988
print(fp.GetNumOnBits())  # 23

# 关于指纹每一位所代表特征的信息、特征的距离矩阵等信息，都可以通过signature factory来查看
print(list(fp.GetOnBits())[:5])  #
print(sigFactory.GetBitDescription(1))  # Acceptor Acceptor |0 1|1 0|

# # 4.3 修改FDef设置
# 如果不需要某个特征，可以直接通过signature factory来跳过某个特征，而不用去修改FDef文件。

# 查看现有药效团（列表）列表：GetFeatureFamilies()

featureFamilies = featFactory.GetFeatureFamilies()
print(featureFamilies)
# ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable',
#  'ZnBinder', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

# 跳过某个药效团（特征家族）：sigFactory.skipFeats
# 每次修改都要初始化以下 ： init()
# 在查看指纹长度 ： GetSigSize()

sigFactory.skipFeats = ['PosIonizable']
sigFactory.Init()
print(sigFactory.GetSigSize())  # 2100

# 重新生成新的指纹：Gen2DFingerprint()
fp2 = Generate.Gen2DFingerprint(mol, sigFactory)
print(fp2.GetNumOnBits())   # 15

# # 4.4 Gobbi  2D 药效团指纹
# Rdkit中还有一种用于生成2D药效团指纹的特征定义方式 ， 根据Gobbi等人的设计实现 ， 在
# rdkit.Chem.Pharm2D.Gobbi_Pharm2D下有一个预定义的signaturefactory,
# 已经包含了这些指纹类型，可以直接调用，操作方法类似

m = Chem.MolFromSmiles('OCC=CC(=O)O')
fp = Generate.Gen2DFingerprint(m, Gobbi_Pharm2D.factory)
fp_num_bit = fp.GetNumOnBits()
print(fp_num_bit)  # 8

print(list(fp.GetOnBits()))  # [23, 30, 150, 154, 157, 185, 28878, 30184]

gobbi_description = Gobbi_Pharm2D.factory.GetBitDescription(157)
print(gobbi_description)   # HA HD |0 3|3 0|

# 本文参考自rdkit(官方文档)[]
