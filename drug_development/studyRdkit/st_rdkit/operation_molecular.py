#! /usr/bin/python
# coding: utf-8
# @Time: 2020-05-27 20:11:44
# @Author: zeoy
# rdkit 操作分子对象

from rdkit import Chem
# 1 获取分子中的原子
smi = 'CC(C)OC(=O)C(C)NP(=O)(OCC1C(C(C(O1)N2C=CC(=O)NC2=O)(C)F)O)OC3=CC=CC=C3'
mol = Chem.MolFromSmiles(smi)
atoms = mol.GetAtoms()
atoms_num = mol.GetNumAtoms()
print(atoms)  # <rdkit.Chem.rdchem._ROAtomSeq object at 0x1050ddc10>
print(atoms[0])  # <rdkit.Chem.rdchem.Atom object at 0x10aa13760>
print('类型=', type(atoms))  # 类型= <class 'rdkit.Chem.rdchem._ROAtomSeq'>
print('类型0=', type(atoms[0]))  # 类型0= <class 'rdkit.Chem.rdchem.Atom'>
print("省略氢的原子数=", atoms_num)

# 2 获取原子的坐标信息
# 前提，导入的原子必须带有坐标信息
print(mol.GetConformer().GetAtomPosition(1)[0])
print(mol.GetConformer().GetAtomPosition(1).x)
print(mol.GetConformer().GetAtomPosition(1).y)
print(mol.GetConformer().GetAtomPosition(1).z)
x, y, z = mol.GetConformer().GetAtomPosition(1)
print(x, y, z)
xyz = list(mol.GetConformer().GetAtomPosition(3))
print(xyz)

# 3访问单个原子的信息
# 对原子进行遍历：mol.GetAtoms()
# 获取原子索引：GetIdx()
# 获取原子序号：GetAtomicNum()
# 获取原子符号：GetSymbol()
# 获取原子连接数（受H是否隐藏影响）：GetDegree()
# 获取原子总连接数（与H是否隐藏无关）：GetTotalDegree()
# 获取原子形式电荷：GetFormalCharge()
# 获取原子杂化方式：GetHybridization()
# 获取原子显式化合价：GetExplicitValence()
# 获取原子隐式化合价：GetImplicitValence()
# 获取原子总的化合价：GetTotalValence()

atom = mol.GetAtomWithIdx(0)
print("标签=", atom.GetSymbol())  # C
print("价电子=", atom.GetExplicitValence())  # 4
print("原子元素周期编号=", atom.GetAtomicNum())  # 6

print("杂化类型=", atom.GetHybridization())  # 返回杂化类型   杂化类型= SP3
print("是否在芳香烃内=", atom.GetIsAromatic())  # 该原子是否在芳香烃内   是否在芳香烃内= False
# 与该原子连接的氢原子个数
print("该原子连接的氢原子个数=", atom.GetTotalNumHs())  # 该原子连接的氢原子个数= 3
neighbors = atom.GetNeighbors()  # 返回该原子的所有邻居原子，以元祖的形式返回
print([x.GetAtomicNum() for x in neighbors])  # [6]

# 3.1 访问所有原子：
print('\t'.join(['id', 'num', 'exp', 'symbol', 'degree', 'charge', 'hybrid']))
for at in atoms:
    print(at.GetIdx(), end='\t')
    print(at.GetAtomicNum(), end='\t')
    print(at.GetExplicitValence(), end='\t')
    print(at.GetSymbol(), end='\t')
    print(at.GetDegree(), end='\t')
    print(at.GetFormalCharge(), end='\t')
    print(at.GetHybridization())

# id	num	exp	symbol	degree	charge	hybrid
# 0	6	1	C	1	0	SP3
# 1	6	3	C	3	0	SP3
# 2	6	1	C	1	0	SP3
# 3	8	2	O	2	0	SP2
# 4	6	4	C	3	0	SP2
# 5	8	2	O	1	0	SP2
# 6	6	3	C	3	0	SP3
# 7	6	1	C	1	0	SP3
# 8	7	2	N	2	0	SP3
# 9	15	5	P	4	0	SP3
# 10	8	2	O	1	0	SP2
# 11	8	2	O	2	0	SP3
# 12	6	2	C	2	0	SP3
# 13	6	3	C	3	0	SP3
# 14	6	3	C	3	0	SP3
# 15	6	4	C	4	0	SP3
# 16	6	3	C	3	0	SP3
# 17	8	2	O	2	0	SP3
# 18	7	3	N	3	0	SP2
# 19	6	3	C	2	0	SP2
# 20	6	3	C	2	0	SP2
# 21	6	4	C	3	0	SP2
# 22	8	2	O	1	0	SP2
# 23	7	3	N	2	0	SP2
# 24	6	4	C	3	0	SP2
# 25	8	2	O	1	0	SP2
# 26	6	1	C	1	0	SP3
# 27	9	1	F	1	0	SP3
# 28	8	1	O	1	0	SP3
# 29	8	2	O	2	0	SP2
# 30	6	4	C	3	0	SP2
# 31	6	3	C	2	0	SP2
# 32	6	3	C	2	0	SP2
# 33	6	3	C	2	0	SP2
# 34	6	3	C	2	0	SP2
# 35	6	3	C	2	0	SP2

# 4 分子中的键操作
bonds = mol.GetBonds()  # 对键进行遍历
print(type(bonds))

# 同样，每一个键也都是对象，可以通过属性和函数来获取键的信息。

# 对键进行遍历：m.GetBonds()
# 获取键的索引：GetIdx()
# 获取键的类型：GetBondType()
# 以数字形式显示键的类型：GetBondTypeAsDouble()
# 是否为芳香键：GetIsAromatic()
# 是否为共轭键：GetIsConjugated()
# 是否在环中：IsInRing()
# 是否在n元环中：IsInRingSize(n)
# 获取起始原子：GetBeginAtom()
# 获取末尾原子：GetEndAtom()

print('\t'.join(['id', 'type', 'double', 'aromic',
                 'conjug', 'ring', 'begin', 'end']))
for bond in bonds:
    print(bond.GetIdx(), end='\t')
    print(bond.GetBondType(), end='\t')
    print(bond.GetBondTypeAsDouble(), end='\t')
    print(bond.GetIsAromatic(), end='\t')
    print(bond.GetIsConjugated(), end='\t')
    print(bond.IsInRing(), end='\t')
    print(bond.GetBeginAtomIdx(), end='\t')
    print(bond.GetEndAtomIdx())

# id	type	double	aromic	conjug	ring	begin	end
# 0	SINGLE	1.0	False	False	False	0	1
# 1	SINGLE	1.0	False	False	False	1	2
# 2	SINGLE	1.0	False	False	False	1	3
# 3	SINGLE	1.0	False	True	False	3	4
# 4	DOUBLE	2.0	False	True	False	4	5
# 5	SINGLE	1.0	False	False	False	4	6
# 6	SINGLE	1.0	False	False	False	6	7
# 7	SINGLE	1.0	False	False	False	6	8
# 8	SINGLE	1.0	False	False	False	8	9
# 9	DOUBLE	2.0	False	False	False	9	10
# 10	SINGLE	1.0	False	False	False	9	11
# 11	SINGLE	1.0	False	False	False	11	12
# 12	SINGLE	1.0	False	False	False	12	13
# 13	SINGLE	1.0	False	False	True	13	14
# 14	SINGLE	1.0	False	False	True	14	15
# 15	SINGLE	1.0	False	False	True	15	16
# 16	SINGLE	1.0	False	False	True	16	17
# 17	SINGLE	1.0	False	False	False	16	18
# 18	AROMATIC	1.5	True	True	True	18	19
# 19	AROMATIC	1.5	True	True	True	19	20
# 20	AROMATIC	1.5	True	True	True	20	21
# 21	DOUBLE	2.0	False	True	False	21	22
# 22	AROMATIC	1.5	True	True	True	21	23
# 23	AROMATIC	1.5	True	True	True	23	24
# 24	DOUBLE	2.0	False	True	False	24	25
# 25	SINGLE	1.0	False	False	False	15	26
# 26	SINGLE	1.0	False	False	False	15	27
# 27	SINGLE	1.0	False	False	False	14	28
# 28	SINGLE	1.0	False	False	False	9	29
# 29	SINGLE	1.0	False	True	False	29	30
# 30	AROMATIC	1.5	True	True	True	30	31
# 31	AROMATIC	1.5	True	True	True	31	32
# 32	AROMATIC	1.5	True	True	True	32	33
# 33	AROMATIC	1.5	True	True	True	33	34
# 34	AROMATIC	1.5	True	True	True	34	35
# 35	SINGLE	1.0	False	False	True	17	13
# 36	AROMATIC	1.5	True	True	True	24	18
# 37	AROMATIC	1.5	True	True	True	35	30

# 4.1 也可以通过索引获取键：
print('通过索引获取键', mol.GetBondWithIdx(3).GetBondType())  # 通过索引获取键 SINGLE


# 5 获取分子中所有的环

# 查看所有最小环(smallest set of smallest rings, SSSR）的信息：GetSymmSSSR()
ssr = Chem.GetSymmSSSR(mol)
num_ring = len(ssr)
print("环的个数=", num_ring)  # 环的个数= 3
for ring in ssr:
    print('ring consisted of atoms id:', list(ring))
# ring consisted of atoms id: [14, 13, 17, 16, 15]
# ring consisted of atoms id: [19, 20, 21, 23, 24, 18]
# ring consisted of atoms id: [31, 32, 33, 34, 35, 30]

# 直接获取环的信息：GetRingInfo()
# 查看一共有几个环：NumRings()
# 查看原子在几个环中：NumAtomRings()
# 查看id为n的原子是否在n1元环中.IsAtomInRingOfSize(n, n1)
# 查看id为n的键是否在n1元环中.IsBondInRingOfSize(n , n1)

ri = mol.GetRingInfo()
print('分子中环的个数=', ri.NumRings())  # 分子中环的个数= 3
print(ri.NumAtomRings(2))  # 0
print(ri.IsAtomInRingOfSize(3, 3))  # False
print(ri.IsBondInRingOfSize(2, 3))  # False
