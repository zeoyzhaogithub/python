#! /usr/bin/python
# -*- coding:utf-8 -*-

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

fss = open('./data/iris.csv')
df = pd.read_csv(fss, index_col=False, header=1) # 以第一行作为索引

fss2 = open('./data/iris2.csv')
df2 = pd.read_csv(fss2, index_col=False) 
def ana():
    '''
    分析数字型数据
    '''
    print(df.tail())   # 展示数据末尾5行

    # 属种是字符串，无法进行统计分析
    print(df.describe())  #相当于对数据集进行概览，会输出该数据集的计数、最大值、最小值等

def xname():
    '''
    分析字符串型数据
    '''
    d10 = df['属种'].value_counts()
    print(d10)
    # virginica（维吉尼亚鸢尾）     51
    # versicolor（变色鸢尾）    50
    # setosa（山鸢尾）        50

def updateIrisSet():
    df.loc[df['属种'] == 'virginica', 'xid'] = 1
    df.loc[df['属种'] == 'setosa', 'xid'] = 2
    df.loc[df['属种'] == 'versicolor', 'xid'] = 3
    df['xid'] = df['xid'].astype(int)
    df.to_csv('./data/iris2.csv', index=False)
    print('succ')

def ana_iris2():
    '''
    分析数字型数据
    '''
    print(df2.tail())   # 展示数据末尾5行

    # 属种是字符串，无法进行统计分析
    print(df2.describe())  #相当于对数据集进行概览，会输出该数据集的计数、最大值、最小值等

def xname_iris2():
    '''
    分析字符串型数据
    '''
    d10 = df2['xname'].value_counts()
    print(d10)
    # virginica（维吉尼亚鸢尾）     51
    # versicolor（变色鸢尾）    50
    # setosa（山鸢尾）        50


def break_data():
    '''
    设置测试数据集
    '''
    # 设置总数据源
    xlst,ysgn = ['x1','x2','x3','x4'], 'xid'
    x,y = df2[xlst],df2[ysgn]
    # print(x.tail())
    # print(y.tail())

    # 生成x_train,x_test,y_train,y_test
    # train_test_split随机划分训练集和测试集
    x_train,x_test,y_train,y_test = train_test_split(x,y,random_state=1)
    x_test.index.name,y_test.index.name = 'xid','xid'

    # 数据类型
    # print('type(x_train)',type(x_train))
    # print('type(x_test)',type(x_test))
    # print('type(y_train)',type(y_train))
    # print('type(y_test)',type(y_test))

    # 保存数据
    fs0 = './data/iris_'
    x_train.to_csv(fs0 + 'xtrain.csv',index=False)
    x_test.to_csv(fs0 + 'xtest.csv',index=False)
    y_train.to_csv(fs0 + 'ytrain.csv',index=False,header=True)
    y_test.to_csv(fs0 + 'ytest.csv',index=False,header=True)  # 保留数据头


def trading_data():
    '''
    用线性回归的方法训练数据
    '''
    # 读取训练数据集
    fs0 = './data/iris_'
    x_train = pd.read_csv(fs0 + 'xtrain.csv',index_col=False)
    y_train = pd.read_csv(fs0 + 'ytrain.csv',index_col=False)
    print(x_train.tail())
    # 对数据进行分析、学习，并建立模型LinearRegression
    mx = LinearRegression()
    mx.fit(x_train.values, y_train.values) # 建模
    # print(mx)
    # print(mx.intercept_)
    # print(mx.coef_)

    # 测试数据
    x_test  = pd.read_csv(fs0 + 'xtest.csv',index_col=False)
    df9 = x_test.copy()
    # 生成预测结果
    y_pred = mx.predict(x_test.values)
    df9['y_predsr'] = y_pred

    # 训练数据
    y_test  = pd.read_csv(fs0 + 'ytest.csv',index_col=False)
    df9['y_test'],df9['y_pred'] = y_test,y_pred
    df9['y_pred'] = round(df9['y_predsr']).astype(int)  # 训练结果
    df9.to_csv('./data/iris_9.csv', index=False)
    print(df9.tail())

    # 检测测试结果

# ana()
# xname()

#updateIrisSet()

# ana_iris2()
# xname_iris2()
# break_data()

trading_data()
