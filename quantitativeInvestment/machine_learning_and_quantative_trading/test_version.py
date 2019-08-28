#! /usr/bin/python
# -*- coding:utf-8 -*-

# 查看库模块版本

import tensorflow as tf      # 深度学习、神经网络开发平台
import tensorlayer as tl
import keras as ks
import nltk


import pandas as pd          # 新一代数据分析工具
import tushare as ts         # 国内股票数据采集模块
import matplotlib as mlp     # 绘图模块
import plotly                # 互动型数据可视化工具


import arrow

# import tflearn
tflearn = tf.contrib.learn # 兼容版本冲突

def queryModelVersion():
    print ('\n#1 tensorflow.ver:', tf.__version__)  #1 tensorflow.ver: 1.14.0
    print ('\n#2 tensorlayer.ver:', tl.__version__) #2 tensorlayer.ver: 2.1.0
    print ('\n#3 keras.ver:', ks.__version__)       #3 keras.ver: 2.2.5
    print ('\n#4 nltk.ver:', nltk.__version__)      #4 nltk.ver: 3.4
    print ('\n#5 pandas.ver:', pd.__version__)      #5 pandas.ver: 0.24.2
    print ('\n#6 tushare.ver:', ts.__version__)     #6 tushare.ver: 1.2.41
    print ('\n#7 matplotlib.ver:', mlp.__version__) #7 matplotlib.ver: 3.0.3
    print ('\n#8 plotly.ver:', plotly.__version__)  #8 plotly.ver: 4.1.0
    print ('\n#9 arrow.ver:', arrow.__version__)    # 9 arrow.ver: 0.14.5
    print ('\n#10 tflearn.ver:')                    #10 tflearn.ver:   无版本信息

queryModelVersion()

