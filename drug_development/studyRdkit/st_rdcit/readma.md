# Rdkit 学习记录

## 简介

Rdkit 是干个什么呢，RDkit 是一款开源化学信息学与机器学习工具包，提供C++ 和python 的API 接口。是著名的开源化学信息学工具之一，基于BSD协议，核心数据结构与算法由C++编写。支持Python2与python3，支持KNIME，支持机器学习方面的分子描述符的产生。

!(文档：)[https://www.rdkit.org/]

## 安装

1.Conda模式 官方建议使用Conda进行安装与管理，Conda可以使用清华的源进行下载，安装完成后，再次更换其安装源。换源的命令参考

```bash
conda create -c rdkit -n my-rdkit-env rdkit
```

激活的方法：

```bash
conda activate my-rdkit-env
```

退出环境

```bash
conda deactivate
```

2：Pycharm模式 Pycharm并不能直接安装RDkit，当使用上一步Conda安装完成后，将python的环境配置修改为conda的目录，即可在pycharm中使用Rdkit。

