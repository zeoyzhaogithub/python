# stockAnalysis
python 股票数据分析

### 开发环境

pandas_datareader 是重要的 pandas 相关包，原来是 `pandas.io.data` 方法，是一个远程获取金融数据的Python工具

    pip install pandas
    pip install pandas-datareader

雅虎财经和谷歌财经的接口变换频繁。如果用 `pip install pandas_datareader`，已经无法得到雅虎财经。

该修正版本的安装方法是

    $ git clone https://github.com/rgkimball/pandas-datareader
    $ cd pandas-datareader
    $ git checkout fix-yahoo
    $ pip install -e