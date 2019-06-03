# searchHotel
### python抓取去哪网当天的酒店信息

. 主要内容
> 环境准备
> selenium 使用
> 数据抓取
---

#### 环境准备
‘’‘

安装selenium
sudo pip install selenium

’‘’
selenium2.x 调用高版本浏览器会出现不兼容问题，调用低版本浏览器正常
selenium3.x 调用浏览器必须下载一个类似不定的文件，比如firefox的geckodriver，chrome的chromedriver
各个浏览器的补丁[下载地址:](http://www.seleniumhq.org/download/)

‘’‘

安装 BeautifulSoup
sudo pip install BeautifulSoup

’‘’

#### selenium 使用
注意事项：

'''

from selenium import webdriver
dr = webdriver.Firefox()

'''

如果运行报错，提示geckodriver(或者其他浏览器对应的补丁)必须在‘PATH’,添加对应的路径到环境变量中，重启，如果还报错，改用下列写法

'''

dr = webdriver.Firefox(execute_path=r"/Users/software/chromedriver.exe")

'''

#### 数据抓取
1. 搜索功能，在搜索框中输入时间地点，点击搜索按钮
2. 获取一页完整数据。由于去哪网一个页面数据分为两次加载，第一次加载15条，这时需要将页面拉到底部，完成第二次数据加载
3. 获取一页完整且经过渲染的HTML文档，使用BeautifulSoup将其中的酒店信息提取出来进行存储
4. 解析完成，点击下一页，继续抽取数据

