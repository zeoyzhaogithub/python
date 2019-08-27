## 第十三章、补充知识Tabix 和SQLite

&emsp;&emsp;在本章中，我们将研究内存不足的方法-围绕存储和处理磁盘上内存之外的数据构建的计算策略。从磁盘读取数据比使用存储器中的数据要慢得多(参见第45页上的“全能Unix管道：速度和美观合一”)，但在许多情况下，这是我们在内存中(例如，将整个数据集加载到R中)或流式方法(例如，使用Unix管道，就像我们在第7章中所做的那样)是不合适的。具体地说，我们将研究两种处理内存中数据的工具：Tabix和SQLite数据库。

## 使用BGZF和Tabix快速访问索引的制表符分隔的文件

### 使用Bgzip压缩Tabix的文件
### 使用Tabix索引文件
### 使用Tabix

## 通过SQLite引入关系数据库

### 在生物信息学中何时使用关系数据库
### 安装SQLite
### 使用命令行界面浏览SQLite数据库

!['suggestion'](../img/suggestion.png)
### 查询数据：全能的SELECT命令
!['suggestion'](../img/suggestion.png)

#### 使用limit限制结果
#### 使用SELECT选择列
#### 使用ORDER BY对行进行排序
#### 使用WHERE筛选行

### SQLite函数
### SQLite聚合函数
#### 使用GROUP BY对行进行分组
### 子查询
### 用joins连接关系数据库
#### 连接关系数据库
#### 内连接
#### 左连接

!['suggestion'](../img/suggestion.png)
### 写入数据库
#### 创建表
#### 在表中写入数据
#### 索引
### 删除数据表和数据库
### 从Python与SQLite交互
#### 连接到SQLite数据库并从Python创建表
#### 将数据从Python加载到表中

### 转储数据库
&emsp;&emsp;最后，让我们谈谈数据库转储。数据库转储是完全复制数据库所必需的所有SQL命令。数据库转储对于备份和复制数据库非常有用。转储在共享数据库方面也很有用，尽管在SQLite中，简单地共享数据库文件更容易(但是对于其他数据库引擎，如MySQL和PostgreSQL，这是不可能的)。SQLite使转储数据库变得非常容易。我们可以使用sqlite3命令行工具：

```sql
$ sqlite3 variants.db ".dump"
PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE  variants(
    id integer primary key,
    chrom test,
    start integer,
    end integer,
    strand text,
    rsid text);
INSERT INTO  "variants" VALUES(1,'chr10',114808901,114808902,'+','rs12255372');
INSERT INTO "variants" VALUES(2,'chr10',114808901,114808902,'+','rs12255372');
INSERT INTO "variants" VALUES(3,'chr10',114808901,114808902,'+','rs12255372');
INSERT INTO "variants" VALUES(4,'chr10',114808901,114808902,'+','rs12255372');
COMMIT;
```
 
&emsp;&emsp;如果希望转储单个表而不是整个数据库，.dump dot命令还会采用可选的表名称参数。我们可以使用数据库转储来创建数据库：
```shell
$ sqlite3 variants.db ".dump" &gt; dump.sql
$ sqlite3 variants-duplicate.db &lt; dump.sql
```

&emsp;&emsp;这一系列命令将variants.db数据库中的所有表转储到SQL文件dump.sql中。然后，将这个SQL文件加载到新的空数据库variantsDuplicate.db中，创建所有表并将所有数据插入到原始variants.db数据库中。
