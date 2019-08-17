## 第六章、生物信息学数据

&emsp;&emsp;到目前为止，我们已经介绍了生物信息学入门的许多基础知识：组织项目目录、中间Unix、使用远程机器以及使用版本控制软件。然而，我们忽略了一个新的生物信息学项目的一个重要组成部分：数据。
      
&emsp;&emsp;数据是任何生物信息学项目的必要条件。我们通过将大量数据提炼到可以从中提取意义的点来进一步理解复杂的生物系统。不幸的是，对于基因组学中常见的大型和复杂数据集，许多简单的小或中型数据集的任务是一个挑战。这些挑战包括：
 
*检索数据。*
        
&emsp;&emsp;无论是下载大型测序数据集还是访问Web应用程序数百次以下载特定文件，在生物信息学中检索数据都需要特殊的工具和技能。

*确保数据完整性。*
        
&emsp;&emsp;跨网络传输大型数据集会带来更多数据损坏的机会，这可能会导致以后的分析不正确。因此，在继续分析之前，我们需要使用工具确保数据的完整性。同样的工具也可以用来验证我们在分析中使用的数据版本是否正确。

*压缩。*
        
&emsp;&emsp;我们在生物信息学中使用的数据足够大，经常需要压缩。因此，处理压缩数据是生物信息学中的一项基本技能

## 检索生物信息学数据
       
&emsp;&emsp;假设您刚刚被告知您的项目的排序已经完成：您有六个通道的Illumina数据要从您的排序中心下载。通过Web浏览器下载如此大量的数据是不可行的：Web浏览器不是为下载如此大的数据集而设计的。此外，您需要将这些测序数据下载到您的服务器，而不是您浏览Internet的本地工作站。要做到这一点，您需要SSH进入数据处理服务器，并使用命令行工具将数据直接下载到这台机器上。我们将在本节中查看其中的一些工具。

## 使用wget和curl下载数据
      
&emsp;&emsp;用于从Web下载数据的两个常用命令行程序是wget和curl。根据您的系统，这些可能尚未安装；您必须使用包管理器(例如，Homebrew或apt-get)安装它们。虽然curl和wget在基本功能上是相似的，但它们的相对优势非常不同，以至于我经常使用这两个

## wget
      
&emsp;&emsp;wget对于从命令行快速下载文件非常有用-例如，从GRCH37(也称为hg19)汇编版本下载人类染色体22：
```shell
$ wget http://hgdownlod.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz--2013-06-30 00:15:45-- http://[...]/goldenPath/hg19/chromosomes/chr22.fa.gz
Resolving hgdownload.soe.ucsc.edu... 128.114.119.163
Connecting to hgdownload.soe.ucsc.edu|128.114.119.163|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 11327826 (11M) [application/x-gzip]
Saving to: ‘chr22.fa.gz’
17% [======&gt; ] 1,989,172 234KB/s eta 66s
```
       
&emsp;&emsp;Wget将此文件下载到您的当前目录，并提供一个有用的进度指示器。请注意，到22号染色体的链接以“http”(超文本传输协议的缩写)开头。Wget还可以处理FTP链接(以“ftp”开头，即文件传输协议的缩写)。通常，对于大文件，FTP比HTTP更可取(UCSC Genome Browser等网站经常推荐使用FTP)。
       
&emsp;&emsp;因为UCSC慷慨地公开提供人类参考基因组，我们不需要任何身份验证就可以访问这个文件。但是，如果您从Lab Information Management System(LIMS)下载数据，则可能需要首先使用用户名和密码进行身份验证。对于简单的HTTP或FTP身份验证，可以使用wget的-user=和-ask-password选项进行身份验证。一些站点使用更复杂的身份验证过程；在这些情况下，您需要与站点管理员联系。
       
&emsp;&emsp;wget的优势之一是它可以递归地下载数据。当使用递归选项(-recursive或-r)运行时，wget还将跟踪和下载链接到的页面，甚至跟踪和下载这些页面上的链接，等等。默认情况下(为了避免下载整个Internet)，wget将自身限制为只跟踪五个深度的链接(但这可以使用选项-level或-l进行自定义)。
       
&emsp;&emsp;递归下载可用于从页面下载特定类型的所有文件。例如，如果另一个实验室有一个包含我们希望下载的许多GTF文件的网页，我们可以使用：
```shell
$ wget --accept "*.gtf" --no-directories --recursive --no-parent \http://genomics.someuniversity.edu/labsite/annotation.html
```
       
&emsp;&emsp;但要小心！Wget的递归下载可以是相当激进的。如果不受约束，wget将在-level设置的最大深度内下载它可以到达的所有内容。在前面的示例中，我们通过两种方式限制wget：使用-no-parent防止wget下载目录结构中较高的页面，以及使用-accept“*.gtf”，它只允许wget下载与此模式匹配的文件名。
       
&emsp;&emsp;使用wget的递归选项时要谨慎；它可能会利用大量网络资源并使远程服务器过载。在某些情况下，如果您下载的太多太快，远程主机可能会阻止您的IP地址。如果您计划从远程主机下载大量数据，最好向网站的sysadmin查询建议的下载限制，这样您的IP就不会被阻止。wget的-limitrate选项可以用来限制wget下载文件的速度。
       
&emsp;&emsp;Wget是一个非常强大的工具；前面的例子仅仅是其功能的表面。请参阅表6-1了解一些常用的选项，或者man wget获取详尽的列表。
       
表6-1、常用的wget选项
 
| option | values | use |
|:-------|:-------|:----|
|--A, --accept| sufx(如“.fast q”)或带有*、？或[和]的模式，可选逗号分隔列表。|仅下载符合此条件的文件。|
|--R, --reject| 和 --accept 选项的功能一样 |不要下载与此匹配的文件；例如，要下载页面上除Excel文件之外的所有文件，请使用-reject“.xls”。|
|-nd，--n- 目录|No value|不要将本地下载的文件与远程文件放在相同的目录层次结构中。|
|-r, --recursive| 
No value |跟踪并下载页面上的链接，默认情况下最大深度为5。|
|-np, --no-parent| No value |不要移动到父目录的上方。|
|--limit-rate |每秒允许的字节数。|节流下载带宽。|
|--user=user|FTP或HTTP用户名。|HTTP或FTP身份验证的用户名。|
|-ask-password|No value |提示输入FTP身份验证的HTTP密码；-password=也可以使用，但是您的密码在shell的历史记录中|
|-O|输出文件名|将文件下载到指定的文件名；如果链接没有信息性名称(例如，http://lims.squencingcenter.com/seqs.html?id=sample_A_03)，则非常有用|

## Curl
      
&emsp;&emsp;curl与wget的用途略有不同。Wget非常适合通过HTTP或FTP下载文件，并使用其递归选项从网页中抓取数据。curl的行为类似，尽管默认情况下会将文件写入标准输出。要像我们使用wget那样下载22号染色体，我们将使用：
```shell
$ curl http://[...]/goldenPath/hg19/chromosomes/chr22.fa.gz &gt; chr22.fa.gz % Total % Received % Xferd Average Speed Time Time Time Current
 Dload Upload Total Spent Left Speed
 14 10.8M 14 1593k 0 0 531k 0 0:00:20 0:00:02 0:00:18 646k
```

&emsp;&emsp;请注意，我必须截断URL，以便这个示例适合页面；URL与我们之前在wget中使用的相同。
        
&emsp;&emsp;如果不希望重定向curl的输出，请使用-O<filename>将结果写入<flename>。如果省略filename参数，curl将使用与远程主机相同的文件名。
        
&emsp;&emsp;curl的优点是它可以使用比wget更多的协议来传输文件，包括SFTP(安全FTP)和SCP(安全拷贝)。curl的一个特别好的特性是，如果启用了-L/-location选项，它可以跟随页面重定向。启用此功能后，curl将下载链接重定向到的最终页面，而不是重定向页面本身。最后，Curl本身也是一个库，这意味着除了命令行curl程序之外，Curl的功能由RCurl和pycurl等软件库包装。

## rsync和安全拷贝(SCP)
        
&emsp;&emsp;Wget和curl适用于从命令行快速下载文件，但不适合某些更繁重的任务。例如，假设一位同事需要项目目录中被Git忽略的所有大型测序数据集(例如，在.gitignore中)。在网络上同步这些整个目录的一个更好的工具是rsync。
        
&emsp;&emsp;对于这些类型的任务，rsync是一个更好的选项，原因有几个。首先，rsync通常更快，因为它只发送文件版本之间的差异(当副本已经存在或部分存在时)，并且它可以在传输过程中压缩文件。其次，rsync有一个存档选项，可以保留链接、修改时间戳、权限、所有权和其他文件属性。这使得rsync成为整个目录的网络备份的最佳选择。rsync还具有许多功能和选项来处理不同的备份场景，例如，如果远程主机上存在文件，该怎么办。
        
&emsp;&emsp;rsync的基本语法是rsync source destination，其中source是要复制的文件或目录的源，destination是要将这些文件复制到的目标。源或目标可以是user@host：/path/to/directory/格式指定的远程主机。
        
&emsp;&emsp;让我们看一个例子，说明如何使用rsync将整个目录复制到另一台机器上。假设您希望将zea_mays/data中的所有项目数据复制到主机192.168.237.42上同事的目录/home/deborah/zea_mays/data。用于复制整个目录的最常见的rsync选项组合是-avz。选项-a启用wrsync的归档模式，-z启用文件传输压缩，-v使rsync的进度更详细，以便您可以看到正在传输的内容。因为我们将通过SSH连接到远程主机，所以我们还需要使用-e ssh。我们的目录复制命令将如下所示：
```shell
$ rsync -avz -e ssh zea_mays/data/ vinceb@[...]:/home/deborah/zea_mays/databuilding file list ... done
zmaysA_R1.fastq
zmaysA_R2.fastq
zmaysB_R1.fastq
zmaysB_R2.fastq
zmaysC_R1.fastq
zmaysC_R2.fastq
sent 2861400 bytes received 42 bytes 107978.94 bytes/sec
total size is 8806085 speedup is 3.08
```
        
&emsp;&emsp;rsync的一个微妙但重要的行为是，当指定rsync中的路径时，尾随斜杠(例如，data/vs data)是有意义的。源路径中的尾部斜杠表示复制源目录的内容，而没有尾部斜杠表示复制整个目录本身。因为我们想复制zea_mays/data/的所有内容到/home/deborah/zea_mays/data我们的例子中，我们使用了一个尾部斜杠。如果远程目标主机上尚不存在data/directory，我们将希望使用zea_mays/data复制它及其内容(例如，省略尾随的斜杠)。
        
&emsp;&emsp;因为rsync只在文件不存在或已更改的情况下传输文件，所以您可以(并且应该)在文件传输后再次运行rsync。这是一种简单的检查，以确保两个目录之间的所有内容都是同步的。在脚本中调用rsync时检查rsync的退出状态也是一个好主意；如果在传输文件时遇到问题，rsync将以非零状态退出。最后，rsync可以使用SSH配置文件中指定的主机别名(参见第57页“使用SSH连接到远程机器”中的第一个提示)。如果通过SSH主机别名连接到主机，则可以省略-e ssh。
        
&emsp;&emsp;有时，我们只需要通过SSH快速复制单个文件-对于Unix的cp就足够了，但需要通过SSH连接工作的任务。rsync可以工作，但有点过火了。安全拷贝(SCP)是这一目的完美之选。Secure Copy的工作方式与cp类似，不同之处在于我们需要同时指定host和path(使用与wget相同的user@host：/path/to/file符号)。例如，我们可以将单个GTF文件传输到192.168.237.42：/home/deborah/zea_mays/data/使用：
```shell
$ scp Zea_mays.AGPv3.20.gtf 192.168.237.42:/home/deborah/zea_mays/data/Zea_mays.AGPv3.20.gtf 100% 55 0.1KB/s 00:00
```

## 数据完整性
      
&emsp;&emsp;我们下载到项目目录中的数据是所有未来分析和结论的起点。虽然看起来不太可能，但在传输大型数据集时，传输过程中数据损坏的风险是一个令人担忧的问题。这些大文件需要很长时间才能传输，这意味着网络连接中断和位丢失的机会更多。除了验证传输已完成而没有错误外，使用校验和显式检查传输数据的完整性也很重要。校验和是非常压缩的数据汇总，以这样一种方式计算，即使只改变了数据的一位，校验和也会不同。
      
&emsp;&emsp;数据完整性检查也有助于跟踪数据版本。在合作项目中，我们的分析可能依赖于同事的中间结果。当这些中间结果发生变化时，所有依赖于这些结果的下游分析都需要重新运行。对于许多中间文件，并不总是清楚哪些数据已经更改，哪些步骤需要重新运行。如果数据改变了，即使是最小的比特，校验和也会不同，所以我们可以用它们来计算数据的版本。校验和还有助于重复性，因为我们可以将特定的分析和结果集链接到由数据的校验和值汇总的数据的精确版本。

## SHA和MD5校验      
      
&emsp;&emsp;用于确保数据完整性的两个最常用的校验和算法是MD5和SHA-1。我们已经在第4章中遇到了SHA-1，因为这是Git用于其提交ID的。MD5是一种较旧的校验和算法，但仍是常用的一种算法。MD5和SHA-1的行为相似，但SHA-1较新，通常首选。然而，MD5更为常见；如果服务器对一组文件进行了预先计算的校验和，则很可能会遇到这种情况。
      
&emsp;&emsp;让我们了解一下使用SHA-1的校验和。我们可以通过以下中的标准将任意字符串传递给程序shasum(在某些系统上，它是sha1sum)：
```shell
$ echo "bioinformatics is fun" | shasumf9b70d0d1b0a55263f1b012adab6abf572e3030b -
$ echo "bioinformatic is fun" | shasum
e7f33eedcfdc9aef8a9b4fec07e58f0cf292aa67 
```
      
&emsp;&emsp;由数字和字母组成的长字符串是SHA-1校验和。校验和以十六进制格式报告，其中每个数字可以是16个字符中的一个：数字0到9以及字母a、b、c、d、e和f。后面的破折号表示这是来自标准输入的SHA-1校验和。请注意，当我们在“生物信息学”中省略“s”并计算SHA-1校验和时，校验和值完全改变。这就是使用校验和的强大之处：当输入的最小部分发生变化时，校验和就会发生变化。校验和是完全确定的，这意味着无论计算校验和的时间或使用的系统如何，校验和值只有在输入不同的情况下才会不同。
      
&emsp;&emsp;我们也可以对文件输入使用校验和(请注意，Csyrichta_Tag-GACT_L008_R1_001.fastq的内容是假示例数据)：
```shell
$ shasum Csyrichta_TAGGACT_L008_R1_001.fastqfea7d7a582cdfb64915d486ca39da9ebf7ef1d83 Csyrichta_TAGGACT_L008_R1_001.fastq
```
&emsp;&emsp;如果我们的测序中心说Csyrichta_Tag-GACT_L008_R1_001.fastq.gz测序文件的校验和是“069bf5894783db241e26f4e44201bd12f2d5aa42”，那么我们本地的SHA校验和是“069bf5894783db241e26f4e44201bd12f2d5aa42。“Feed7d7a582cdfb64915d486ca39da9ebf7ef1d83”，我们知道我们的文件与原始版本有所不同。
      
&emsp;&emsp;当下载许多文件时，逐个检查每个校验和可能会变得相当单调乏味。程序shasum有一个方便的解决方案-它可以针对包含文件校验和的文件进行创建和验证。我们可以为data/目录中的所有FASTQ文件创建一个SHA-1校验和文件，如下所示：
```shell
$ shasum data/*fastq &gt; fastq_checksums.sha$ cat fastq_checksums.sha
524d9a057c51b1[...]d8b1cbe2eaf92c96a9 data/Csyrichta_TAGGACT_L008_R1_001.fastq
d2940f444f00c7[...]4f9c9314ab7e1a1b16 data/Csyrichta_TAGGACT_L008_R1_002.fastq
623a4ca571d572[...]1ec51b9ecd53d3aef6 data/Csyrichta_TAGGACT_L008_R1_003.fastq
f0b3a4302daf7a[...]7bf1628dfcb07535bb data/Csyrichta_TAGGACT_L008_R1_004.fastq
53e2410863c36a[...]4c4c219966dd9a2fe5 data/Csyrichta_TAGGACT_L008_R1_005.fastq
e4d0ccf541e90c[...]5db75a3bef8c88ede7 data/Csyrichta_TAGGACT_L008_R1_006.fastq
```
       
&emsp;&emsp;然后，我们可以使用shasum的check选项(-c)来验证这些文件是否与原始版本匹配：
```shell
$ shasum -c fastq_checksums.shadata/Csyrichta_TAGGACT_L008_R1_001.fastq: OK
data/Csyrichta_TAGGACT_L008_R1_002.fastq: OK
data/Csyrichta_TAGGACT_L008_R1_003.fastq: OK
data/Csyrichta_TAGGACT_L008_R1_004.fastq: OK
data/Csyrichta_TAGGACT_L008_R1_005.fastq: OK
data/Csyrichta_TAGGACT_L008_R1_006.fastq: FAILED
shasum: WARNING: 1 computed checksum did NOT match
```
        
&emsp;&emsp;如果文件的校验和不一致，shasum将显示哪个文件未通过验证并以非零错误状态退出。
        
&emsp;&emsp;程序md5sum(或OSX上的md5)计算md5散列，其操作类似于shasum。但是，请注意，在OSX上，md5命令没有-c选项，因此您需要为该选项安装GNU版本。此外，一些服务器使用过时的校验和实现，例如sum或chsum。我们如何使用这些旧的命令行校验和程序类似于使用shasum和md5sum。
        
&emsp;&emsp;最后，您可能想知道如何通过40个字符长的SHA-1校验和来总结所有文件。他们不能-只有1640种可能的不同校验和。然而，1640是一个巨大的数字，因此校验和冲突的概率非常小。为了检查数据完整性，发生冲突的风险可以忽略不计。

## 查看数据之间的差异
       
&emsp;&emsp;虽然校验和是检查文件是否不同的一种很好的方法，但它们不能告诉我们文件的不同之处。一种方法是使用Unix工具diff计算两个文件之间的dif。UNIX的diff逐行工作，并输出文件之间不同的块(称为大块)(类似于我们在第4章中看到的Git的git diff命令)。
       
&emsp;&emsp;假设您注意到协作者正在使用与您正在使用的文件版本不同的文件版本。她的版本是gene-2.bed，而你的版本是gene-1.bed(这些文件在GitHub上，如果你愿意跟随的话)。因为下游结果依赖于这个数据集，所以需要检查文件是否确实不同。在比较SHA-1校验和之后，您发现文件不完全相同。在使用协作者的版本重新运行分析之前，您想要检查文件是否存在显著差异。我们可以通过计算gene-1.bed和gene-2.bed之间的差异来做到这一点：
```shell
$ diff -u gene-1.bed gene-2.bed--- gene-1.bed 2014-02-22 12:53:14.000000000 -0800  (1)
+++ gene-2.bed 2015-03-10 01:55:01.000000000 -0700
@@ -1,22 +1,19@@   (2)
 1 6206197 6206270 GENE00000025907
 1 6223599 6223745 GENE00000025907   (3)
 1 6227940 6228049 GENE00000025907
+1 6222341 6228319 GENE00000025907   (4)
 1 6229959 6230073 GENE00000025907
-1 6230003 6230005 GENE00000025907   (5)
 1 6233961 6234087 GENE00000025907
 1 6234229 6234311 GENE00000025907
 1 6206227 6206270 GENE00000025907
 1 6227940 6228049 GENE00000025907
 1 6229959 6230073 GENE00000025907
-1 6230003 6230073 GENE00000025907  (6)
+1 6230133 6230191 GENE00000025907
 1 6233961 6234087 GENE00000025907
 1 6234229 6234399 GENE00000025907
 1 6238262 6238384 GENE00000025907
-1 6214645 6214957 GENE00000025907
 1 6227940 6228049 GENE00000025907
 1 6229959 6230073 GENE00000025907
-1 6230003 6230073 GENE00000025907
 1 6233961 6234087 GENE00000025907
 1 6234229 6234399 GENE00000025907
-1 6238262 6238464 GENE00000025907
 1 6239952 6240378 GENE00000025907
```      
&emsp;&emsp;选项-u告诉diff以Unifed dif格式输出，这是一种与git diff使用的格式几乎相同的格式。我选择使用统一的Diffs而不是diff的默认diff格式，因为统一的Diffs提供了更多的上下文。

&emsp;&emsp;统一差异被分解为两个文件之间不同的区块。让我们逐步了解此格式的关键部分：
       
(1)这两行是统一diff的表头。原始文件gene-1.bed前缀为-，修改后的文件gene-2.bed前缀为+。这两行中的日期和时间是这些文件的修改时间

(2)这一行表示更改后的大块的开始。@@和@@之间的整数对分别表示块在原始文件(-1，22)和修改文件(+1，19)中的开始位置以及长度。
  
(3)diff中以空格开头的行表示修改后的文件的行没有更改。
  
(4)diff中以+开头的行表示已将行添加到修改后的文件中。
  
(5)同样，-表示在修改的文件中删除的行。
  
(6)相邻行删除和添加表示此行在修改后的文件中发生了更改。
      
&emsp;&emsp;diff文件一开始看起来非常神秘，但随着时间的推移，您会逐渐熟悉它们。diff的输出也可以重定向到文件，这会创建一个补丁文件。修补程序文件充当有关如何更新纯文本文件的说明，进行diff文件中包含的更改。Unix工具补丁可以将更改应用于需要打补丁的文件。补丁在软件开发中比生物信息学更常用，因此我们不会详细介绍它们。最后，需要注意的是，在大型文件上计算Diffs的计算成本可能很高，因此在对大型数据集运行diff时要谨慎。

## 压缩数据和处理压缩数据
        
&emsp;&emsp;数据压缩，压缩数据的过程，使其占用更少的空间(在磁盘驱动器上，在内存中，或通过网络传输)，是现代生物信息学中不可缺少的技术。例如，来自最近的Illumina HiSeq运行的序列在使用Gzip压缩时占用21，408，674，240字节，这是低于20G字节的比特。未压缩，这个文件是惊人的63，203，414，514字节(约58G字节)。这个FASTQ文件有1.5亿个200bp的读数，这是人类基因组的10倍覆盖，拟南芥基因组的190倍覆盖，或者六倍体小麦基因组的2倍覆盖。该数据的压缩比(未压缩大小/压缩大小)约为2.95，这意味着显著节省了约66%的空间。你自己的生物信息学项目可能会包含更多的数据，特别是随着测序成本持续下降，有可能将基因组测序到更高的深度，在表达研究中包括更多的生物复制或时间点，或者在基因分型研究中对更多的个体进行测序。
        
&emsp;&emsp;在大多数情况下，数据可以在整个处理和分析过程中保持压缩在磁盘上。大多数编写良好的生物信息学工具都可以压缩数据作为输入，而不需要我们首先将其解压缩到磁盘。使用管道和重定向(第3章中介绍)，我们可以流式传输压缩数据并将压缩文件直接写入磁盘。此外，常见的Unix工具(如cat、grep和less)都有处理压缩数据的变体，而Python的gzip模块允许我们从Python中读取和写入压缩数据。因此，虽然在生物信息学中处理大型数据集可能具有挑战性，但使用unix和软件库中的压缩工具使我们的生活变得容易多了。
## gzip
        
&emsp;&emsp;Unix上使用的两个最常用的压缩系统是gzip和bzip2。两者都有各自的优点：gzip压缩和解压缩数据的速度比bzip2快，但bzip2具有更高的压缩率(前面提到的FASTQ文件在使用bzip2压缩时只有大约16 GB)。通常，gzip在生物信息学中用于压缩最大的文件，而bzip2更常见于长期数据归档。我们将主要关注gzip，但是bzip2的工具的行为与gzip非常相似。
        
&emsp;&emsp;命令行工具gzip允许您以几种不同的方式压缩文件。首先，gzip可以压缩标准输入的结果。这是非常有用的，因为我们可以直接压缩来自另一个生物信息学程序的标准输出的结果。例如，假设我们有一个从FASTQ文件中删除低质量碱基的程序，称为trimmer(这是一个虚构的程序)。我们的trimmer程序可以在本地处理gzip压缩的输入文件，但将未压缩的修剪后的FASTQ结果写入标准输出。使用gzip，我们可以在写入磁盘之前就地压缩trimmer的输出：
```shell
$ trimmer in.fastq.gz | gzip &gt; out.fastq.gz
```
        
&emsp;&emsp;gzip从标准输入中获取输入，对其进行压缩，并将压缩后的输出写入标准输出。
        
&emsp;&emsp;gzip还可以就地压缩磁盘上的文件。如果我们的in.fastq.gz文件没有压缩，我们可以按如下方式压缩它：
```shell
$ lsin.fastq
$ gzip in.fastq
$ ls
in.fastq.gz
```
       
&emsp;&emsp;gzip将在适当的位置压缩此文件，用压缩文件替换原始的未压缩版本(将扩展名.gz附加到原始文件名)。类似地，我们可以使用gunzip命令将文件解压缩到适当的位置：
```shell
$ gunzip in.fastq.gz$ ls
in.fastq
```
      
&emsp;&emsp;请注意，这将用解压缩版本替换我们的in.fastq.gz文件(也删除了.gz后缀)。gzip和gunzip都可以将它们的结果输出到标准输出(而不是在适当的地方更改文件)。可以使用-c选项启用此功能：
```shell
$ gzip -c in.fastq &gt; in.fastq.gz$ gunzip -c in.fastq.gz &gt; duplicate_in.fastq
```
       
&emsp;&emsp;gzip压缩算法的一个很好的特性是，您可以将gzip压缩输出直接连接到现有的gzip文件中。例如，如果我们想压缩in2.fastq文件并将其附加到压缩的in.fastq.gz文件中，我们就不必首先解压缩in.fastq.gz，连接两个文件，然后压缩连接的文件。相反，我们可以执行以下操作：

```shell
$ lsin.fastq.gz in2.fastq
$ gzip -c in2.fastq &gt;&gt; in.fastq.gz
```
      
&emsp;&emsp;重要的是，请注意，我们使用的重定向运算符是>>如果我们使用>，我们会将in2.fast q的压缩版本覆盖为in.fast q.gz(而不是附加到它)。使用重定向时一定要谨慎，并确保您使用的是适当的运算符(并保留文件备份！)。通过将所有内容压缩在一起(例如，使用cat in.fastq|gzip>|in.fast q.gz)，您可能会获得稍好的压缩比，但是附加到现有gzip文件的便利性是有用的。另外，请注意，gzip不会分离这些压缩文件：压缩在一起的文件是连接在一起的。如果需要将多个单独的文件压缩到单个归档文件中，请使用tar实用程序(有关详细信息，请参阅man tar的示例部分)。

## 使用压缩文件
      
&emsp;&emsp;也许gzip(和bzip2)的最大优势是许多常见的Unix和生物信息学工具可以直接处理压缩文件。例如，我们可以使用grep的gzip文件类似物zgrep来搜索压缩文件。同样，cat具有zcat(在某些系统如OS X上，这是gzcat)，diff具有zdiff，less具有zless。如果程序无法处理压缩输入，您可以使用zcat和管道输出直接到另一个程序的标准输入。
      
&emsp;&emsp;这些处理压缩输入的程序的行为与标准的对应程序完全相同。例如，grep中的所有可用选项都在zgrep中可用：
```shell
$ zgrep --color -i -n "AGATAGAT" Csyrichta_TAGGACT_L008_R1_001.fastq.gz2706: ACTTCGGAGAGCCCATATATACACACTAAGATAGATAGCGTTAGCTAATGTAGATAGATT
```
     
&emsp;&emsp;在处理gzip文件时可能会有轻微的性能损失，因为您的CPU必须首先解压缩输入。通常，zgrep、zless和zcat等z工具的便利性和节省的磁盘空间超过了任何潜在的性能命中。

## 案例研究：重复下载数据
       
&emsp;&emsp;可重现地下载数据可能具有欺骗性的复杂性。我们通常通过互联网从远程服务器下载基因组资源，如序列和注释文件，这在未来可能会发生变化。此外，序列和注释数据的新版本可能会发布，因此我们必须记录关于如何获取数据的所有内容，以实现完全的可重复性。作为这一点的演示，让我们通过一个案例研究。我们将为老鼠(Mus Musculus)下载一些基因组和序列资源，并记录我们是如何获得它们的。
       
&emsp;&emsp;对于这个例子，我们将下载GRCm38小鼠参考基因组和附带的注释。请注意，此案例研究涉及下载大型文件，因此您可能不想跟随这些示例。通过基因组参考联盟协调小鼠、人类(Homo Sapiens)和斑马鱼(Danio Rerio)基因组的释放。GRCm38中的“GRC”前缀指的是基因组参考联盟。我们可以使用wget从Ensembl(该联盟的成员)下载GRCm38。对于这个和本节中的其他示例，我不得不截断URL，使它们适合图书的页面宽度；如果您跟随，请参阅本章GitHub上的Readme.md以获得复制和粘贴的完整链接。
```shell
$ wget ftp://ftp.ensembl.org/[...]/Mus_msculus.GRCm38.74.dna.toplevel.fa.gz
```
          
&emsp;&emsp;Ensembl的网站提供了许多有机体的参考基因组、注释、变异数据和其他有用文件的链接。此FTP链接来自于导航[http://www.ensembl.org](http://www.ensembl.org)点击鼠标项目页面，然后点击“下载DNA序列”链接。如果我们要记录如何下载此文件，我们的Markdown Readme.md可能包括如下内容：
```shell
Mouse (*Mus musculus*) reference genome version GRCm38 (Ensemblrelease 74) was downloaded on Sat Feb 22 21:24:42 PST 2014, using:

$ wget ftp://ftp.ensembl.org/[...]/Mus_msculus.GRCm38.74.dna.toplevel.fa.gz
```
        
&emsp;&emsp;我们可能需要查看这些文件包含的染色体、脚手架和重叠群，以此作为健全性检查。此文件是一个gzip压缩的Fasta文件，因此我们可以通过浏览正则表达式“^>”快速查看所有序列标头，该表达式匹配以>(Fasta标头)开头的所有行。我们可以使用zgrep程序提取这个gzip文件上的FASTA标头：
```shell
$ zgrep "^>" Mus_musculus.GRCm38.74.dna.toplevel.fa.gz | less
```
       
&emsp;&emsp;Ensembl还在父目录中提供了一个名为CHECKSUMS的校验和文件。此校验和文件包含使用较旧的Unix工具SUM计算的校验和。我们可以使用SUM程序将我们的校验和值与校验和中的值进行比较：
```shell
$ wget ftp://ftp.ensembl.org/pub/release-74/fasta/mus_msculus/dna/CHECKSUMS$ sum Mus_musculus.GRCm38.74.dna.toplevel.fa.gz
53504 793314
```
      
&emsp;&emsp;校验和53504与条目mus_musculus.GRCm38.74.dna.toplevel.fa.gz的校验和文件中的条目一致。我还喜欢将所有重要数据的SHA-1总和包括在我的数据README.md文件中，以便未来的合作者可以验证他们的数据文件与我使用的数据文件完全相同。让我们使用shasum计算SHA-1总和：
```shell
$ shasum Mus_musculus.GRCm38.74.dna.toplevel.fa.gz01c868e22a981[...]c2154c20ae7899c5f Mus_musculus.GRCm38.74.dna.toplevel.fa.gz
```
       
&emsp;&emsp;然后，我们可以将此SHA-1总和复制并粘贴到我们的readme.md中。接下来，我们可以从Ensembl下载附带的GTF和该目录的校验和文件：
```shell
$ wget ftp://ftp.ensembl.org/[...]/Mus_msculus.GRCm38.74.gtf.gz$ wget ftp://ftp.ensembl.org/[...]/CHECKSUMS
```
       
&emsp;&emsp;同样，让我们确保我们的校验和与校验和文件中的校验和匹配，并对此文件运行shasum以获得我们自己的文档：
```shell
$ sum Mus_musculus.GRCm38.74.gtf.gz00985 15074
$ shasum cf5bb5f8bda2803410bb04b708bff59cb575e379 Mus_musculus.GRCm38.74.gtf.gz
```
       
&emsp;&emsp;再一次，我们将SHA-1复制到我们的readme.md中。到目前为止，我们的readme.md可能如下所示：
```shell
## Genome and Annotation DataMouse (*Mus musculus*) reference genome version GRCm38 (Ensembl
release 74) was downloaded on Sat Feb 22 21:24:42 PST 2014, using:

wget ftp://ftp.ensembl.org/[...]/Mus_musculus.dna.toplevel.fa.gz

Gene annotation data (also Ensembl release 74) was downloaded from Ensembl on
Sat Feb 22 23:30:27 PST 2014, using:

wget ftp://ftp.ensembl.org/[...]/Mus_msculus.GRCm38.74.gtf.gz
```

&emsp;&emsp;虽然这不是很多文档，但这比不记录数据是如何获取的要好得多。如本例所示，只需很少的努力就可以正确跟踪进入您的项目的数据，从而确保可重复性。记录您的工作最重要的一步是您的工作是一致的，并使其成为一种习惯。