# 第三部分 实践：生物信息学数据技能

## 第七章、Unix 数据工具

>We often forget how science and engineering function. Ideas come from previous
exploration more often than from lightning strokes.
—John W. Tukey

&emsp;&emsp;在第3章中，我们学习了Unix shell的基础知识：使用流、重定向输出、管道和处理进程。这些核心概念不仅允许我们使用shell来运行命令行生物信息学工具，而且可以利用Unix作为处理生物信息学数据的模块化工作环境。在本章中，我们将看到如何将Unix shell与命令行数据工具结合起来，以快速浏览和操作数据。

## UNIX数据工具和Unix单行方法：Pearls编程经验

&emsp;&emsp;了解如何在生物信息学中使用Unix数据工具不仅仅是学习每个工具做什么，它还涉及掌握将工具连接在一起的实践-从Unix管道创建程序。通过将数据工具与管道连接在一起，我们可以构建分析、操作和汇总数据的程序。UNIX管道可以在shell脚本中开发，也可以作为“one-liners”(通过将Unix工具与shell上的管道直接连接而构建的微型程序)来开发。无论是在脚本中还是作为一个线人，从小的、模块化的工具构建更复杂的程序利用了Unix的设计和哲学(在“为什么我们在生物信息学中使用Unix？模块化和Unix哲学“(参见第37页)。构建程序的流水线方法在Unix(和生物信息学)中是一个很好的传统，因为它是解决问题的快速方法，非常强大，并且可以适应各种问题。简单Unix管道的力量的一个例证来自两位杰出的计算机科学家之间著名的交流：Donald Knuth和Doug McIlroy(回想一下第3章，McIlroy发明了Unix管道)。

&emsp;&emsp;1986年，在ACM杂志的通讯专栏中有一篇文章“Pearls编程”，专栏作家Jon Bentley让计算机科学家Donald Knuth编写了一个简单的程序，将k个最常用的单词按降序计数并打印在一个文件中。选择Knuth编写这个程序是为了演示文字编程，这是Knuth开创的一种编程方法。识字程序以文本文档的形式编写，解释如何解决编程问题(用简单的英语)，代码散布在整个文档中。然后，可以使用识字的编程工具将文档中的代码“缠绕”在文档之外(这种方法可能会被熟悉R的编织或Sweave的读者所识别-两者都是这个概念的现代后裔)。Knuth的识字程序有七页长，也针对这个特定的编程问题进行了高度定制；例如，Knuth实现了一个自定义的数据结构，用于计算英语单词的数量。Bentley然后要求McIlroy批评Knuth的七页长的解决方案。McIlroy赞扬了Knuth的识字编程和新颖的数据结构，但总体上不同意他的工程方法。McIlroy回复了一个6行的Unix脚本，解决了同样的编程问题：

```shell
tr -cs A-Za-z '\n' | (1)tr A-Z a-z |  (2)
sort |   (3)
uniq -c |  (4)
sort -rn |  (5)
sed ${1}q   (6)
```

&emsp;&emsp;您现在不用担心不能完全理解这写脚本代码(我们将在本章学习这些工具)，McIlroy的基本方法是：

(1) 将所有非字母字符(-c取第一个参数的补码)转换为换行符，并在翻译后将所有相邻字符挤在一起(-s)。这将为整个输入流创建单字线。

(2) 将所有大写字母转换成小写。

(3) 排序输入，将相同的单词放在连续的行上

(4) 删除所有重复的连续行，只保留一行，并对出现次数进行计数(-c)。

(5) 按反向(-r)数字顺序(-n)排序。

(6) 打印脚本(${1})的第一个参数提供的前k行，然后退出。

&emsp;&emsp;shell(假设k在这里是10)。但是，我不得不在这里添加一个换行符，这样代码就不会扩展到页边距之外:

```shell
$ cat input.txt \| tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 10q
```

&emsp;&emsp;McIlroy的脚本无疑比Knuth的程序实现得快得多，而且工作得也一样好(而且可以说更好，因为Knuth的解决方案中有一些小bug)。此外，他的解决方案是建立在可重用的Unix数据工具(或如他所称的“Unix订书机”)上，而不是“从头开始单片编程”，使用McIlroy的措辞。这种方法的速度和力量就是为什么它是生物信息学工作的核心部分。

## 何时使用Unix管道方法以及如何安全地使用它

&emsp;&emsp;虽然McIlroy的例子很有吸引力，但Unix的一行程序方法并不适用于所有问题。许多生物信息学任务都可以通过一个定制的、有良好文档记录的脚本更好地完成，这更类似于Knuth在“Programming Pearls”中的程序。知道何时使用快速而简单的工程解决方案(如Unix管道)以及何时求助于编写文档化的Python或R脚本需要经验。与生物信息学中的大多数任务一样，选择一个合适的方法，则项目已经成功了一半。

&emsp;&emsp;直接进入命令行的UNIX管道作为一个快速、低级的数据操作工具包，用于探索数据、在格式之间转换数据以及检查数据是否存在潜在问题。在这种情况下，我们并不是在寻找彻底的、理论上支离破碎的答案-我们通常只是想快速的获取要我们的数据生成的图片。我们愿意牺牲一个文档化的实现来解决特定的问题，而不是倾向于从模块化的Unix工具构建一个快速的粗略图景。正如麦克罗伊在他的回复中解释的那样：

>The simple pipeline … will suffice to get answers right now, not next week or next month. It could well be enough to finish the job. But even for a production project … it would make a handsome down payment, useful for testing the value of the answers and for smoking out follow-on questions.
—Doug McIlroy (my emphasis)

&emsp;&emsp;生物信息学中的许多任务都具有这种性质：我们希望得到快速的答案，并继续推进我们的项目。我们可以编写一个定制的脚本，但是对于简单的任务来说，这样做可能有些过分，花费的时间会比必要的多。正如我们将在本章后面看到的，构建Unix管道的速度很快：我们可以直接在shell中迭代组装和测试Unix管道。

&emsp;&emsp;对于更大、更复杂的任务，通常更可取的是用Python(如果工作涉及大量数据分析，则为R)这样的语言编写自定义脚本。虽然shell方法(无论是一行程序还是shell脚本)都很有用，但这些方法不允许在检查输入数据、结构化程序、使用数据结构、代码文档以及添加Assert语句和测试等方不如如Python和R灵活。这些语言还具有更好的工具，用于逐步记录更大数据量的分析，就像R的knitr(在第254页的“使用Knitr和Rmarkdown的再生性”中介绍)和IPython。相比之下，冗长的Unix管道可能比自定义脚本脆弱且不够健壮。

&emsp;&emsp;那么，在必须使用Unix管道的情况下，我们可以采取哪些步骤来确保它们是可复制的？正如第1章中提到的，产生结果的每一个步骤都要有文档记录，这一点很重要。因为Unix单行程序是直接在shell中输入的，所以特别容易忘记哪个单行程序产生了哪个版本的输出。记住记录一行程序需要额外的勤奋(而且经常被忽略，特别是在生物信息学工作中)。将管道存储在脚本中是一种很好的方法-脚本不仅可以用作对数据执行哪些步骤的文档，而且还允许重新运行管道，并可以将其检入Git存储库。我们将在第12章中更详细地研究脚本。

## 使用Unix工具检查和操作文本数据

&emsp;&emsp;在本章中，我们的重点是学习如何使用核心Unix工具来操作和探索纯文本格式的数据。生物信息学中的许多数据都是由字符分隔的简单表格或纯文本文件。生物信息学中使用的最常见的表格纯文本文件格式是制表符分隔的。这并非偶然：大多数Unix工具(如cut和awk)在默认情况下将制表符视为分隔符。生物信息学逐渐倾向于制表符分隔的格式，因为使用Unix工具处理这些文件非常方便。制表符分隔的文件格式也易于使用Python和Perl等脚本语言进行解析，并且易于加载到R中。

<hr>

**表格纯文本数据格式**

&emsp;&emsp;表格纯文本数据格式在计算中被广泛使用。基本格式非常简单：每一行(也称为记录)保存在自己的行上，每一列(也称为字段)由一些分隔符分隔。您将遇到三种风格：制表符分隔、逗号分隔和变量空格分隔。

&emsp;&emsp;在这三种格式中，制表符分隔是生物信息学中最常用的格式。BED、GTF/GFF、SAM、表格BLAST输出和VCF等文件格式都是制表符分隔文件的示例。制表符分隔文件的列由单个制表符分隔(其转义码为\t)。常见的约定(但不是标准)是在制表符分隔的文件的前几行包含元数据。这些元数据行以#开头，以区别于表格数据记录。由于制表符分隔的文件使用制表符分隔列，因此不允许使用制表符。

&emsp;&emsp;逗号分隔值(CSV)是另一种常用格式。CSV与制表符类似，不同之处在于分隔符是逗号字符。虽然在生物信息学中并不常见，但以CSV格式存储的数据可能包含逗号(这会干扰解析它的能力)。有些变体就是不允许这样做，而另一些变体则在可能包含逗号的条目周围使用引号。不幸的是，没有标准的CSV格式定义如何处理这个问题和CSV的许多其他问题-尽管在[RFC 4180](https://tools.ietf.org/html/rfc4180)中给出了一些指导原则。

&emsp;&emsp;最后，还有空格分隔的格式。一些顽固的生物信息学程序使用可变数量的空格来分隔列。一般来说，与空格分隔格式相比，制表符分隔格式和CSV是更好的选择，因为经常会遇到包含空格的数据。

&emsp;&emsp;尽管表格数据格式很简单，但有一个主要的常见问题：行是如何分离的。Linux和OS X使用单个换行符(使用转义代码\n)来分隔行，而Windows使用回车符和换行符的DOS样式的行分隔符(\r\n)。CSV文件通常也使用这种DOS样式，因为这是在CSV规范RFC-4180(实际上是松散地遵循)中指定的。偶尔，您可能也会遇到仅由回车分隔的文件。
​​
<hr>

&emsp;&emsp;在本章中，我们将使用非常简单的基因组特征格式：BED(三列)和GTF文件。这些文件格式以制表符分隔的格式存储特征的位置，例如基因、外显子和变体。不要太担心这些格式的细节；我们将在第9章中更详细地讨论这两种格式。本章的目标主要是开发使用Unix数据工具自由操作纯文本文件或流的技能。我们将分别学习每个工具，并累积工作到更高级的管道和程序。

### 用头和尾检查数据

&emsp;&emsp;我们在生物信息学中遇到的许多文件太长了，无法用cat-running cat对一个文件进行检查，一百万行长的文件很快就会用文本滚动来填充你的外壳，滚动速度太快，无法理解。更好的选择是用head查看文件的头部。在这里，让我们看一下文mus_musculus.GRCm38.75_chr1.bed：

```shell
$ head Mus_musculus.GRCm38.75_chr1.bed1 3054233 3054733
1 3054233 3054733
1 3054233 3054733
1 3102016 3102125
1 3102016 3102125
1 3102016 3102125
1 3205901 3671498
1 3205901 3216344
1 3213609 3216344
1 3205901 3207317
```

&emsp;&emsp;我们还可以通过-n参数控制使用head看到的行数：

```shell
$ head -n 3 Mus_musculus.GRCm38.75_chr1.bed1 3054233 3054733
1 3054233 3054733
1 3054233 3054733
```

&emsp;&emsp;head对于快速检查文件很有用。head-n3允许您快速检查文件，以查看是否存在列标头、有多少列、正在使用什么分隔符、一些示例行等等。

&emsp;&emsp;HEAD有一个相关的命令，用于查看文件的末尾或尾部。尾部的工作方式与head相同：

```shell
$ tail -n 3 Mus_musculus.GRCm38.75_chr1.bed1 195240910 195241007
1 195240910 195241007
1 195240910 195241007
```
       
&emsp;&emsp;我们也可以使用tail来删除文件的头。通常，-n参数指定要包含文件的最后几行，但如果给-n一个数字x，前面有一个+号(例如，+x)，尾部将从第X行开始。因此，为了删除标题，我们从第二行开始，使用-n+2。在这里，我们将使用命令seq生成一个包含3个数字的文件，并剪切第一行：
```shell
$ seq 3 &gt; nums.txt$ cat nums.txt
1
2
3
$ tail -n +2 nums.txt
2
3
```
        
&emsp;&emsp;有时查看文件的开始和结束都是有用的-例如，如果我们有一个排序的BED文件，并且我们想要查看第一个特征和最后一个特征的位置。我们可以使用数据科学家(和前生物信息学家)塞斯·布朗(Seth Brown)的技巧来做到这一点：
       
&emsp;&emsp;这是一个有用的技巧，但是打字有点长。为了方便使用，我们可以在您的shell配置文件中创建一个快捷方式，即~/.bashrc或~/.Profle：
```shell
# inspect the first and last 3 lines of a filei() { (head -n 2; tail -n 2) &lt; "$1" | column -t}
```

&emsp;&emsp;然后，在shell配置文件上运行source，或者启动一个新的终端会话并确保其正常工作。然后我们可以使用i(表示Inspect)作为普通命令：
```shell
$ i Mus_musculus.GRCm38.75_chr1.bed1 3054233 3054733
1 3054233 3054733
1 195240910 195241007
1 195240910 195241007
```
     
&emsp;&emsp;head对于查看Unix管道产生的数据也很有用。例如，假设我们想要grep包含字符串gene_id“ENSMUSG00000025907”的行的mus_musculus.GRCm38.75_chr1.gtf文件(因为我们的GTF结构良好，可以安全地假设这些都是属于这个基因的特征-但情况可能并不总是这样！)。我们将使用grep的结果作为流水线中下一个程序的标准输入，但是首先我们想检查grep的标准输出，看看是否一切看起来都是正确的。我们可以通过管道直接将grep中的标准输出到head查看：
```shell
$ grep 'gene_id "ENSMUSG00000025907"' Mus_musculus.GRCm38.75_chr1.gtf | head -n 11 protein_coding gene 6206197 6276648 [...] gene_id "ENSMUSG00000025907" [...]
```
&emsp;&emsp;注意，为了清楚起见，我省略了这个GTF的整行，因为它相当长。

&emsp;&emsp;打印完数据的前几行以确保管道正常工作后，head进程退出。这是一个重要的功能，有助于确保您的管道不会继续处理不必要地数据。当head退出时，shell会捕捉到这一点并停止整个管道，包括grep进程。在幕后，shell向管道中的其他程序发送一个名为SIGPIPE的信号-非常类似于按Control-c时发送的信号(该信号是SIGINT)。当构建处理大量数据的复杂管道时，这是极其重要的。这意味着在像这样的管道中：
```shell
$ grep "some_string" huge_file.txt | program1 | program2 | head -n 5
```
       
&emsp;&emsp;grep不会继续搜索ge_fle.txt，并且program 1和program 2在head输出5行并退出后不会继续处理输入。虽然HEAD很好地说明了管道的这一特性，但SIGPIPE适用于所有程序(除非程序明确捕捉并忽略此符号-这是一种可能性，但我们在生物信息学程序中不会遇到这种可能性)。

### less
       
&emsp;&emsp;less也是一个有用的程序，用于检查文件和管道的输出。less是一个终端寻呼机，一个允许我们在终端中查看大量文本的程序。通常，如果我们将长文件复制到屏幕上，文本会在瞬间闪过-允许我们查看和滚动长文件，并一次标准输出一个屏幕。其他应用程序可以调用默认的终端分页程序来处理大量输出的显示；这就是git log显示整个Git存储库的提交历史的方式。您可能会遇到另一个常见的，但较老的终端寻呼机称为More，但less具有更多功能，并且通常是首选的(less这个名称是对“Less is More”的玩法)。
       
&emsp;&emsp;less运行起来更像是一个应用程序而不是命令：一旦我们启动less，它将保持打开直到我们退出它。让我们回顾一个例子-在本章的GitHub存储库中的目录中，有一个名为contentinated.fastq的文件。让我们用更少的内容来看这个：
```shell
$ less contaminated.fastq
```
      
&emsp;&emsp;将会使用less命令在您的终端中打开程序，并显示一个充满序列的FASTQ文件。首先，如果您需要退出Less，请按q。在任何时候，您都可以使用h调出Less的所有命令的帮助页面。
       
&emsp;&emsp;在less中移动很简单：按空格键向下一页，按b向上一页。您可以使用j和k一次上下一行(这些键与编辑器Vim用于上下移动的键相同)。要返回到文件顶部，请输入g；要返回到文件底部，请按G。Table 7-1列出了最常用的less命令。稍后我们将讨论当less通过管道从另一个程序获取输入时，这是如何工作的。
      
表7-1。常用的less命令

|快捷键|作用|
|:-----|:---|
|空格|下一页|
|b|上一页|
|g|第一页|
|G|最后一页|
|j|向下一行|
|k|向上一行|
|/<pattern>|向上搜索字符<pattern>|
|？<pattern>|向下搜索字符<pattern>|
|n|向下重复上次搜索|
|N|向上重复上次搜索|
      
&emsp;&emsp;less最有用的特性之一是，它允许您搜索文本并突出显示匹配项。视觉上突出显示匹配可能是查找数据中潜在问题的一种非常有用的方法。例如，让我们使用LEASH来快速了解continated.fastq文件中是否存在3'适配器污染物。在本例中，我们将查找AGATCGGAAGAGCACGTCTGAACTCCAGTCAC(来自Illumina TruSeq®kit1的已知适配器)。我们的目标不是进行详尽的测试或移除这些适配器-我们只想利用30秒，看看是否有任何迹象表明可能存在污染。

&emsp;&emsp;搜索整个字符串帮助不大，原因如下：

1.	很可能只有一部分适配器在序列中。
2.	序列中有一些不匹配是常见的，使得精确匹配无效(特别是因为碱基调用准确度通常在Illumina测序读数的3'末端下降)

&emsp;&emsp;为了解决这个问题，让我们搜索前11个碱基AGATCGGAAGA。首先，我们在Less中打开contentinated.fast q，然后按/并输入AGATCGG。结果如图7-1所示，它通过了眼睛间的测试-结果正好击中了你的眼睛之间。注意在序列读数末尾的匹配位置上的倾斜(我们预期在那里被污染)，以及在匹配后碱基的高度相似性。虽然只是一个快速的目视检查，这是相当丰富的信息。
       
&emsp;&emsp;LESS在调试我们的命令行管道时也非常有用。Unix管道的一大优点是可以在任何时候轻松调试-只需通过管道将您要调试的命令的输出减少到较少，然后删除之后的所有内容。当您运行管道时，Less将捕获最后一个命令的输出并暂停，以便您可以检查它。
       
&emsp;&emsp;在迭代构建管道时，Less也是至关重要的-这是构建管道的最佳方式。假设我们有一个假想的管道，它涉及三个程序，Step1，Step2和Step3。我们完成的管道看起来像step1 input.txt|Step2|step3>output.txt。然而，我们希望将其分成几个部分，首先运行step1 input.txt并检查其输出，然后添加step3并检查输出，依此类推。这样做的自然方法是使用更少的：
```shell
$ step1 input.txt | less # inspect output in less
$ step1 input.txt | step2 | less
$ step1 input.txt | step2 | step3 | less
```
![7-1](../img/7-1.png)
图7-1。使用less搜索以“AGATCGG”开头的污染适配器序列；注意匹配后的核苷酸都非常相似

&emsp;&emsp;管道的一个有用的行为是，当LESS具有全屏幕数据时，将暂停将输出管道到LESS的程序的执行。这是由于管道在管道已满时阻止程序写入管道的方式。当您将程序的输出通过管道传送到Less并对其进行检查时，Less会停止从管道读取输入。很快，管道就会变满，并阻止将数据放入管道的程序继续。其结果是，我们可以在复杂的管道处理大数据之后抛出更少的东西，并且不用担心浪费计算能力-管道将会阻塞，我们可以根据需要花费尽可能多的时间来检查输出。

### 包含wc、ls和awk的纯文本数据摘要信息
      
&emsp;&emsp;除了使用head、tail或less查看文件外，我们可能还需要关于纯文本数据文件的其他摘要信息，如行数或列数。对于像制表符分隔的和CSV文件这样的纯文本数据格式，行数通常是行数。我们可以使用程序wc(用于字数)检索这一点：
```shell      
$ wc Mus_musculus.GRCm38.75_chr1.bed
 81226 243678 1698545 Mus_musculus.GRCm38.75_chr1.bed
``` 
 
&emsp;&emsp;默认情况下，wc输出所提供文件的字数、行数和字符数。它还可以处理许多文件：
```shell
 $ wc Mus_musculus.GRCm38.75_chr1.bed Mus_musculus.GRCm38.75_chr1.gtf
 81226 243678 1698545 Mus_musculus.GRCm38.75_chr1.bed
 81231 2385570 26607149 Mus_musculus.GRCm38.75_chr1.gtf
 162457 2629248 28305694 total
```
 
&emsp;&emsp;通常，我们只关心行数。我们可以使用选项-l只返回行数：
```shell
 $ wc -l Mus_musculus.GRCm38.75_chr1.bed
 81226 Mus_musculus.GRCm38.75_chr1.bed
 ```
         
&emsp;&emsp;对于这个染色体1鼠标注释，您可能已经注意到BED文件和GTF文件之间的差异。到底怎么回事？使用HEAD，我们可以检查Mus_musculus.GRCm38.75_chr1.gtf文件，并看到前几行是注释：
```shell
 $ head -n 5 Mus_musculus.GRCm38.75_chr1.gtf
!genome-build GRCm38.p2
!genome-version GRCm38
!genome-date 2012-01
!genome-build-accession NCBI:GCA_000001635.4
!genebuild-last-updated 2013-09
```
      
&emsp;&emsp;我们看到wc-l的五行差异是由于这个标题。使用散列标记(#)作为元数据的注释字段是一种常见的约定；这是我们在使用Unix数据工具时需要考虑的一个约定。
        
&emsp;&emsp;另一个我们通常想要的关于文件的信息是它的大小。要做到这一点，最简单的方法是与我们的老朋友ls一起使用-l选项：
```shell   
$ ls -l Mus_musculus.GRCm38.75_chr1.bed
-rw-r--r-- 1 vinceb staff 1698545 Jul 14 22:40 Mus_musculus.GRCm38.75_chr1.bed
```
     
&emsp;&emsp;在第四列(创建数据之前的列)中，ls-l报告文件大小(以字节为单位)。如果我们希望使用可读懂的大小，我们可以使用ls-lh：

```shell    
 $ ls -lh Mus_musculus.GRCm38.75_chr1.bed
-rw-r--r-- 1 vinceb staff 1.6M Jul 14 22:40 Mus_musculus.GRCm38.75_chr1.bed
```
       
&emsp;&emsp;这里，“M”表示兆字节；如果文件的大小是千兆字节，ls-lh将输出以千兆字节为单位的结果“G”。

​​![warning](../img/warning.png)数据格式和假设

    尽管wc-l是确定纯文本数据文件(例如TSV或CSV文件)中有多少行的快速方法，但它会假设您的数据格式良好。例如，假设脚本写入数据输出如下：

    $ cat some_data.bed
    1 3054233 3054733
    1 3054233 3054733
    1 3054233 3054733

    $ wc -l data.txt
    5 data.txt
         
    这里有一个微妙的问题：虽然只有三行数据，但有五行数据。这两个额外的行是文件末尾的空换行符。因此，尽管wc-l是计算文件中行数的一种快速而简单的方法，但它不是检查文件中有多少行数据的最健壮的方法。尽管如此，wc-l在大多数情况下仍然可以很好地工作，当我们只需要大致了解有多少行时。如果我们希望排除只有空格(空格、制表符或换行符)的行，可以使用grep：
 
    $ grep -c "[^ \\n\\t]" some_data.bed
    3

    稍后我们将讨论更多关于grep的内容。

&emsp;&emsp;关于一个文件，我们还经常需要另外一点信息：它包含多少列。我们总是可以用head-n1手动计算第一行的列数，但更简单的方法是使用awk。AWK是一种简单的小型编程语言，擅长处理TSV和CSV文件等文本数据。我们将在第157页的“使用awk进行文本处理”中更详细地介绍awk作为一种语言，但是让我们使用awk单行代码返回一个文件包含多少个字段：
```shell
 $ awk -F "\t" '{print NF; exit}' Mus_musculus.GRCm38.75_chr1.bed
3
```
       
&emsp;&emsp;AWK是为表格纯文本数据处理而设计的，因此有一个内置变量NF设置为当前数据集的字段数。这个简单的awk oneliner只打印Mus_musculus.GRCm38.75_chr1.bed文件第一行的字段数，然后退出。默认情况下，awk将空白(制表符和空格)视为字段分隔符，但我们可以通过设置awk的-F参数将其更改为仅制表符(因为我们正在使用的BED和GTF格式的示例都是制表符分隔的)。
       
&emsp;&emsp;查找mus_musculus.GRCm38.75_chr1.gtf中有多少列有点棘手。请记住，我们的mus_musculus.GRCm38.75_chr1.gtf文件前面有一系列注释：以散列符号(#)开头的五行，其中包含有用的元数据，如基因组构建、版本、日期和登录号。因为这个文件的第一行是注释，我们的awk技巧不起作用-它不报告数据列数，而是返回第一个注释的列数。要查看有多少列数据，我们需要首先删除注释，然后将结果传递给awk一行程序。一种方法是使用我们之前看到的尾部技巧：

```shell       
$ tail -n +5 Mus_musculus.GRCm38.75_chr1.gtf | head -n 1  (1)
!genebuild-last-updated 2013-09
$ tail -n +6 Mus_musculus.GRCm38.75_chr1.gtf | head       (2)
1 pseudogene gene 3054233 3054733 . + . [...]
$ tail -n +6 Mus_musculus.GRCm38.75_chr1.gtf | awk -F "\t" '{print NF; exit}'   (3)
```
(1)	 使用带-n+5参数的tail(注意前面的加号)，我们可以砍掉一些行。在将这些结果传送到awk之前，我们通过管道将其传送到head以检查我们所拥有的。实际上，我们发现我们犯了一个错误：返回的第一行是最后一条注释-我们需要再删掉一行。

(2)	将-n参数递增到6，并使用head检查结果，我们得到了想要的结果：标准输出流的第一行是mus_musculus.GRCm38.75_chr1.gtf GTF文件的第一行。

(3)	现在，我们可以通过管道将这些数据传输到awk单行代码中，以获取此文件中的列数。

&emsp;&emsp;虽然使用tail删除文件开头的注释标头块确实有效，但它并不是非常优雅，而且作为工程解决方案也有不足之处。随着你对计算越来越熟悉，你会认识到这样的解决方案是脆弱的。虽然我们已经设计了一个可以实现我们想要的解决方案，但它会在其他文件上起作用吗？它的健壮性吗？这样做好吗？这三个问题的答案都是否定的。认识到解决方案过于脆弱是发展Unix数据技能的重要组成部分。
       
&emsp;&emsp;使用tail-n+6从文件中删除注释的标题行的缺点是，此解决方案必须针对特定文件进行定制。这不是一个通用的解决方案，而从文件中删除注释行是一个通用的问题。使用tail需要计算出需要删除多少行，然后在Unix管道中对这个值进行硬编码。在这里，更好的解决方案是简单地排除与注释行模式匹配的所有行。使用程序grep(我们将在第140页的“全能grep”中详细讨论)，我们可以轻松地删除以“#”开头的行：

```shell
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | head -n 3
1 pseudogene gene 3054233 3054733 . + . [...]
1 unprocessed_pseudogene transcript 3054233 3054733 . + . [...]
1 unprocessed_pseudogene exon 3054233 3054733 . + . [...]
```
     
&emsp;&emsp;这个解决方案更快、更简单(因为我们不必计算有多少注释的标题行)，除了不那么脆弱和更健壮之外。总体而言，它是一个更好的工程解决方案-健壮性的最佳平衡，可通用化，并且能够快速实现。这些是您在使用Unix数据工具时应该寻找的解决方案类型：它们可以完成工作，既不过度设计，也不太脆弱。     

### 使用cut和Columns处理列数据
 
&emsp;&emsp;在处理纯文本表格数据格式(如制表符分隔的文件和CSV文件)时，我们通常需要从原始文件或流中提取特定的列。例如，假设我们只想提取Mus_musculus.GRCm38.75_chr1.bed文件的开始位置(第二列)。最简单的方法是使用CUT。这个程序从文本文件中剪切出指定的列(也称为字段)。默认情况下，CUT将制表符视为分隔符，因此要提取我们使用的第二列：
```shell
$ cut -f 2 Mus_musculus.GRCm38.75_chr1.bed | head -n 3
3054233
3054233
3054233
```
     
&emsp;&emsp;f参数是我们指定要保留哪些列的方式。参数-f还允许我们指定列的范围(例如-f 3-8)和列集(例如-f 3，5，8)。请注意，不可能使用CUT对列进行重新排序(例如，不幸的是，-f6，5，4，3将不起作用)。要对列进行重新排序，需要使用awk，稍后将对此进行讨论。
       
&emsp;&emsp;使用cut，我们可以将Mus_musculus.GRCm38.75_chr1.gtf的GTF转换为三列制表符分隔的基因组范围文件(例如染色体、起始和结束位置)。我们将使用前面介绍的grep命令剪切元数据行，然后使用cut提取第一列、第四列和第五列(染色体、开始、结束)：

```shell
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f1,4,5 | head -n 3
1 3054233 3054733
1 3054233 3054733
1 3054233 3054733
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f1,4,5 &gt; test.txt
```
      
&emsp;&emsp;请注意，尽管我们的基因组位置的三列文件看起来像BED格式的文件，但这并不是由于基因组范围格式的细微差异引起的。我们将在第9章中学习更多关于这方面的知识。
      
&emsp;&emsp;cut还允许我们指定列分隔符字符。因此，如果我们遇到一个包含染色体名称、起始位置和结束位置的CSV文件，我们也可以从中选择列：
```shell
$ head -n 3 Mus_musculus.GRCm38.75_chr1_bed.csv
1,3054233,3054733
1,3054233,3054733
1,3054233,3054733
$ cut -d, -f2,3 Mus_musculus.GRCm38.75_chr1_bed.csv | head -n 3
3054233,3054733
3054233,3054733
3054233,3054733
```

### 使用列设置表格数据的格式
       
&emsp;&emsp;正如您可能已经注意到的，在处理制表符分隔的文件时，不能很好的区分哪些元素属于特定的列。例如：
```shell
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f1-8 | head -n3
1 pseudogene gene 3054233 3054733 . + .
1 unprocessed_pseudogene transcript 3054233 3054733 . + .
1 unprocessed_pseudogene exon 3054233 3054733 . + .
```
       
&emsp;&emsp;虽然制表符在纯文本数据文件中是一个很好的分隔符，但我们的可变宽度数据导致列堆叠得不好。在unix：program column-t中对此进行了修复(-t选项告诉column将数据视为表)。column-t生成更易于阅读的整洁的列：
```shell
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f 1-8 | column -t
 | head -n 3
1 pseudogene gene 3054233 3054733 . + .
1 unprocessed_pseudogene transcript 3054233 3054733 . + .
1 unprocessed_pseudogene exon 3054233 3054733 . + .
```
         
&emsp;&emsp;请注意，您应该只使用column nt-t来可视化终端中的数据，而不是重新格式化数据以写入文件。制表符分隔的数据比由可变空格分隔的数据更可取，因为程序更容易解析。
         
&emsp;&emsp;与cut类似，Column的默认分隔符是制表符(\t)。我们可以使用-s选项指定不同的分隔符。因此，如果我们想更容易地可视化mus_musculus.GRCm38.75_chr1_bed.csv文件的列，我们可以使用：
```shell
$ column -s"," -t Mus_musculus.GRCm38.75_chr1_bed.csv | head -n 3
1 3054233 3054733
1 3054233 3054733
1 3054233 3054733   
```

&emsp;&emsp;column说明了一个关于我们应该如何处理数据的重要观点：没有理由以牺牲程序的可读性为代价来使数据格式具有吸引力。这与建议“为人类编写代码，为计算机编写数据”(第11页的“为人类编写代码，为计算机编写数据”)有关。尽管单字符分隔的列(如CSV或制表符分隔)可能难以让人阅读，但请考虑以下几点：  

1.	 它们可以立即与几乎所有Unix工具一起工作。   
2.	 使用column-t很容易将它们转换为可读的格式。

&emsp;&emsp;一般而言，使计算机可读的数据对人类具有吸引力，比使计算机可读的人类友好格式的数据更容易。不幸的是，在生物信息学中，将人类可读性置于计算机可读性之上的格式的数据仍然挥之不去。

### 全能的grep

&emsp;&emsp;在前面，我们已经看到grep是一个有用的工具，用于提取匹配(或不匹配)模式的文件行。grep-v允许我们以比tail更健壮的方式排除GTF文件的标题行。但是正如我们将在本节中看到的，这仅仅是grep功能的表面；grep是最强大的Unix数据工具之一。

&emsp;&emsp;首先，重要的是要提到grep是快速的。真的很快。如果需要在文件中找到模式(固定字符串或正则表达式)，grep将比用Python编写的任何内容都要快。图7-2显示了在文件中查找精确匹配行的四种方法的运行时：grep、sed、awk和一个简单的自定义Python脚本。如您所见，grep在这些基准中占主导地位：它比最快的替代方案Python快五倍。然而，这是有点不公平的比较：grep是快速的，因为它被调优为完成一项非常好的任务：查找与模式匹配的文件行。这个基准中包含的其他程序更通用，但是在这个特定任务的效率方面要付出代价。这证明了一个观点：如果计算速度是我们最优先考虑的事情(在很多情况下，它并不像我们想象的那么重要)，调整为执行某些任务的Unix工具实际上通常是最快的实现。

![7-2](../img/7-2.png)
图7-2、在玉米基因组中搜索精确字符串“AGATGCATG”所需时间的基准

&emsp;&emsp;虽然我们已经看到grep以前在本书中使用过，但让我们简要回顾一下它的基本用法。grep需要两个参数：模式(要搜索的字符串或基本正则表达式)和用于搜索它的文件(或多个文件)。作为一个非常简单的示例，让我们使用grep在文件mus_muscu-lus.GRCm38.75_chr1_genes.txt(包含染色体1上所有蛋白质编码基因的所有Ensembl基因标识符和基因名称)中查找基因“olfr418-ps1”：
```shell
$ grep "Olfr418-ps1" Mus_musculus.GRCm38.75_chr1_genes.txt
ENSMUSG00000049605 Olfr418-ps1
```

&emsp;&emsp;模式周围的引号不是必需的，但是使用引号是最安全的，这样我们的shell就不会试图解释任何符号。grep返回任何与模式匹配的行，即使是仅部分匹配的行：

```shell
$ grep Olfr Mus_musculus.GRCm38.75_chr1_genes.txt | head -n 5
ENSMUSG00000067064 Olfr1416
ENSMUSG00000057464 Olfr1415
ENSMUSG00000042849 Olfr1414
ENSMUSG00000058904 Olfr1413
ENSMUSG00000046300 Olfr1412
```

&emsp;&emsp;使用grep时一个有用的选项是-color=auto。此选项启用端子颜色，因此图案的匹配部分在端子中着色。
​​
![warning](../img/warning.png)GNU、BSD和grep的异同

    到目前为止，我们已经掩盖了一个非常重要的细节：Unix工具有不同的实现。grep、cut和sort等工具来自两种风格之一：BSD utils和GNU coreutil.这两种实现都包含我们在本章中使用的所有标准Unix工具，但它们的功能可能略有不同。BSD的工具可以在MAX OS X和其他Berkeley Software Distribution派生的操作系统(如FreeBSD)上找到。GNU的coreutils是Linux系统上的标准工具集。了解您正在使用的实现是很重要的(通过阅读手册页可以很容易地看出这一点)。如果你使用的是MacOSX并且想要使用GNU coreutils，你可以通过Homebrew和BREW install coreutils安装这些。每个程序都将安装前缀“g”(例如，cut将别名为gcut)，以便不干扰系统的默认工具。

    与BSD的utils不同，GNU的coreutils仍在积极开发中。GNU的coreutils也比BSD的utils具有更多的特性和扩展，我们在本章中将使用其中的一些功能和扩展。通常，我建议您使用GNU的coreutils而不是BSD utils，因为文档更全面，GNU扩展是有帮助的(有时是必要的)。在本章中，我将指出特定功能何时依赖于GNU版本。

&emsp;&emsp;前面，我们看到了如何使用grep只返回与指定的模式-这就是我们从GTF文件中排除注释行的方式。我们使用的选项是-v，用于反转。例如，假设您想要一个包含“olfr”的所有基因的列表，“olfr1413”除外。使用-v和链接一起调用带有管道的grep，我们可以使用：
```shell
$ grep Olfr Mus_musculus.GRCm38.75_chr1_genes.txt | grep -v Olfr1413
```

&emsp;&emsp;但要小心！这个可能会出什么问题？部分匹配可能会在这里咬我们：虽然我们想排除“Olfr1413”，但这个命令也会排除像“Olfr1413a”和“Olfr14130”这样的基因。但是我们可以通过使用-w来绕过这个问题，它匹配整个单词(由空格包围)。让我们用一个更简单的玩具例子来看看这是如何工作的：
```shell
$ cat example.txt
bio
bioinfo
bioinformatics
computational biology
$ grep -v bioinfo example.txt
bio
computational biology
$ grep -v -w bioinfo example.txt
bio
bioinformatics
computational biology
```

&emsp;&emsp;通过将匹配限制为单词，我们使用了一种更具限制性的模式。通常，我们的模式应该总是尽可能的限制性，以避免由部分匹配引起的无意匹配。

&emsp;&emsp;当我们需要通过肉眼检查结果时，grep的默认输出通常不能为我们提供足够的匹配上下文；只有匹配行被打印到标准输出。有三个有用的选项可以绕过这个上下文(-B)，上下文：之后(-A)，以及上下文前后(-C)。这些参数中的每一个都需要提供多少行上下文：

```shell
$ grep -B1 "AGATCGG" contam.fastq | head -n 6      (1)
@DJB775P1:248:D0MDGACXX:7:1202:12362:49613
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA
--
@DJB775P1:248:D0MDGACXX:7:1202:12782:49716
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
--
$ grep -A2 "AGATCGG" contam.fastq | head -n 6    (2)
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA
+
JJJJJIIJJJJJJHIHHHGHFFFFFFCEEEEEDBD?DDDDDDBDDDABDDCA
--
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
+
```
(1)在(-B)匹配行之前打印一行上下文。

(2)在(-A)匹配行之后打印两行上下文。

&emsp;&emsp;grep还支持一种称为POSIX基本正则表达式(BRE)的正则表达式。如果您熟悉Perl或Python中的正则表达式，您会注意到grep的正则表达式没有这些语言中的正则表达式那么强大。尽管如此，对于许多简单的应用程序来说，它们仍然工作得很好。例如，如果我们想要找到“Olfr1413”和“Olfr1411”的Ensembl基因标识符，我们可以使用：
```shell
$ grep "Olfr141[13]" Mus_musculus.GRCm38.75_chr1_genes.txt
ENSMUSG00000058904 Olfr1413
ENSMUSG00000062497 Olfr1411
```
         
&emsp;&emsp;在这里，我们在这两个基因名称之间使用一个共享前缀，并允许最后一个字符是'1'或'3'。然而，如果我们有更多不同的模式需要搜索，这种方法就不那么有用了。例如，构建匹配'Olfr218'和“Olfr1416”的BRE模式将是复杂的并且容易出错。对于这样的任务，使用grep对POSIX扩展正则表达式(ERE)的支持要容易得多。grep允许我们使用-E选项打开ere(在许多系统上，它的别名是egrep)。ERES允许我们使用交替(正则表达式行话，用于匹配几种可能的模式之一)来匹配“Olfr218”或“Olfr1416”。该方法的分隔符为竖线(|)：

```shell        
$ grep -E "(Olfr1413|Olfr1411)" Mus_musculus.GRCm38.75_chr1_genes.txt
ENSMUSG00000058904 Olfr1413
ENSMUSG00000062497 Olfr1411
```

&emsp;&emsp;我们现在只是触及了bre和ere的表面；我们没有空间在这里深入介绍这两种正则表达式风格(有关正则表达式的一些资源，请参见第16页上的“本书做出的假设”)。重要的是，你要认识到其中的差异，并知道必要的术语，以便在需要时寻求进一步的帮助。

&emsp;&emsp;grep有一个选项来计算与模式匹配的行数：-c。例如，假设我们想快速查看有多少基因以“olfr”开头：

```shell
$ grep -c "\tOlfr" Mus_musculus.GRCm38.75_chr1_genes.txt
27
```

&emsp;&emsp;或者，我们可以通过管道将匹配的行连接到wc-l：
```shell
$ grep "\tOlfr" Mus_musculus.GRCm38.75_chr1_genes.txt | wc -l
 27
```
        
&emsp;&emsp;对匹配的行数进行计数非常有用-尤其是对于行表示行的纯文本数据，并且计算与模式匹配的行数可用于对数据中的出现情况进行计数。例如，假设我们想知道Mus_musculus.GRCm38.75_chr1.gtf 文件中有多少snRNA。在此GTF文件的最后一列中，snRNA被注释为gene_biotype“snRNA”。计算这些特征的简单方法是：
```shell
 $ grep -c 'gene_biotype "snRNA"' Mus_musculus.GRCm38.75_chr1.gtf
315
```

&emsp;&emsp;注意这里我们是如何使用单引号来指定我们的模式的，因为我们的模式包括双引号字符(“)。

&emsp;&emsp;目前，grep正在输出整个匹配行。事实上，这就是grep如此快速的原因之一：一旦找到匹配，它就不会再费心搜索行的其余部分，只需将其发送到标准输出即可。然而，有时使用grep只提取模式的匹配部分是有用的。我们可以使用-o：
```shell
$ grep -o "Olfr.*" Mus_musculus.GRCm38.75_chr1_genes.txt | head -n 3
Olfr1416
Olfr1415
Olfr1414
```
&emsp;&emsp;或者，假设我们想从Mus_musculus.GRCm38.75_chr1.gtf文件的最后一列提取“gene_id”字段的所有值。使用-o很容易完成这项任务：
```shell
$ grep -E -o 'gene_id "\w+"' Mus_musculus.GRCm38.75_chr1.gtf | head -n 5
gene_id "ENSMUSG00000090025"
gene_id "ENSMUSG00000090025"
gene_id "ENSMUSG00000090025"
gene_id "ENSMUSG00000064842"
gene_id "ENSMUSG00000064842"
```

&emsp;&emsp;这里，我们使用扩展正则表达式来捕获字段中的所有基因名称。但是，正如您所看到的，存在大量冗余：我们的GTF文件具有多个特征(转录本、外显子、起始密码子等)。都有相同的基因名称。作为后面几节内容的一种尝试，Example 7-1展示了我们如何快速地将grep的这个凌乱的输出转换成一个独特的、排序的基因名称列表。

示例7-1、使用Unix数据工具清理一组基因名称。
```shell
$ grep -E -o 'gene_id "(\w+)"' Mus_musculus.GRCm38.75_chr1.gtf | \
 cut -f2 -d" " | \
 sed 's/"//g' | \
 sort | \
 uniq > mm_gene_id.txt
``` 
 
&emsp;&emsp;尽管看起来很复杂，但编写这篇文章只花了不到一分钟(还有其他可能的解决方案，省略cut，或者只使用awk)。这个文件的长度(根据wc-l)是2，027行-这是我们在Ensembl的BioMart数据库界面上点击获得相同信息时得到的相同的数字。在本章的其余部分中，我们将学习这些工具，以便您可以在工作中使用这种类型的快速管道。

### 解码纯文本数据：十六进制转储

&emsp;&emsp;在生物信息学中，我们处理的纯文本数据通常以ASCII编码。ASCII是一种字符编码方案，它使用7位来表示128个不同的值，包括字母(大写和小写)、数字和特殊的不可见字符。虽然ASCII只使用7位，但现在的计算机使用8位字节(表示8位的单位)来存储ASCII字符。有关ASCII的更多信息可通过man ascii在您的终端中获得。因为纯文本数据使用字符来编码信息，所以我们的编码方案很重要。当使用纯文本文件时，98%的时间不必担心ASCII的细节和文件的编码方式。然而，编码时2%的时间确实很重要-通常是在一个不可见的非ASCII字符输入数据时-它可能会导致重大的头痛。在本节中，我们将介绍在低级别检查文本数据以解决这些类型的问题的基础知识。如果您现在想跳过这一节，请将其添加到书签中，以防您在某个时候遇到这个问题。

&emsp;&emsp;首先，使用程序文件查看文件的编码，它从文件的内容推断编码是什么。例如，我们看到本章中使用的许多示例文件都是ASCII编码的：
```shell
$file Mus_musculus.GRCm38.75_chr1.bed Mus_musculus.GRCm38.75_chr1.gtf。
mus_musculus.GRCm38.75_chr1.bed：ASCII text
mus_musculus.GRCm38.75_chr1.gtf：ASCII text, with very long lines
```

&emsp;&emsp;某些文件将具有非ASCII编码方案，并且可能包含特殊字符。最常见的字符编码方案是UTF-8，它是ASCII的超集，但允许使用特殊字符。例如，本章的GitHub目录中包含的utf8.txt是一个UTF-8文件，从file的输出可以看出：

```shell
$ file utf8.txt
utf8.txt: UTF-8 Unicode English text
```

&emsp;&emsp;因为UTF-8是ASCII的超集，如果我们删除这个文件中的特殊字符并保存它，file将返回这个文件是ASCII编码的。

&emsp;&emsp;从Ensembl、NCBI和UCSC的Genome浏览器等数据源下载的大多数文件都不会有特殊字符，并且将采用ASCII编码(也就是没有这些特殊字符的UTF-8)。通常，我遇到的问题是来自人类生成的数据，通过复制和粘贴来自其他来源的数据可能会导致无意的特殊字符。例如，GitHub存储库中本章目录中的improper.fa文件在第一次检查时看起来就像常规的FASTA文件：

```shell
$ cat improper.fa
>good-sequence
AGCTAGCTACTAGCAGCTACTACGAGCATCTACGGCGCGATCTACG
>bad-sequence
GATCAGGCGACATCGAGCTATCACTACGAGCGAGΑGATCAGCTATT
```

&emsp;&emsp;然而，使用Bioawk找到这些序列的反向补码(暂时不要担心这个程序的细节-我们将在后面介绍它)会导致奇怪的结果：
```shell
$ bioawk -cfastx '{print revcomp($seq)}' improper.fa
CGTAGATCGCGCCGTAGATGCTCGTAGTAGCTGCTAGTAGCTAGCT
AATAGCTGATC
```

&emsp;&emsp;到底怎么回事？我们的第二个序列中有一个非ASCII字符：
```shell
$ file improper.fa
improper.fa: UTF-8 Unicode text
```
&emsp;&emsp;使用十六进制转储程序，我们可以确定是哪个字母导致了这个问题。十六进制转储程序返回每个字符的十六进制值。使用-c选项时，还会打印字符：
```shell
$ hexdump -c improper.fa
```

&emsp;&emsp;正如我们可以看到的，第二个序列中“CGAGCGAG”后面的字符显然不是ASCII字符。查看非ASCII字符的另一种方法是使用grep。这个命令有点棘手(它搜索十六进制范围之外的字符)，但是它是这样一个特定的用例，没有什么理由深入解释它：

```shell
$ LC_CTYPE=C grep --color='auto' -P "[\x80-\xFF]" improper.fa
GATCAGGCGACATCGAGCTATCACTACGAGCGAG[m�GATCAGCTATT
```

&emsp;&emsp;请注意，这不适用于BSD grep(Mac OS X附带的版本)。另一个有用的grep选项是-n，它将行号添加到每个匹配的行。在我的系统上，我的shell配置文件中有以下方便的行别名为nonascii(通常是~/.bashrc或~/.Profle)：
```shell
$ alias nonascii="LC_CTYPE=C grep --color='auto' -n -P '[\x80-\xFF]'"
```

&emsp;&emsp;总体而言，file、hexdump和grep命令在某些情况下非常有用，在这种情况下，您可能会怀疑文件的编码可能是原因(即使在准备本书的测试数据期间也会发生这种情况！)。这对于人工管理的数据尤其常见；总是要小心在没有检查的情况下将这些文件传递到分析管道中。

### 使用排序对纯文本数据进行排序

&emsp;&emsp;在生物信息学中，我们经常需要对纯文本数据进行排序处理，为什么需要对数据进行排序有一下两个最常见的原因：

·当对排序的数据执行某些操作时，效率要高得多。

·使用Unix sort|uniq成语对数据进行排序是查找所有唯一行的先决条件。

&emsp;&emsp;我们将在下一节中更多地讨论sort|uniq相关内容。这里我们重点讨论如何使用sort对数据进行排序。

&emsp;&emsp;首先，与CUT一样，SORT设计用于处理带有列的纯文本数据。在不带任何参数的情况下运行sort只需按行按字母数字方式对文件进行排序：

```shell
$ cat example.bed
chr1 26 39
chr1 32 47
chr3 11 28
Inspecting and Manipulating Text Data with Unix Tools | 147
chr1 40 49
chr3 16 27
chr1 9 28
chr2 35 54
chr1 10 19
$ sort example.bed
chr1 10 19
chr1 26 39
chr1 32 47
chr1 40 49
chr1 9 28
chr2 35 54
chr3 11 28
chr3 16 27
```
&emsp;&emsp;因为染色体是第一列，所以按线排序可以有效地将染色体组合在一起，因为这些是排序顺序中的“ties”。正如我们将看到的，分组数据对排序非常有用。

!['suggestion'](../img/suggestion.png)将不同的分隔符与排序一起使用
    默认情况下，排序将空白字符(如制表符或空格)视为字段分隔符。如果您的文件使用其他分隔符(例如CSV文件的逗号)，则可以使用-t(例如，-t“，”)指定字段分隔符。

&emsp;&emsp;但是，使用按行按字母数字排序的默认排序方式不能正确处理表格数据。我们需要两个新功能：

.按特定列进行排序。

·告诉sort某些列是数字值(而不是字母数字文本；有关差异的示例，请参阅第2章的提示“前导零和排序”)

&emsp;&emsp;sort有一个简单的语法来做到这一点。让我们看看如何按染色体(第一列)和起始位置(第二列)对example.bed进行排序：
```shell
$ sort -k1,1 -k2,2n example.bed
chr1 9 28
chr1 10 19
chr1 26 39
chr1 32 47
chr1 40 49
chr2 35 54
chr3 11 28
chr3 16 27
```

&emsp;&emsp;在这里，我们指定按-k参数的列排序(及其顺序)。在技术术语中，-k指定了排序键及其顺序。每个-k参数都将一系列列作为start，end，因此要按单个列排序，我们使用start。在前面的示例中，我们首先按第一列(染色体)排序，因为第一个-k参数是-k1，1。仅按第一列排序会导致具有相同染色体的行中存在许多关联(例如，“chr1”和“chr3”)。添加带有不同列的第二个-k参数告诉sort如何打破这些联系。在我们的示例中，-k2，2n告诉sort按第二列(开始位置)排序，将此列视为数值数据(因为-k2，2n中有n)。

&emsp;&emsp;最终结果是行按染色体分组，并按起始位置排序。然后我们可以将排序的标准输出流重定向到文件：
```shell
$ sort -k1,1 -k2,2n example.bed &gt; example_sorted.bed
```

&emsp;&emsp;如果需要对所有列进行数字排序，可以使用参数-n，而不是使用类似于-k2，2n的语法指定哪些特定列是数字列。

&emsp;&emsp;理解-k参数语法非常重要，我们将逐步介绍另一个示例。

&emsp;&emsp;Mus_musculus.GRCm38.75_ChR1_Random.gtf文件是带有置换行(没有元数据标头)的Mus_musculus.GRCm38.75_chr1.gtf文件。假设我们想再次按染色体分组行，并按位置排序。因为这是一个GTF文件，所以第一列是染色体，第四列是起始位置。因此，要对此文件进行排序，我们将使用：
```shell
$ sort -k1,1 -k4,4n Mus_musculus.GRCm38.75_chr1_random.gtf &gt; \
 Mus_musculus.GRCm38.75_chr1_sorted.gtf
 ```

​​![warning](../img/warning.png)排序稳定性

    关于排序有一个棘手的技术细节值得注意：排序稳定性。为了理解稳定的排序，我们需要返回并考虑如何处理具有相同排序关键字的行。如果根据我们指定的所有排序关键字，两行完全相同，那么它们在排序时是不可区分的和等价的。当行相等时，排序将根据整个行对它们进行排序，这是将它们按某种顺序排列的最后手段。这意味着即使这两行根据排序关键字是相同的，它们的排序顺序也可能与它们在原始文件中的出现顺序不同。此行为使排序成为不稳定的排序。

    如果我们不希望sort根据我们的排序键更改相等的行的顺序，我们可以指定-s选项。-s关闭这种最后的排序，从而使排序成为一种稳定的排序算法。
    
    排序可以是计算密集型的。与Unix工具不同，Unix工具一次只对一行进行操作，排序必须比较多行才能对文件进行排序。如果您有一个怀疑已经排序过的文件，那么验证它确实已经排序比重新排序要快速。我们可以使用-c：检查文件是否根据-k参数排序。

### 在uniq中查找唯一值
### 连接 
### 使用AWK进行文本处理
​​
!['suggestion'](../img/suggestion.png)

### Bioawk：生物格式的Awk
### 使用SED进行流编辑

!['suggestion'](../img/suggestion.png)

## 高级Shell技巧
### 子Shell
### 命名管道和进程替换

## Unix哲学再探

&emsp;&emsp; 在本章中，Unix哲学-配备了Unix管道的力量-允许我们使用一组丰富的Unix工具快速缝合小程序。Unix管道工作流不仅构建速度快，易于调试，而且通用性强，而且它们通常也是计算效率最高的解决方案。这证明了Unix令人难以置信的设计，我们处理现代生物信息学的很多方式都是由全能的Unix管道驱动的，这是一项40多年前由肯·汤普森(Ken Thompson)在“一个发烧的夜晚”(Doug McIlroy所描述的那样)发明的技术。