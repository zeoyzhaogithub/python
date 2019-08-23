## 第十章、使用序列数据

>One of the core issues of Bioinformatics is dealing with a profusion of (often poorly
defined or ambiguous) file formats. Some ad hoc simple human readable formats have
over time attained the status of de facto standards.
—Peter Cock et al. (2010)

>Good programmers know what to write. Great ones know what to rewrite (and reuse).
— Te Cathedral and the Bazaar Eric
S. Raymond

&emsp;&emsp;核苷酸(和蛋白质)序列以两种在生物信息学中广泛使用的纯文本格式存储：FASTA和FASTQ-分别发音为fast-ah(或fast-A)和fast-q。我们将在本节中讨论每种格式及其限制，然后查看一些用于处理这些格式的数据的工具。这是一个简短的章节，但有一个重要的教训：在使用特殊的生物信息学格式时，要小心常见的陷阱。像文件格式这样的小细节上的简单错误可能会消耗不成比例的时间和精力来发现和修复，所以在早期注意这些细节。

## FASTA格式

&emsp;&emsp;FASTA格式源自由William R.Pearson和David J.Lipman创建的FASTA对齐套件。FASTA格式用于存储不需要每个碱基对质量分数的任何类型的序列数据。这包括参考基因组文件、蛋白质序列、编码DNA序列(CDS)、转录本序列等。FASTA也可以用于存储多个比对数据，但我们不会在这里讨论这种格式的专用变体。我们在前面的章节中已经遇到了FASTA格式，但在本节中，我们将更详细地介绍这种格式，看看常见的陷阱，并介绍一些使用这种格式的工具。

&emsp;&emsp;FASTA文件由序列条目组成，每个条目包含两个部分：描述和序列数据。描述行以大于符号(>；)开头，并包含序列标识符和其他(可选)信息。序列数据从描述之后的下一行开始，并一直持续到有另一描述行(以>；开头的行)或文件结束。本章GitHub目录中的EGFR_fank.fasta文件是一个Fasta文件示例：
```shell
$ head -10 egfr_flank.fasta
>ENSMUSG00000020122|ENSMUST00000138518
CCCTCCTATCATGCTGTCAGTGTATCTCTAAATAGCACTCTCAACCCCCGTGAACTTGGT
TATTAAAAACATGCCCAAAGTCTGGGAGCCAGGGCTGCAGGGAAATACCACAGCCTCAGT
TCATCAAAACAGTTCATTGCCCAAAATGTTCTCAGCTGCAGCTTTCATGAGGTAACTCCA
GGGCCCACCTGTTCTCTGGT
>ENSMUSG00000020122|ENSMUST00000125984
GAGTCAGGTTGAAGCTGCCCTGAACACTACAGAGAAGAGAGGCCTTGGTGTCCTGTTGTC
TCCAGAACCCCAATATGTCTTGTGAAGGGCACACAACCCCTCAAAGGGGTGTCACTTCTT
CTGATCACTTTTGTTACTGTTTACTAACTGATCCTATGAATCACTGTGTCTTCTCAGAGG
CCGTGAACCACGTCTGCAAT
```

&emsp;&emsp;FASTA格式的简单性和灵活性带来了一个不幸的缺点：FASTA格式是一种定义松散的特别格式(不幸的是，它在生物信息学中非常常见)。因此，您可能会遇到FASTA格式的变体，这些变体可能会导致细微的错误，除非您的程序对这些变体很健壮。这就是为什么通常使用现有的FASTA/FASTQ解析库而不是实现您自己的解析库的原因；现有的库已经经过了开源社区的审查(稍后将对此进行更多介绍)。

&emsp;&emsp;关于FASTA格式最麻烦的是，描述中没有针对标识符格式的通用规范。例如，以下FASTA描述是否应引用相同的条目？
```shell
>ENSMUSG00000020122|ENSMUST00000138518
> ENSMUSG00000020122|ENSMUST00000125984
>ENSMUSG00000020122|ENSMUST00000125984|epidermal growth factor receptor
>ENSMUSG00000020122|ENSMUST00000125984|Egfr
>ENSMUSG00000020122|ENSMUST00000125984|11|ENSFM00410000138465
```
&emsp;&emsp;如果没有标识符的标准方案，我们就不能使用简单的精确匹配来检查标识符是否与FASTA条目标题行匹配。相反，我们需要依赖于FASTA描述和我们的标识符之间的模糊匹配。这可能很快就会变得非常混乱：我们的模式应该有多大的允许性？我们是否冒着将错误的序列与过于允许的正则表达式匹配的风险？从根本上说，模糊匹配是一种脆弱的策略。

&emsp;&emsp;幸运的是，这个问题有一个更好的解决方案(也非常简单)：不是依赖于后点模糊匹配来纠正不一致的命名，而是从严格的命名约定开始并保持一致。然后，通过几次健全的检查来运行来自外部来源的任何数据，以确保它遵循您的格式。这些检查不需要复杂(检查重复的名称，手动检查一些条目，检查>；和标识符之间的错误空格，检查不同文件之间的名称重叠，等等)。

&emsp;&emsp;如果需要整理外部数据，请始终保留原始文件，并编写一个脚本，将更正后的版本写入新文件。这样，脚本就可以很容易地在您收到的原始数据集的任何新版本上重新运行(但是您仍然需要检查一切-不要盲目地信任数据！)。

&emsp;&emsp;一种常见的命名约定是在第一个空格处将描述行分成两部分：标识符和注释。此格式的序列如下所示：
```shell
>gene_00284728 length=231;type=dna
GAGAACTGATTCTGTTACCGCAGGGCATTCGGATGTGCTAAGGTAGTAATCCATTATAAGTAACATGCGCGGAATATCCG
GAGGTCATAGTCGTAATGCATAATTATTCCCTCCCTCAGAAGGACTCCCTTGCGAGACGCCAATACCAAAGACTTTCGTA
GCTGGAACGATTGGACGGCCCAACCGGGGGGAGTCGGCTATACGTCTGATTGCTACGCCTGGACTTCTCTT
```
&emsp;&emsp;这里的gene_00284728是标识符，长度=231type=dna是注释。另外，ID应该是唯一的。虽然肯定不是标准，但将第一个空格之前的所有内容视为标识符、将第一个空格之后的所有内容视为非必要内容的惯例在生物信息学程序中很常见(例如，BEDtools、Samtools和BWA都这样做)。有了这个约定，通过标识符查找特定的序列是很容易的-我们将在本章的末尾看到如何使用索引的FASTA文件高效地做到这一点。

## FASTQ格式

&emsp;&emsp;FASTQ格式通过包括序列中每个碱基的数字质量分数来扩展FASTA。FASTQ格式被广泛用于存储高吞吐量排序数据，该数据与指示每个基本呼叫的置信度的每个基本质量分数一起报告。不幸的是，像FASTA一样，FASTQ也有变种和缺陷，这些变体和缺陷会使看似简单的格式在使用时变得令人沮丧。

FASTQ格式如下所示：
```shell
@DJB775P1:248:D0MDGACXX:7:1202:12362:49613 (1)
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA (2)
+(3)
JJJJJIIJJJJJJHIHHHGHFFFFFFCEEEEEDBD?DDDDDDBDDDABDDCA (4)
@DJB775P1:248:D0MDGACXX:7:1202:12782:49716
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
+
IIIIIIIIIIIIIIIHHHHHHFFFFFFEECCCCBCECCCCCCCCCCCCCCCC
```
(1)描述行，以@开头。这包含记录标识符和其他信息。

(2)序列数据，可以在一行或多行上。

(3)序列行后面以+开头的行表示序列的结束。在较旧的FASTQ文件中，在这里重复描述行是很常见的，但这是多余的，并导致不必要的大FASTQ文件。

(4)质量数据，也可以在一条或多条线上，但必须与序列长度相同。使用我们稍后讨论的方案(第344页的“基本质量”)，用ASCII字符对每个数字基本质量进行编码。

&emsp;&emsp;与FASTA一样，FASTQ文件中的一个常见约定是按第一个空格将描述行拆分为两个部分：记录标识符和注释。

&emsp;&emsp;FASTQ很难正确解析。一个常见的陷阱是将以@开头的每一行都视为描述行。但是，@也是有效的质量字符。FASTQ序列和质量行可以换行到下一行，因此以@开头的行有可能是质量行，而不是标题行。因此，编写总是将以@开头的行作为标题行的解析器可能会导致脆弱和不正确的解析。然而，我们可以使用质量分数字符数必须等于序列字符数这一事实来可靠地解析此格式-这就是后面介绍的readfq解析器的工作原理。
​​
<hr>

**计数FASTA/FASTQ条目的输入和输出**

&emsp;&emsp;作为纯文本格式，使用Unix工具可以轻松地使用FASTQ和FASTA。常见的命令行生物信息学习惯用法是：
```shell
$ grep -c "^>" egfr_flank.fasta
5
```
&emsp;&emsp;如第47页的“Pipes in Action：Creating Simple Programs with grep and Pipes”中所示，您必须将>字符引起来，以防止shell将其解释为重定向运算符(并覆盖您的FASTA文件！)。这是一种计算FASTA文件数量的安全方法，因为虽然格式是松散定义的，但每个序列都有一个online描述，并且只有这些行以>开头。

&emsp;&emsp;我们可能会尝试对FASTQ文件使用类似的方法，使用@而不是>；
```shell
$ grep -c "^@" untreated1_chr4.fq
208779
```

&emsp;&emsp;它告诉我们untreed1_chr4.fq有208，779个条目。但是通过仔细阅读unTreated1_chr4.fq，您会注意到每个FASTQ条目占用四行，但总行数是：
```shell
$ wc -l untreated1_chr4.fq
 817420 untreated1_chr4.fq
```
&emsp;&emsp;和817，420/4=204，355，这与grep-c给我们的结果完全不同！发生了什么？请记住，@是有效的质量字符，质量行可以此字符开头。可以使用grep“^@”untreated1_chr4.fq|less查看这方面的示例。

&emsp;&emsp;如果您绝对肯定，您的FASTQ文件每个序列条目使用四行，您可以通过使用wc-l估计行数并除以四来估计序列的数量。如果不确定某些FASTQ条目是否跨多行包装，一种更健壮的对序列进行计数的方法是使用Bioawk：

 $ bioawk -cfastx 'END{print NR}' untreated1_chr4.fq
204355
​​

## 核苷酸编码

&emsp;&emsp;在介绍了基本的FASTA/FASTQ格式之后，让我们看一下这些格式中编码核苷酸和基础质量分数的标准。显然，编码核苷酸很简单：A，T，C，G代表核苷酸腺嘌呤，胸腺嘧啶，胞嘧啶和鸟嘌呤。小写碱基通常用于指示掩蔽重复序列或低复杂性序列(由RepeatMasker和Tandem Repeats Finder等程序)。重复和低复杂性序列也可能是硬掩蔽的，其中核苷酸被N(或有时是X)替换。

&emsp;&emsp;简并(或二义性)核苷酸编码用于表示两个或多个碱基，例如，N用于表示任何碱基。国际纯粹和应用化学联合会(IUPAC)有一套标准化的核苷酸，既有明确的，也有歧义的(见表10-1)。

表10-1。IUPAC核苷酸编码
```shell
todo
```

&emsp;&emsp;一些生物信息学程序可能会以不同的方式处理含糊的核苷酸。例如，BWA读取比对工具将参考基因组中不明确的核苷酸字符转换为随机碱基(Li和Durbin，2009)，但使用随机种子集，因此两次重新生成比对索引不会导致两个不同的版本。

## 基本素质

&emsp;&emsp;FASTQ条目的每个序列基数在质量线中都有相应的数字质量分数。每个基本质量分数都被编码为单个ASCII字符。质量行看起来像一串随机字符，就像这里的第四行：

## 示例：检查和修整低质量碱基

&emsp;&emsp;注意在前面的例子中碱基的精确度是如何下降的；这是Illumina测序的特征误差分布。本质上，碱基不正确的概率越远(朝向3‘末端)，我们就会在由这种测序技术产生的读数中增加。这可能会对下游分析产生深远的影响！当使用排序数据时，您应该始终

·注意测序技术的误差分布和限制(例如，它是否受GC内容的影响)。

·考虑这可能会如何影响您的分析。
所有这些都是实验特有的，并且需要相当多的计划。

&emsp;&emsp;我们的Python基本准确度列表作为一个学习工具非常有用，可以用来了解如何将质量转换为概率，但它对我们理解数百万序列的质量概况没有多大帮助。从这个意义上说，一张图片胜过千言万语-而且有软件可以帮助我们在阅读中看到不同碱基之间的质量分布。最流行的是Java程序FastQC，它易于运行并输出有用的图形和质量指标。如果您更喜欢在R中工作，您可以使用名为qrqc的Bioconductor包(由您真正的编写)。我们将在示例中使用qrqc，以便我们可以自己修补如何可视化这些数据。

&emsp;&emsp;让我们首先安装本例所需的所有程序。首先，使用以下命令在R中安装qrqc：
```shell
>library(BiocInstaller)
> biocLite('qrqc')
```

&emsp;&emsp;接下来，让我们安装两个程序，这两个程序将允许我们修剪低质量的基础：镰刀和seqtk。seqtk是Heng Li编写的通用序列工具包，其中包含一个子命令，用于修剪序列末尾的低质量碱基(除了许多其他有用的函数)。镰刀和seqtk都可以轻松安装在带有Homebrew的Mac OS X上(例如，使用BREW安装seqtk和BREW安装镰刀)。

!['suggestion'](../img/suggestion.png)