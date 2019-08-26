## 第十二章、生物信息学shell脚本，编写管道和并行化任务脚本

&emsp;&emsp;我一直等到这本书的倒数第二章来分享一个令人遗憾的事实：日常的生物信息学工作往往涉及大量乏味的数据处理。生物信息学家经常需要对几十个(有时甚至数百个)文件运行一系列命令。因此，生物信息学的很大一部分是将各种处理步骤修补到一个管道中，然后将此管道重复应用于许多文件。这不是令人兴奋的科学工作，但在处理更令人兴奋的分析之前，这是一个必要的障碍。

&emsp;&emsp;虽然编写shell脚本是生物信息学家的日常负担，但重要的是要编写出健壮和可重现的shell脚本。shell脚本必须对数据处理过程中可能发生的问题具有健壮性。当我们直接在shell中对数据执行一系列命令时，我们通常会清楚地看到是否出现了错误-当输出文件应该包含数据或程序退出并出现错误时，输出文件是空的。但是，当我们通过处理管道运行数据时，我们牺牲了对每个步骤输出的仔细关注，以获得自动处理大量文件的能力。问题是，错误不仅可能仍然会发生，而且更有可能发生，因为我们正在对更多数据文件和使用更多步骤进行自动化处理。由于这些原因，构建健壮的管道是至关重要的。

&emsp;&emsp;同样，管道在可再现性方面也发挥着重要作用。精心设计的管道可以完美地记录数据是如何被处理的。在最好的情况下，个人可以下载您的处理脚本和数据，并轻松复制您的确切步骤。然而，不幸的是，很容易编写模糊或草率的管道，这会阻碍可重复性。我们将看到一些原则，可以帮助您避免这些错误，从而导致更多的可重复性项目。

&emsp;&emsp;在本章中，我们将学习构建健壮和可复制管道的基本工具和技能。我们将了解如何编写可重新运行的Bash shell脚本，使用find和xargs自动执行文件处理任务，并行运行管道，并看到一个简单的makefile。请注意，本章所涵盖的技术子集并不针对特定的集群或高性能计算(HPC)体系结构-这些是在任何机器上都能很好工作的通用Unix解决方案。对于特定于您的HPC系统的并行化技术，您将需要参考其文档。

## 碱基Bash脚本

>We’ve found that duct tape is not a perfect solution for anything. But with a little
ingenuity, in a pinch, it’s an adequate solution for just about everything.
— Mythbusters’ Jamie Hyneman

&emsp;&emsp;Bash，我们在整本书中交互使用的shell，也是一种成熟的脚本语言。像本书中介绍的许多其他工具一样，在生物信息学中有效使用Bash脚本的诀窍是知道何时使用它们，何时不使用它们。与Python不同，Bash不是一种通用语言。Bash被明确设计为使运行和连接命令行程序尽可能简单(shell的一个好特性！)。由于这些原因，Bash经常担当生物信息学的管道磁带语言(也称为胶水语言)的角色，因为它被用来将许多命令捆绑在一起，形成一个有凝聚力的工作流程。

&emsp;&emsp;在深入研究如何在Bash中创建管道之前，重要的是要注意，Python可能是一种更适合通常重用或高级管道的语言。Python是一种比Bash更现代、功能齐全的脚本语言。与Python相比，Bash缺少几个对数据处理脚本有用的优秀特性：更好的数值类型支持，有用的数据结构，更好的字符串处理，精炼的选项解析，大量库的可用性，以及帮助构建程序的强大功能。然而，与Bash相比，从Python脚本调用命令行程序(称为调出或shellout)会有更多的开销。尽管Bash缺少Python的一些特性，但Bash通常是最好和最快的“管道胶带”解决方案(我们在生物信息学中经常需要)。

## 编写和运行健壮的Bash脚本

&emsp;&emsp;生物信息学中的大多数Bash脚本只是组织成可重新运行的脚本的命令，添加了一些铃声和口哨来检查文件是否存在，并确保任何错误都会导致脚本中止。这些类型的Bash脚本编写起来非常简单：您已经学习了重要的shell特性，如管道、重定向和后台进程，它们在Bash脚本中扮演着重要的角色。在本节中，我们将介绍编写和执行Bash脚本的基础知识，特别关注如何创建健壮的Bash脚本。

## 健壮的Bash报头

&emsp;&emsp;按照惯例，Bash脚本的扩展名为.sh。您可以在您最喜欢的文本编辑器中创建它们(并且像样的文本编辑器将支持Bash脚本语法突出显示)。任何时候编写Bash脚本时，都应该使用以下Bash脚本头，它设置了一些Bash选项，可以生成更健壮的脚本(在GitHub上本章目录下的template.sh文件中也有这个头的副本)：

```shell
!/bin/bash    (1)
set -e           (2)
set -u           (3)
set -o pipefail (4)
```
(1)这称为shebang，它指示用于执行此脚本的解释器的路径。这只在将脚本作为程序运行时才是必要的(稍后将对此进行更多介绍)。无论您计划如何运行Bash脚本，最好包含一个shebang行。
(2)默认情况下，包含失败命令(以非零退出状态退出)的shell脚本不会导致整个shell脚本退出-shell脚本只会继续到下一行。这不是一个可取的行为；我们总是希望错误是响亮的和值得注意的。set-e通过在任何命令以非零退出状态退出时终止脚本来防止这种情况。但是请注意，set-e具有复杂的规则来适应非零退出状态指示“false”而不是失败的情况。例如，如果test-d file.txt的参数不是目录，则它将返回非零的退出状态，但在此上下文中，这并不意味着表示错误。由于这个原因，set-e忽略if条件中的非零状态(我们将在后面讨论)。此外，set-e会忽略Unix管道中除最后一个以外的所有退出状态-这与set-o pipefailure有关，我们稍后将对此进行讨论。

(3)set-u修复了Bash脚本的另一个不幸的默认行为：任何包含对未设置变量名的引用的命令仍将运行。作为这可能导致的可怕示例，请考虑：rm-rf$TEMP_DIR/*。如果没有设置shell变量$TEMP_DIR，Bash仍然会用它的值(没有什么)代替它。最终结果是RM-RF/*！你可以亲眼看到这个：
```shell
$ echo "rm $NOTSET/blah"
rm /blah
```
set-u通过在变量的值未设置时中止脚本来防止此类错误。

(4)	如前所述，如果遇到非零退出状态，set-e将导致脚本中止，但有一些例外。例如，如果在Unix管道中运行的程序未成功退出；除非该程序是管道中的最后一个程序，否则即使使用set-e，这也不会导致脚本中止。包括set-o管道失败将防止这种不良行为-任何在管道中返回非零退出状态的程序都将导致整个管道返回非零状态。在启用set-e的情况下，这将导致脚本中止。

!['note'](../img/note.png)本章Bash脚本示例中的健壮Bash标头

    为了清晰和节省空间，我将在本章的Bash脚本中省略这个标题，但您应该始终在自己的工作中使用它

    这三个选项是针对具有静默错误和不安全行为的Bash脚本的第一层保护。不幸的是，Bash是一种脆弱的语言，我们需要注意其他一些奇怪的地方，以便在生物信息学中安全地使用它。我们将在学习更多关于这门语言的知识时看到这些。

## 运行Bash脚本

&emsp;&emsp;可以通过以下两种方式之一运行Bash脚本：直接使用bash程序(例如，bash script.sh)，或将脚本作为程序调用(./script.sh)。就我们的目的而言，没有技术上的理由选择一种方法而不是另一种方法。实际上，养成使用./script.sh运行收到的脚本的习惯是明智的，因为他们可能使用/bin/bash以外的解释器(例如zsh、csh等)。但是，虽然我们可以使用bash script.sh运行任何脚本(只要它具有读取权限)，但将脚本作为可执行文件调用需要它具有可执行权限。我们可以使用以下方法设置这些参数：
```shell
$ chmod u+x script.sh
```

&emsp;&emsp;这将为拥有文件(U)的用户添加可执行权限(+x)。然后，可以使用./script.sh运行脚本。

## 变量和命令参数

&emsp;&emsp;Bash变量在健壮、可重复的Bash脚本中扮演着极其重要的角色。处理具有大量设置的管道，这些设置应存储在变量中(例如，存储结果的目录、命令的参数值、输入文件等)。将这些设置存储在文件顶部定义的变量中，可使调整设置和重新运行管道变得更加容易。使用变量来存储设置意味着您只需要更改一个值-您为变量分配的值，而不是必须在脚本中更改大量的硬编码值。Bash还将命令行参数读入变量，因此您需要熟悉访问变量的值才能使用命令行参数。
&emsp;&emsp;与其他编程语言不同，Bash的变量没有数据类型。将Bash的变量视为字符串是有帮助的(但根据上下文的不同可能会有不同的行为)。我们可以使用创建变量并为其赋值(请注意，在设置Bash变量时，空格很重要-请勿使用等号周围的空格！)：

```shell
results_dir="results/"
```

&emsp;&emsp;要访问变量的值，我们在变量名称前使用美元符号(例如，$results_dir)。您可以在Bash脚本中进行试验，也可以直接在命令行上进行试验：
```shell
$ results_dir="results/"
$ echo $results_dir
results/
```
&emsp;&emsp;正如上一节中提到的，您应该始终设置-u，以便在没有设置变量的情况下强制Bash脚本退出。
     
&emsp;&emsp;尽管使用美元符号语法访问变量的值是可行的，但它有一个缺点：在某些情况下，不清楚变量名在哪里结束，相邻字符串从哪里开始。例如，假设Bash脚本的一部分为样本的对齐数据创建了一个目录，名为<sample>_aln/，其中<sample>被样本的名称替换。这将如下所示：
```shell
sample="CNTRL01A"
mkdir $sample_aln/
```

&emsp;&emsp;虽然这个代码块的目的是创建一个名为CNTRL01A_ALN/的目录，但这实际上会失败，因为Bash将尝试检索名为$SAMPLE_ALN的变量的值。要防止出现这种情况，请将变量名放在大括号中：
```shell
sample=“CNTRL01A”
mkdir${sample}_aln/
```
&emsp;&emsp;现在，将创建名为CNTRL01A_ALN/的目录。虽然这解决了Bash将SAMPLE_ALN解释为变量名的直接问题，但为了使其更加健壮，我们还应该采取另一个步骤：引用变量。这可以防止命令解释变量可能包含的任何空格或其他特殊字符。我们的最后一个命令如下所示：
```shell
sample=“CNTRL01A”
mkdir“${sample}_aln/”
```

### 命令行参数

&emsp;&emsp;现在让我们看一下Bash如何处理命令行参数(分配给值$1、$2、$3等)。变量$0存储脚本的名称。我们可以通过一个简单的示例脚本看到这一点：
```shell
!/bin/bash
Basic Bash Scripting | 399
echo "script name: $0"
echo "first arg: $1"
echo "second arg: $2"
echo "third arg: $3"
```
&emsp;&emsp;运行此文件将打印分配给$0、$1、$2和$3的参数：
```shell
$ bash args.sh arg1 arg2 arg3
script name: args.sh
first arg: arg1
second arg: arg2
third arg: arg3
```
&emsp;&emsp;Bash将命令行参数的数量分配给变量$#(这不会将脚本名称$0作为参数计算)。这对于用户友好的消息很有用(这使用了Bash if Conditional，我们将在下一节中更深入地讨论)：
```shell
!/bin/bash
if [ "$#" -lt 3 ] # are there less than 3 arguments?
then
 echo "error: too few arguments, you provided $#, 3 required"
 echo "usage: script.sh arg1 arg2 arg3"
 exit 1
fi
echo "script name: $0"
echo "first arg: $1"
echo "second arg: $2"
echo "third arg: $3"
```

&emsp;&emsp;使用太少的参数运行此命令会产生错误(并导致进程以非零的退出状态退出-如果您对退出状态的含义感到生疏，请参阅第52页的“退出状态：如何以编程方式告诉您的命令是否工作”)：
```shell
$ ./script.sh some_arg
error: too few arguments, you provided 1, 3 required
usage: script.sh arg1 arg2 arg3
```
&emsp;&emsp;可以使用unix工具getopt进行更复杂的选项和参数解析。这超出了本书的范围，但是getopt的手动输入非常彻底。然而，如果您发现您的脚本需要大量或复杂的选项，那么使用Python而不是Bash可能会更容易。Python的argparse模块比getopt更易于使用。

!['suggestion'](../img/suggestion.png)再现性与环境变量

    一些生物信息学家利用环境变量来使用命令export来存储设置，但通常这会降低脚本的可移植性和可重复性。相反，所有重要的设置都应该作为变量存储在脚本中，而不是作为外部环境变量。这样，脚本是独立的和可重现的。

&emsp;&emsp;在Bash脚本中创建的变量仅在运行该脚本的Bash进程期间可用。例如，运行一个创建变量的脚本，其中某些_var=3将不会在您当前的shell中创建一些_var，因为该脚本在一个完全独立的shell进程中运行。

### Bash脚本中的条件：if语句

&emsp;&emsp;与其他脚本语言一样，Bash支持标准的if条件语句。使Bash有点独特的是，命令的退出状态提供了true和false(请记住：与其他语言相反，0代表true/successful，而其他任何东西都是false/failure)。基本语法为：
```bash
if [commands]    (1)
then
 [if-statements]   (2)
else
 [else-statements] (3)
fi
```

(1)[Commands]是任何命令、命令集、管道或测试条件(稍后将看到)的占位符。如果这些命令的退出状态为0，则之后继续执行该块；否则，在其他情况下继续执行该块。Then关键字可以放在相同的行上，就像if一样，但是随后需要一个分号：if[Commands]；Then。

(2)如果[Commands]的计算结果为TRUE(0)，则[IF-Statements]是执行的所有语句的占位符。

(3)如果[Commands]的计算结果为False(1)，则[Else-Statements]是执行的所有语句的占位符。else块是可选的。

&emsp;&emsp;与Python等语言相比，Bash的if条件语句可能看起来有点奇怪，但请记住：Bash主要是为了将其他命令缝合在一起而设计的。这是Bash在编写管道时相对于Python的优势：Bash允许您的脚本直接使用命令行程序，而不需要任何调用程序的开销。虽然用Bash编写复杂的程序可能是不愉快的，但是编写简单的程序是非常容易的，因为Unix工具和Bash很好地协调。例如，假设我们希望仅当文件包含某个字符串时才运行一组命令。因为grep仅在与文件中的模式匹配时返回0，否则返回1，所以我们可以使用：
```bash
!/bin/bash
if grep "pattern" some_file.txt &gt; /dev/null   (1)
then
 # commands to run if "pattern" is found
 echo "found 'pattern' in 'some_file.txt"
fi
```
(1)	这个grep命令是我们的条件语句。重定向是整理此脚本的输出，以便grep的输出重定向到/dev/null，而不是脚本的标准输出。

&emsp;&emsp;if条件中的一组命令可以使用到目前为止我们已经掌握的Unix的所有特性。例如，使用&(逻辑AND)和|(逻辑OR)等逻辑运算符链接命令：
```bash
!/bin/bash
if grep "pattern" file_1.txt &gt; /dev/null &&
 grep "pattern" file_2.txt &gt; /dev/null
then
 echo "found 'pattern' in 'file_1.txt' and in 'file_2.txt'"
fi
```
&emsp;&emsp;我们还可以使用！：否定程序的退出状态。
```bash
!/bin/bash
if ! grep "pattern" some_file.txt &gt; /dev/null
then
 echo "did not find 'pattern' in 'some_file.txt"
fi
```

&emsp;&emsp;最后，可以在if条件语句中使用管道。然而，请注意，行为取决于set-o pipefailure。如果设置了pipefailure，则条件语句中管道中的任何非零退出状态都将导致继续执行，跳过if-Statements部分(如果存在，则转到else块)。但是，如果未设置pipefailure，则只考虑最后一个命令的退出状态。而不是试图记住所有这些规则，只需使用前面提供的健壮头部-pipefailure是更明智的默认设置。

&emsp;&emsp;理解Bash的if语句所需的最后一个组件是test命令。与其他程序一样，test退出时要么为0，要么为1。然而，test的退出状态指示通过其参数指定的测试的返回值，而不是退出成功或错误。test支持许多标准比较运算符(两个字符串是否相等，两个整数是否相等，一个整数是否大于或等于另一个整数，等等)。Bash不能依赖熟悉的语法，如>表示“大于”，因为这是用于重定向的：相反，test有自己的语法(有关完整列表，请参见表12-1)。通过直接在命令行(使用echo“$？”打印退出状态)：
```bash
test "ATG" = "ATG" ; echo "$?"
0
$ test "ATG" = "atg" ; echo "$?"
1
$ test 3 -lt 1 ; echo "$?"
1
$ test 3 -le 3 ; echo "$?"
0
```
表12-1。字符串和整数比较运算符
```python
 todo
```

&emsp;&emsp;在实践中，您要检查的最常见的条件不是某个整数是否小于另一个整数，而是检查文件或目录是否存在以及您是否可以写入它们。test支持许多与文件和目录相关的测试操作(在生物信息学中最有用的几个在表12-2中)。让我们看几个基本的命令行示例：

```shell
$ test -d some_directory ; echo $? # is this a directory?
0
$ test -f some_file.txt ; echo $? # is this a file?
0
$ test -r some_file.txt ; echo $? $ is this file readable?
Basic Bash Scripting | 403
0
$ test -w some_file.txt ; echo $? $ is this file writable?
1
```
表12-2。文件和目录测试表达式

``` python
todo
```

&emsp;&emsp;将test与if语句结合使用很简单；test是一个命令，因此我们可以使用：
```bash
if test -f some_file.txt
then
 [...]
fi
```
&emsp;&emsp;然而，Bash为test语句提供了一种更简单的语法替代：[-f some_file.txt]。注意括号周围和括号内的空格-这些是必需的。这使得涉及比较的IF语句更加简单：
```bash
if [ -f some_file.txt ]
then
 [...]
fi
```
&emsp;&emsp;当使用这种语法时，我们可以用-a作为逻辑，-o作为逻辑OR，！作为否定，并将括号用于分组语句。我们熟悉的&和|运算符在测试中不起作用，因为这些是shell运算符。例如，假设我们想要确保脚本有足够的参数，并且输入文件是可读的：
```bash
!/bin/bash
set -e
set -u
set -o pipefail
if [ "$#" -ne 1 -o ! -r "$1" ]
then
 echo "usage: script.sh file_in.txt"
 exit 1
fi
```

&emsp;&emsp;如前所述，我们引用变量(特别是来自人工输入的变量)；这是一种很好的做法，可以防止出现具有特殊字符的问题。

&emsp;&emsp;当与-a或-e链接在一起时，test的语法使用短路评估。这意味着test只会根据需要计算尽可能多的表达式，以确定整个语句是真还是假。在本例中，如果没有提供确切的一个参数，test将不会检查文件参数$1是否可读(第一个条件为true)。这两个表达式与逻辑OR组合在一起，逻辑OR只需要一个表达式为TRUE，整个条件就为TRUE。

## 使用for循环和Globbing处理Bash文件

&emsp;&emsp;在生物信息学中，我们的大多数数据被分成多个文件(例如，不同的处理、复制、基因型、物种等)。任何处理管道的核心都是将相同的工作流应用于这些文件中的每个文件的方法，注意跟踪样本名。用Bash的for循环遍历文件是实现这一点的最简单的方法。这是生物信息学处理管道的一个非常重要的部分，我们将在本章的下一节中介绍其他有用的工具和方法。

&emsp;&emsp;创建管道以处理一组文件有三个基本部分：

·选择要应用命令的文件。

·在数据上循环并应用命令。

·跟踪创建的任何输出文件的名称。

&emsp;&emsp;有不同的计算技巧来完成这些任务中的每一个。让我们首先看一下选择要应用命令的文件的简单方法。

&emsp;&emsp;有两种常见的方法来选择要应用生物信息学工作流的文件：以包含有关样本的信息(其样本名称、文件路径等)的文本文件开始的方法，以及使用某些标准选择目录中的文件的方法。任何一种方法都是好的-它主要归结为对给定任务最有效的方法。我们将首先研究一种从示例名称开始的方法，然后返回到如何查找文件。假设您有一个名为samples.txt的文件，它告诉您关于原始数据的基本信息：示例名称、读取对和文件的位置。下面是一个示例(也在本章GitHub上的目录中)：
```bash
$ cat samples.txt
zmaysA R1 seq/zmaysA_R1.fastq
zmaysA R2 seq/zmaysA_R2.fastq
zmaysB R1 seq/zmaysB_R1.fastq
zmaysB R2 seq/zmaysB_R2.fastq
zmaysC R1 seq/zmaysC_R1.fastq
zmaysC R2 seq/zmaysC_R2.fastq
```

&emsp;&emsp;第一列提供示例名称，第二列包含读取对，最后一列包含此示例/读取对组合的FASTQ文件的路径。前两列称为元数据(关于数据的数据)，这对于将样本信息与其物理文件相关联至关重要。请注意，元数据也在文件名本身中，这很有用，因为它允许我们在需要时从文件名中提取元数据。

&emsp;&emsp;有了这个samples.txt文件，创建管道的第一步就完成了：关于要处理的文件的所有信息，包括它们的路径，都是可用的。第二步和第三步是对这些数据进行循环，并以一种保持样本平直的方式进行。我们如何做到这一点，取决于具体的任务。如果您的命令接受单个文件并返回单个文件，则解决方案很简单：文件是我们正在处理的单元。我们只需对每个文件进行循环，并使用该文件名称的修改版本进行输出。

&emsp;&emsp;让我们来看一个例子：假设我们想要遍历每个文件，收集每个文件的质量统计信息(使用假想的程序FASTQ_STAT)，并将这些信息保存到输出文件中。每个输出文件都应该有一个基于输入文件的名称，所以如果摘要文件表明有错误，我们就可以知道哪个文件受到了影响。这里有很多小部分，所以我们将一步地逐步完成这一步，一次一块地学习Bash数组、Basename和其他几个shell技巧，在此过程中，我们将学习一些关于Bash数组、basename和其他几个shell技巧的知识。

&emsp;&emsp;首先，我们将文件名加载到Bash数组中，然后可以对其进行循环。可以使用以下方法手动创建Bash阵列：
```shell
$ sample_names=(zmaysA zmaysB zmaysC)
```
&emsp;&emsp;并且可以使用提取特定元素(注意Bash数组是0索引的)：
```shell
$ echo ${sample_names[0]}
zmaysA
$ echo ${sample_names[2]}
zmaysC
```
&emsp;&emsp;所有元素都是用看起来很神秘的${sample_files[@]}提取的：
```shell
$ echo ${sample_names[@]}
zmaysA zmaysB zmaysC
```

&emsp;&emsp;还可以使用以下内容访问数组中有多少个元素(以及每个元素的索引)：
```shell
$ echo ${#sample_names[@]}
3
$ echo ${!sample_names[@]}
0 1 2
```

&emsp;&emsp;但是手工创建Bash数组非常繁琐并且容易出错，特别是因为在sample.txt文件中已经有了我们的文件名。Bash的美妙之处在于我们可以使用命令替换(在第54页的“命令替换”中讨论)来构造Bash数组(尽管这可能是危险的；请参见下面的警告)。因为我们想要循环遍历每个文件，所以我们需要使用cut-f 3从samples.txt中提取第三列。在shell中演示这一点：
```shell
$ sample_files=($(cut -f 3 samples.txt))
$ echo ${sample_files[@]}
seq/zmaysA_R1.fastq seq/zmaysA_R2.fastq seq/zmaysB_R1.fastq
seq/zmaysB_R2.fastq seq/zmaysC_R1.fastq seq/zmaysC_R2.fastq
```
&emsp;&emsp;同样，这只在您可以对您的文件名做出强有力的假设时才有效-即它们只包含字母数字字符、(_)和(-)！如果空格、制表符、换行符或*等特殊字符最终出现在文件名中，则会破坏这种方法。
!['warning'](../img/warning.png)内部字段分隔符、单词拆分和文件名

    当使用sample_files=($(cut-f3 samples.txt)通过命令替换创建Bash数组时，Bash使用Word Split通过内部字段分隔符(IFS)中的字符将字段拆分为数组元素。内部字段分隔符存储在Bash变量IFS中，默认情况下包括空格、制表符和换行符。您可以使用以下命令检查IFS的值：

    ```shell
    $ printf %q "$IFS"
    $' \t\n'
    ```
    请注意，空格包含在IFS中(第一个字符)。当文件名包含空格时，这可能会引入问题，因为Bash将在空格字符上拆分，将文件名分成几个部分。同样，避免问题的最佳方法是不在文件名中使用空格、制表符、换行符或特殊字符(例如，*)-仅使用字母数字字符、(-)和(_)。本节中教授的技术假定文件已正确命名，并且对不正确命名的文件没有健壮性。如果文件名中存在空格，则可以将IFS的值设置为仅使用制表符和换行符；有关此主题的其他详细信息，请参阅本章GitHub上的自述文件。

&emsp;&emsp;有了Bash数组中的文件名，我们几乎可以循环访问它们了。最后一个组件是剥离每个文件名的路径和扩展名，留给我们可以用来创建输出文件名的最基本的文件名。Unix程序基本名称从文件名中剥离路径：
```shell
$ basename seqs/zmaysA_R1.fastq
zmaysA_R1.fastq
```
&emsp;&emsp;Basename还可以从文件名中剥离作为第二个参数提供的后缀(例如，扩展名)(或者使用参数-s)：
```shell
$ basename seqs/zmaysA_R1.fastq .fastq
zmaysA_R1
$ basename -s .fastq seqs/zmaysA_R1.fastq
zmaysA_R1
```

&emsp;&emsp;我们使用basename返回每个文件名的基本部分，然后使用它为FASTQ_STAT的结果创建输出文件名。

&emsp;&emsp;现在，所有的部分都准备好构建我们的处理脚本：
```shell
#!/bin/bash
set -e
set -u
set -o pipefail
#specify the input samples file, where the third
#column is the path to each sample FASTQ file
sample_info=samples.txt

#create a Bash array from the third column of $sample_info
sample_files=($(cut -f 3 "$sample_info"))    (1)
for fastq_file in ${sample_files[@]}     (2)
do
 # strip .fastq from each FASTQ file, and add suffix
 # "-stats.txt" to create an output filename for each FASTQ file
 results_file="$(basename $fastq_file .fastq)-stats.txt"        (3)
 # run fastq_stat on a file, writing results to the filename we've
 # above
 fastq_stat $fastq_file &gt; stats/$results_file     (4)
done
```
(1)此行使用命令替换创建包含所有FASTQ文件的Bash数组。请注意，这将使用变量$sample_info中包含的文件名，如果管道要在不同的样本集上运行，则稍后可以轻松更改该文件名。

(2)接下来，我们使用Bash for循环遍历每个示例文件名。表达式${sample_files[@]}返回Bash数组中的所有项。

(3)这一重要的行将使用输入文件的文件名创建输出文件名。命令替换$(basename $fastq_file.fast q)采用存储在$fastq_file中的迭代的当前文件名，并去掉.fast扩展名。剩下的是示例名称中标识原始文件的部分，后缀-stats.txt被添加到该文件中。

(4)最后，运行命令FASTQ_STAT，使用当前迭代的文件名作为输入，并将结果写入stats/$results_file。

&emsp;&emsp;非那样做不行。更精炼的脚本可能会添加一些额外的功能，例如，如果FASTQ文件不存在，则使用if语句提供友好的错误，或者调用echo来报告当前正在处理的样本。

&emsp;&emsp;这个脚本很容易编写，因为我们的处理步骤将单个文件作为输入，并创建单个文件作为输出。在这种情况下，只需为每个文件名添加一个后缀就足以使我们的示例保持直截了当。然而，许多生物信息学管道将两个或多个输入文件组合成单个输出文件。对齐成对端读是一个主要的例子：大多数对齐器获取两个输入FASTQ文件并返回一个输出对齐文件。当编写脚本来对齐成对的末端读取时，我们不能像之前那样循环遍历每个文件。相反，每个样本，而不是每个文件，是处理单元。我们的对齐步骤为每个示例获取两个FASTQ文件，并将其转换为此示例的单个对齐文件。因此，我们的循环必须迭代唯一的样本名称，并且我们使用这些样本名称重新创建用于对齐的输入FASTQ文件。这在一个例子中会更清楚；假设我们使用校准器BWA，我们的基因组参考命名为zmays_AGPv3.20.fa：
```shell
#!/bin/bash
set -e
set -u
set -o pipefail
#specify the input samples file, where the third
#column is the path to each sample FASTQ file
sample_info=samples.txt
#our reference
reference=zmays_AGPv3.20.fa
#create a Bash array from the first column, which are
#sample names. Because there are duplicate sample names
#(one for each read pair), we call uniq
sample_names=($(cut -f 1 "$sample_info" | uniq))    (1)
for sample in ${sample_names[@]}
do
 # create an output file from the sample name
 results_file="${sample}.sam"                                               (2)
 bwa mem $reference ${sample}_R1.fastq ${sample}_R2.fastq \    (3)
 &gt; $results_file
done
```
(1)这与前面的示例非常相似，除了现在我们使用CUT来获取第一列(对应于样本名称)，并且(最重要的)将这些样本名称通过管道传输到uniq，以便删除相同样本名称的重复项。这是必要的，因为我们的第一列将每个样本名称重复两次，每个配对结束文件一次。

(2)与前面一样，我们为正在迭代的当前示例创建输出文件名。在本例中，所需的只是存储在$sample中的示例名称。

(3)我们对bwa的调用提供了引用，并为这个示例提供了两个成对结束的FASTQ文件作为输入。注意我们如何为给定的样本名称重新创建两个输入FASTQ文件，因为这些文件在样本之间的命名是一致的。在实践中，这对于大量的生物信息学数据是可能的，这些数据通常来自于一致命名文件的机器。最后，bwa的输出被重定向到$results_file。为了清楚起见，我在这个命令中省略了引号变量，但您可能希望添加这一点。

&emsp;&emsp;最后，在某些情况下，直接遍历文件可能比处理包含samples.txt这样的样本信息的文件更容易。做到这一点最简单(也是最安全)的方法是使用Bash的通配符来遍历要循环的文件(回想一下我们在第26页的“组织数据以自动执行文件处理任务”中介绍了globbing)。它的语法非常简单：
```shell
#!/bin/bash
set -e
set -u
set -o pipefail
for file in *.fastq
do
 echo "$file: " $(bioawk -c fastx 'END {print NR}' $file)
done
```

&emsp;&emsp;这个简单的脚本使用Bioawk来计算一个文件中有多少个FASTQ记录，对于当前目录中以.fastq结尾的每个文件。

&emsp;&emsp;Bash的循环是将命令应用于大量文件的便捷方式，但也有一些缺点。首先，与Unix工具find(我们将在下一节中看到)相比，globbing不是选择某些文件的非常强大的方法。其次，对于简单的操作，Bash的循环语法很长，并且有点过时。最后，没有一种简单的方法可以约束所使用的子进程的数量来并行化Bash循环。我们将在下一节中看到一个强大的文件处理Unix习惯用法，它可以更好地处理一些Bash脚本可能不是最优的任务。

### 使用find和xargs自动化文件处理

&emsp;&emsp;在本节中，我们将学习一种使用Unix find指定符合某些条件的文件的更强大的方法。我们还将了解如何将find打印的文件传递给另一个名为xargs的工具，以创建功能强大的基于Unix的处理工作流。

&emsp;&emsp;使用find和xargs

&emsp;&emsp;首先，让我们看一下find和xargs解决的一些常见shell问题。假设您有一个名为process_fq的程序，它通过标准进入进程接受多个文件名。如果要对所有后缀为.fq的文件运行此程序，可以运行：
```shell
$ ls *.fq
treatment-01.fq treatment 02.fq treatment-03.fq
$ ls *.fq | process_fq
```

&emsp;&emsp;您的shell将这个通配符扩展到当前目录中的所有匹配文件，ls打印这些文件名。不幸的是，这会导致常见的复杂性，使得ls和通配符成为脆弱的解决方案。假设您的目录包含一个名为Treatment 02.fq的文件名。在本例中，ls返回Treatment 02.fq以及其他文件。但是，由于文件由空格分隔，并且此文件包含空格，process_fq会将处理02.fq解释为两个单独的文件，名为Treatment和02.fq。这个问题以不同的方式周期性地出现，在编写文件处理管道时需要注意。请注意，在参数中使用file globbing时不会发生这种情况-如果process_fq采用多个文件作为参数，则shell可以正确处理：

```shell
$ process_fq *.fq
```

&emsp;&emsp;在这里，您的shell自动转义文件名Treatment 02.fq中的空格，因此process_fq将正确接收参数Treatment-01.fq，Treatment 02.fq，Treatment-03.fq。那么为什么不使用这种方法呢？遗憾的是，可以指定为参数的文件数量是有限制的。在处理大量文件时，不太可能遇到这个限制。例如，假设您有一个tmp/目录，其中包含数千个临时文件，您希望在重新运行脚本之前删除这些临时文件。您可以尝试rm tmp/*，但会遇到问题：

```shell
$ rm tmp/*
/bin/rm: cannot execute [Argument list too long]
```

&emsp;&emsp;新的生物信息学家经常会遇到这两个问题(就我个人而言，不同的同事至少每隔一个月就会问我如何解决这些问题)。这两个问题的解决方案都是通过find和xargs来实现的，我们将在下面的章节中看到这一点。

### 使用find查找文件

&emsp;&emsp;与ls不同，find是递归的(它将在子目录和子目录的子目录中搜索匹配的文件，等等)。如果您的项目目录嵌套很深，并且您希望在整个项目中搜索文件，则“查找”非常有用。事实上，在目录上运行find(没有其他参数)可以快速查看项目目录的结构。同样，使用我们在第26页的“Organizing Data to Automate File Processing Tasks”中创建的zmays-SNPs/TOILY目录：
```shell
$ find zmays-snps
zmays-snps
zmays-snps/analysis
zmays-snps/data
zmays-snps/data/seqs
zmays-snps/data/seqs/zmaysA_R1.fastq
zmays-snps/data/seqs/zmaysA_R2.fastq
zmays-snps/data/seqs/zmaysB_R1.fastq
zmays-snps/data/seqs/zmaysB_R2.fastq
zmays-snps/data/seqs/zmaysC_R1.fastq
zmays-snps/data/seqs/zmaysC_R2.fastq
zmays-snps/scripts
```

&emsp;&emsp;find的递归搜索可以被限制为只搜索几个带有参数-maxdeep的深度目录。例如，要仅在当前目录中搜索，请使用-maxepth 1；要在当前目录及其子目录内(但不在那些子目录内)进行搜索，请使用-maxepth 2。

&emsp;&emsp;find的基本语法是find path表达式。在这里，path指定find是在哪个目录中搜索文件(如果您当前在这个目录中，它就是find)。find语法的表达式部分是find的真正力量所在。表达式是我们描述我们想要find返回的文件的方式。表达式由谓词构建，谓词由逻辑AND和OR运算符链接在一起。find仅在表达式计算为true时返回文件。通过表达式，find可以基于创建时间或文件权限等条件匹配文件，以及这些条件的高级组合，例如“查找上周之后创建的所有具有只读权限的文件”。

&emsp;&emsp;为了了解简单谓词是如何工作的，让我们看看如何使用find通过-name谓词按文件名匹配文件。前面我们对ls使用了未加引号的通配符，shell将其扩展为所有匹配的文件名。使用find时，我们引用模式(就像我们对grep所做的那样)，以避免shell解释*这样的字符。例如，假设我们想要查找所有匹配模式“zmaysB*FASTQ”的文件(例如，样本“B”中的FASTQ文件，两个读取对)以传递到管道。我们将使用Example 12-1中所示的命令：

示例12-1。通过名称匹配查找。
```shell
$ find data/seqs -name "zmaysB*fastq"
data/seqs/zmaysB_R1.fastq
data/seqs/zmaysB_R2.fastq
```
&emsp;&emsp;正如我们所预期的那样，这给出了与ls zmaysB*FASTQ类似的结果。主要区别是find报告结果由换行符分隔，默认情况下，find是递归的。

### Find‘s表达式

&emsp;&emsp;find的表达式允许您使用简单的语法缩小特定文件的范围。在前面的示例中，find命令(示例12-1)将返回同样匹配模式“zmaysB*FASTQ”的目录。因为我们只想返回FASTQ文件(而不是具有匹配名称的目录)，所以我们可能希望使用-type选项来限制我们的结果：
```shell
$ find data/seqs -name "zmaysB*fastq" -type f
data/seqs/zmaysB_R1.fastq
data/seqs/zmaysB_R2.fastq
```
&emsp;&emsp;您可以搜索许多不同的类型(例如，文件、目录、命名管道、链接等)，但最常用的是f表示文件，d表示目录，l表示链接。

&emsp;&emsp;默认情况下，Find使用逻辑AND连接表达式的不同部分。在这种情况下，find命令返回名称与“zmaysB*FASTQ”匹配并且是文件(类型为“f”)的结果。Find还允许使用不同的运算符显式连接表达式的不同部分。前面的命令相当于：
```shell
$ find data/seqs -name "zmaysB*fastq" -and -type f
data/seqs/zmaysB_R1.fastq
data/seqs/zmaysB_R2.fastq
```

&emsp;&emsp;我们可能还需要来自示例A或C的所有FASTQ文件。在这种情况下，我们希望将表达式与另一个运算符-or链接起来(有关完整列表，请参见表12-3)：
```shell
$ find data/seqs -name "zmaysA*fastq" -or -name "zmaysC*fastq" -type f
data/seqs/zmaysA_R1.fastq
data/seqs/zmaysA_R2.fastq
data/seqs/zmaysC_R1.fastq
data/seqs/zmaysC_R2.fastq
```
表12-3、常见的fnd表达式和运算符
```shell
todo
```

&emsp;&emsp;选择这些文件的相同方法是使用否定：
```shell
$ find seqs -type f "!" -name "zmaysC*fastq"
seqs/zmaysA_R1.fastq
seqs/zmaysA_R2.fastq
seqs/zmaysB_R1.fastq
seqs/zmaysB_R2.fastq
```

&emsp;&emsp;让我们来看一个更高级的例子。假设一个乱七八糟的合作者决定在seqs/中创建一个名为zmaysB_r1-temp.fastq的文件。您会注意到这个文件，因为现在您的find命令正在匹配它(我们仍然在zmays/data目录中)：
```shell
$ touch seqs/zmaysB_R1-temp.fastq
$ find seqs -type f "!" -name "zmaysC*fastq"
seqs/zmaysB_R1-temp.fastq
seqs/zmaysA_R1.fastq
seqs/zmaysA_R2.fastq
seqs/zmaysB_R1.fastq
seqs/zmaysB_R2.fastq
```
&emsp;&emsp;您不想删除或重命名他的文件，因为您的合作者可能需要该文件和/或依赖具有该特定名称的文件。所以，处理这个问题的最好方法似乎是改变你的find命令，然后和你的合作者谈谈这个神秘的文件。幸运的是，find允许这种高级文件查询：
```shell
$ find seqs -type f "!" -name "zmaysC*fastq" -and "!" -name "*-temp*"
seqs/zmaysA_R1.fastq
seqs/zmaysA_R2.fastq
seqs/zmaysB_R1.fastq
seqs/zmaysB_R2.fastq
```

&emsp;&emsp;请注意，应该将find的运算符(如！、(和)用引号引起来，以避免shell解释这些运算符。

### find‘s-exec：运行find的结果

&emsp;&emsp;虽然find对于定位文件很有用，但它在生物信息学中的真正优势是作为一种工具，以编程方式访问您想要运行命令的某些文件。在上一节中，我们看到了find的表达式如何允许您选择符合特定条件的不同文件。在本节中，我们将了解find如何允许您使用find的-exec选项对每个文件find return运行命令。
&emsp;&emsp;让我们看一个简单的例子来理解-exec是如何工作的。从上一个示例继续，假设一个乱七八糟的协作者创建了许多临时文件。让我们模拟一下(在zmays-snps/data/seqs目录中)：
```shell
$ touch zmays{A,C}_R{1,2}-temp.fastq
$ ls
zmaysA_R1-temp.fastq zmaysB_R1-temp.fastq zmaysC_R1.fastq
zmaysA_R1.fastq zmaysB_R1.fastq zmaysC_R2-temp.fastq
zmaysA_R2-temp.fastq zmaysB_R2.fastq zmaysC_R2.fastq
zmaysA_R2.fastq zmaysC_R1-temp.fastq
```
&emsp;&emsp;假设您的协作者允许您删除所有这些临时文件。删除这些文件的一种方法是使用rm*-temp.fast。然而，在充满重要数据文件的目录中使用通配符的RM太危险了。如果您意外地在*和-temp.fast q之间留了一个空格，通配符*将匹配当前目录中的所有文件，并将它们传递给RM，从而导致目录中的所有内容被意外删除。使用find的-exec是删除这些文件的更安全的方法。

&emsp;&emsp;find‘s-exec的工作方式是将与find的表达式匹配的每个文件传递给-exec指定的命令。使用-exec时，有必要在命令的末尾使用分号来表示命令的结束。例如，让我们使用find-exec和rm-i删除这些临时文件。rm的-i强制rm是交互式的，提示您确认要删除文件。我们的查找和删除命令是：
```shell
$ find . -name "*-temp.fastq" -exec rm -i {} \;
remove ./zmaysA_R1-temp.fastq? y
remove ./zmaysA_R2-temp.fastq? y
remove ./zmaysB_R1-temp.fastq? y
remove ./zmaysC_R1-temp.fastq? y
remove ./zmaysC_R2-temp.fastq? y
```

&emsp;&emsp;在一行中，我们能够务实地识别并对匹配特定模式的文件执行命令。我们的命令是RM，但也可以很容易地成为生物信息学程序。使用这种方法允许您对目录中任意数量的文件调用程序。使用find和-exec时，像用程序处理包含100，000个文本文件的目录这样令人畏惧的任务非常简单。

!['suggestion'](../img/suggestion.png)使用find-exec删除文件

    使用find-exec删除文件是一种常见的操作，find也有一个-delete选项，您可以使用它来代替-exec-rm{}(但它不会是交互式的，不像使用-i的rm)。

&emsp;&emsp;使用-exec时，始终首先编写表达式，并检查返回的文件是否是要应用命令的文件。然后，在find命令返回正确的文件子集之后，在-exec语句中添加。find-exec最适合于快速、简单的任务(如删除文件、更改权限等)。对于较大的任务，xargs(我们将在下一节中看到)是更好的选择。

### xargs：一个Unix工具

&emsp;&emsp;如果有一个Unix工具向我介绍了Unix令人难以置信的原始力量，那就是xargs。xargs允许我们获取从standard in传递给它的输入，并将此输入用作另一个程序的参数，这允许我们从通过standard in接收的值以编程方式构建命令(通过这种方式，它有点类似于R的do.call()。将find与xargs一起使用非常类似于使用-exec的find，但是具有一些额外的优势，这些优势使xargs成为更好的选择来执行更大的生物信息学文件处理任务。
&emsp;&emsp;让我们重新创建凌乱的临时文件目录示例(从zmayssnps/data/seqs目录)：

```shell
$ touch zmays{A,C}_R{1,2}-temp.fastq # create the test files
```

&emsp;&emsp;xargs的工作方式是从标准中获取输入，并通过空格、制表符和换行符将其拆分为参数。然后，将这些参数传递给提供的命令。例如，为了使用rm模拟find-exec的行为，我们将xargs与rm一起使用：

```shell
$ find . -name "*-temp.fastq"
./zmaysA_R1-temp.fastq
./zmaysA_R2-temp.fastq
./zmaysC_R1-temp.fastq
./zmaysC_R2-temp.fastq
$ find . -name "*-temp.fastq" | xargs rm
```
!['warning'](../img/warning.png)使用find和xargs进行安全操作

    find和xargs有一个重要的问题：文件名中的空格可能会破坏东西，因为xargs将空格视为参数分隔符。这将导致像Treatment 02.fq这样的文件名被解释为两个独立的参数，Treatment和02.fq。find和xargs开发人员创建了一个聪明的解决方案：两者都允许使用空字节作为分隔符的选项。下面是如何使用空字节分隔符运行find和xargs的示例：

    ```shell
    $ find . -name "samples [AB].txt" -print0 | xargs -0 rm
    ```
    
    除了这种预防措施，简单地不使用包含空格或其他奇怪字符的文件名也是明智的。带有破折号或下划线的简单字母数字名称最好。为了简化示例，我将省略-print0和-0，但在实践中应该始终使用它们。

&emsp;&emsp;本质上，xargs将find的输出拆分成参数，并运行：
```shell
$ rm ./zmaysA_R1-temp.fastq ./zmaysA_R2-temp.fastq \
 ./zmaysC_R1-temp.fastq ./zmaysC_R2-temp.fastq
 ```
&emsp;&emsp;xargs将通过standard in接收的所有参数传递给提供的程序(本例中为rm)。这对于rm、touch、mkdir等采用多个参数的程序非常有效。但是，其他程序一次只接受一个参数。我们可以使用xargs的-n参数设置传递给每个命令调用的参数数量。例如，我们可以使用以下命令调用RM四次(每次在一个文件上)：
```shell
$ find . -name "*-temp.fastq" | xargs -n 1 rm
```

&emsp;&emsp;xargs的一大好处是，它将指定要操作的文件(Find)的进程与对这些文件应用命令的过程(通过xargs)分开。如果我们想在对此列表中的所有文件运行rm之前检查一个长的文件列表find return，我们可以使用：
```shell
$ find . -name "*-temp.fastq" &gt; files-to-delete.txt
$ cat files-to-delete.txt
./zmaysA_R1-temp.fastq
./zmaysA_R2-temp.fastq
./zmaysC_R1-temp.fastq
./zmaysC_R2-temp.fastq
$ cat files-to-delete.txt | xargs rm
```

&emsp;&emsp;另一个常见的技巧是使用xargs构建写入简单Bash脚本的命令。例如，我们可以在rm上调用echo，而不是直接运行rm，然后允许xargs在这个命令之后放置参数(请记住，xargs的行为非常简单：它只是在您提供的命令之后放置参数)。例如：
```shell
$ find . -name "*-temp.fastq" | xargs -n 1 echo "rm -i" &gt; delete-temp.sh
$ cat delete-temp.sh
rm -i ./zmaysA_R1-temp.fastq
rm -i ./zmaysA_R2-temp.fastq
rm -i ./zmaysC_R1-temp.fastq
rm -i ./zmaysC_R2-temp.fastq
```
&emsp;&emsp;以这种方式分解任务允许我们检查使用xargs构建的命令(因为xargs运行的命令是echo，它只是打印所有内容)。然后，我们可以使用以下命令运行这个简单的脚本：
```shell
$ bash delete-temp.sh
remove ./zmaysA_R1-temp.fastq? y
remove ./zmaysA_R2-temp.fastq? y
remove ./zmaysC_R1-temp.fastq? y
remove ./zmaysC_R2-temp.fastq? y
```

### 使用带有替换字符串的xargs将命令应用于文件

&emsp;&emsp;到目前为止，xargs纯粹通过在您提供的命令的末尾添加参数来构建命令。然而，有些程序通过选项接受参数，比如program-in file.txt-out-file out.txt；其他程序有很多位置参数，比如program arg1 arg2。xargs的-i选项通过用单个参数替换占位符字符串的所有实例，允许将参数更细粒度地放置到命令中。按照惯例，我们与-i一起使用的占位符字符串是{}。

&emsp;&emsp;让我们来看一个例子。假设假想的程序FASTQ_STAT通过选项-in获取输入文件，收集FASTQ统计信息，然后将摘要写入由-out选项指定的文件。正如我们的Bash循环示例(第405页上的“使用For Loops和globbing处理带有Bash的文件”)一样，我们希望输出文件名与输入文件名成对，并具有相应的名称。我们可以使用find、xargs和basename来解决这个问题。第一步是使用find获取要处理的文件，然后使用xargs和basename提取示例名称。basename允许我们通过参数-s删除后缀：
```shell
$ find . -name "*.fastq" | xargs basename -s ".fastq"
zmaysA_R1
zmaysA_R2
zmaysB_R1
zmaysB_R2
zmaysC_R1
zmaysC_R2
```

&emsp;&emsp;然后，我们想要运行命令fastq_stat-in file.fastq-out../summary/file.txt，但是用文件的基本名称替换file。为此，我们将前面创建的示例名称通过管道传输到另一个运行FASTQ_STAT的xargs命令：
```shell
$ find . -name "*.fastq" | xargs basename -s ".fastq" | \
 xargs -I{} fastq_stat --in {}.fastq --out ../summaries/{}.txt
```
!['warning'](../img/warning.png)BSD和GNU xargs

    不幸的是，-i的行为在BSD xargs(这是OS X使用的)和GNU xargs之间是不同的。默认情况下，BSD xargs最多只能替换-I指定的字符串的五个实例，除非使用-R设置更多实例。一般来说，最好使用GNU xargs。如果你在Mac上，你可以用Homebrew安装GNU Coreutils。为了防止与系统的xargs(BSD版本)发生冲突，Homebrew在其版本前加上g，因此xargs的GNU版本名称是gxargs
 
&emsp;&emsp;将xargs与basename结合使用是一种强大的习惯用法，用于将命令应用于许多文件，以便跟踪特定输入文件创建的输出文件。虽然我们可以通过其他方式(例如，通过Bash for loops或自定义Python脚本)实现这一点，但xargs允许非常快速和增量的命令构建。然而，正如我们将在下一节中看到的，xargs与for loop相比还有一个非常大的优势：它允许在预先指定数量的进程上进行并行化。总体而言，可能需要一些练习才能让这些xargs技巧掌握在你的手指下，但它们将为你几十年的生物信息学工作提供很好的服务。

### xargs和并行化

&emsp;&emsp;xargs令人难以置信的强大特性是它可以并行启动有限数量的进程。我强调数量有限，因为这是xargs相对于Bash的for循环的优势之一。我们可以使用Bash for循环启动许多个后台进程，在具有多个处理器的系统上，这些进程将并行运行(取决于其他正在运行的任务)：
```bash 
for filename in *.fastq; do
    program "$filename" &
done
```
&emsp;&emsp;但是这会启动很多后台进程，在*.fastq中有文件！在共享服务器上，这肯定不是良好的计算礼仪，即使您是此服务器的唯一用户，这也可能导致瓶颈，因为所有进程都开始从磁盘读取和写入磁盘。因此，当并行运行多个进程时，我们希望明确限制同时运行的进程数量。Xargs允许我们使用选项-P<；num>；进行此操作，其中<；num>；是同时运行的进程数。

&emsp;&emsp;让我们看一个简单的例子-并行运行我们想象的程序FASTQ_STAT，最多使用六个进程。我们通过在第二个xargs调用中添加-P6来实现这一点(并行化basename命令没有意义，因为这将非常快)：
```shell
$ find . -name "*.fastq" | xargs basename -s ".fastq" | \
 xargs -P 6 -I{} fastq_stat --in {}.fastq --out ../summaries/{}.txt
 ```

&emsp;&emsp;通常，FASTQ_STAT可以是任何程序，甚至可以是每个示例执行许多任务的shell脚本。关键的一点是，我们提供了程序或脚本运行示例名称所需的所有信息，这就是替换字符串{}的内容。

!['suggestion'](../img/suggestion.png)xargs、管道和重定向

    初学者经常遇到的一个绊脚石是试图使用管道和使用xargs进行重定向。这不起作用，因为读取xargs命令的shell进程会将管道和重定向解释为如何处理xarg的输出，而不是作为xargs运行的命令的一部分。绕过这个限制的最简单和最干净的技巧是创建一个小的Bash脚本，其中包含处理单个样本的命令，并让xargs在许多并行的Bash进程中运行这个脚本。例如：

    ```bash
    #!/bin/bash

    set -e
    set -u
    set -o pipefail
    sample_name=$(basename -s ".fastq" "$1")
    some_program ${sample_name}.fastq | another_program &gt;
    ${sample_name}-results.txt
    ```
    然后运行:
    ```bsh
    $ find . -name "*.fastq" | xargs -n 1 -P 4 bash script.sh
    ```
    其中-n 1强制xargs一次处理一个输入参数。这可以通过指定使用-P运行多少个进程来轻松并行化。

&emsp;&emsp;无可否认，一些功能强大的xargs工作流的代价是复杂性。如果您发现自己主要使用xargs来并行化任务，或者您正在编写使用basename的复杂xargs命令，那么学习GNU parallel可能是值得的。GNU parallel扩展和增强了xargs的功能，并修复了xargs的几个限制。例如，GNU parallel可以处理命令中的重定向，并且有一个快捷方式({/.})来提取不带baseName的基本文件名。这允许使用非常短、功能强大的命令：
```bash
$ find . -name "*.fastq" | parallel --max-procs=6 'program {/.} &gt; {/.}-out.txt'
```
&emsp;&emsp;GNU parallel有许多其他选项和特性。如果您发现自己经常使用xargs来处理复杂的工作流，我建议您学习更多关于GNU parallel的知识。GNU Parallel网站有很多例子和详细的教程。

### Make和Makefles：管道的另一种选择

&emsp;&emsp;尽管本章主要关注从Bash构建管道，并使用find和xargs将命令应用于某些文件，但我不能忽视快速介绍另一个用于创建生物信息管道的非常强大的工具。这个工具是make，它解释makefles(用它们自己的makefile语言编写)。make旨在编译软件，这是一个复杂的过程，因为编译文件的每个步骤都需要确保每个依赖项都已编译或可用。与SQL(我们在第13章中讨论)一样，makefile语言是声明性的-与Bash脚本不同，makefile不会从上到下运行。相反，makefile被构建为一组规则，其中每个规则包含三个部分：目标、先决条件和配方。每个菜谱都是一组用于构建目标的命令，目标是一个文件。先决条件指定配方需要哪些文件来构建目标文件(依赖项)。make的惊人独创性在于，该程序从先决条件和目标中找出了如何使用所有规则为您构建文件。让我们看一个简单的示例-我们想要编写一个简单的管道，它从Internet下载一个文件并创建它的摘要：
```shell
FASTA_FILE_LINK=http://bit.ly/egfr_flank (1)
.PHONY: all clean      (2)
all: egfr_comp.txt        (3)
egfr_flank.fa:              (4)
 curl -L $(FASTA_FILE_LINK) &gt; egfr_flank.fa
egfr_comp.txt: egfr_flank.fa                               (5)
 seqtk comp egfr_flank.fa &gt; egfr_comp.txt
clean:   (6)
 rm -f egfr_comp.txt egfr_flank.fa
```
(1)我们在Makefile中定义变量，就像在Bash脚本中一样。我们保持这个链接到这个FASTA文件在顶部，所以它是值得注意的，可以很容易地调整。

(2)这个makefile中的目标ALL和CLEAN不是文件的名称，而只是我们可以引用的目标的名称。我们通过在特殊目标.PHONY中将它们指定为先决条件来指示这些目标不是文件。

(3)All是用于构建此makefile要构建的所有文件的目标的常规名称。在这里，这个简单的示例Makefile的最终目标是从Internet下载一个Fasta文件并在其上运行seqtk comp，返回这个Fasta文件的序列组成。我们写入序列合成的最后一个文件是EGFR_comp.txt，所以这是all目标的先决条件。

(4)此规则创建文件EGFR_fank.fa。在这个规则中没有先决条件，因为创建egfr_fank.fa不需要本地文件(因为我们正在下载这个文件)。我们的配方使用curl下载存储在变量FASTA_FILE_LINK中的链接。由于这是一个缩短的链接，我们使用curl的-L选项来跟踪重定向。最后，请注意，要引用makefile中的变量值，我们使用语法$(Variable)。

(5)此规则创建文件EGFR_fank.fa。此规则创建包含排序组成数据EGFR_comp.txt的文件。因为我们需要FASTA文件EGFR_fank.fa来创建此文件，所以我们将EGFR_fank.fa指定为先决条件。配方在先决条件上运行seqtk comp，并将输出重定向到目标文件egfr_comp.txt。

(6)最后，有一个名为clean的目标是很常见的，它包含一个用于清理这个makefile生成的所有文件的食谱。这允许我们运行make clean并将目录返回到运行makefile之前的状态。

&emsp;&emsp;我们使用make命令运行makefile。对于前面的makefile，我们使用make all运行它，其中all参数指定make应该以all目标开始。然后，程序make将首先在当前目录中搜索名为Makefle的文件，加载它，然后从目标ALL开始。这将如下所示：

&emsp;&emsp;make的一个特别强大的特性是，它只在目标不存在或其先决条件被修改时才生成目标。这是非常强大的：如果你有一个长而复杂的makefile，并且你修改了一个文件，make只会重新运行依赖于这个修改过的fle的目标的配方(假设你完全指定了依赖项)。注意如果我们再次运行make all会发生什么：
```shell
$ make all
make: Nothing to be done for `all'.
```
&emsp;&emsp;因为已经创建了所有目标，并且没有修改输入文件，所以没有什么可做的。现在，看看如果我们使用touch更改egfr_fank.fa文件的修改时间会发生什么：
```shell
$ touch egfr_flank.fa
$ make all
seqtk comp egfr_flank.fa &gt; egfr_comp.txt
```
&emsp;&emsp;因为EGFR_fank.fa是创建EGFR_comp.txt文件的先决条件，所以请使用最新版本的EGFR_fank.txt重新运行此规则以更新EGFR_comp.txt。

&emsp;&emsp;最后，我们可以删除使用CLEAN目标创建的所有文件：
```shell
$ make clean
rm -f egfr_comp.txt egfr_flank.fa
```
&emsp;&emsp;在这个例子中，我们只是触及make的功能的表面；关于这种语言的完整教程超出了本书的范围。不幸的是，像Bash一样，make的语法对于一些更高级的任务来说可能非常神秘和复杂。此外，因为makefile是以声明性的方式编写的(并以非线性方式执行)，所以调试makefile可能非常棘手。尽管如此，make仍然是一个有用的工具，如果您需要简单任务和工作流的选项，您应该知道它是一个有用的工具。有关更多信息，请参见GNU make文档。