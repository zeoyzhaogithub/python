# 第二部分 先决条件：生物信息学项目入门必备技能
## 第二章、建立和管理生物信息学项目

&emsp;&emsp;就像一个组织良好的实验室会让科学家的生活变得更容易一样，一个组织良好且有完整文档记录的项目会让一个生物信息学家的工作变得更容易。无论您正在处理的是什么特定的项目，您的项目目录都应该以一致和可理解的方式进行布局。清晰的项目组织使你和合作者更容易弄清楚每件事的确切位置和内容。此外，当文件分布规范并且命名明确时，自动执行任务要容易得多。例如，使用脚本处理存储在单独的FASTA文件中的300个基因序列是微不足道的，如果这些文件被组织在单个目录中并且名称一致。
    
&emsp;&emsp;每个生物信息学项目都是从一个空的项目目录开始的，所以本书从项目组织的一章开始是合适的。在本章中，我们将研究组织您的生物信息学项目目录的一些最佳实践，以及如何使用纯文本Markdown文件对您的工作进行数字文档记录。我们还将了解为什么项目目录组织不仅仅是关于整洁，而是对于跨大量文件自动执行任务的方式至关重要(我们在生物信息学中经常这样做)。

### 项目目录及其结构

&emsp;&emsp;创建组织良好的目录结构是可复制生物信息学项目的基础。实际过程非常简单：布置一个项目只需要用mkdir创建几个目录，用touch创建空的readme文件(命令我们将在后面进行更深入地研究)。但从长远来看，这种简单的初始规划是很有益的。对于大型项目，研究人员可能会花费数年时间在此目录结构中工作。
    
&emsp;&emsp;其他研究人员也注意到了良好的项目组织的重要性(Noble 2009)。虽然最终您将开发一个适合您的项目组织方案，但我们将在本章从我在工作中使用的一个方案开始(与Noble的方案类似)。
    
&emsp;&emsp;项目中使用的所有文件和目录都应该存在于一个具有明确名称的单个项目目录中。在项目过程中，您将积累数据文件、笔记、脚本等-如果这些文件分散在您的硬盘驱动器上(或者更糟，分布在许多计算机的硬盘驱动器上)，那么跟踪所有内容将是一场噩梦。更糟糕的是，这样一个混乱的项目以后会让你的研究几乎不可能重现。
    
&emsp;&emsp;将所有文件保存在单个目录中将极大地简化您和您的合作者的工作，并促进项目的可复现性(我们将在第5章讨论如何与Git协作处理项目目录)。假设您使用SNP调用玉米(Zea Mays)。您的第一步将是选择一个简短的、适当的项目名称，并创建一些基本目录：

``` shell
$ mkdir zmays-snps
$ cd zmays-snps
$ mkdir data
$ mkdir data/seqs scripts analysis
$ ls -l
total 0
drwxr-xr-x 2 vinceb staff 68 Apr 15 01:10 analysis
drwxr-xr-x 3 vinceb staff 102 Apr 15 01:10 data
drwxr-xr-x 2 vinceb staff 68 Apr 15 01:10 scripts
```

&emsp;&emsp;这是一个合理的项目布局方案。在这里，data/包含所有原始和中间数据。正如我们将看到的，数据处理步骤在这个data/目录中被视为单独的子项目。我将通用项目范围的脚本保存在scripts/目录中。如果脚本包含许多文件(例如，多个Python模块)，它们应该驻留在自己的子目录中。在开发这些脚本时(当它们生成测试输出文件时)，将脚本隔离在它们自己的子目录中也可以保持项目目录的整洁。
    
&emsp;&emsp;生物信息学项目包含许多较小的分析-例如，分析原始序列的质量、比对输出和最终数据，这些数据将为论文生成图表。我更喜欢将这些保存在单独的analysis/目录中，因为它允许合作者查看这些高级分析，而不必更深入地挖掘子项目目录。

![suggestion](../../img/suggestion.png) **名字里有什么？**

&emsp;&emsp;在计算机上命名文件和目录比你想象的更重要。在从基于图形用户界面(GUI)的操作系统过渡到Unix命令行的过程中，许多人带来了在文件和目录名称中使用空格的坏习惯。这在基于Unix的环境中是不合适的，因为空格用于分隔命令中的参数。例如，假设您从GUI(例如，通过OS X的Finder)创建名为raw sequences的目录，然后尝试使用以下命令删除它及其内容：

``` shell
$ rm -rf raw sequences
```

&emsp;&emsp;如果你幸运的话，会看到警示语‘不存在名称为raw和sequences的目录或文件’。这里发生了什么事？空格问题：您的shell将这个rm命令解释为“删除raw和sequences文件/目录”，而不是“删除名称为raw sequences的单个文件或目录”。
    
&emsp;&emsp;如果你不够幸运，存在一个名称为raw或sequences的文件或目录，这个rm命令将会删除它。可以通过使用引号(例如rm -r "raw sequences")对其进行转义，但更好的做法是不要在文件或目录名称中使用空格。最好在文件和目录名称中只使用字母、数字、下划线和破折号。

&emsp;&emsp;尽管Unix不需要文件扩展名，但在文件名中包含扩展名有助于指示每个文件的类型。例如，名为osativa-genes.fasta的文件表明这是FASTA格式的序列文件。相反，名为osativa-gene的文件可以是基因模型的文件、关于这些Oryza sativa基因来自哪里的注释或序列数据。当有疑问时，当涉及到文件名、文档和编写代码时，显式总是比隐式更好。
    
&emsp;&emsp;脚本和分析通常需要引用项目层次结构中的其他文件(例如数据)。这可能需要引用目录层次结构中的父目录(例如，使用...)。在这些情况下，一定要始终使用相对路径(例如../data/stats/qual.txt)而不是绝对路径(例如/home/vinceb/projects/zmays-snps/data/stats/qual.txt).只要您的内部项目目录结构保持不变，这些相对路径将始终有效。相反，绝对路径依赖于项目目录级以上的特定用户帐户和目录结构细节(不好)。使用绝对路径会使您的工作在协作者之间的可移植性降低，并降低可重复性。
    
&emsp;&emsp;我这里的项目目录方案绝不是唯一的方案。我与生物信息学家合作过，他们的项目使用了完全不同的方案，而且做得非常出色。然而，不管组织方案如何，一个好的生物信息学家总是会广泛地记录一切，并使用计算机可以解析的清晰的文件名，这两点我们将在稍后讨论。

### 项目文档
    
&emsp;&emsp;除了有一个组织良好的目录结构外，你的生物信息学项目也应该有很好的文档记录。糟糕的文档可能导致不可重复性和严重错误。生物信息学工作中潜伏着大量的复杂性：复杂的工作流程，多个文件，无数的程序参数，以及不同的软件版本。防止这种复杂性导致问题的最好方法是详细记录各个步骤的所有内容。当您需要返回并重新运行分析，编写关于论文步骤的详细方法，或在目录中查找某些数据的来源时，详细的文档会让你的工作变得更轻松。那么你到底应该记录什么呢？以下是一些建议：

*记录您的方法和工作流程*
    
&emsp;&emsp;这应该包括通过运行shell生成数据或中间结果的完整命令行(复制和粘贴)。即使你在软件中使用默认值，也一定要写下这些值；程序的较高版本可能使用不同的默认值。脚本自然地记录了所有的步骤和参数(我们将在第12章中讨论这个主题)，但是一定要记录用于运行这个脚本的所有命令行及其参数。通常，任何在您的工作中使用的产生结果的命令都需要记录在某个地方。

*记录项目目录中所有数据的来源*

&emsp;&emsp;您需要跟踪数据是从哪里下载的，是谁提供给您的，以及任何其他相关信息。“数据”不仅仅是指项目的实验数据-它是程序用来创建输出的任何数据。这包括您的合作者通过各自的分析、基因注释轨迹、参考基因组等发送给您的文件。记录这些关于您的数据或元数据的重要数据是至关重要的。例如，如果您下载了一组基因区域，请记录网站的URL。这似乎是一个明显的建议，但是我无数次遇到了一个分析步骤，因为有人忘记记录数据的来源，这一步骤很难被复制。

*记录下载数据的时间*
    
&emsp;&emsp;重要的是要包括下载数据的时间，因为外部数据源(如网站或服务器)将来可能会发生变化。例如，直接从数据库下载数据的脚本如果在更新外部数据库后重新运行，可能会产生不同的结果。因此，记录数据何时进入存储库是很重要的。

*记录数据版本信息*
    
&emsp;&emsp;许多数据库都有明确的发布号、版本号或名称(例如，拟南芥的TAIR10版本的基因组注释，或线虫的Wormbase版本WS231)。重要的是在文档中记录所有版本信息，包括次要版本号。

*描述您是如何下载数据的*
     
&emsp;&emsp;例如，您是否使用MySQL下载了一组基因？还是UCSC基因组浏览器？这些详细信息对跟踪问题非常有用，例如当协作者之间的数据不同时。

*记录您运行的软件的版本*
      
&emsp;&emsp;这可能看起来并不重要，但请记住第6页“可重现研究”中的例子，我和我的同事将不一致的结果追溯到正在更新的单个软件上。这些细节很重要。好的生物信息学软件通常有一个命令行选项来返回当前版本信息。使用版本控制系统(如Git)管理的软件对每个版本都有明确的标识符，可以用来记录您运行的精确版本(我们将在第5章中了解更多)。如果没有可用的版本信息，发布日期、软件链接和下载日期就足够了。

&emsp;&emsp;所有这些信息最好存储在 README文件中。可以直接从命令行轻松读取、搜索和编辑纯文本，使其成为便携和可访问的README文件的最佳选择。它还可以在所有计算机系统上使用，这意味着您可以直接在服务器或计算机集群上工作时记录您的步骤。纯文本也缺乏复杂的格式，这可能会在将命令从文档复制并粘贴回命令行时产生问题。最好避免像Microsoft Word for README文档这样的格式，因为这些格式对生物信息学中常见的Unix系统的可移植性较差。
      
&emsp;&emsp;你应该把你的README文件放在哪里？一种好的方法是在项目的每个主目录中保存README文件。这些README文件不一定需要很长，但它们至少应该解释这个目录中的内容以及它是如何到达那里的。即使是这个小小的努力也可以节省人们探索您的项目目录的大量时间，并防止混淆。这个人可以是你的顾问或合作者，可能是在你搬到另一个实验室后试图重现你的工作的同事，甚至可能是你六个月后完全忘记自己做过的事情时的自己(每个人都会发生这种情况！)。
       
&emsp;&emsp;例如，data/README文件将包含data/目录中有关数据文件的元数据。即使你认为你可以记住关于你的数据的所有相关信息，仅仅把它放在自述文件中要容易得多(而且合作者不需要通过电子邮件来询问你有什么文件或者它们在哪里)。让我们使用touch创建一些空的README文件。touch创建一个文件(如果该文件尚不存在)，或者更新文件的修改时间。我们可以将其用于前一个目的，以创建空模板文件来布局我们的项目结构：

```shell
$ touch README data/README
```

&emsp;&emsp;根据刚刚讨论的文档指南，此data/READM文件将包括您在data/中下载数据的位置、何时下载以及如何下载。当我们在第6章中了解更多关于数据的内容时，我们将看到一个如何在项目目录中下载和正确记录数据的研究案例示例(第120页上的“案例研究：可重现下载数据”)。
       
&emsp;&emsp;通过记录这些信息，我们开始记录关于我们的实验和分析的一切，使其可重现。请记住，随着项目的增长和数据文件的积累，记录数据文件的这些信息对于你梳理自己的项目是非常用益的。

### 使用目录将项目划分为多个子项目

&emsp;&emsp;生物信息学项目涉及许多子项目和子分析。例如，在通过生物信息学工具(如校准器或装配器)运行之前，应对原始实验数据的质量进行评估，并删除低质量数据(我们在第346页的“示例：检查和修剪低质量碱基”中将会看到这方面的一个例子)。
       
&emsp;&emsp;创建目录以便从逻辑上分隔子项目(例如，排序数据质量改进，对齐，分析对齐结果等)可以简化复杂的项目并帮助保持文件的组织结构。它还有助于降低使用有漏洞的脚本意外破坏文件的风险，因为子目录有助于隔离事故。将一个项目分解为几个子项目，并将它们保存在单独的子目录中，这也使得记录您的工作更容易；每个READM文件都与其所在的目录有关。最终，您将做到让这个系统为您的项目组织工作；其要点是：利用目录来帮助保持组织结构。

### 组织数据以实现文件处理任务自动化
&emsp;&emsp;因为自动化文件处理任务是生物信息学不可或缺的一部分，所以组织我们的项目以促进这一点是至关重要的。将数据组织到子目录中并使用清晰一致的文件命名方案是势在必行的-这两种做法都允许我们以编程方式引用文件，这是自动化任务的第一步。以编程方式做某事意味着通过代码而不是手动完成，使用一种可以毫不费力地扩展到多个对象(例如文件)的方法。以编程方式引用多个文件比键入所有文件更容易、更安全(因为它不太容易出错)。

![suggestion](../../img/suggestion.png)**shell扩展提示**

&emsp;&emsp;生物信息学家、软件工程师和系统管理员花费大量时间在终端上打字。不足为奇的是，这些人收集技巧，使这个过程尽可能有效。当你花更多的时间在shell上时，你会发现花一点时间学习这些技巧可以为你节省很多时间。
        
&emsp;&emsp;一个有用的技巧是shell扩展。shell扩展是指当您的shell(例如，Bash，很可能是您正在使用的shell)为您展开文本时，您不必键入它。如果您曾经键入cd~转到您的home目录，那么您已经使用过shell展开-是您的shell将代字号(~)展开为您的home目录的完整路径(如 /users/vinceb/)。您的shell也会将通配符如星号(*)扩展到所有匹配的文件。

&emsp;&emsp;一种称为支架扩展的shell扩展可用于使用单个命令快速创建简单的zmays-snps/项目目录结构。大括号扩展通过扩展大括号内的逗号分隔值来创建字符串。通过一个简单的例子可以更容易的理解这一点：

```shell
$ echo dog-{gone,bowl,bark}
dog-gone dog-bowl dog-bark
```
&emsp;&emsp;使用相同的策略，我们可以使用下面的命令创建zmays-snps/项目目录：

```shell
$mkdir -p zmays-snps/{data/seqs,scripts,analysis}
```
&emsp;&emsp;用这个命令行生成的zmays-snps项目布局与我们在第21页的“项目目录和目录结构”中的四个单独步骤构建的项目布局：analysis/、data/seqs和scripts/相同。因为mkdir接受多个参数(为每个参数创建一个目录)，所以这会创建三个子目录(并且省去了您键入“zmays-snps/”三次的麻烦)。注意，我们需要使用mkdir的-p标签，它告诉mkdir创建它需要的任何必要的子目录(在我们的例子中，data/来创建data/seqs/)。
      
&emsp;&emsp;我们将通过一个玩具示例来说明这一点，并在此过程中学习一些重要的shell通配符技巧。在本例中，将数据文件组织到具有一致文件名的单个目录中，使我们做好了迭代所有数据的准备，无论是本例中使用的四个示例文件，还是实际项目中的40，000个文件。这样想：还记得当你发现你可以用鼠标光标选择很多文件的时候吗？使用这个技巧，你可以像移动六个文件一样轻松地移动60个文件。您也可以选择某些文件类型(如照片)，并将它们全部附加到一封电子邮件中，只需一次移动即可。通过使用一致的文件命名和目录组织，您可以使用Unix shell和其他编程语言以编程方式完成相同的操作。我们将看到一个Unix示例，使用shell通配符在许多文件中自动执行任务。稍后在第12章中，我们将看到更高级的批量文件操作策略。
      
&emsp;&emsp;让我们创建一些假的空数据文件，看看一致的名称如何帮助以编程方式处理文件。假设我们有三个玉米样品，“A”、“B”和“C”，以及每种样品的成对末端测序数据：

```shell     
cd data$ touch seqs/zmays{A,B,C}_R{1,2}.fastq
$ ls seqs/
zmaysA_R1.fastq zmaysB_R1.fastq zmaysC_R1.fastq
zmaysA_R2.fastq zmaysB_R2.fastq zmaysC_R2.fastq
```
     
&emsp;&emsp;在此文件命名方案中，每个文件名的两个变量部分表示样本名称(zmaysA，zmaysB，zmaysC)和读取对(r1和r2)。假设我们想要以编程方式检索具有示例名称zmaysB的所有文件(无论读取对是什么)，而不是必须手动指定每个文件。为此，我们可以使用Unix shell通配符，星号字符(*)：

```shell    
$ ls seqs/zmaysB*zmaysB_R1.fastq zmaysB_R2.fastq
```
      
&emsp;&emsp;通配符扩展到所有匹配的文件或目录名称(此过程称为全局匹配)。在前面的示例中，shell将表达式zmaysB*扩展为zmaysB_R1.fastq和zmaysB_R2.fastq，因为这两个文件以zmaysB开头。如果这个目录包含数百个zmaysB文件，那么所有文件都可以很容易地引用和使用shell通配符进行处理。

​​
![warning](../../img/warning.png)**通配符和“参数列表太长”**
       
&emsp;&emsp;OSX和Linux系统对可以提供给命令行的参数数量有限制(更严格地说，限制参数的总长度)。当使用与数万个文件匹配的通配符时，我们有时会达到这个限制。当这种情况发生时，您将看到一条“参数列表太长”的错误消息，表明您的参数长度已经达到了上限。幸运的是，有一个解决这个问题的聪明方法(参见411页上的“使用find和xargs”获得解决方案)。

       
&emsp;&emsp;通常，尽可能严谨的使用通配符。这可以防止意外匹配。例如，如果一个粗心的同事在这个目录中创建了一个名为zmaysB-interest-snps-ound.xls的Excel文件，这会意外地匹配通配符表达式zmaysB*。如果需要处理所有zmaysB FASTQ文件，使用zmaysB*匹配它们会包含此Excel文件并导致异常情况。这就是为什么尽可能严谨的使用通配符。使用zmaysB*fastq或zmaysB_R？.fastq(？仅匹配单个字符)，代替zmaysB*。
       
&emsp;&emsp;还有其他一些简单的shell通配符，它们在以编程方式访问文件时使用非常方便。假设一个合作者告诉你C样本序列的质量很差，所以当C被重新排序时，你将不得不只处理A和B样本。在收到新的样本之前，您不想删除zmaysC_R1.fastq和zmaysC_R2.fast，因此在此期间您希望忽略这些文件。发明通配符的人预见到了这样的问题，所以他们创建了shell通配符，允许您匹配特定的字符或字符范围。例如，我们可以用[UVWXY]或[U-Y]匹配字符U、V、W、X和Y(两者都是等效的)。回到我们的示例，我们可以使用以下任一方法排除C示例：

```shell       
$ ls zmays[AB]_R1.fastqzmaysA_R1.fastq zmaysB_R1.fastq
$ ls zmays[A-B]_R1.fastq
zmaysA_R1.fastq zmaysB_R1.fastq
```
      
&emsp;&emsp;使用A和B之间的范围并不是真正必要的，但是如果我们有样本A到I，那么使用像zmays[C-I]_R1.fastq这样的范围将比输入zmays[CDEFGHI]_R1.fastq更好。有一个非常重要的警告：范围只可以对字符范围进行操作，而不能操作像13到30这样的数字范围。这意味着像snps_[10-13].txt这样的通配符将不匹配文件snps_10.txt、snps_11.txt、snps_12.txt和snps_13.txt。
      
&emsp;&emsp;然而，shell确实提供了对数值范围的扩展解决方案-通过我们前面看到的大括号扩展。在我们看到这个快捷方式之前，请注意，虽然通配符匹配和花括号扩展看起来可能表现类似，但它们略有不同。通配符仅扩展到与其匹配的现有文件，而大括号扩展始终扩展，而不管是否存在相应的文件或目录。如果我们知道文件snps_10.txt到snps_13.txt确实存在，我们可以将它们与大括号扩展序列表达式(如snps_{10..13}.txt)进行匹配。这将扩展为整数序列10到13(但请记住，这些文件是否存在不会通过大括号展开进行检查)。表2-1列出了常见的Unix通配符。

表2-1、通用Unix文件名通配符

| 通配符  | 匹配的内容 |
|:-------|:---------|
|*|零个或多个字符(但忽略以句点开头的隐藏文件)。|
|？|一个字符(也忽略隐藏文件)。|
|[A-Z]  |提供的字母数字范围内的任何字符(在本例中为A和Z之间的任何字符)；这适用于任何字母数字字符范围(例如，[0-9]匹配0到9之间的任何字符)。|

&emsp;&emsp;到目前为止，您应该开始看到shell通配符的效用：它们允许我们轻松地处理多个文件。由于大量日常生物信息学工作涉及文件处理，以编程方式访问文件使我们的工作变得更容易，并消除了因错误键入文件名或忘记样本而造成的错误。然而，只有当我们的文件名一致时，我们才能以编程方式使用通配符(或R或Python中的其他方法)访问文件。虽然通配符很强大，但如果文件命名不一致，它们是无用的。例如，处理名称为zmays sampleA-1.fastq、zmays_sampleA-2.fastq、sampleB1.fastq、sample-B2.fastq的文件子集是不必要的，因为这些文件名不一致。不幸的是，命名不一致的现象在整个生物学中普遍存在，并且这是世界各地生物信息学家的祸害。总的来说，生物信息学家可能已经浪费了数千个小时来对抗他人对文件、基因和代码的糟糕命名方案。

![suggestion](../../img/suggestion.png)**前导零和排序**

&emsp;&emsp;另一个有用的技巧是在命名文件时使用前导零(例如，file-0021.txt而不是file-21.txt)。这很有用，因为按照字典顺序对文件进行排序(就像ls一样)会导致正确的排序。例如，如果我们有文件名，例如gene-1.txt，gene-2.txt，…，gene-14.txt，按词典顺序对这些进行排序将产生：

```shell
$ ls -l
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:24 genes-1.txt
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:24 genes-11.txt
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:24 genes-12.txt
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:24 genes-13.txt
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:24 genes-14.txt
[...]
```

&emsp;&emsp;但是如果我们使用前导零(例如gene-001.txt，gene-002.txt，…，gene-014.txt)，文件将按正确的顺序排序：

```shell       
$ ls -l-rw-r--r-- 1 vinceb staff 0 Feb 21 21:23 genes-001.txt
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:23 genes-002.txt
[...]
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:23 genes-013.txt
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:23 genes-014.txt
```

&emsp;&emsp;使用前导零不仅在命名文件名时有用；这也是命名基因、转录本等的最好方法。像Ensembl这样的项目在命名他们的基因时使用这种命名方案(例如，ENSG00000164256)。
      
&emsp;&emsp;除了简化对文件的处理之外，一致的文件命名也是健壮的生物信息学中经常被忽略的组成部分。不好的样本命名方案很容易导致切换样本。当您或合作者认为您使用的是正确的数据时，选择不当的文件名也会导致严重的错误，但实际上它已经过时或错误的文件。我保证，在过去十年发表的所有论文中，由于文件命名问题，至少有一些，很可能更多的论文含有错误的结果。

### 用Markdown记录项目文档
      
&emsp;&emsp;记录一个项目笔记非常重要，其中包含有关计算工作年表的详细信息，所采取的步骤，关于您做出决定的原因的信息，当然还有复制您的工作所需的所有相关信息。一些科学家在手写的笔记本上这样做，其他的在微软的Word文档中。与README文件一样，生物信息学家通常喜欢将项目笔记本保存在简单的纯文本中，因为这些笔记本可以从命令行读取、搜索和编辑，也可以通过网络连接到服务器。纯文本也是一种面向未来的格式：20世纪60年代编写的纯文本文件今天仍然可读，而文字处理程序中只有10年历史的文件可能很难打开和编辑或不可能打开和编辑。此外，纯文本项目笔记本也可以置于版本控制之下，我们将在第5章中讨论这一点。
      
&emsp;&emsp;虽然纯文本很容易在文本编辑器中编写，但对于不熟悉命令行的协作者来说，阅读起来可能会很不方便。一种称为Markdown的轻量级标记语言是一种纯文本格式，易于阅读，并且可以轻松地合并到键入的笔记中，并且还可以转换为HTML或PDF格式。
      
&emsp;&emsp;Markdown源自纯文本电子邮件中使用的简单格式约定。早在HTML进入电子邮件之前，电子邮件就用简单的标记来修饰强调、列表和文本块。随着时间的推移，这成为事实上的纯文本电子邮件格式化方案。这个方案非常直观：文本两侧的下划线或星号表示强调，列表只是以破折号开头的文本行。
      
&emsp;&emsp;Markdown只是纯文本，这意味着它是可移植的，并且将存在编辑和读取它的程序。任何用旧版本的文字处理器写过笔记或论文的人都可能熟悉尝试共享或更新过时的专有格式的麻烦。由于这些原因，Markdown成为一种简单而优雅的笔记本格式。

### Markdown格式基础知识
       
&emsp;&emsp;Markdown的格式化功能与生物信息学笔记本的所有需求相匹配：文本可以分解为层次结构部分，代码块和内联代码都有语法，并且很容易嵌入链接和图像。虽然Markdown格式非常简单，但有几个不同的变体。在我们的示例中，我们将使用John Gruber发明的原始Markdown格式。John Gruber的完整的[降价语法规范](https://daringfireball.net/projects/markdown/syntax)可以在他的网站上找到。下面是说明格式的基本Markdown文档：

```shell      
# *Zea Mays* SNP Calling
We sequenced three lines of *zea mays*, using paired-end
sequencing. This sequencing was done by our sequencing core and we
received the data on 2013-05-10. Each variety should have **two**
sequences files, with suffixes `_R1.fastq` and `_R2.fastq`, indicating
which member of the pair it is.

## Sequencing Files

 All raw FASTQ sequences are in `data/seqs/`:
 $ find data/seqs -name "*.fastq"
 data/seqs/zmaysA_R1.fastq
 data/seqs/zmaysA_R2.fastq
 data/seqs/zmaysB_R1.fastq
 data/seqs/zmaysB_R2.fastq
 data/seqs/zmaysC_R1.fastq
 data/seqs/zmaysC_R2.fastq
 
## Quality Control Steps

After the sequencing data was received, our first stage of analysis
was to ensure the sequences were high quality. We ran each of the
three lines' two paired-end FASTQ files through a quality diagnostic
and control pipeline. Our planned pipeline is:
1. Create base quality diagnostic graphs.
2. Check reads for adapter sequences.
3. Trim adapter sequences.
4. Trim poor quality bases.

Recommended trimming programs:
 - Trimmomatic
 - Scythe
```

&emsp;&emsp;下面这些内容显示了由PandocHTML5呈现的这个示例Markdown笔记本。Markdown成为实验室笔记本的一种很好的格式是，它在未呈现的纯文本中的阅读就像在呈现的HTML中一样容易。接下来，让我们看一下本例中使用的Markdown语法(参见表2-2以获得参考)。
      Markdown展示内容 //todo
      
表2-2、最小Markdown内联语法
  
|Markdown语法|结果|
|:----|:--------|
|*强调*|强调|
|*粗体*|粗体|
|`inline code`|inline code|
|<http://websitecom/link>| Hyperlink to http://website.com/link|
|[link text](http://websitecom/link)|Hyperlink to http://website.com/link|
|![alt text](path/to/figure.png)|Image, with alternative text “alt text”|

&emsp;&emsp;像标头、列表和代码块这样的块元素很简单。可以使用不同数量的哈希值(#)指定不同级别的标头。例如：

```shell
# Header level 1
## Header level 2
### Header level 3
```

&emsp;&emsp;Markdown支持最多六个级别的标题。还可以对最多两个级别的标头使用替代语法：

```shell 
Header level 1
==============
Header level 2
--------------
```

&emsp;&emsp;有序列表和无序列表都很容易。对于无序列表，使用破折号、星号或加号作为列表元素标记。例如，使用破折号：

```shell
- Francis Crick
- James D. Watson
- Rosalind Franklin
```

&emsp;&emsp;要对列表进行排序，只需使用数字(即1.、2.、3.等)。顺序并不重要，因为HTML会自动为您递增这些值。然而，清晰地给要点编号仍然是一个好主意，这样纯文本版本就可以阅读了。
代码块也很简单-只需在每个代码行之前添加四个空格或一个制表符：

```shell
我运行了以下命令：
$find seqs/-name“*.fast”
```

&emsp;&emsp;如果要在列表项中放置代码块，请将这八个空格或两个制表符设置为：

```shell
1. 我使用以下命令搜索所有FASTQ文件: find seqs/ -name "*.fastq"
2. 最后，VCF文件具有:
 find vcf/ -name "*.vcf"
 ```
      
&emsp;&emsp;我们在这里介绍的内容应该足以让你开始用数字方式记录你的生物信息学工作。此外，还有Markdown的扩展，例如[MultiMarkdown](https://fletcherpenney.net/multimarkdown/)和[GitHub](https://help.github.com/en/categories/writing-on-github)风格的Markdown。这些变体添加了功能(例如，MultiMarkdown添加表、脚注、LaTeX数学支持等)，并更改了一些默认渲染选项。如果您喜欢基于GUI的编辑器，也有许多专门的Markdown编辑应用程序(有关一些建议，请参阅本章GitHub上的自述文件)。

### 使用Pandoc将Markdown呈现为HTML
     
&emsp;&emsp;我们将使用Pandoc，一个流行的文档转换器，将Markdown文档呈现为有效的HTML。然后，这些HTML文件可以与合作者共享或托管在网站上。有关如何在您的系统上安装Pandoc的说明，请参阅[Pandoc](https://pandoc.org/installing.html)安装页面。
     
&emsp;&emsp;Pandoc可以在各种不同的标记和输出格式之间进行转换。使用Pandoc非常简单-要将Markdown转换为HTML，请使用-from mark down和-to html选项，并提供输入文件作为最后一个参数：
 
```shell
$ pandoc --from markdown --to html notebook.md &gt; output.html
```
      
&emsp;&emsp;默认情况下，Pandoc将输出写入标准输出，我们可以将其重定向到文件(我们将在第3章了解更多关于标准输出和重定向的信息)。我们还可以使用-output output.html指定输出文件。最后，请注意，Pandoc可以在多种格式之间进行转换，而不仅仅是Markdown到HTML。我在本章关于GitHub的自述文件中包含了更多的例子，包括如何从HTML转换为PDF。