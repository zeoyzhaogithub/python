# 化学信息学-smarts规则速查表

@Time: 2020-06-17 11:07:11
@Author:
化学信息学-smarts规则速查表

## 一、SMARTS简介

### 1.1 什么是SMARTS

SMARTS（SMiles ARbitrary Target Specification）是一种用于描述分子模式和属性的语言。SMILES所有的符号和属性在SMARTS中同样适用，因此它也是SMILES的延伸。此外，SMARTS还包括了逻辑操作符和额外的分子描述符，后文会一一介绍。

### 1.2 SMARTS能做什么

SMARTS可以从不同程度来概括和表示结构模式。举个例子：
甲烷的SMILES可以用"C"或"[CH4]“表示。
而”[CH4]"在SMARTS中，高度特异地表示与甲烷一致的结构，即只能匹配一个带有4个氢的脂肪族碳原子。
而"C"在SMARTS中特异程度较低，可以表示带有任意数量氢的脂肪族碳原子，比如乙烷、乙烯、环戊烷。

## 二、原子属性

SMARTS | 匹配结构 | 说明
-------|---------|-------
[+1]|带有一个正电荷的原子|SMILES对电荷、氢、同位素、键、手性等描述方式，在SMARTS中都可以兼容。一个"+“表示”+1"，两个"++“表示”+2"|
[a]|带有芳香性的原子|"a"表示任何带有芳香性质的原子|
[A]|带有脂肪族属性的原子|"A"表示任何带有脂肪族性质的原子|
[  # 6]|原子序数为6的原子（c或C）|"#<number>"表示序数为<number>的原子，无论是脂肪族还是芳香族|
[R2]|在两个环中的原子|"R<number>"表示在 < number > 元环中的原子，默认{R}为在任何环中的原子|
[r5]|在五元环中的原子|"r<number>"表示 < number > 元环中的原子|
[v4]|四价原子|"v<number>“表示任何含有键的数量为<number>的原子。另外”=“表示双键，”#"表示三键|
[X2]|与两个原子链接的原子|"X<number>"表示任何与 < number > 个原子相连的原子（包括氢原子）|
[H]|氢原子|一个氢原子（通常也叫一个显式氢，explicit hydrogen）具有一些特殊的性质[H+], [2H], [H][H]等。[H+]和[2H]含义相似|
[H1]|与一个氢相连的原子|"H<number>"表示任何与 < number > 个氢（显式或隐式氢）相连的原子。[*H]表示没有氢相连的原子|
[D3]|三个显式键的原子|"D<number>"表示任何与 < number > 个重原子相连的原子（与氢相连的不算）
*| 任何原子| "*"表示通配原子，匹配任意重原子（非氢原子）|

## 三、键属性

SMARTS|匹配结构|说明
------|--------|--------
CC|两个由单键相连接的脂肪碳|所有SMILES的键的属性在SMARTS中都可以使用，包括隐式单键、显式单键（-）、双键（=）、三键（  # ）、芳香键（:）
[  # 6]~[#6]|两个由任意键相连的碳|"~"表示通配键
[  # 6]@[#6]|两个在同一个环中相连的碳|"@"表示在同一个环中
[F /?[  # 6]=C/Cl]|氟原子通过"/"（"/“指定了顺反异构构型）或未指明的键与碳原子相连（比如可以匹配到"F/C=C/Cl"或"FC=C/Cl”，不能匹配到"F\C=C/Cl"）|“?“表示"或不确定”，还可以和手性描述符”@"一起使用|

## 四、逻辑操作符

SMARTS|匹配结构|说明
------|-------|---------
!c|非芳香的碳|"!“表示"非”
[N,  # 8]|匹配脂肪族氮或匹配一个氧|“,“表示"或”，优先级高于"与”（";"），低于另一个"与"（"&"）
[  # 7,C&+O,+1] or [#7,C+O,+1]|氮原子或中性脂肪碳原子或带一个正电荷的原子|"&“表示"与”（优先级高），是默认的逻辑操作符，可以省略
[  # 7,C;+0,+1]|氮或脂肪碳，且不带或带一个正电荷|";“也表示"与”，但优先级低

## 五、递归SMARTS

SMARTS|匹配结构|说明
------|-------|---------
[$(*O); $(*CC)]|一个与脂肪氧相连的原子，或一个连接有两个脂肪碳的原子|"$<SMARTS>"表示匹配周围具有某种结构的原子
[$([CX3]=[OX1]), $([CX3+]-[OX1-])]|与一个羰基相连的原子，或相对合理的结构|
[$([A]aaO); $([A]aaaN)]|芳环上位于氧的邻位，氮的间位的原子

## 六、组合匹配

SMARTS|匹配结构|说明
------|-------|---------
[  # 8].[#8]|匹配两个氧（例如O=O, OCCO, O.CCO）|"."表示无需连接
([  # 8].[#8])|在同一个结构中匹配两个氧（例如O=O, OCCO，无法匹配O.CCO）|可以在SMARTS外加圆括号，表示括号内的结构需要在同一组分中出现|
([  # 8]).([#8])|在不同的结构中匹配两个氧（例如O.CCO，无法匹配O=O, OCCO）|可以使用多个圆括号，表示需要在不同的组分中进行匹配

## 七、反应SMARTS

SMARTS|匹配结构|说明
------|-------|---------
[  # 6]=,:[#6]|由一个双键或芳香键连接的碳|分子SMARTS（没有">"符号）可以对任意反应组分（反应物、试剂或产物）中进行匹配|
>>[  # 6]=,:[#6]|产物中碳由一个双键或芳香键连接|反应SMARTS（带有">"符号）不能用于分子的匹配
[C:1] >> [C:1]|加入":"和数字表示对原子进行标记和映射|原子如果做了标记，就会与对应的被标记的原子匹配，不会匹配到没有标记的原子上
[C:1] >> C|反应碳|标记的原子无法完全配对时，将会忽略|
[C:1][C:1] >> [C:1]|多次标记的反应碳|
参考DAYLIGHT的SMARTS介绍。
英文版在这里[英文版在这里]('https://www.daylight.com/dayhtml_tutorials/languages/smarts/index.html)。

原文件在这里