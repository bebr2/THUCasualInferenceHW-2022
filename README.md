#	THU-因果推断导论大作业

## 简介

本仓库是2022年秋季学期清华大学《因果推断导论》课程中，本人的大作业实现。

### 研究过程

本文主要研究废除死刑对凶杀案数量的影响，研究对象是2010年和2014年的美国各州。本人首先在[Kaggle平台](www.kaggle.com)和[美国人口普查局网站](www.census.gov)收集了美国各州的凶杀案数据，以及各州的经济、人口、种族、教育、医疗等数据，同时在维基百科收集了废死州和未废死州的分布。由于数据缺失等原因，最后选取了2010年和2014年的美国各州为个体。

简单线性回归发现，废死 ~ 凶杀案数量几乎没有任何关系，但同时观察数据分布可以看出废死和非废死州之间的协变量非常不均衡。

由于数据是观察性数据，选取Rubin框架进行研究，并作出一些假设。接着用逻辑回归的方法平衡倾向得分，使用trimming的方法对得到的倾向得分进行处理，接着采用无放回匹配和有放回匹配的方法分别探究，同时也研究了是否加入协变量修正的影响。

最后，对于上述方法，得到的4个平均因果作用都非常接近于0，方差也在合理范围内，因此本文得出的结论是废死对凶杀案数量没有影响。

研究细节均在正文展示。

### 免责声明

本文结论对相关学术领域无任何正面或负面的结果，本文无意传达任何地域、种族歧视或其他不良内容。仅供因果推断的入门学习者的实践参考。

关于本文的研究局限性，已经在正文文末提出，这里不再赘述。

## 文件结构

pdf文件是正文，`Rcode.md`文件是正文相关结论的R代码实现（不含输出结果），`data`文件夹下包含raw data和processed data，以及处理代码。



