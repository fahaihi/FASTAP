# 软件：基于Python的生物信息数据处理与序列比对软件
# 当前版本：V【2.1】
# IDE：JetBrains PyCharm Community Edition 2019.1.3 x64
# 下载地址：
# 联系方式：adairmillersh@gmail.com
# 基础使用说明：
![](https://github.com/fahaihi/FASTAP/blob/master/1.JPG)
   一：模块简介
   1.1：FASTA文件格式处理
            FASTA格式是一种基于文本用于表示核酸序列或多肽序列的格式。其中核酸
       或氨基酸均以单个字母来表示，且允许在序列前添加序列名及注释。如下为格式
       信息举例：
       -----------------------------------------------------------------------
       >X02662.1 E. coli gap gene for GAPDH (glyceraldehyde-3-phosphate dehydr
       ogenase)GATCAAACAGTGATATACGCCGTCACGCTTGTTATGCAGTAAACGACCCGTAAATGGCGGCTC
       TGTCCCA.......
       -----------------------------------------------------------------------
       转换文件为TXT格式，文件路径在当前软件运行路径。


   1.2：FASTQ文件格式处理
            astQ格式是序列格式中常见的一种，FASTQ格式的序列一般都包含有四行，
       第一行由'@'开始，后面跟着序列的描述信息，这点跟FASTA格式是一样的。第二行
       是序列。第三行由'+'开始，后面也可以跟着序列的描述信息。第四行是第二行序
       列的质量评价（quality values，注：应该是测序的质量评价），字符数跟第二行
       的序列是相等的。例如在NCBI看到的FASTQ格式如下：
       -----------------------------------------------------------------------
       @HWUSI-EAS100R:6:73:941:1973#0/1
       GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTT
       +HWUSI-EAS100R:6:73:941:1973#0/1
       !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC6
       -----------------------------------------------------------------------
       转换文件为TXT格式，文件路径在当前软件运行路径。


   1.3：简单序列比对
              简单序列比对就是最简单的匹配算法，有叫做天真匹配，其匹配速度和
       Boyer-Moore以及KMP相比，相差甚远。但是由于其简单，易于实现的特点，可以
       用在比对体量比较小的场合，如下展示了简单匹配函数。
       -----------------------------------------------------------------------
       def naive(p,t):
          occurr ents = []        #创建一个跟踪p和t发生匹配的列表
          for i in range(len(t) - len(p) + 1):
            match = True
          for j in range(len(p)):
            if not t[i+j] == p[j]:
                match = False
                break #因为存在一个不匹配说明此匹配项一定不存在
          if match:
            occurrents.append(i)
       return occurrents
       -----------------------------------------------------------------------


   1.4：Boyer-Mooer序列比对
           假设将主串中自位置i起往左的一个子串与模式进行从右到左的匹配过程中，
      若发现不匹配，则下次应从主串的i + dist(si)位置开始重新进行新一轮的匹配，
      其效果相 当于把模式和主串向右滑过一段距离distance（si），即跳过distance
      （si）个字符而无需进行比较。如下展示了其匹配函数。
       -----------------------------------------------------------------------
       def BoyerMoore(pattern, m, text, n):
          bmBc = PreBmBc(pattern, m)
          bmGs = PreBmGs(pattern, m)
          j = 0
          flag = 0
          results = []
          while j <= n - m:
          for i in reversed(range(m)):

            if pattern[i] == text[i + j]:
                if i <= 0:
                    flag += 1
                    results.append(i+j)#" "*j + pattern + " "*(n - m - j) + "  "
                    + str(j) + "  " + str(flag))
                    j += bmGs[0]
            else:
                j += max(bmBc[ord(text[i + j])] - m + 1 + i, bmGs[i])
                break
       return results

       -----------------------------------------------------------------------

   二：Interface
    本软件预留的功能接口，其功能暂未开放，用户可以暂时忽略相关模块。

