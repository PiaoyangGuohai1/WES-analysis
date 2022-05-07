# fastp软件详细介绍

### 基本介绍

对于下机的FASTQ数据需要进行质控和预处理，以保证下游分析输入的数据都是干净可靠的。通常我们都是使用FASTQC等软件进行质控，使用cutadapt软件去除接头，使用Trimmomatic等软件进行剪裁，然后使用一些自已开发的脚本进行过滤。这一过程可能需要使用多个软件，相当繁琐，而且速度较慢，这些软件大多又不支持多线程，遇到较大的FASTQ文件，处理起来可真是让人等得心急如焚。



fastp是一款数据质控过滤软件，作者是陈实富，来自深圳海普洛斯公司。他们将这款工具开源免费使用，这一点是非常值得称赞的。

1、fastp可以实现处理数据的一次性处理，包括过滤低质量，过滤adapter，截取reads，split分割大文件等操作

2、支持长reads，也就是不仅仅适用与illumina测序平台，还可以处理Pacbio和Ion torrent的测序数据

3、直接输出质控和统计报告，包括json格式和html格式；

4、使用c++写的，执行效率非常高。



其github地址：https://github.com/OpenGene/fastp



### 使用方法

```shell
## 单末端测序数据，非压缩格式
fastp -i in.fq -o out.fq

## 双末端测序数据，gzip压缩格式
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

## 使用案例： 
fastp -i reads.1.fq.gz -I reads.2.fq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z 4 -q 20 -u 30 -n 
```



### 主要参数

```
## I/O 相关
-i, --in1    输入read1文件名  
-o, --out1   输出read1文件名
-I, --in2    输入read2文件名  
-O, --out2   输出read2文件名，软件默认是根据扩展名识别压缩文件，所以输出文件需要加上*.gz扩展名；  
-6, --phred64  指定质量体系是phred64。目前主流测序数据都采用phred33，如果从NCBI下载以前hiseq 2000以及之前的数据，可能是Phred 64质量体系。  
-z, --compression 输出压缩格式。给定一个数字1-9，调整压缩比率和效率的平衡；  

## adapter相关选项
-A 关闭adapter trimming，默认软件会切出adapter，如果设置-A，则关闭这个功能；  
-a 给定一个adapter序列文件；

## 全局裁剪选项
-f, --trim_front1  裁剪read1前多少个碱基，默认0；
-F, --trim_front2  裁剪read2前多少个碱基，默认0，如果没有指定，将保持与read1相同设置；
-t, --trim_tail1   裁剪read1末尾多少个碱基，默认0；
-T, --trim_tail2   裁剪read2末尾多少个碱基，默认0，如果没有指定，将保持与read1相同设置；
-b, --max_len1     如果read1比max_len1长，就从末尾截取一段使read1与max_len1等长，默认0，代表没有限制；
-B, --max_len2     如果read2比max_len2长，就从末尾截取一段使read2与max_len2等长，默认0，代表没有限制，如果没有指定，就保持与read1相同设置；

## polyG尾裁剪，针对NextSeq/NovaSeq数据
-g, --trim_poly_g  截取polyG尾，对于Illumina NextSeq/NovaSeq测序数据时默认自动操作的；
    --poly_g_min_len 检测read末尾的polyG的长度，默认10；
-G, --disable_trim_poly_g 取消-g的功能；

## polyX tail trimming
-x, --trim_poly_x  截取3'末端polyX
    --poly_x_min_len    检测read末尾的polyX的长度，默认10；

## 通过质量值对每条read进行裁剪
-5, --cut_front  从read的5'端至末尾移动窗口，去除窗口中平均质量值小于'<'阈值的碱基；
-3, --cut_tail   从read的3'端值至开头移动窗口，去除窗口中平均质量值小于'<'阈值的碱基；
-r, --cut_right  从read的开头到末尾移动窗口，如果某一窗口的平均质量值小于阈值，去除窗口中的碱基及其右侧部分，并停止；
-W, --cut_window_size  滑动窗口过滤，这个类似于计算kmer，1~1000, 默认是4个碱基；
-M -W选择的窗口中，碱基平均质量值，范围1~36，默认是Q20，如果这个区域窗口平均低于20，则认为是一个低质量区域，处理掉；

## 质量过滤选项
-Q 控制是否去除低质量，默认自动去除，设置-Q关闭；
-q 设置低质量的标准，默认是15，也就是质量值小于15认为是低质量碱基，一般我们设置20，常说的Q20；
-u 低质量碱基所占百分比，并不是包含低质量碱基就把一条reads丢掉，而是设置一定的比例，默认40代表40%，也就是150bpreads，包含60个以上低质量的碱基就丢掉，只要有一条reads不满足条件就成对丢掉；
-n 过滤N碱基过多的reads，如果N碱基含量大于n，这条read/pair将被舍弃，默认5；

## 长度过滤选项
-L 关闭reads长度过滤选项；
-l 接一个长度值，小于这个长度reads被丢掉，默认是15，这个在处理非illumina测序数据时很有用。

## 低复杂度过滤
-y, --low_complexity_filter    使用低复杂度过滤，这里低复杂度的定义是与其下一个碱基不同的碱基比例(base[i] != base[i+1]).
-Y, --complexity_threshold    低复杂度的阈值(0~100)，默认30；

## 根据indexes过滤reads--删除可能的污染
--filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
--filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
--filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])

## base correction by overlap analysis options
-c 是对overlap的区域进行纠错，所以只适用于pairend reads。

## UMI processing 分子标签处理
-U, --umi                        enable unique molecular identifier (UMI) preprocessing
  --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
  --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])
  --umi_prefix                   if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
  --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])

# overrepresented sequence analysis
-p, --overrepresentation_analysis    enable overrepresented sequence analysis.
-P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])

# reporting options
-j, --json        输出json格式报告文件名(string [=fastp.json])
-h, --html        输出html 格式报告文件名，可以用浏览器直接查看(string [=fastp.html])
-R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])

# threading options
-w, --thread     使用线程数，默认是2(int [=2])

# 控制split选项，有时候单条reads文件太大，可以分割为多份分别比对，在合并bam结果，这样可以提高效率。
-s, --split      切割数目(2~999)，默认是0，不分割
-S, --split_by_lines  split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
-d, --split_prefix_digits    输出前缀位数，默认是4，0001,0002这种命名，如果设置为3，就是001,002这种；

# help
-?, --help       输出帮助信息
```



### 基本功能简介

- 过滤
  **fastp可以对低质量序列，较多N的序列，该功能默认是启用的**，但可以使用-Q参数关闭。使用-q参数来指定合格的phred质量值，比如-q 15表示质量值大于等于Q15的即为合格，然后使用-u参数来指定最多可以有多少百分比的质量不合格碱基。比如-q 15 -u 40表示一个read最多只能有40%的碱基的质量值低于Q15，否则会被扔掉。使用-n可以限定一个read中最多能有多少个N。
  **fastp还默认启用了read长度过滤**，但可以使用-L参数关闭。使用-l参数指定最低要求一个read有多长，比如-l 30表示低于30个碱基的read会被扔掉。这个功能可以用于实现常用的discard模式，以保证所有输出的序列都一样长。
  在fastp的HTML报告中，最头上的Summary表格很清楚地显示了过滤的统计信息，
- 接头处理
  接头（adapter）污染的处理是FASTQ文件预处理中很重要的一步。**fastp默认启用了接头处理**，但是可以使用-A命令来关掉。fastp可以**自动化地查找接头序列并进行剪裁，也就是说你可以不输入任何的接头序列**，fastp全自动搞定了！**对于SE数据，你还是可以-a参数来输入你的接头，而对于PE数据则完全没有必要**，fastp基于PE数据的overlap分析可以更准确地查找接头，去得更干净，而且对于一些接头本身就有碱基不匹配情况处理得更好。fastp对于接头去除会有一个汇总的报告。
- 滑窗质量剪裁
  很多时候，一个read的低质量序列都是集中在read的末端，也有少部分是在read的开头。fastp支持像Trimmomatic那样对滑动窗口中的碱基计算平均质量值，然后将不符合的滑窗直接剪裁掉。使用-5参数开启在5’端，也就是read的开头的剪裁，使用-3参数开启在3’端，也就是read的末尾的剪裁。使用-W参数指定滑动窗大小，默认是4，使用-M参数指定要求的平均质量值，默认是20，也就是Q20。
- PE数据的碱基校正
  fastp支持对PE数据的每一对read进行分析，查找它们的overlap区间，然后对于overlap区间中不一致的碱基，如果发现其中一个质量非常高，而另一个非常低，则可以将非常低质量的碱基改为相应的非常高质量值的碱基值。**该校正功能默认没有开启**使用-c参数可以启用，**对于一些对噪声容忍度低的应用，比如液体活检，建议开启。**
- 全局剪裁
  fastp可以对所有read在头部和尾部进行统一剪裁，该功能在去除一些测序质量不好的cycle比较有用，比如151*2的PE测序中，最后一个cycle通常质量是非常低的，需要剪裁掉。使用-f和-t分别指定read1的头部和尾部的剪裁，使用-F和-T分别指定read2的头部和尾部的剪裁。
- polyG剪裁
  对于两色发光法的Illumina设备（NextSeq/NovaSeq），因为在没有光信号情况下base calling的结果会返回G，所以在序列的尾端可能会出现较多的polyG，需要被去除。fastp会自动化地识别NextSeq/NovaSeq的数据，然后进行polyG识别和剪裁。如果你想强制开启该功能，可以指定-g参数，如果想强制关闭该功能，则可以指定-G参数。
- 分子标签UMI处理
  UMI在处理ctDNA类似的超低频突变检测应用中是十分有用的，为了更好地对带UMI的FASTQ文件进行预处理，fastp也很好地支持了UMI预处理功能。该功能默认没有启用，需要使用-U参数开启，另外需要使用--umi_loc来指定UMI所在的位置，它可以是（index1、index2、read1、read2、per_index、per_read ）中的一种，分别表示UMI是在index位置上，还是在插入片段中。如果指定了是在插入序列中，还需要使用--umi_len参数来指定UMI所占的碱基长度。
- 输出文件切分
  很多时候我们需要对输出的FASTQ进行切分，分成大小均匀的多个文件，这样可以使用比对软件并行地比对，提高并行处理的速度。fastp软件也提供了相应的功能，并且支持了两种模式，分别是使用参数-s指定切分后文件的个数，或者使用-S参数指定每个切分后文件的行数。



### 质控报告解读

fastp的报告在单一文件中同时包含了过滤前和过滤后的统计结果，**如果是PE数据，则同时包含了read1和read2的统计结果**。fastp会生成HTML的报告和JSON格式的报告。HTML报告的默认文件名是fastp.html，但是可以通过-h参数修改，JSON报告的默认文件名是fastp.json，但是可以通过-j参数修改。而且fastp报告还有一个标题，默认是fastp report，这个也可以通过-R参数修改为你想要的标题。JSON格式的报告是优化过的，人机皆可读，适合进阶的用户使用程序解析，而这里我们重点关注HTML格式的报告。

- **质量分布曲线图**
  我们第一关注的当然是质量，所以**fastp提供了质量分布曲线**，即每一个cycle的平均质量值，而且fastp同时提供了A/T/C/G四种不同碱基的平均质量，以及总的平均质量图，从图中我们可以看到，一共有5条曲线，分别是A/T/C/G和mean。而且HTML报告中的每一个项目和分项目都是可以点击进行隐藏和展开的。
- **碱基含量分布曲线**
  和质量分布曲线类似，碱基含量分布曲线也是按照每一个cycle来的，显示了每一个位置的碱基含量。从图中可以看到，fastp同时显示了A/T/C/G/N/GC的每一个位置的比例和总的比例。而且如果你觉得头部那里比较乱看不清的话，可以用鼠标拉一个框，它就放大了。
- **KMER统计表格**
  fastp对5个碱基长度的所有组合的出现次数进行了统计，然后把它放在了一张表格中，表格的每一个元素为深背景白字，**背景越深，则表示重复次数越多**。这样，一眼望去，就可以发现有哪一些异常的信息。比如，从KMER表格中，我们可以发现，GGGGG的颜色特别深，从鼠标移上去之后显示的信息中我们可以发现它的出现次数是平均次数的12.8倍，这是不正常的，因为GGGGG的正常倍数应该在1倍左右。幸好我们有fastp，它可以过滤掉这种polyG，让数值较多地回归正常。
- **过表达序列(overrepresented sequence)**
  fastp提供了overrepresetned sequence的分析，而且不但提供了这些overrepresented sequence的序列个数和占比，还提供了他们在测序cycles中的分布情况，这有利于分析各种问题。



参考：https://www.cnblogs.com/dataanaly/p/13197991.html
