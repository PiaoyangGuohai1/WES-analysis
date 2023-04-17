# WES分析全流程

## 0、整体文件结构

创建文件夹

```shell
mkdir 00_ref # 参考数据文件夹
mkdir 01_fastp # 质控及质控报告
mkdir 02_bwa_out # 比对
mkdir 03_markdup # 标记重复
mkdir 04_bqsr # 重新校准碱基质量值
mkdir 05_gvcf # gvcf文件
mkdir 06_joint_genotype # joint-genotype
mkdir 07_VQSR # 变异质控
mkdir 08_annovar # 位点注释
mkdir tmp # 临时文件存放，gatk过程中会用到
```

- 00_ref
  - bed
    - `interval.list`：存放测序时外显子位置信息的bed文件，一般从测序公司获得。
  - gatk_call_vcf
  - ucsc-human-hg38
- 01_fastp
- 02_bwa_out
- 03_markdup
- 04_bqsr
- 05_gvcf
- 06_joint_genotype

从rawdata处获取所有样本名：

```shell
# 指定fasta文件路径
path=~/Rawdata/AAD_WGS_20220127/GHDLB20070348 # rawdata位置
files=$(ls $path | grep "\\_1.fq.gz")
for filename in $files
do
   echo $filename | awk -F"_" '{print $1;exit}'>> ./sample_id.list
done
```



以下流程全部以单样本分析为案例，多样本分析查看snakemake流程。

## 1、质控分析

### 质控目的：

过滤掉低质量、太短以及太多N的reads（注释：N表示为测到的位点） 

（fastp不需要提供接头序列就可以自动识别并去除接头）

### fastp软件介绍：

fastp一款超快速全功能的FASTQ文件自动化质控+过滤+校正+预处理软件

- 拥有FASTQC + cutadapt + Trimmomatic加起来的功能；
- 支持多线程，因为它使用C++开发，速度快；
- 对overlap区域碱基质量值自动校正；
- 软件的使用非常简单，默认情况下只需要指定输入和输出文件，就可以很好地工作。

### 基础代码：

```shell
fastp -i in.R1.fq.gz -o out.R1.fq.gz -I in.R2.fq.gz -O out.R2.fq.gz
```



### 实操代码：

```shell
path=/share/home/longxinyang/Rawdata/AAD_WGS_20220127/GHDLB20070348 # 原始数据位置
fastp \
-i ${path}/AD2021001_1.fq.gz -o 01_fastp/AD2021001_1.clean.fq.gz \
-I ${path}/AD2021001_2.fq.gz -O 01_fastp/AD2021001_2.clean.fq.gz \
-z 4 -f 5 -t 5 -F 5 -T 5 -5 -W 5 -M 20 -Q -l 50 -c -w 10 \
-j 01_fastp/AD2021001.fastp.json \
-h 01_fastp/AD2021001.fastp.html
```



### 参数说明：

[[[fastp质控报告解读]]](https://github.com/PiaoyangGuohai1/WES-analysis/blob/main/fastp%E8%B4%A8%E6%8E%A7%E6%8A%A5%E5%91%8A%E8%A7%A3%E8%AF%BB.md)

### 质控结果查看：

[fastp质控报告解读](https://github.com/PiaoyangGuohai1/WES-analysis/blob/main/fastp%E8%B4%A8%E6%8E%A7%E6%8A%A5%E5%91%8A%E8%A7%A3%E8%AF%BB.md)

## 2、比对，alignment



**序列比对**

NGS测序下来的短序列（read）存储于FASTQ文件里，reads随机打断，不同read之间的前后顺序关系就已经全部丢失，需要先把这一大堆的短序列，一个个去跟该物种的参考基因组比较，找到每一条read在参考基因组上的位置，然后按顺序排列好，这个过程就称为测序数据的比对。



### 第一步：下载基因组，构建索引

参考基因组的选择：[GRCh38/hg38](http://www.bio-info-trainee.com/1469.html)

选择使用来自[UCSC的ref](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest)

```shell
cd 00_ref/ucsc-human-hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz & # 下载参考基因组
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/md5sum.txt & # 下载md5校验文件
md5sum hg38.fa.gz # 检查数据完整性
gunzip hg38.fa.gz # 解压参考基因组
bwa index hg38.fa # 构建bwa比对所需的参考基因组的index数据库
samtools faidx hg38.fa # 创建fasta序列格式索引
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict # 生成参考基因组的dict文件
```



### 第二步：比对

将样本测序数据reads与人类参考基因组进行比对，并将输出文件转化为bam格式，可有效节省磁盘空间。

使用`bwa men`完成数据比对，`bwa men`对于任何长度大于40bp小玉200bp的read都非常有效。

```shell
mkdir 02_bwa_out
bwa mem -t 4 -M -R '@RG\tID:laneID\tPL:illumina\tSM:AD2021001' \
00_ref/ucsc-human-hg38/hg38.fa \
01_fastp/AD2021001_1.clean.fq.gz \
01_fastp/AD2021001_2.clean.fq.gz \
| samtools view -b -@ 4 - > 02_bwa_out/AD2021001.align.bam
```

**参数解析：**

* bwa mem：比对

  - -t，线程数

  - -M，将shorter split hits标记为次优，以兼容Picard markDuplicates软件

  - -R，接的是Read Group的字符串信息，以@RG开头，以\t分离各选项。这是用于将比对的read进行分组，对于后续对比对数据进行错误率分析和Mark duplicate时非常重要。
    - ID，这是Read Group的分组ID，一般设置为测序的**lane ID**（不同lane之间的测序过程认为是独立的），下机数据中我们都能看到这个信息的，一般都是包含在fastq的文件名中；
    - PL，指的是所用的测序平台，**这个信息不要随便写！**在GATK中，PL只允许被设置为（不区分大小写）：
      - ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，PACBIO，IONTORRENT，CAPILLARY，HELICOS
      - UNKNOWN：如果实在不知道，那么必须设置为UNKNOWN
    - SM，**样本ID，非常重要**，有时候我们测序的数据比较多的时候，那么可能会分成多个不同的lane分布测出来，这个时候SM名字就是可以用于区分这些样本；
    - LB，测序文库的名字，这个重要性稍微低一些，主要也是为了协助区分不同的group而存在。文库名字一般可以在下机的fq文件名中找到，如果上面的lane ID足够用于区分的话，也可以不用设置LB。




- samtools view，由于bwa mem输出文件为**SAM，占用大量的空间**。因此此处通过管道符"|"将bwa mem的输出作为samtools的输入，**生成BAM格式**的文件。
  - -b，输出格式为BAM文件
  - -@，线程数，一般设置为4，不可设置太高，反而影响服务器性能-，表示管道符引流的数据所存在的位置




### 第三步：排序

sort排序：BWA比对后输出的BAM文件是没顺序的，比对后得到的结果文件中，每一条记录之间位置的先后顺序是乱的，我们后续去重复等步骤都需要在比对记录按照顺序从小到大排序下来才能进行。

```shell
samtools sort \
-@ 4 -m 20G -O bam \
-o 02_bwa_out/sorted/AD2021001.sorted.bam 02_bwa_out/AD2021001.align.bam
```

- -@，线程数
- -m，内存数（注意，此处为给每个线程的内存数，不宜设置过大）
- -O，输出格式为bam文件
- -o，输出文件名
- 最后跟上：输入文件名



通过对比结果发现，**未排序的数据（5G）是排序后数据的2倍**



**同时进行比对+排序，结果与分开做一致：**

```shell
bwa mem -t 4 -M -R '@RG\tID:AD2021001\tPL:illumina\tSM:AD2021001' \
00_ref/ucsc-human-hg38/hg38.fa \
01_fastp/AD2021001_1.clean.fq.gz \
01_fastp/AD2021001_2.clean.fq.gz \
| samtools sort -@ 4 -m 20G -O bam \
-o 02_bwa_out/AD2021001.sorted.bam -
```



### 第四步：比对结果统计【可选项】

```shell
samtools view -h in.bam
samtools index in.bam  # 生成in.bam的索引文件in.bam.bai
samtools view in.bam chr22            # 跳转到chr22染色体
samtools view in.bam chr22:16050103   # 跳转到chr22:16050103位置
samtools view in.bam chr22:16050103-16050103  # 只查看该位置
```



## 3、标记重复序列：mark duplication

### 标记重复

```shell
gatk MarkDuplicates \
-I 02_bwa_out/AD2021001.sorted.bam \
-O 03_markdup/AD2021001.sorted.markdup.bam \
-M 03_markdup/AD2021001.markdup_matrics.txt
```

参数解析：

- -I，输入文件（排序好的bam文件）
- -O，输出标记好重复的bam文件
- -M，输出重复的统计信息

注意Picard标记的重复序列只是在BAM文件中的[FLAG信息](https://broadinstitute.github.io/picard/explain-flags.html)中标记出来（标记值为1024）而不删除，因此这些重复序列依然会被留在文件中，只是我可以在变异检测的时候识别到它们，并进行忽略。



### 查看被标记的reads【可选项】

```shell
samtools view -f 1024 03_markdup/AD2021001.markdup.bam | less
```



### 针对低map质量的reads进行过滤

可以先抽取一部分bam文件在IGV中进行查看，可以看到大量由于map质量低（MAPQ值）导致的空reads。

参考：[MAPQ](https://www.jianshu.com/p/9c87bba244d8)和[高通量比对的质量MAPQ](https://zhuanlan.zhihu.com/p/35495052)

过滤后，构建索引。

```
samtools view -bSq 10 03_markdup/AD2021001.sorted.markdup.bam > 03_markdup/AD2021001.map.filtered.bam
samtools index 03_markdup/AD2021001.map.filtered.bam
```



### 为BAM文件创建索引

作用：能够让我们可以随机访问这个文件中的任意位置，而且**后面的步骤也要求这个BAM文件一定要有索引**，完成后得到bam.bai文件

```shell
samtools index 03_markdup/AD2021001.sorted.markdup.bam
```



### 使用IGV查看bam文件

[肿瘤外显子数据处理系列教程（番外篇）bam文件载入igv可视化](https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247488497&idx=2&sn=ccee8550494e4aef673969857dfe7dc5&scene=21#wechat_redirect)



## 4、变异检测（GATK）

需要提前了解[GATK背景知识](https://github.com/PiaoyangGuohai1/WES-analysis/blob/main/GATK%E6%B5%81%E7%A8%8B%E5%AD%A6%E4%B9%A0%E8%83%8C%E6%99%AF%E7%9F%A5%E8%AF%86.md)

此处主要遵循：GATK call variant（snp/indel）流程

### 第一步：数据库下载

前往[GATK数据库](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)下载数据

下载完后将数据放入在GATK_call_vcf内（需要注意，有些vcf.gz的文件，在windows下载后，后缀变成了vcf，此时需要**手动修改文件名**） 从上述数据库中可以下载vcf文件及其tbi文件。除此之外，我们还需要为每个后续需要用到的文件构建索引（index）

```shell
cd ./00_ref/gatk_call_vcf
# 先解压（4个文件）
gunzip resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# 后构建索引（4个文件）
gatk IndexFeatureFile -I resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf
gatk IndexFeatureFile -I resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf
gatk IndexFeatureFile -I resources_broad_hg38_v0_hapmap_3.3.hg38.vcf
gatk IndexFeatureFile -I resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf
```

所有数据库文件如下图所示：

![image-20220506114513973](https://markdown-1300560293.cos.ap-guangzhou.myqcloud.com/markdown/202205061757274.png)

### 第二步：重新校准碱基质量值

假设原始数据中就存在着一些由于测序仪器产生的系统性误差，那么变异位点识别过程中找到的variant，就会存在大量的**假阳性**。这一步主要目的是**调整原始碱基的质量分数**。**变异检测极度依赖测序碱基质量值**，因为这个质量值是衡量我们测序出来的这个碱基到底有多正确的重要指标。它来自于测序图像数据的base calling，因此，基本上是由测序仪和测序系统来决定的，计算出来的碱基质量值未必与真实结果统一。**BQSR （Base Quality Score Recalibration）**这个步骤主要是通过机器学习的方法构建测序碱基的错误率模型，然后对这些碱基的质量值进行相应的调整。



这里包含了**两个步骤**：

- 第一步，BaseRecalibrator，这里计算出了所有需要进行重校正的read和特征值，然后把这些信息输出为一份**校准表文件**（wes.recal_data.table）
- 第二步，ApplyBQSR，这一步利用第一步得到的校准表文件（wes.recal_data.table）重新调整原来BAM文件中的碱基质量值，并使用这个新的质量值重新输出一份**新的BAM文件**。

```shell
# 指定参考数据文件位置
hg38_vcf=00_ref/gatk_call_vcf
hg38_ref=00_ref/ucsc-human-hg38/hg38.fa
bed=00_ref/bed/interval.list
```

#### step1：BaseRecalibrator

利用已有的snp数据库，建立相关性模型，**产生重校准表**(recalibration table)，输入已知的多态性位点数据库，用于屏蔽那些不需要重校准的部分。

```shell
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" BaseRecalibrator \
-R $hg38_ref \
-I 03_markdup/AD2021001.map.filtered.bam \
--known-sites $hg38_vcf/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
--known-sites $hg38_vcf/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
--known-sites $hg38_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-L $bed \
-O 04_bqsr/AD2021001.recal_data.table
```



#### step2：ApplyBQSR

根据这个模型对原始碱基进行调整，只会调整非已知SNP区域。该过程中结束后会生成bam文件及其bai索引文件。

```shell
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" ApplyBQSR \
-R $hg38_ref -L $bed \
-bqsr 04_bqsr/AD2021001.recal_data.table \ # 上一步产生的重校准表
-I 03_markdup/AD2021001.map.filtered.bam \
-O 04_bqsr/AD2021001.sorted.markdup.BQSR.bam
```

参数解析：

-L参数，--intervals，是指将分析限制在特定的间隔内，一般用于WES或者靶向测序时使用。

[什么时候需要使用这个参数？](https://gatk.broadinstitute.org/hc/en-us/articles/360035889551?id=4133)

在外显子组分析和其他靶向测序中，由于其数据不涵盖整个参考基因组，需要提供intervals lists，不仅可以**排除靶外的噪声**，还能够**加速分析**，启用并行。但在全基因组是，不需要严格限制intervals，除非你希望屏蔽掉已知的不可靠或混乱区域（如着丝粒）。如何获得或构建规范的bed文件？

经过检验，跑单样本时，若设置-L所花费时间为17.2min，若不设置所花费时间为27.9min



### 第三步：HaplotypeCaller进行变异检测

#### HaplotypeCaller简介

HaplotypeCaller，简称HC。能过通过对活跃区域（也就是<u>与参考基因组不同处较多的区域</u>）**局部重组装**，同时寻找SNP和INDEL。



##### HaplotypeCaller的核心操作就是四步：

1. 寻找活跃区域，就是和参考基因组不同部分较多的区域
2. 通过对该区域进行局部重组装，确定单倍型（**haplotypes**）。就是这一步可以省去indel realignment
3. 在给定的read数据下，计算单倍型的可能性。
4. 分配样本的基因型。

![img](https://markdown-1300560293.cos.ap-guangzhou.myqcloud.com/markdown/202205100847335.webp)

#### HaplotvpeCaller操作方法

HaplotvpeCaller的应用有两种做法，**差别只在于要不要在中间生成一个gVCF**：

1. 直接进行HaplotypeCaller，这适合于单样本，或者那种固定样本数量的情况，也就是只执行一次HaplotypeCaller。如果增加一个样本就得重新运行这个HaplotypeCaller，而这个时候算法需要重新去读取所有人的BAM文件，浪费大量时间；
2. 每个样本先**各自生成gVCF**，然后再进行**群体joint-genotype**, gVCF全称是genome VCF，是每个样本用于变异检测的中间文件，格式类似于VCF, 它把joint-genotype过程中所需的所有信息都记录在这里面，文件无论是大小还是数据量都远远小于原来的BAM文件。当新增样本时，只需要重新执行joint-genotype即可。



**我们在此选择第二种方法**，对于每个样本gvcf的合并，同样存在两种方法：

* `CombineGVCFs`：将多个样本的gvcf文件直接合并为一个gvcf文件，这是一个总的gvcf文件。
  * `CombineGVCFs`的主要优点是能够一次组合多个区间，而无需构建 GenomicsDB。
  * CombineGVCFs 比 GenomicsDBImport 慢，因此建议仅在要合并的样本很少时使用 CombineGVCFs。

* `GenomicsDBImport`：将多个样本的gvcf文件生成一个工作空间，这是一个GenomicsDB工作目录。
  * 适用于>1000的样本量。




#### step1：生成中间文件gvcf

```shell
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" HaplotypeCaller \
-ERC GVCF \ # 标记生成中间文件gvcf
-R $hg38_ref \
-L $bed \
-I 04_bqsr/AD2021001.sorted.markdup.BQSR.bam \
-D $hg38_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-O 05_gvcf/AD2021001.gvcf
```

获得的gvcf文件如下所示：

```shell
grep -v "#" AD2021001.gvcf | less
chr1	23980414	.	C	<NON_REF>	.	.	END=23980416	GT#:DP:GQ:MIN_DP:PL	0/0:35:81:33:0,81,1215
chr1	23980417	.	G	<NON_REF>	.	.	END=23980419	GT:DP:GQ:MIN_DP:PL	0/0:33:75:33:0,75,1125
chr1	23980420	.	A	<NON_REF>	.	.	END=23980420	GT:DP:GQ:MIN_DP:PL	0/0:31:69:31:0,69,1035
chr1	23980421	.	A	<NON_REF>	.	.	END=23980423	GT:DP:GQ:MIN_DP:PL	0/0:31:57:31:0,57,855
chr1	23980424	.	C	<NON_REF>	.	.	END=23980429	GT:DP:GQ:MIN_DP:PL	0/0:27:54:25:0,54,810
chr1	23980430	.	C	<NON_REF>	.	.	END=23980434	GT:DP:GQ:MIN_DP:PL	0/0:19:51:18:0,51,765
chr1	23980435	.	G	<NON_REF>	.	.	END=23980436	GT:DP:GQ:MIN_DP:PL	0/0:19:54:19:0,54,810
chr1	23980437	.	T	<NON_REF>	.	.	END=23980438	GT:DP:GQ:MIN_DP:PL	0/0:19:48:19:0,48,720
chr1	23980439	.	T	<NON_REF>	.	.	END=23980439	GT:DP:GQ:MIN_DP:PL	0/0:19:42:19:0,42,630
chr1	23980440	.	T	<NON_REF>	.	.	END=23980441	GT:DP:GQ:MIN_DP:PL	0/0:18:39:18:0,39,585
chr1	23980442	.	C	<NON_REF>	.	.	END=23980443	GT:DP:GQ:MIN_DP:PL	0/0:18:36:18:0,36,540
chr1	23980444	.	C	<NON_REF>	.	.	END=23980457	GT:DP:GQ:MIN_DP:PL	0/0:13:30:10:0,30,315
chr1	23980458	.	G	<NON_REF>	.	.	END=23980460	GT:DP:GQ:MIN_DP:PL	0/0:10:27:10:0,27,405
chr1	23980461	.	T	<NON_REF>	.	.	END=23980461	GT:DP:GQ:MIN_DP:PL	0/0:8:24:8:0,24,254
chr1	23980462	.	C	<NON_REF>	.	.	END=23980462	GT:DP:GQ:MIN_DP:PL	0/0:10:27:10:0,27,405
chr1	23980463	.	C	<NON_REF>	.	.	END=23980467	GT:DP:GQ:MIN_DP:PL	0/0:10:24:10:0,24,360
chr1	23980468	.	G	<NON_REF>	.	.	END=23980472	GT:DP:GQ:MIN_DP:PL	0/0:9:18:9:0,18,270
chr1	23980473	.	T	<NON_REF>	.	.	END=23980489	GT:DP:GQ:MIN_DP:PL	0/0:6:15:5:0,15,176
```



**到目前为止，所有的样本操作（质控、比对、排序、标记重复序列、BQSR和gVCF）均在单样本模式下进行，所有流程将使用snakemake-wes构建标准流程。**



#### step2：多个gvcf文件整合

样本较少时，通过`-V`参数手动输入样本

```shell
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" CombineGVCFs \
-R $hg38_ref \
-V 05_gvcf/AD2021001.gvcf \
-V 05_gvcf/AD2021002.gvcf \
-L $bed \
-O 06_joint_genotype/AD1to2joint.gvcf
```



样本较多时，使用`-V`参数读取列表

```shell
gatk --java-options "-Xmx100G -Djava.io.tmpdir=./tmp" CombineGVCFs \
-R $hg38_ref $(for i in $(ls 05_gvcf/*.gvcf);do echo "--variant $i";done) \
-L $bed \
-O 06_joint_genotype/combined.gvcf
```

注意内存使用量，若设置过小，会报错 `java.lang.OutOfMemoryError: GC overhead limit exceeded`



由于该过程太慢，因此选用分染色体进行循环用以加速。

```shell
# 将bed文件按照不同染色体进行切分
for i in $(cat 00_ref/bed/interval.list | awk -F: '{print $1}' | uniq)
do
cat 00_ref/bed/interval.list | grep "$i:" > 06_joint_test/small_bed_$i.list
done
# 按照bed文件大小将较小的bed文件进行merge，来均衡并行数。
cat 06_joint_test/small_bed_chr4.list 06_joint_test/small_bed_chr5.list > 06_joint_test/small_bed_chr_4+5.list
cat 06_joint_test/small_bed_chr6.list 06_joint_test/small_bed_chr7.list > 06_joint_test/small_bed_chr_6+7.list
cat 06_joint_test/small_bed_chr8.list 06_joint_test/small_bed_chr9.list > 06_joint_test/small_bed_chr_8+9.list
cat 06_joint_test/small_bed_chr13.list 06_joint_test/small_bed_chr16.list > 06_joint_test/small_bed_chr_13+16.list
cat 06_joint_test/small_bed_chr14.list 06_joint_test/small_bed_chr15.list > 06_joint_test/small_bed_chr_14+15.list
cat 06_joint_test/small_bed_chr17.list 06_joint_test/small_bed_chr18.list > 06_joint_test/small_bed_chr_17+18.list
cat 06_joint_test/small_bed_chr19.list 06_joint_test/small_bed_chr20.list 06_joint_test/small_bed_chr21.list > 06_joint_test/small_bed_chr_19+20+21.list
cat 06_joint_test/small_bed_chr22.list 06_joint_test/small_bed_chrX.list 06_joint_test/small_bed_chrY.list > 06_joint_test/small_bed_chr_22+X+Y.list

for i in chr1 chr2 chr3 chr_4+5 chr_6+7 chr_8+9 chr10 chr11 chr12 chr_13+16 chr_14+15 chr_17+18 chr_19+20+21 chr_22+X+Y
do
nohup gatk --java-options "-Xmx100G" CombineGVCFs -R 00_ref/ucsc-human-hg38/hg38.fa -L 06_joint_test/small_bed_$i.list $(for j in $(ls 05_gvcf/*.gvcf);do echo "--variant $j";done) \
-O 06_joint_genotype/combined_$i.gvcf > 06_joint_genotype/output_$i.log 2>&1 &
done
```



#### step3：joint genotyping

在上述的两个snp call过程中，只是标记了哪些位点上存在snp，但没有进行基因分型。

```shell
gatk --java-options "-Xmx200G -Djava.io.tmpdir=./tmp" GenotypeGVCFs \
-R $hg38_ref \
-V 06_joint/combined.gvcf \
-L $bed \
-G StandardAnnotation \
--only-output-calls-starting-in-intervals \
--use-new-qual-calculator \
-O 06_joint_genotype/WES_AD_variants.vcf
```



同样的，如果要分割染色体进行genotyping，使用以下代码，然后使用gatk merge进行合并。

```shell
for i in chr1 chr2 chr3 chr_4+5 chr_6+7 chr_8+9 chr10 chr11 chr12 chr_13+16 chr_14+15 chr_17+18 chr_19+20+21 chr_22+X+Y
do
nohup gatk --java-options "-Xmx200G -Djava.io.tmpdir=./tmp" GenotypeGVCFs \
-R 00_ref/ucsc-human-hg38/hg38.fa \
-L 06_joint_test/small_bed_$i.list \
-V 06_joint_test/combined_$i.gvcf \
-G StandardAnnotation \
--only-output-calls-starting-in-intervals \
--use-new-qual-calculator \
-O 06_joint_test/WES_AD_variants_$i.vcf >> 06_joint_test/output_$i.log 2>&1 &
done
gatk MergeVcfs \
$(for i in $(ls 06_joint_test/WES_AD_variants_*.vcf);do echo "-I $i";done) \
-O 06_joint_test/WES_AD_variants.vcf
```



## 5、变异质控，VQSR

参考：

[GATK4.0和全基因组数据分析实践（下）](https://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798455&idx=1&sn=67a7407980a57ce138948eb46992b603&chksm=83c1d52bb4b65c3dde31df94e9686654bf616166c7311b531213ebf0010f67a32ce827e677b1&scene=21#wechat_redirect)

[如何正确设置GATK VQSR的模型训练参数](https://zhuanlan.zhihu.com/p/40823886)



**质控的含义和目的是指通过一定的标准，最大可能地剔除假阳性的结果，并尽可能地保留最多的正确数据**。

在GATK HaplotypeCaller之后，**首选的质控方案是GATK VQSR** ，[原理](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612?id=39)是利用自身数据和已知变异位点集的overlap，通过GMM模型构建一个分类器来对变异数据进行打分，从而评估每个位点的可信度。



具体来说：

1. GATK认为VQSR比根据各种annotations进行hard-filtering过滤要好，减少了人为阈值的局限性，避免了一刀切的弊端，从而在sensitivity和specificity之间达到一定的平衡；
2. VQSR根据机器学习算法从highly validated变异位点数据集（每个位点的annotation profile，一般使用5-8个annotation）中获取到good variants/bad variants；
3. 根据上述的位点从我们自己数据集中挑选出一个变异子集（probably true positives）来建模训练，获得一个可识别good variants的模型；bad variants的模型也是如此获得；
4. 然后根据上述获得的模型，对自己数据集的每个变异位点进行一个总的打分；
5. 最后根据**设定的sensitivity阈值**对变异位点进行过滤。



### 两个前提条件：

1. 需要一个**精心准备的已知变异集**，它将作为训练质控模型的真集。对于人可以使用Hapmap、OMNI，1000G和dbsnp等这些国际性项目的数据作为高质量的已知变异集。
   - 导致很多非人的物种无法使用BQSR

2. 要求新检测的结果中**有足够多的变异**，不然VQSR在进行模型训练的时候会因为可用的变异位点数目不足而无法进行。
   - 对于样本数较少的WES或者小panel测序，由于最后的变异位点不够，也无法使用VQSR。



不满足VQSR时的处理：

- GATK对于样本不达标的WES，建议用1000 Genomes Project中的数据代替（将1000G的bam文件和WES样本的bam生成g.vcf文件，再一起做joint call）
- 不满足这两个条件就只能采用硬过滤进行质控了（通过人为设定一个或者若干个指标阈值，如QUAL，然后把所有不满足阈值的变异位点采用一刀切掉的方法），硬过滤将调用gatk SelectVariants，VariantFiltration， MergeVcfs等功能（具体查看上面的参考）。




### 质控（VQSR的2个步骤）

* VariantRecalibrator：即上述原理的1-4，对每个变异位点打分注释VQSLOD value，从而生成一个recalibration文件以及一个xxx.plots.R.pdf(Gaussian mixture model plots)。

* ApplyRecalibration：对应5步骤，根据recalibration文件生成recalibrated VCF文件，并且根据过滤参数进行过滤（标记PASS）。



**对于设置sensitivity阈值：**

* 所以当你tranche sensitivity（--truth-sensitivity-filter-level）设置为99.9%时，则表示将VQSLOD值高于整体99.9%的变异位点标记为PASS，表示通过过滤阈值，剩下的位点则被认为是假阳性的了；
* 这个参数可以根据你的具体需求来设定，看你是需要more specific or more sensitive（tranche从90到100，specific随之降低，而sensitive随之升高）；或者生成多个tranche的结果，从中挑选你满意的阈值（可以看结果文件SNP中的xxx.tranches，而indel是没这个文件的；人全基因组Ti/Tv（转换和颠换的比例）的值一般在2.0-2.1...）



**对于annotation参数：**

GATK的VQSR用的annotations指标有以下几种（VQSR的文档中还提到了InbreedingCoeff，这个是在大于10个样本中才使用的一个指标；如果样本过少或者是closely related samples(such as a family)的话，建议剔除；而且QD不适合用于外显子测序的数据的VQSR）：

- Coverage（DP）
- QualByDepth（QD）
  - QD是变异质量值(Quality)除以覆盖深度(Depth)得到的比值。这里的变异质量值就是VCF中QUAL的值——用来衡量变异的可靠程度，这里的覆盖深度是这个位点上所有含有变异碱基的样本的覆盖深度之和，通俗一点说，就是这个值可以通过累加每个含变异的样本（GT为非0/0的样本）的覆盖深度（VCF中每个样本里面的DP）而得到。

- FisherStrand (FS)
  - FS是一个通过Fisher检验的p-value转换而来的值，它要描述的是测序或者比对时对于只含有变异的read以及只含有参考序列碱基的read是否存在着明显的正负链特异性（Strand bias，或者说是差异性）。这个差异反应了测序过程不够随机，或者是比对算法在基因组的某些区域存在一定的选择偏向。如果测序过程是随机的，比对是没问题的，那么不管read是否含有变异，以及是否来自基因组的正链或者负链，只要是真实的它们就都应该是比较均匀的，也就是说，不会出现链特异的比对结果，FS应该接近于零。

- StrandOddsRatio (SOR)
  - SOR原理和FS类似，但是用了不同的统计检验方法计算差异程度，用的是symmetric odds ratio test，它更适合于高覆盖度的数据。

- RMSMappingQuality (MQ)
  - MQ这个值是所有比对至该位点上的read的比对质量值的均方根(先平方、再平均、然后开方)。它和平均值相比更能够准确地描述比对质量值的离散程度。

- MappingQualityRankSumTest (MQRankSum)
  - 对杂合位点进行的不同碱基之间比对质量的曼惠特尼秩和检验结果，通过ref和alt碱基的比对质量差异来评估位点的可信度。

- ReadPosRankSumTest (ReadPosRankSum)
  - 对杂合位点进行的秩和检验，看不同的碱基是否倾向于出现在reads上的特定位置(例如接近reads的起始或者终止)。




对于GATK的VQSR，**snp和indel过滤是分开的**，两者参数整体上差不多，但有点略微的区别

#### SNP的VQSR的过滤

一般的参数可以如下：

**生成校准表**

```shell
# rescource文件需要先index，已经在GATK流程中的第一步中完成。
gatk VariantRecalibrator -R $hg38_ref -V 06_joint/WES_AD_variants.vcf \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hg38_vcf/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf \
--resource:omini,known=false,training=true,truth=false,prior=12.0 $hg38_vcf/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 $hg38_vcf/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $hg38_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-an DP -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum --mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
--output 07_BQSR/WES_AD.snp.recal \
--tranches-file 07_BQSR/WES_AD.snp.tranches \
--rscript-file 07_BQSR/WES_AD.snp.plots.R
```

参数：

* `-tranche`：默认是输出[100,99.9,99.0,90.0]4个tranche阈值的统计结果，如果想看其他阈值的结果，需要自行加上
* `-an QD`：不适用于WES，此处没有添加
* `-resource`参数：
  * 训练集名字，这个名字是可以随便改动的，但是为了便于交流，一般还是默认按照数据集的名字来设置(如上面的例子)；
  * known：该数据是否作为已知变异数据集，用于对变异数据的标注；
  * training：该数据是否作为模型训练的数据集，用于训练VQSR模型；
  * truth：该数据是否作为验证模型训练的真集数据，这个数据同时还是VQSR训练bad model时自动进行参数选择的重要数据；
  * prior：该数据集在VQSR模型训练中的权重，或者叫Prior likelihood(这里转化为Phred-scale，比如20代表的是0.99)。

* **SNP的VQSR过滤，选用的resource datasets为：**
  - HapMap，
    - 来自国际人类单倍体型图计划，这个数据集包含了大量家系数据，并且有非常严格的质控和严密的实验验证，因此它的准确性是目前公认最高的。
    - hapmap_3.3.hg38.vcf.gz，`truth=true`表示VQSR将这个数据集中的变异位点作为真位点true sites，`training=true`表示VQSR将true sites用于训练recalibration model，并赋予这些变异位点prior likelihood值为`Q15 (96.84%)`
  - Omni，
    - 源自Illumina的Omni基因型芯片，大概2.5百万个位点，它也是一个高可信的变异结果。
    - 1000G_omni2.5.hg38.vcf.gz，`truth=true`，`training=false`（文档中写着是true，参数建议中写着的是false。。。我就按照参数上的来了），`Q12 (93.69%)`
  - 1000G，
    - 是千人基因组计划(1000 genomes project)质控后的变异数据，一般不作为模型验证的真集数据。
    - 1000G_phase1.snps.high_confidence.hg38.vcf.gz，`truth=false`表示VQSR考虑到在1000G数据集中的不仅包含了true variants还有false positives，`training=true`，`Q10 (90%)`
  - dbSNP，
    - 是一个绝对不可以作为训练集位点的数据——收集的数据实际都是研究者们发表了相关文章提交上来的变异，这些变异很多是没做过严格验证的，很多甚至还是假的，在没被反复验证之前，是不可信的，dbSNP的唯一作用就是用于标注我们的变异集中哪些是已经在其它研究中出现过的。
    - dbsnp_146.hg38.vcf.gz，`truth=false`表示VQSR未将dbSNP数据集中的位点作为可信数据集，`training=false`表示不用于训练数据集，`known=true`表示stratify output metrics such as Ti/Tv ratio by whether variants are present in dbsnp or not，`Q2 (36.90%)`



**应用VQSR**

```shell
gatk ApplyVQSR  -R $hg38_ref -V 06_joint/WES_AD_variants.vcf \
-O 07_BQSR/WES_AD.snp.VQSR.vcf \
--truth-sensitivity-filter-level 99.5 \
--tranches-file 07_BQSR/WES_AD.snp.tranches \
--recal-file 07_BQSR/WES_AD.snp.recal \
--mode SNP
```

参数:

`--truth-sensitivity-filter-level`：与上一步的tranche参数对应，上一步为标记，这一步为过滤。可以根据自身需求来设定，此处设定为99.5。



注意，該步驟會使用上述生成的R脚本，如果error與R脚本有關，則需要檢查R的環境，最好在conda裏重裝一個R環境。

```shell
conda install -c conda-forge r-base=4.1.3
conda install -c r r-ggplot2
```



通过阈值标准的SNP将被标记为“PASS”，此处进行提取，未标记PASS的位点不纳入下游分析

```shell
cat 07_BQSR/WES_AD.snp.VQSR.vcf | grep "PASS" > 07_BQSR/WES_AD.snp.filtered.vcf
cat 07_BQSR/WES_AD.snp.VQSR.vcf | wc -l # 统计全部位点数目
## 877950
cat 07_BQSR/WES_AD.snp.filtered.vcf | wc -l # 统计通过质控的位点数目
## 789954
```



#### Indel的VQSR的过滤

其实INDEL由于置信度较低，在后续分析中将不予考虑。

**生成校准表**

```shell
gatk VariantRecalibrator  -R $hg38_ref -V 06_joint/WES_AD_variants.vcf \
--max-gaussians 4 \
--resource:mills,known=false,training=true,truth=true,prior=12.0 $hg38_vcf/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $hg38_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum --mode INDEL \
--output 07_BQSR/WES_AD.indel.recal \
--tranches-file 07_BQSR/WES_AD.indel.tranches \
--rscript-file 07_BQSR/WES_AD.indel.plots.R
```

参数：

- `--max-gaussians 4`用于设定Gaussians（clusters of variants that have similar properties）的数目，即减少聚类的组数，从而使得每个组的变异位点数目达到要求
- 未设置`MQ`过滤
- **INDEL的VQSR过滤，选用的resource datasets为：**
  - Mills，
    - 对于Indel来说能正在算得上真集的并不多，Mills_and_1000G_gold_standard.indels.hg38.vcf算是其中一个，并被专门做过验证。
    - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz，`truth=true`，`training=true`，`Q12 (93.69%)`
  - dbSNP，dbsnp_146.hg38.vcf.gz，`truth=false`，`training=false`，`known=true`，`Q2 (36.90%)`
    - 和上面SNP模式下的作用是一样的。不过假如这一步对Indel进行VQSR的VCF数据是顺着上面SNP VQSR后下来的话，那么这个dbSNP的参数可以省略。



**应用VQSR**

```shell
gatk ApplyVQSR  -R $hg38_ref -V 06_joint/WES_AD_variants.vcf \
-O 07_BQSR/WES_AD.indel.VQSR.vcf \
--truth-sensitivity-filter-level 99.0 \
--tranches-file 07_BQSR/WES_AD.indel.tranches \
--recal-file 07_BQSR/WES_AD.indel.recal \
--mode INDEL
```





#### 合并变异

由于INDEL置信度低，在后续分析中不予考虑。该步骤不会执行。

其实就是简单的将上面的vcf文件进行叠加，完成后每个位点有SNP和INDEL两行信息。

```
gatk MergeVcfs \
-I 07_BQSR/WES_AD.snp.VQSR.vcf \
-I 07_BQSR/WES_AD.indel.VQSR.vcf \
-O 07_BQSR/WES_AD.all.VQSR.vcf
```





### 硬过滤（不采用）

* `gatk SelectVariants `：该工具可以根据各种标准选择变量的子集，以便于进行某些分析。例如比较case与control，提取满足特定要求的变异或非变异位点，或排除一些意外结果，等等。 
* `gatk VariantFiltration`：基于INFO和/或FORMAT注释的过滤变量调用。  

```shell
# 使用SelectVariants，选出SNP
gatk SelectVariants -select-type SNP -V 07_BQSR/WES_AD.snp.VQSR.vcf -O 07_BQSR/WES_AD.snp.select.vcf
# 为SNP作过滤
gatk VariantFiltration -V 07_BQSR/WES_AD.snp.select.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "PASS" \
-O 07_BQSR/WES_AD.snp.filter.vcf

# 使用SelectVariants，选出Indel
gatk SelectVariants -select-type INDEL -V wes.raw.vcf -O wes.indel.vcf
# 为Indel作过滤
gatk VariantFiltration -V wes.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "PASS" -O wes.indel.filter.vcf
```



## 6、位点注释

主要使用的工具：`ANNOVAR`



### ANNOVAR下载数据库（数据库类型和解释也前往[官网](https://annovar.openbioinformatics.org/en/latest/user-guide/download/#annovar-main-package)进行查看下载）

注意有些数据库（如dbnsfp30a和avsnp147）十分大，使用命令行下载会failed。建议复制好网址，在windows上下载。

```shell
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/ # 来自NCBI的注释
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar knownGene humandb/ # 来自UCSC的注释
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
# -buildver 表示version
# -downdb 下载数据库的指令
# -webfrom annovar 从annovar提供的镜像下载，不加此参数将寻找数据库本身的源
# humandb/ 存放于humandb/目录下
# 注意有些gz压缩的数据库需要解压
gunzip humandb/hg38_avsnp147.txt.gz
gunzip humandb/hg38_dbnsfp30a.txt.gz
```



### ANNOVAR输入准备

1. 方案一：直接使用vcf文件。`table_annovar.pl`程序可以直接使用vcf程序，需要添加参数`-vcfinput`
2. 方案二：先转换成`aviput`格式，最后进行注释



#### 方案一

ANNOVAR使用.avinput格式，如以上代码所示，该格式每列以tab分割，最重要的地方为前5列，分别是

1. 染色体(Chromosome)
2. 起始位置(Start)
3. 结束位置(End)
4. 参考等位基因(Reference Allele)
5. 替代等位基因(Alternative Allele)
6. 剩下为注释部分（可选）。

实操：

```shell
perl ~/software/annovar/annovar/convert2annovar.pl -format vcf4old 07_BQSR/WES_AD.snp.filtered.vcf \
-includeinfo  -comment -outfile 08_annovar/WES_snp.avinput
## -format vcf4 指定格式为vcf，在多样本时设置为vcf4old，或增加-allsample
## -outfile 输出
## -includeinfo 会保留vcf文件中的所有信息(可选)
## -comment 会保留vcf文件的头部注释信息(以#开头的行,可选)
## -allsample 转换格式时vcf中的每一个样本会单独生成一个待注释的vcf文件(可选)
## -withzyg 输出杂合性，质量，read覆盖度等信息
## -withfreq  输出等位基因频率
## -comment 添加VCF的header信息
cat 08_annovar/WES_snp.avinput | less -S
## chr1    13273   13273   G       C       hom     43379.72        17354   41.08   6.46
## chr1    13284   13284   G       A       hom     3010.60 12783   42.78   11.67
## chr1    13303   13303   G       A       hom     883.60  12586   41.82   3.40
## chr1    13334   13334   G       A       hom     2613.61 12548   42.95   10.37
## chr1    13372   13372   G       C       hom     494.61  12259   41.70   3.96
## chr1    13417   13417   -       GAGA    unknown 17598.40        9401    37.30   13.53
## chr1    13418   13418   G       A       hom     245.77  9036    34.33   22.34
## chr1    13504   13504   G       A       hom     62.47   1492    24.28   2.40
humandb=~/software/annovar/annovar/humandb
perl ~/software/annovar/annovar/table_annovar.pl 08_annovar/WES_snp.avinput \
$humandb --buildver hg38 --thread 12 \
-out 08_annovar/snp_anno \
-remove -protocol refGene,knownGene,cytoBand,exac03,avsnp147,dbnsfp30a \
-operation g,g,r,f,f,f -nastring . 
## -protocol: 注释的数据库
## -buildver： 基因组版本
## -nastring：缺省值用 . 填充
## --otherinfo： 输入文件中第5例后面的info信息也进行输出
## --operation：数据库的类型，需要与前面的-protocol顺序严格对应，且逗号分割。g为基因注释类型 ，r为区域注释类型，f为过滤注释类型，gx表示基于基因的交叉引用注释（与-xref参数有关）。
## 数据库数据哪种类型可以前往官网查看说明书，具有详细说明。
## --gff3dbfile：注释使用的 gff 文件
## -vcfinput：输入为vcf格式的文件，输出也为vcf格式
## --thread：线程数
## --xreffile：为基于基因的注释指定一个交叉引用文件，以便最终输出包括基因的额外列
```

其他选择：

annovar其实提供了好多个gene-based数据库下载，refGene是NCBI提供的，还有ensembl、UCSC提供的ensGene、knowGene，这三个数据库是有所差别的，至于怎么选择，建议就和使用的参考基因组来源一致即可，这里使用了来自NCBI的refGene和UCSC的knownGene进行基于基因的注释。

-csvout 表示最后输出.csv文件，使用这个参数得到的是csv格式的文件，可以直接用excel打开，但是我没有使用，输出是tab分隔的txt文件，方便后续我写脚本处理。为什么“,”分隔的csv文件不适合处理呢？那是因为refgene的注释结果也是逗号分隔的，在分析的时候就会出现问题。

```
grep -v "##" out.avinput | cut -f1-9 --complement > gt.txt
paste hg19_multianno.txt gt.txt > merge.anno.txt
```



#### 方案二

直接使用vcf作为输入文件。

```shell
perl ~/software/annovar/annovar/table_annovar.pl 07_BQSR/WES_AD.snp.VQSR.vcf \
$humandb -buildver hg38 -out 08_annovar/snp_anno \
-remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
-operation gx,r,f,f,f -nastring . -vcfinput -thread 6
```





