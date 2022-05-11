# snakemake-wes代码

## yaml配置文件

```yaml
samples:
  AD2021001:
  AD2021002:
```



```
# 获得所有文件名
path=~/Rawdata/AAD_WGS_20220127/GHDLB20070348 # rawdata位置
files=$(ls $path | grep "\\_1.fq.gz")
for filename in $files
do
   echo $filename | awk -F"_" '{print $1;exit}'>> ./sample_id.list
done
```



```
# 生成配置文件
echo "samples:" > config.yaml
for line in $(cat sample_id.list)
do
echo "  ${line}:" >> config.yaml
done
```



## snakemake讲解

`rule all`是每个Snakefile文件中必有的一个`rule`，比较特殊，只需要一个`input`，用来定义流程最终输出的结果。

除了rule all，其余必须有output。

注意：如果你的流程有不同的分支，最终会生成多个需要的结果，那么这些结果都需要在这里定义。

```python
rule all:
    input:
        expand("05_gvcf/{sample}.gvcf", sample=config["samples"])
```

`rule fastp`等为具体步骤，主要有三个部分组成：

* input：输入文件，用`，`隔开
* output：输出文件，用`，`隔开；没有被shell定义的输出可以不写，但写了也没关系，也方便软件识别流程顺序
* shell：正常的shell脚本
* 其他：
  * log：日志文件
  * threads：定义线程数
  * params：各类参数

**snakemake规则：**

- Snakemake可以自动确定不同规则的输入输出的依赖关系，根据**时间戳**来判断文件是否需要重新生成
- Snakemake以{sample}.fa形式进行**文件名通配**，用{wildcards.sample}获取sample的实际文件名
- Snakemake用expand()生成多个文件名，本质是Python的列表推导式
- Snakemake可以在规则外直接写Python代码，在规则内的run里也可以写Python代码。
- Snakefile的第一个规则通常是rule all，根据all里的文件决定执行哪些rule。如上面的例子，注释掉all里的input则不执行第二条rule
- 在output中的结果文件可以是未存在目录中的文件，这时会自动创建不存在的目录（不需要事先建文件夹，这个功能实在是方便）

## snakemake代码



```python
###########################################################################################################
# 2022-04-29
# Long Xinyang
# wes pipeline with snakemake
# 对于并行计算使用8线程
###########################################################################################################

configfile: "config.yaml"
hg38_vcf = "/share/home/longxinyang/bio/2022_AAD_WES/shell/00_ref/gatk_call_vcf"
hg38_ref = "/share/home/longxinyang/bio/2022_AAD_WES/shell/00_ref/ucsc-human-hg38/hg38.fa"
interval = "/share/home/longxinyang/bio/2022_AAD_WES/shell/00_ref/bed/interval.list"


rule all:
    input:
        "07_BQSR/WES_AD.all.VQSR.vcf"


# 质控：Quality Control
rule fastp:
    input:
        "/share/home/longxinyang/Rawdata/AAD_WGS_20220127/GHDLB20070348/{sample}_1.fq.gz",
        "/share/home/longxinyang/Rawdata/AAD_WGS_20220127/GHDLB20070348/{sample}_2.fq.gz"
    output:
        "01_fastp/{sample}_1.clean.fq.gz",
        "01_fastp/{sample}_2.clean.fq.gz",
        "01_fastp/{sample}.fastp.json",
        "01_fastp/{sample}.fastp.html"
    # log:
    #     "01_fastp/{sample}.log"
    threads: 8
    shell:
        """fastp \
-i {input[0]} -o {output[0]} \
-I {input[1]} -O {output[1]} \
-z 4 -f 5 -t 5 -F 5 -T 5 -5 -W 5 -M 20 -Q -l 50 -c -w {threads} \
-j {output[2]} \
-h {output[3]}"""

# 比对：alignment
rule bwa_mem:
    input:
        "01_fastp/{sample}_1.clean.fq.gz",
        "01_fastp/{sample}_2.clean.fq.gz"
    output:
        "02_bwa_out/{sample}.sorted.bam"
    threads: 8
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    shell:
        "bwa mem -t {threads} -M -R '{params.rg}' {hg38_ref} {input[0]} {input[1]} | samtools sort -@ {threads} -m 20G -O bam -o {output} - "


# 标记重复序列（MarkDuplicates）
rule MarkDuplicates:
    input:
        "02_bwa_out/{sample}.sorted.bam"
    output:
        "03_markdup/{sample}.sorted.markdup.bam",
        "03_markdup/{sample}.markdup_matrics.txt",
        "03_markdup/{sample}.sorted.markdup.bam.bai"
    shell:
        "gatk MarkDuplicates -I {input} -O {output[0]} -M {output[1]}\n"
        "samtools index {output[0]}"

# 碱基质量分数重校准（BaseRecalibrator）
rule BaseRecalibrator:
    input:
        "03_markdup/{sample}.sorted.markdup.bam",
        "03_markdup/{sample}.sorted.markdup.bam.bai"
    output:
        "04_bqsr/{sample}.recal_data.table"
    shell:
        """gatk --java-options '-Xmx30G -Djava.io.tmpdir=./tmp' BaseRecalibrator \
-R {hg38_ref} --input {input[0]} \
--known-sites {hg38_vcf}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
--known-sites {hg38_vcf}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
--known-sites {hg38_vcf}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-L {interval} -O {output}"""
        
# 应用碱基质量分数重校准（ApplyBQSR）
rule ApplyBQSR:
    input:
        "03_markdup/{sample}.sorted.markdup.bam",
        "04_bqsr/{sample}.recal_data.table"
    output:
        "04_bqsr/{sample}.sorted.markdup.BQSR.bam"
    shell:
        "gatk --java-options '-Xmx30G -Djava.io.tmpdir=./tmp' ApplyBQSR -R {hg38_ref} -I {input[0]} -bqsr {input[1]} -L {interval} -O {output}"

rule HaplotypeCaller:
    input:
        "04_bqsr/{sample}.sorted.markdup.BQSR.bam"
    output:
        "05_gvcf/{sample}.gvcf"
    shell:
        """gatk --java-options '-Xmx30G -Djava.io.tmpdir=./tmp' HaplotypeCaller \
-ERC GVCF -R {hg38_ref} -I {input} -D {hg38_vcf}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf -L {interval} -O {output}"""


rule CombineGVCFs:
    input:
        expand("05_gvcf/{sample}.gvcf", sample=config["samples"])
    output:
        "06_joint_genotype/combined.gvcf"
    params:
        " -V ".join(expand("05_gvcf/{sample}.gvcf", sample=config["samples"]))
    shell:
        """gatk --java-options "-Xmx100G -Djava.io.tmpdir=./tmp" CombineGVCFs \
-R {hg38_ref} --variant {params} -L {interval} -O {output}"""

rule GenotypeGVCFs:
    input:
        "06_joint_genotype/combined.gvcf"
    output:
        "06_joint_genotype/WES_AD_variants.vcf"
    shell:
        """gatk --java-options "-Xmx200G -Djava.io.tmpdir=./tmp" GenotypeGVCFs \
-R {hg38_ref} -V {input} -L {interval} -G StandardAnnotation \
--only-output-calls-starting-in-intervals --use-new-qual-calculator \
-O {output}"""

rule VQSR_snp:
    input:
        "06_joint_genotype/WES_AD_variants.vcf"
    output:
        "07_BQSR/WES_AD.snp.recal",
        "07_BQSR/WES_AD.snp.tranches",
        "07_BQSR/WES_AD.snp.plots.R",
        "07_BQSR/WES_AD.snp.VQSR.vcf"

    shell:
        """gatk VariantRecalibrator \
-R {hg38_ref} -V {input} \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hg38_vcf}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf \
--resource:omini,known=false,training=true,truth=false,prior=12.0 {hg38_vcf}/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 {hg38_vcf}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {hg38_vcf}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-an DP -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum --mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
--output {output[0]} \
--tranches-file {output[1]} \
--rscript-file {output[2]}"""
        """gatk ApplyVQSR \
-R {hg38_ref} -V {input} \
-O {output[3]} \
--truth-sensitivity-filter-level 99.5 \
--tranches-file {output[1]} \
--recal-file {output[0]} \
--mode SNP"""

rule VQSR_indel:
    input:
        "06_joint_genotype/WES_AD_variants.vcf"
    output:
        "07_BQSR/WES_AD.indel.recal",
        "07_BQSR/WES_AD.indel.tranches",
        "07_BQSR/WES_AD.indel.plots.R",
        "07_BQSR/WES_AD.indel.VQSR.vcf"

    shell:
        """gatk VariantRecalibrator \
-R {hg38_ref} -V {input} \
--max-gaussians 4 \
--resource:mills,known=false,training=true,truth=true,prior=12.0 {hg38_vcf}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {hg38_vcf}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum --mode INDEL \
--output {output[0]} \
--tranches-file {output[1]} \
--rscript-file {output[2]}"""
        """gatk ApplyVQSR \
-R {hg38_ref} -V {input} -O {output[3]} \
--truth-sensitivity-filter-level 99.0 \
--tranches-file {output[1]} \
--recal-file {output[0]} \
--mode SNP"""

rule VQSR_merge:
    input:
        "07_BQSR/WES_AD.snp.VQSR.vcf",
        "07_BQSR/WES_AD.indel.VQSR.vcf"
    output:
        "07_BQSR/WES_AD.all.VQSR.vcf"
    shell:
        "gatk MergeVcfs -I {input[0]} -I {input[1]} -O {output}"



# 试运行
# module load anaconda/3-2021.11 
# source activate
# conda activate wes
# snakemake -s snakemake_wes.py -n -p
# 绘制流程图
# snakemake -s snakemake_wes.py --dag | dot -Tpdf > test_dag.pdf
```







## 提交集群

对于snakemake来说，不能通过像shell提交循环那种方式进行任务投递。由于snakemake不允许有多个snakemake在同一路径下运行。



对于snakemake来说，支持lsf集群的使用。

参考:[Snakemake LSF profile](https://github.com/Snakemake-Profiles/lsf)



### 安装

使用Cookiecutter部署配置文件

```shell
pip install --user cookiecutter
# or
conda install -c conda-forge cookiecutter
```



### 在集群上下载profile文件



```shell
# create configuration directory that snakemake searches for profiles
profile_dir="${HOME}/.config/snakemake"
mkdir -p "$profile_dir"
# use cookiecutter to create the profile in the config directory
template="gh:Snakemake-Profiles/lsf"
cookiecutter --output-dir "$profile_dir" "$template"
```





### 配置profile文件

#### LSF_UNIT_FOR_LIMITS

**Default**: `KB`
**Valid options:** `KB`, `MB`, `GB`, `TB`, `PB`, `EB`, `ZB`

重要！通过以下代码进行查看

```shell
grep '^LSF_UNIT_FOR_LIMITS' ${LSF_ENVDIR}/lsf.conf
# LSF_UNIT_FOR_LIMITS=MB
```

如果您的集群上这个值是MB，那么当使用-M 1000设置内存限制时，这个值就被认为是MB。因为snake允许你用参数resources: mem_mb来设置rule的内存，所以这个配置文件需要知道在提交时是否需要转换成其他单位。



#### UNKWN_behaviour

**Default**: `wait`
**Valid options**: `wait`, `kill`

当LSF返回的作业状态是UNKWN时，是等待考虑再次运行还是kill掉？



#### ZOMBI_behaviour

**Default**: `ignore`
**Valid options**: `ignore`, `kill`

当LSF返回的作业状态是ZOMBI时，是忽略掉还是kill掉？（任务都会被认为失败）



#### latency_wait

**Default:** `5`

等同于设置 `snakemake` 中的 `--latency-wait/--output-wait/-w` 参数。作业完成后，等待输出文件的时间（seconds）。

```
  --latency-wait SECONDS, --output-wait SECONDS, -w SECONDS
                        Wait given seconds if an output file of a job is not
                        present after the job finished. This helps if your
                        filesystem suffers from latency (default 5).
```



#### use_conda

**Default**: `False`
**Valid options:** `False`, `True`

相当于 `snakemake` 中的 `--use-conda` 参数。是否识别conda指令？

```
  --use-conda           If defined in the rule, run job in a conda
                        environment. If this flag is not set, the conda
                        directive is ignored.
```



#### use_singularity

**Default**: `False`
**Valid options:** `False`, `True`

相当于 `snakemake` 中的 `--use-singularity` 参数。是否在`singularity`容器中运行作业。

```
  --use-singularity     If defined in the rule, run job within a singularity
                        container. If this flag is not set, the singularity
                        directive is ignored.
```



#### restart_times

**Default**: `0`

相当于 `snakemake` 中的 `--restart-times` 参数。对于失败作业，重新运行的次数。

```
  --restart-times RESTART_TIMES
                        Number of times to restart failing jobs (defaults to 0).
```



#### print_shell_commands

**Default**: `False`
**Valid options:** `False`, `True`

相当于 `snakemake` 中的 `--printshellcmds/-p` 参数。是否打印执行的shell脚本。

```
  --printshellcmds, -p  Print out the shell commands that will be executed.
```



#### jobs

**Default**: `500`

相当于 `snakemake` 中的 `--cores/--jobs/-j` 参数。

* 最多使用N个核进行并行。
* 在集群的背景下，-j表示同时提交给集群的作业数量。

```
  --cores [N], --jobs [N], -j [N]
                        Use at most N cores in parallel. If N is omitted or
                        'all', the limit is set to the number of available
                        cores.
```



#### default_mem_mb

**Default**: `1024`

设置默认内存（单位MB）。



#### default_cluster_logdir

**Default**: `"logs/cluster"`

* 日志文件的目录。
* 默认值是一个相对路径，相对于该pipeline的work directory。如果不存在，会自动创建。
* logs会按照不同的`rule`被分配到不同的子目录中。
* 其标准输出：`logs/cluster/foo/sample=a,ext=fq/jobid<jobid>-<uuid>.out`



#### default_queue

**Default**: None

默认提交的队列名称。如果不设置，将会使用cluster默认队列。相当于`bsub`命令中的参数`-q`



#### default_project

**Default**: None

默认提交的默认项目名。如果不设置，将使用cluster默认项目名。相当于`bsub`命令中的参数`-P`



#### max_status_checks_per_second

**Default**: `10`

相当于 `snakemake` 中的 `--max-status-checks-per-second` 参数。每秒作业状态检查的最大数目，默认为10。

```
  --max-status-checks-per-second MAX_STATUS_CHECKS_PER_SECOND
                        Maximal number of job status checks per second,
                        default is 10, fractions allowed.
```



#### max_jobs_per_second

**Default**: `10`

相当于 `snakemake` 中的 `--max-jobs-checks-per-second` 参数。每秒cluster/drmaa任务最大数，默认为10。

```
  --max-jobs-per-second MAX_JOBS_PER_SECOND
                        Maximal number of cluster/drmaa jobs per second,
                        default is 10, fractions allowed.
```



#### profile_name

**Default**: `lsf`

该配置文件的名称。位于`$HOME/.config/snakemake/<profile_name>`

使用时传递给`--profile`参数

```shell
snakemake --profile <profile_name>
```





### 为每个rule提供单独的资源配置



#### rule特异性的标准集群资源设置

可以在snakemake文件中的每条rule中指定以下资源配置：

* `threads`: <INT>：指定job需要的线程数。
* resources:
  * `mem_mb = <INT>`:该rule所需的内存，单位为MB。如果未指定，将使用profile中设置的值。
  * `time_min: <INT>`：该rule所需的运行时限制，以分钟为单位。

**注意:这些设置将覆盖配置文件默认值。**



#### rule特异性的非标准集群资源设置

在snakemake中cluster configuration被弃用后，为每个rule进行资源配置成为snakemake-profile专属功能。



**具体实现：**

会在snakemake file之外建立一个cluster config文件（JSON或YAML格式），如下图所示，其包含与snakemake中的rule名称一样的对象，从而达到为每个rule都进行单独的资源设置。

![image-20220509150032289](https://raw.githubusercontent.com/PiaoyangGuohai1/Typora-image/main/202205091500350.png)



这些`rule`的配置必须放在一个名为`lsf.yaml`的文件中，且必须位于pipeline的working directory。

**注意:这些设置只对本`profile`文件有效，不保证在非lsf集群系统上有效。**



在`yaml`文件的设置中:

* 所有设置都以规则名称作为`key`；
* cluter配置以字符串(标量)或列表(序列)的形式给出；
* 在`__default__`中指定适用于所有`rule`的配置；
* 在各`rule`中的配置会覆盖掉默认配置。



我的配置文件，使用josn格式。该文件放置于`profile_dir`，通过指定`--cluster-config`可以成功运行。预计放在working directory也可以识别。

```json
{
    "__default__" :
    {
        "queue"     : "fat2",
        "nCPUs"     : "2",
        "memory"    : 35000,
        "resources" : "\"select[mem>30720] rusage[mem=30720] span[hosts=4]\"",
        "name"      : "JOBNAME.{rule}.{wildcards}",
        "output"    : "logs/cluster/{rule}.{wildcards}.out",
        "error"     : "logs/cluster/{rule}.{wildcards}.err"
    },


    "fastp" :
    {
        "queue"     : "mpi41",
        "memory"    : 5120,
        "nCPUs"     : "4"
    },

    "bwa_mem" :
    {
        "queue"     : "mpi21",
        "nCPUs"     : "8"
    },

    "MarkDuplicates" :
    {
        "queue"     : "fat2"
    },

    "BaseRecalibrator" :
    {
        "queue"     : "mpi41",
        "memory"    : 5120
    },

    "ApplyBQSR" :
    {
        "queue"     : "fat3"
    },

    "HaplotypeCaller" :
    {
        "queue"     : "fat2"
    },

    "CombineGVCFs" :
    {
        "queue"     : "fat2",
        "nCPUs"     : "10",
        "memory"    : 100000
    },

    "GenotypeGVCFs" :
    {
        "queue"     : "fat2",
        "nCPUs"     : "10",
        "memory"    : 100000
    },

    "VQSR_snp" :
    {
        "queue"     : "fat2",
        "nCPUs"     : "5",
        "memory"    : 10000
    },

    "VQSR_indel" :
    {
        "queue"     : "fat2",
        "nCPUs"     : "5",
        "memory"    : 10000
    },

    "VQSR_merge" :
    {
        "queue"     : "fat2",
        "nCPUs"     : "5",
        "memory"    : 10000
    }
}
```



#### Quote-escaping

注意：有些LSF命令需要多级引用转义。



### 提交任务



在服务器上进行试运行

```
module load anaconda/3-2021.11 
source activate
conda activate wes
# 试运行，得到流程，检查代码
snakemake -s snakemake_wes.py -n -p
# 绘制流程图，检查rule顺序
snakemake -s snakemake_wes.py --dag | dot -Tpdf > test_dag.pdf
```





```shell
nohup snakemake -s snakemake_wes.py --profile lsf --rerun-incomplete --cluster-config cluster.json 2>&1 &
```

参数：

* --profile lsf：提供profile配置信息的文件夹
* --rerun-incomplete：如果存在不完整的文件，重新运行
* --cluster-config：指定cluster配置文件，由json或yaml存储。



运行任务时，若需要中断所有任务，重新运行：

```
bjobs | awk '{print $1}' > job.list
for i in `cat job.list`; do bkill $i; done
snakemake -s snakemake_wes.py --unlock
```

