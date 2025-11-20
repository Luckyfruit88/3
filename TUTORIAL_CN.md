# ChIP-seq 分析流程教程

## 概述

本教程提供使用 Nextflow 运行完整 ChIP-seq 分析流程的详细步骤说明。该流程处理单端 ChIP-seq 数据，从原始测序数据到质量控制、比对、峰调用和注释。

## 流程概览

该流程执行以下步骤：
1. 使用 FastQC 进行质量控制
2. 使用 Trimmomatic 修剪接头序列
3. 使用 Bowtie2 进行序列比对
4. BAM 文件处理（排序、索引、统计）
5. 生成 BigWig 覆盖度轨迹
6. 样本相关性分析
7. 使用 HOMER 进行峰调用
8. 识别可重复的峰
9. 过滤黑名单区域
10. 峰注释
11. 基因信号可视化（跨基因体）
12. Motif 富集分析

## 环境要求

### 软件需求

流程使用 Singularity 容器，您需要：
- Nextflow（版本 20.0 或更高）
- Singularity（用于容器执行）
- 计算集群访问权限（可选但推荐）

### 输入数据

流程需要：
- 样本表 CSV 文件，包含列：`name`、`path`
- 参考基因组 FASTA 文件
- GTF 注释文件
- 接头序列文件
- ENCODE 黑名单 BED 文件
- UCSC 基因 BED 文件（TSS/TTS 位置）

### 分析环境

用于结果可视化和解释：
- Conda 环境（在 `environment.yml` 中指定）
- Jupyter notebook，包含 Python 3.10、pandas、matplotlib、seaborn
- 交互式分析笔记本（`analysis.ipynb`）

## 快速开始

### 1. 基础测试运行（桩模式）

在不运行实际计算的情况下测试工作流逻辑：

```bash
nextflow run main.nf -stub-run -profile local
```

这将：
- 验证所有进程定义
- 测试通道逻辑
- 创建空的桩输出文件
- 在几秒钟内完成

### 2. 开发运行（子采样数据）

使用子采样数据进行更快的测试：

```bash
nextflow run main.nf -profile singularity,local
```

默认使用 `params.subsampled_samplesheet`，指向较小的 FASTQ 文件。

### 3. 生产运行（完整数据）

使用完整数据集进行完整分析：

```bash
nextflow run main.nf -profile singularity,cluster
```

## 配置说明

### 样本表格式

流程需要以下格式的 CSV 文件：

```csv
name,path
INPUT_rep1,/path/to/INPUT_rep1.fastq.gz
INPUT_rep2,/path/to/INPUT_rep2.fastq.gz
IP_rep1,/path/to/IP_rep1.fastq.gz
IP_rep2,/path/to/IP_rep2.fastq.gz
```

**重要的命名约定：**
- 样本必须以 `INPUT` 或 `IP` 开头
- 必须包含重复标识符：`rep1`、`rep2` 等
- 格式：`{类型}_rep{编号}` 或 `{类型}_rep{编号}_subset`

### 在数据集之间切换

要从子采样切换到完整数据，编辑 `main.nf` 第 25 行：

```groovy
# 从这个改为：
Channel.fromPath(params.subsampled_samplesheet)

# 改为这个：
Channel.fromPath(params.samplesheet)
```

### 关键参数

编辑 `nextflow.config` 进行自定义：

```groovy
params {
    // 参考文件
    genome = "/path/to/genome.fa"
    gtf = "/path/to/annotation.gtf"
    adapter_fa = "/path/to/adapters.fa"
    blacklist = "/path/to/blacklist.bed"
    
    // 样本表
    samplesheet = "$projectDir/full_samplesheet.csv"
    subsampled_samplesheet = "$projectDir/subsampled_samplesheet.csv"
    
    // 输出目录
    outdir = "$projectDir/results/"
}
```

## 执行配置文件

### 本地配置

在单台机器上运行：

```bash
nextflow run main.nf -profile singularity,local
```

- 使用本地执行器
- 每个进程使用单个 CPU
- 适合测试

### 集群配置

在 SGE 集群上运行：

```bash
nextflow run main.nf -profile singularity,cluster
```

- 使用 SGE 执行器
- 根据进程标签分配资源：
  - `process_single`: 1 CPU
  - `process_low`: 4 CPUs
  - `process_medium`: 8 CPUs
  - `process_high`: 16 CPUs

## 流程工作流详情

### 阶段 1：质量控制和预处理

```
原始 FASTQ → FastQC (QC) → Trimmomatic (修剪接头) → 修剪后的 FASTQ
```

**输出：**
- `results/fastqc/` - FastQC HTML 和 ZIP 报告
- `results/trimmomatic/` - 修剪后的 FASTQ 文件和日志

### 阶段 2：比对

```
修剪后的 FASTQ → Bowtie2 (比对) → SAM → Samtools (转换) → BAM
                       ↑
                Bowtie2 索引（只构建一次）
```

**输出：**
- `refs/bowtie2_index/` - 基因组索引文件
- `results/bowtie2_align/` - BAM 文件和比对日志

### 阶段 3：BAM 处理

```
BAM → Samtools Sort → 排序的 BAM → Samtools Index → 索引的 BAM
                           ↓
                    Samtools Flagstat (统计)
```

**输出：**
- `results/samtools_sort/` - 排序的 BAM 文件
- `results/samtools_idx/` - BAM 索引文件 (.bai)
- `results/samtools_flagstat/` - 比对统计

### 阶段 4：覆盖度轨迹和相关性

```
索引的 BAM → BamCoverage → BigWig 文件 → MultiBigwigSummary → 矩阵
                                                                    ↓
                                                            PlotCorrelation → 热图
```

**输出：**
- `results/bamcoverage/` - BigWig 覆盖度轨迹 (.bw)
- `results/multibwsummary/` - 相关性矩阵（NPZ 格式）
- `results/plotcorrelation/` - 相关性热图（PNG）和矩阵（TSV）

### 阶段 5：峰调用

```
索引的 BAM → makeTagDirectory → 标签目录
                                       ↓
                     IP + INPUT 按重复配对
                                       ↓
                              findPeaks (HOMER)
                                       ↓
                              峰 TXT 文件 → pos2bed → BED 文件
```

**输出：**
- `results/homer_tagdir/` - HOMER 标签目录
- `results/homer_findpeaks/` - HOMER 格式的峰文件 (.txt)
- `results/homer_pos2bed/` - BED 格式的峰文件 (.bed)

### 阶段 6：可重复峰和过滤

```
rep1_peaks.bed ──┐
                 ├─→ Bedtools Intersect → 可重复峰
rep2_peaks.bed ──┘                              ↓
                                    Bedtools Remove (黑名单)
                                                ↓
                                    过滤后的可重复峰
```

**输出：**
- `results/bedtools_intersect/` - 可重复峰 (repr_peaks.bed)
- `results/bedtools_remove/` - 过滤后的峰 (repr_peaks_filtered.bed)

### 阶段 7：峰注释

```
过滤后的峰 + 基因组 + GTF → annotatePeaks → 注释的峰
```

**输出：**
- `results/homer_annotate/` - 注释的峰 (annotated_peaks.txt)

### 阶段 8：基因信号可视化

```
BigWig（仅 IP）+ 基因 BED → computeMatrix → 矩阵 → plotProfile → 信号图
```

**输出：**
- `results/computematrix/` - IP 样本的读数覆盖度矩阵
- `results/plotprofile/` - 基因信号分布图（PNG）

**注意：** 仅处理 IP 样本进行基因信号可视化，因为 INPUT 样本代表背景噪音。

### 阶段 9：Motif 富集分析

```
过滤后的峰 + 基因组 → findMotifsGenome → Motif 结果
```

**输出：**
- `results/homer_motifs/motifs/` - HOMER motif 分析结果
  - `homerResults.html` - 主要结果页面
  - `knownResults.html` - 已知 motif 匹配
  - Motif 标识和统计

### 阶段 10：QC 汇总

```
FastQC 报告 ──┐
Trimmomatic 日志 ├─→ MultiQC → 综合 QC 报告
Bowtie2 日志    ─┤
Flagstat 文件 ──┘
```

**输出：**
- `results/multiqc/` - MultiQC HTML 报告和数据目录

## 监控流程执行

### 查看进度

Nextflow 提供实时进度监控：

```bash
nextflow run main.nf -profile singularity,cluster
```

输出显示：
```
executor >  sge (12)
[ab/cd1234] process > FASTQC (1)           [100%] 4 of 4 ✔
[ef/gh5678] process > TRIM (2)             [100%] 4 of 4 ✔
[ij/kl9012] process > BOWTIE2_BUILD        [100%] 1 of 1 ✔
...
```

### 恢复失败的运行

如果流程失败，从最后一个成功的步骤恢复：

```bash
nextflow run main.nf -profile singularity,cluster -resume
```

Nextflow 使用缓存来跳过已完成的进程。

## 输出结构

```
project-3-Luckyfruit88/
├── main.nf                      # 主工作流文件
├── nextflow.config              # 配置文件
├── modules/                     # 进程定义
├── results/                     # 流程输出
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bowtie2_align/
│   ├── samtools_sort/
│   ├── samtools_idx/
│   ├── samtools_flagstat/
│   ├── bamcoverage/
│   ├── multibwsummary/
│   ├── plotcorrelation/
│   ├── homer_tagdir/
│   ├── homer_findpeaks/
│   ├── homer_pos2bed/
│   ├── bedtools_intersect/
│   ├── bedtools_remove/
│   ├── homer_annotate/
│   └── multiqc/
├── work/                        # 临时工作目录
└── .nextflow/                   # Nextflow 元数据
```

## 故障排除

### 常见问题

**1. 容器拉取错误**
```bash
# 预先拉取容器
singularity pull docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
```

**2. 内存错误**
- 在 `nextflow.config` 中增加内存分配
- 使用子采样数据进行测试

**3. 文件路径错误**
- 在配置中使用绝对路径
- 验证文件权限

**4. 样本配对问题**
- 检查样本命名约定
- 确保 `IP_rep1` 与 `INPUT_rep1` 配对

### 调试

**检查进程日志：**
```bash
# 查找失败的进程目录
ls -lt work/
cd work/ab/cd1234...

# 查看日志文件
cat .command.log
cat .command.err
```

**验证样本表：**
```bash
# 测试 CSV 解析
nextflow run main.nf -stub-run -profile local
```

## 高级用法

### 自定义资源分配

在 `nextflow.config` 中编辑进程标签：

```groovy
process {
    withLabel: process_high {
        cpus = 32
        memory = '64 GB'
    }
}
```

### 自定义相关性方法

在 `modules/deeptools_plotcorrelation/main.nf` 中更改相关性方法：

```groovy
plotCorrelation -in ${matrix} \
    --corMethod pearson \  // 从 spearman 更改
    --whatToPlot heatmap \
    --plotFile correlation_plot.png
```

### 自定义峰调用参数

在 `modules/homer_findpeaks/main.nf` 中修改 HOMER findPeaks：

```groovy
findPeaks ${ip_tagdir} \
    -style factor \
    -i ${input_tagdir} \
    -F 3 \              // 富集倍数阈值
    -P 0.001 \          // P 值阈值
    -o ${rep}_peaks.txt
```

## 性能提示

1. **使用集群配置**进行生产运行
2. **启用恢复**使用 `-resume` 标志
3. **监控工作目录**大小并定期清理
4. **使用子采样数据**进行开发
5. **根据数据大小调整进程标签**

## 预期运行时间

### 子采样数据（本地配置）
- 总计：约 30-60 分钟
- 比对：每个样本约 5-10 分钟
- 峰调用：每个重复约 2-5 分钟

### 完整数据（集群配置）
- 总计：约 4-8 小时
- 比对：每个样本约 30-60 分钟
- 峰调用：每个重复约 15-30 分钟

## 结果解读

### 需要查看的关键输出

1. **MultiQC 报告** (`results/multiqc/multiqc_report.html`)
   - 整体数据质量
   - 比对率
   - 样本统计

2. **相关性热图** (`results/plotcorrelation/correlation_plot.png`)
   - 样本相似性
   - 预期：IP 样本聚在一起，INPUT 样本聚在一起

3. **注释的峰** (`results/homer_annotate/annotated_peaks.txt`)
   - 峰的基因组位置
   - 最近的基因
   - 峰特征

4. **最终峰集** (`results/bedtools_remove/repr_peaks_filtered.bed`)
   - 重复之间的可重复峰
   - 已过滤黑名单
   - 准备进行下游分析

5. **基因信号分布图** (`results/plotprofile/IP_rep*_signal_coverage.png`)
   - 跨基因体的信号分布
   - TSS、基因体和 TTS 的结合模式
   - 显示 2kb 侧翼区域

6. **Motif 富集** (`results/homer_motifs/motifs/homerResults.html`)
   - 富集的序列 motif
   - 已知 motif 匹配
   - 潜在的共同因子

### 使用分析笔记本

流程包含一个全面的 Jupyter 笔记本用于结果可视化：

```bash
# 设置 conda 环境
conda env create -f environment.yml
conda activate chipseq-analysis

# 启动 Jupyter notebook
jupyter notebook analysis.ipynb
```

**笔记本章节：**
- 质量控制摘要与 MultiQC 集成
- 相关性分析解释（Spearman 理由）
- 峰调用统计和可视化
- 可重复峰策略说明
- 峰注释分析
- 基因信号分布图显示
- Motif 富集结果
- 总结和结论

## 引用

如果您使用此流程，请引用以下工具：

- Nextflow: Di Tommaso et al., Nature Biotechnology 2017
- FastQC: Andrews S., Babraham Bioinformatics 2010
- Trimmomatic: Bolger et al., Bioinformatics 2014
- Bowtie2: Langmead and Salzberg, Nature Methods 2012
- SAMtools: Li et al., Bioinformatics 2009
- deepTools: Ramírez et al., Nucleic Acids Research 2016
- HOMER: Heinz et al., Molecular Cell 2010
- BEDTools: Quinlan and Hall, Bioinformatics 2010

## 支持

如有流程问题：
1. 查看上面的故障排除部分
2. 查看 `work/` 目录中的进程日志
3. 查阅各个工具文档
4. 验证配置参数

## 许可证

此流程作为 BF528 课程作业的一部分提供用于教育目的。
