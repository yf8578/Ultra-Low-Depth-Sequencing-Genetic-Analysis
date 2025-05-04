version 1.0

workflow ultra_low_depth_seq_pipeline {
  input {
    # 样本信息
    String sample_id
    File fastq_file
    
    # 参考文件数组
    Array[File] bwa_hg38_files  # BWA所需的参考基因组和索引文件
    Array[File] gatk_hg38_files # GATK所需的参考基因组和索引文件
    Array[File] gatk_bundle_files # GATK bundle文件（包含dbSNP等）
    
    # 主要参考文件名（用于实际操作）
    String bwa_reference_filename = "hg38.fa.masked"  # BWA使用的参考基因组文件名
    String gatk_reference_filename = "hg38.fasta"     # GATK使用的参考基因组文件名
    String dbsnp_filename = "dbsnp_146.hg38.withDict.vcf.gz"   # dbSNP文件名
    
    # 配置参数
    String platform = "IONTORRENT"
    Int fastp_threads = 8
    Int bwa_threads = 16
    Int samtools_threads = 8
    Int gatk_java_mem = 4
  }

  # 步骤1: fastp质控 - 单端模式
  call fastp_task {
    input:
      sample_id = sample_id,
      fastq = fastq_file,
      threads = fastp_threads
  }

  # 步骤2: BWA比对和处理
  call bwa_alignment_task {
    input:
      sample_id = sample_id,
      clean_fastq = fastp_task.clean_fastq,
      bwa_reference_files = bwa_hg38_files,
      bwa_reference_filename = bwa_reference_filename,
      platform = platform,
      threads = bwa_threads,
      samtools_threads = samtools_threads
  }

  # 步骤3: GATK碱基质量校正
  call gatk_bqsr_task {
    input:
      sample_id = sample_id,
      input_bam = bwa_alignment_task.sorted_rmdup_bam,
      input_bai = bwa_alignment_task.sorted_rmdup_bai,
      gatk_reference_files = gatk_hg38_files,
      gatk_bundle_files = gatk_bundle_files,
      reference_filename = gatk_reference_filename,
      dbsnp_filename = dbsnp_filename,
      java_mem = gatk_java_mem
  }

  # 步骤4: BAM统计和深度计算
  call bam_stats_task {
    input:
      sample_id = sample_id,
      bqsr_bam = gatk_bqsr_task.bqsr_bam,
      bqsr_bai = gatk_bqsr_task.bqsr_bai
  }

  # 步骤5: 结果整理
  call organize_results_task {
    input:
      sample_id = sample_id,
      bqsr_bam = gatk_bqsr_task.bqsr_bam,
      bqsr_bai = gatk_bqsr_task.bqsr_bai,
      bamstats = bam_stats_task.bamstats,
      depth_info = bam_stats_task.depth_info,
      fastp_html = fastp_task.html_report,
      fastp_json = fastp_task.json_report
  }

  # 定义工作流输出
  output {
    File final_bam = organize_results_task.results_bam
    File final_bai = organize_results_task.results_bai
    File final_bamstats = organize_results_task.results_bamstats
    File final_depth_info = organize_results_task.results_depth_info
    File results_archive = organize_results_task.results_archive
  }
}

# 定义fastp质控任务 - 单端模式
task fastp_task {
  input {
    String sample_id
    File fastq
    Int threads = 6
  }

  command <<<
    set -e
    
    # 创建临时输出目录
    mkdir -p temp_data
    
    # 运行fastp - 单端模式，自动检测接头
    /opt/conda/envs/tools/bin/fastp \
      -i ~{fastq} \
      -o temp_data/~{sample_id}.clean.fq.gz \
      --qualified_quality_phred=5 \
      --unqualified_percent_limit=50 \
      --n_base_limit=10 \
      --disable_trim_poly_g \
      --thread=~{threads} \
      -j temp_data/~{sample_id}.json \
      -h temp_data/~{sample_id}.html \
      -R "~{sample_id}"
  >>>

  output {
    File clean_fastq = "temp_data/~{sample_id}.clean.fq.gz"
    File html_report = "temp_data/~{sample_id}.html"
    File json_report = "temp_data/~{sample_id}.json"
  }

  runtime {
    docker_url: "stereonote_hpc/zhangyifan1_c1b2a81d272d442fa7c97d8a049b24f4_private:latest"
    req_cpu: threads
    req_memory: "4Gi"
  }
}

# 定义BWA比对任务 - 适用于单端数据，复制所有BWA参考文件
task bwa_alignment_task {
  input {
    String sample_id
    File clean_fastq
    Array[File] bwa_reference_files
    String bwa_reference_filename
    String platform
    Int threads = 4
    Int samtools_threads = 8
  }

  command <<<
    set -e
    
    # 创建临时目录
    mkdir -p temp_data
    mkdir -p bwa_reference
    
    # 复制所有BWA参考文件到工作目录
    echo "复制BWA参考文件..."
    for file in ~{sep=' ' bwa_reference_files}; do
      cp -v $file bwa_reference/
    done
    
    # 显示复制后bwa_reference目录中的文件列表
echo "=====================================查看bwa_reference目录内容==========================================="
ls -la bwa_reference/
echo "====================================================================================================="
     
    # 确定参考基因组文件路径
    reference_path="bwa_reference/~{bwa_reference_filename}"
    echo "使用参考基因组: $reference_path"
    
    # 检查参考文件是否存在
    if [ ! -f "$reference_path" ]; then
      echo "错误: 无法找到参考基因组文件: $reference_path"
      ls -la bwa_reference/
      exit 1
    fi
    
    # 确认参考基因组已索引，如果未索引则创建索引
    if [ ! -f ${reference_path}.bwt ]; then
      echo "参考基因组索引不存在，正在创建..."
      /opt/conda/envs/tools/bin/bwa index -a bwtsw ${reference_path}
    fi
    
    # BWA比对 - 单端比对
    echo "执行BWA比对..."
    /opt/conda/envs/tools/bin/bwa aln \
      -e 10 -t ~{threads} -i 5 -q 0 \
      ${reference_path} ~{clean_fastq} > temp_data/~{sample_id}.sai
    
    # 生成SAM文件并转换为BAM - 使用samse处理单端数据
    echo "生成BAM文件..."
    /opt/conda/envs/tools/bin/bwa samse \
      -r "@RG\tID:~{sample_id}\tPL:~{platform}\tSM:~{sample_id}" \
      ${reference_path} temp_data/~{sample_id}.sai ~{clean_fastq} | \
      /opt/conda/envs/tools/bin/samtools view -h -Sb - > temp_data/~{sample_id}.bam
    
    # 排序BAM文件
    echo "排序BAM文件..."
    /opt/conda/envs/tools/bin/samtools sort \
      -@ ~{samtools_threads} \
      -O bam \
      -o temp_data/~{sample_id}.sorted.bam \
      temp_data/~{sample_id}.bam
    
    # 去除重复
    echo "去除重复..."
    /opt/conda/envs/tools/bin/samtools rmdup \
      temp_data/~{sample_id}.sorted.bam \
      temp_data/~{sample_id}.sorted.rmdup.bam
    
    # 索引BAM文件
    echo "索引BAM文件..."
    /opt/conda/envs/tools/bin/samtools index \
      temp_data/~{sample_id}.sorted.rmdup.bam
      
    # 创建完成标记
    touch temp_data/bwa_sort_rmdup.finish
  >>>

  output {
    File sorted_rmdup_bam = "temp_data/~{sample_id}.sorted.rmdup.bam"
    File sorted_rmdup_bai = "temp_data/~{sample_id}.sorted.rmdup.bam.bai"
    File finish_flag = "temp_data/bwa_sort_rmdup.finish"
  }

  runtime {
    docker_url: "stereonote_hpc/zhangyifan1_c1b2a81d272d442fa7c97d8a049b24f4_private:latest"
    req_cpu: threads + samtools_threads
    req_memory: "8Gi"
  }
}

# 定义GATK碱基质量校正任务 - 复制GATK参考文件和bundle文件
task gatk_bqsr_task {
  input {
    String sample_id
    File input_bam
    File input_bai
    Array[File] gatk_reference_files
    Array[File] gatk_bundle_files
    String reference_filename
    String dbsnp_filename
    Int java_mem = 4
  }

  command <<<
    set -e
    
    # 创建临时目录
    mkdir -p temp_data
    mkdir -p gatk_reference
    mkdir -p gatk_bundle
    
    # 复制GATK参考文件到工作目录
    echo "复制GATK参考文件..."
    for file in ~{sep=' ' gatk_reference_files}; do
      cp -v $file gatk_reference/
    done
    
    # 显示复制后gatk_reference目录中的文件列表
echo "=====================================查看gatk_reference目录内容==========================================="
ls -la gatk_reference/
echo "====================================================================================================="
    
    
    
    # 复制GATK bundle文件到工作目录
    echo "复制GATK bundle文件..."
    for file in ~{sep=' ' gatk_bundle_files}; do
      cp -v $file gatk_bundle/
    done
    
    # 显示复制后gatk_bundle目录中的文件列表
echo "=====================================查看gatk_bundle目录内容============================================="
ls -la gatk_bundle/
echo "====================================================================================================="
    
    
    # 确定参考基因组和dbSNP文件路径
    reference_path="gatk_reference/~{reference_filename}"
    dbsnp_path="gatk_bundle/~{dbsnp_filename}"
    
    echo "使用GATK参考基因组: $reference_path"
    echo "使用dbSNP文件: $dbsnp_path"
    
    # 检查参考文件是否存在
    if [ ! -f "$reference_path" ]; then
      echo "错误: 无法找到GATK参考基因组文件: $reference_path"
      ls -la gatk_reference/
      exit 1
    fi
    
    # 检查dbSNP文件是否存在
    if [ ! -f "$dbsnp_path" ]; then
      echo "错误: 无法找到dbSNP文件: $dbsnp_path"
      ls -la gatk_bundle/
      # 尝试查找替代的dbSNP文件
      alt_dbsnp=$(find gatk_bundle/ -name "*.vcf.gz" | head -1)
      if [ ! -z "$alt_dbsnp" ]; then
        echo "使用替代的dbSNP文件: $alt_dbsnp"
        dbsnp_path=$alt_dbsnp
      else
        exit 1
      fi
    fi
    
    # 检查参考基因组索引
    if [ ! -f ${reference_path}.fai ]; then
      echo "创建参考基因组索引..."
      /opt/conda/envs/tools/bin/samtools faidx ${reference_path}
    fi
    
    # 检查参考基因组字典
    ref_dict=$(echo "${reference_path}" | sed 's/\.[^.]*$/\.dict/')
    if [ ! -f "$ref_dict" ]; then
      echo "创建参考基因组字典..."
      /opt/conda/envs/tools/bin/gatk CreateSequenceDictionary -R ${reference_path}
    fi
    
    # 检查dbSNP索引
    if [ ! -f ${dbsnp_path}.tbi ] && [ ! -f ${dbsnp_path}.idx ]; then
      echo "为dbSNP创建索引..."
      /opt/conda/envs/tools/bin/gatk IndexFeatureFile -I ${dbsnp_path}
    fi
    
    # 运行BaseRecalibrator
echo "执行BaseRecalibrator..."
# 首先设置JAVA_HOME环境变量
export JAVA_HOME=/opt/conda/envs/tools
export PATH=$JAVA_HOME/bin:$PATH

# 然后运行GATK命令
/opt/conda/envs/tools/bin/gatk --java-options "-Xmx~{java_mem}g -XX:+UseParallelGC" \
  BaseRecalibrator \
  -R ${reference_path} \
  -I ~{input_bam} \
  --known-sites ${dbsnp_path} \
  -O temp_data/~{sample_id}.recal_data.table

# 创建完成标记
touch temp_data/baseRecalibrator.finish

# 运行ApplyBQSR
echo "执行ApplyBQSR..."
/opt/conda/envs/tools/bin/gatk --java-options "-Xmx~{java_mem}g -XX:+UseParallelGC" \
  ApplyBQSR \
  -R ${reference_path} \
  --bqsr-recal-file temp_data/~{sample_id}.recal_data.table \
  -I ~{input_bam} \
  -O temp_data/~{sample_id}.sorted.rmdup.BQSR.bam

# 创建完成标记  
touch temp_data/PrintReads.finish

# 索引BQSR处理后的BAM
echo "索引BQSR BAM文件..."
/opt/conda/envs/tools/bin/samtools index \
  temp_data/~{sample_id}.sorted.rmdup.BQSR.bam
  
# 创建完成标记
touch temp_data/bam_index.finish
  >>>

  output {
    File bqsr_bam = "temp_data/~{sample_id}.sorted.rmdup.BQSR.bam"
    File bqsr_bai = "temp_data/~{sample_id}.sorted.rmdup.BQSR.bam.bai"
    File recal_table = "temp_data/~{sample_id}.recal_data.table"
    File baserecal_finish = "temp_data/baseRecalibrator.finish"
    File printreads_finish = "temp_data/PrintReads.finish"
    File bamindex_finish = "temp_data/bam_index.finish"
  }

  runtime {
    docker_url: "stereonote_hpc/zhangyifan1_c1b2a81d272d442fa7c97d8a049b24f4_private:latest"
    req_cpu: 4
    req_memory: "~{java_mem + 4}Gi"
  }
}

# 定义BAM统计和深度计算任务
task bam_stats_task {
  input {
    String sample_id
    File bqsr_bam
    File bqsr_bai
  }

  command <<<
    set -e
    
    # 直接在命令中定义
    genome_size=3095693981
    
    # 生成BAM统计信息
    echo "生成BAM统计信息..."
    /opt/conda/envs/tools/bin/samtools stats ~{bqsr_bam} > ~{sample_id}.sorted.rmdup.BQSR.bamstats
    
    # 创建完成标记
    touch bamstats.finish
    
    # 提取比对碱基数并计算深度
    echo "计算测序深度..."
    bases_mapped=$(/opt/conda/envs/tools/bin/samtools stats ~{bqsr_bam} | awk '/bases mapped \(cigar\)/ {print $4}')
    
    # 如果无法提取，尝试备选方法
    if ! [[ "$bases_mapped" =~ ^[0-9]+$ ]]; then
      bases_mapped=$(/opt/conda/envs/tools/bin/samtools stats ~{bqsr_bam} | grep -o "bases mapped (cigar):[[:space:]]*[0-9]*" | grep -o "[0-9]*$")
    fi
    
    # 如果仍然失败，使用默认值0
    if ! [[ "$bases_mapped" =~ ^[0-9]+$ ]]; then
      bases_mapped=0
      echo "警告: 无法提取有效碱基数，使用默认值0"
    fi
    
    # 计算深度 - 使用Shell变量语法
    depth=$(awk -v bases="$bases_mapped" -v genome="$genome_size" 'BEGIN {printf "%.6f", bases/genome}')
    depth_percentage=$(awk -v d="$depth" 'BEGIN {printf "%.4f", d*100}')
    
    echo "样本 ~{sample_id} 的测序深度: ${depth}X (${depth_percentage}%)"
    echo "比对碱基数: $bases_mapped"
    echo "参考基因组大小: $genome_size bp"
    
    # 保存深度信息
    echo "~{sample_id},${depth},${bases_mapped},${depth_percentage}%" > ~{sample_id}.depth_info.csv
  >>>

  output {
    File bamstats = "~{sample_id}.sorted.rmdup.BQSR.bamstats"
    File depth_info = "~{sample_id}.depth_info.csv"
    File bamstats_finish = "bamstats.finish"
  }

  runtime {
    docker_url: "stereonote_hpc/zhangyifan1_c1b2a81d272d442fa7c97d8a049b24f4_private:latest"
    req_cpu: 2
    req_memory: "4Gi"
  }
}

# 定义结果整理任务
task organize_results_task {
  input {
    String sample_id
    File bqsr_bam
    File bqsr_bai
    File bamstats
    File depth_info
    File fastp_html
    File fastp_json
  }

  command <<<
    set -e
    
    # 创建结果目录
    mkdir -p ~{sample_id}_results
    
    # 复制所有结果文件
    cp ~{bqsr_bam} ~{sample_id}_results/~{sample_id}.sorted.rmdup.BQSR.bam
    cp ~{bqsr_bai} ~{sample_id}_results/~{sample_id}.sorted.rmdup.BQSR.bam.bai
    cp ~{bamstats} ~{sample_id}_results/~{sample_id}.sorted.rmdup.BQSR.bamstats
    cp ~{depth_info} ~{sample_id}_results/~{sample_id}.depth_info.csv
    cp ~{fastp_html} ~{sample_id}_results/~{sample_id}.html
    cp ~{fastp_json} ~{sample_id}_results/~{sample_id}.json
    
    # 创建摘要文件
    echo "样本ID: ~{sample_id}" > ~{sample_id}_results/summary.txt
    echo "处理完成时间: $(date)" >> ~{sample_id}_results/summary.txt
    echo "测序深度: $(cat ~{depth_info} | cut -d',' -f2)X" >> ~{sample_id}_results/summary.txt
    
    # 打包结果
    tar -czvf ~{sample_id}_results.tar.gz ~{sample_id}_results
    
    # 创建完成标记
    touch move2final.finish
  >>>

  output {
    File results_archive = "~{sample_id}_results.tar.gz"
    File results_bam = "~{sample_id}_results/~{sample_id}.sorted.rmdup.BQSR.bam"
    File results_bai = "~{sample_id}_results/~{sample_id}.sorted.rmdup.BQSR.bam.bai"
    File results_bamstats = "~{sample_id}_results/~{sample_id}.sorted.rmdup.BQSR.bamstats"
    File results_depth_info = "~{sample_id}_results/~{sample_id}.depth_info.csv"
    File move2final_finish = "move2final.finish"
  }

  runtime {
    docker_url: "stereonote_hpc/zhangyifan1_c1b2a81d272d442fa7c97d8a049b24f4_private:latest"
    req_cpu: 2
    req_memory: "4Gi"
  }
}