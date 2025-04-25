"""
Author: zhangyifan1
Date: 2025-04-23 22:12:17
LastEditors: zhangyifan1 zhangyifan1@genomics.cn
LastEditTime: 2025-04-24 00:05:09
FilePath: //Ultra-Low-Depth-Sequencing-Genetic-Analysis//01_bam_QC.py
Description:

"""

# !/usr/bin/env python3
"""
BAM文件基础信息统计脚本
"""

import os
import sys
import subprocess
import re
import logging
from typing import Dict, Any, List

# 设置日志
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def create_output_dir(output_dir: str) -> bool:
    """创建输出目录，如果不存在"""
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logger.info(f"创建输出目录: {output_dir}")
        return True
    except Exception as e:
        logger.error(f"创建输出目录失败: {e}")
        return False


def check_file_exists(filepath: str) -> bool:
    """检查文件是否存在并且是BAM文件"""
    if not os.path.exists(filepath):
        logger.error(f"文件不存在: {filepath}")
        return False

    if not filepath.lower().endswith(".bam"):
        logger.error(f"文件不是BAM格式: {filepath}")
        return False

    return True


def get_bam_stats(filepath: str) -> Dict[str, Any]:
    """
    获取BAM文件的基础统计信息

    参数:
        filepath: BAM文件路径
    返回:
        统计信息字典
    """
    stats = {
        "文件路径": filepath,
        "文件大小(MB)": round(os.path.getsize(filepath) / (1024 * 1024), 2),
    }

    try:
        # 验证文件是否是有效的BAM文件
        try:
            cmd = ["samtools", "quickcheck", filepath]
            subprocess.run(
                cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            logger.error(f"无效的BAM文件: {filepath}")
            stats["错误"] = "无效的BAM文件"
            return stats

        # 获取总读数
        cmd = ["samtools", "view", "-c", filepath]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stats["总读数"] = int(result.stdout.strip())

        # 获取比对率等信息
        cmd = ["samtools", "flagstat", filepath]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # 解析flagstat结果
        for line in result.stdout.splitlines():
            if "mapped (" in line:
                match = re.search(r"(\d+\.\d+)%", line)
                if match:
                    stats["比对率"] = f"{match.group(1)}%"

            # 提取duplicate信息
            if "duplicates" in line:
                match = re.search(r"(\d+) \+ \d+ duplicates", line)
                if match:
                    stats["duplicate数"] = int(match.group(1))
                    # 计算duplicate率
                    if stats["总读数"] > 0:
                        dup_rate = (stats["duplicate数"] / stats["总读数"]) * 100
                        stats["duplicate率"] = f"{dup_rate:.2f}%"

        # 获取平均深度
        try:
            cmd = ["samtools", "coverage", filepath]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            lines = result.stdout.strip().split("\n")

            if len(lines) > 1:
                # 提取平均深度列
                headers = lines[0].split("\t")
                mean_depth_idx = (
                    headers.index("meandepth") if "meandepth" in headers else -1
                )

                if mean_depth_idx >= 0:
                    depths = []
                    for line in lines[1:]:
                        cols = line.split("\t")
                        if len(cols) > mean_depth_idx:
                            depths.append(float(cols[mean_depth_idx]))

                    if depths:
                        stats["平均深度"] = round(sum(depths) / len(depths), 2)
        except Exception as e:
            logger.warning(f"获取覆盖度信息失败: {e}")

        # 计算GC含量和低质量reads数
        try:
            # 采样一部分reads分析GC含量和质量
            sample_size = min(10000, stats["总读数"])
            cmd = ["samtools", "view", filepath]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            lines = result.stdout.splitlines()[:sample_size]

            total_bases = 0
            gc_bases = 0
            low_qual_reads = 0
            qual_threshold = 20  # 质量阈值，通常Q20被认为是低质量阈值

            for line in lines:
                fields = line.split("\t")
                if len(fields) >= 10:
                    seq = fields[9]
                    qual = fields[10]

                    # 计算GC含量
                    total_bases += len(seq)
                    gc_bases += seq.count("G") + seq.count("C")

                    # 计算低质量reads
                    qual_values = [ord(q) - 33 for q in qual]  # 将ASCII质量值转换为数值
                    if sum(qual_values) / len(qual_values) < qual_threshold:
                        low_qual_reads += 1

            if total_bases > 0:
                stats["GC含量"] = f"{(gc_bases / total_bases) * 100:.2f}%"

            if len(lines) > 0:
                stats["低质量reads比例"] = f"{(low_qual_reads / len(lines)) * 100:.2f}%"
                stats["低质量reads数"] = int(
                    (low_qual_reads / len(lines)) * stats["总读数"]
                )

        except Exception as e:
            logger.warning(f"计算GC含量和质量信息失败: {e}")

    except Exception as e:
        logger.error(f"获取BAM统计信息时出错: {e}")
        stats["错误"] = str(e)

    return stats


def process_bam_file_list(file_list_path: str, output_dir: str = "bam_stats") -> None:
    """
    处理BAM文件列表，统计信息

    参数:
        file_list_path: 包含BAM文件路径列表的文件
        output_dir: 输出目录名称
    """
    # 创建输出目录
    if not create_output_dir(output_dir):
        return

    bam_files = []
    results = []

    try:
        # 读取文件列表
        if os.path.isfile(file_list_path):
            with open(file_list_path, "r") as f:
                bam_files = [
                    line.strip()
                    for line in f
                    if line.strip() and line.strip().lower().endswith(".bam")
                ]
        else:
            # 如果提供的是单个BAM文件路径
            if file_list_path.lower().endswith(".bam") and os.path.isfile(
                file_list_path
            ):
                bam_files = [file_list_path]
            else:
                logger.error(f"无效的输入: {file_list_path}")
                return

        if not bam_files:
            logger.error("没有找到任何BAM文件")
            return

        print(f"读取了 {len(bam_files)} 个BAM文件路径")

        # 处理每个BAM文件
        for i, filepath in enumerate(bam_files, 1):
            if check_file_exists(filepath):
                print(f"处理BAM文件 {i}/{len(bam_files)}: {filepath}")

                # 获取统计信息
                stats = get_bam_stats(filepath)

                # 将结果保存到单独的文件
                file_basename = os.path.basename(filepath)
                stats_file = os.path.join(output_dir, f"{file_basename}_stats.txt")
                with open(stats_file, "w") as f:
                    for key, value in stats.items():
                        f.write(f"{key}: {value}\n")

                print(f"  统计信息已保存到: {stats_file}")
                results.append(stats)
            else:
                print(f"跳过无效文件: {filepath}")

        # 如果没有有效结果，退出
        if not results:
            logger.error("没有生成任何统计结果")
            return

        # 输出汇总结果到文本文件
        summary_file = os.path.join(output_dir, "summary_stats.txt")
        with open(summary_file, "w") as f:
            for i, stat in enumerate(results, 1):
                f.write(f"=== BAM文件 {i}/{len(results)} ===\n")
                for key, value in stat.items():
                    f.write(f"{key}: {value}\n")
                f.write("\n")

        print(f"\n汇总统计结果已保存到 {summary_file}")

        # 创建CSV格式摘要
        csv_file = os.path.join(output_dir, "bam_stats_summary.csv")
        with open(csv_file, "w") as f:
            # 确定所有可能的字段
            all_fields = set()
            for stat in results:
                all_fields.update(stat.keys())

            # 按照特定顺序排列关键字段
            priority_fields = [
                "文件路径",
                "总读数",
                "平均深度",
                "比对率",
                "GC含量",
                "duplicate数",
                "duplicate率",
                "低质量reads数",
                "低质量reads比例",
            ]

            # 确保关键字段在前面，其他字段按字母顺序排列
            headers = [field for field in priority_fields if field in all_fields]
            other_fields = sorted(list(all_fields - set(headers)))
            headers.extend(other_fields)

            # 写入表头
            f.write(",".join(f'"{h}"' for h in headers) + "\n")

            # 写入每个文件的数据
            for stat in results:
                row = []
                for field in headers:
                    value = stat.get(field, "")
                    # 处理可能含有逗号的字段
                    if isinstance(value, str) and ("," in value or '"' in value):
                        value = f'"{value.replace("\"", "\"\"")}"'
                    row.append(str(value))
                f.write(",".join(row) + "\n")

        print(f"CSV格式摘要已保存到 {csv_file}")

        # 打印简短摘要
        print("\nBAM文件统计结果摘要:")
        for stat in results:
            path = stat.get("文件路径", "Unknown")
            reads = stat.get("总读数", "N/A")
            align_rate = stat.get("比对率", "N/A")
            depth = stat.get("平均深度", "N/A")
            dup_rate = stat.get("duplicate率", "N/A")
            gc = stat.get("GC含量", "N/A")
            print(
                f"{os.path.basename(path)}: 读数: {reads}, 比对率: {align_rate}, "
                f"平均深度: {depth}, Duplicate率: {dup_rate}, GC含量: {gc}"
            )

    except Exception as e:
        logger.error(f"处理BAM文件列表时出错: {e}")


def main():
    """主函数"""
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("用法: python bam_stats.py <BAM文件列表路径或单个BAM文件> [输出目录]")
        print("BAM文件列表应包含每行一个BAM文件的绝对路径")
        return

    file_input = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) == 3 else "bam_stats"

    process_bam_file_list(file_input, output_dir)


if __name__ == "__main__":
    main()
