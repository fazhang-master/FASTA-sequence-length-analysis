#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from datetime import datetime

# ===================== 中文显示配置 =====================
mpl.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun']  # 中文字体优先级列表
mpl.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常

def count_long_sequences(fasta_file, min_length=428):
    """
    统计FASTA文件中长度大于指定阈值的序列
    参数:
        fasta_file: FASTA文件路径
        min_length: 最小长度阈值（默认428 aa）
    返回:
        long_count: 长序列数量
        total_count: 总序列数量
        long_ids: 长序列ID列表
        all_lengths: 所有序列长度列表
    """
    long_count = total_count = 0
    long_ids, all_lengths = [], []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_count += 1
        cleaned_seq = str(record.seq).replace("-", "")
        seq_length = len(cleaned_seq)
        all_lengths.append(seq_length)
        
        if seq_length > min_length:
            long_count += 1
            long_ids.append(record.id)
    
    return long_count, total_count, long_ids, all_lengths

def set_chinese_font(font_path=None):
    """动态设置中文字体路径（跨平台兼容）"""
    import matplotlib.font_manager as fm
    
    # 1. 优先使用用户指定路径
    if font_path and os.path.isfile(font_path):
        return fm.FontProperties(fname=font_path)
    
    # 2. 自动搜索常见系统字体路径
    system_fonts = {
        'Windows': 'C:/Windows/Fonts/msyh.ttc',  # 微软雅黑
        'Linux': '/usr/share/fonts/truetype/msttcorefonts/msyh.ttf',
        'Darwin': '/System/Library/Fonts/PingFang.ttc'  # macOS
    }
    candidate = system_fonts.get(sys.platform, '')
    if os.path.isfile(candidate):
        return fm.FontProperties(fname=candidate)
    
    # 3. 回退到Matplotlib内置字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']  # 支持中文的基本字体
    return fm.FontProperties()

def plot_length_distribution(all_lengths, min_length=428, output_file="sequence_length_plot.png", font_prop=None):
    """
    绘制专业的序列长度分布图
    参数:
        all_lengths: 序列长度列表
        min_length: 最小长度阈值
        output_file: 图表输出路径
    """
    # 创建画布
    plt.figure(figsize=(14, 6), dpi=120)
    
    # ===== 直方图 =====
    plt.subplot(1, 2, 1)
    n, bins, patches = plt.hist(
        all_lengths, bins=50, 
        color='#4C72B0', edgecolor='black', 
        alpha=0.8, density=False,
    )
    plt.axvline(
        x=min_length, color='#C44E52', 
        linestyle='--', linewidth=2.5,
        label=f'{min_length} aa'
    )
    
    plt.title('序列长度分布直方图', fontsize=14, fontweight='bold', fontproperties=font_prop)
    plt.xlabel('序列长度 (aa)', fontsize=12, fontproperties=font_prop)
    plt.ylabel('序列数量', fontsize=12, fontproperties=font_prop)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    
    # ===== 箱线图 =====
    plt.subplot(1, 2, 2)
    box = plt.boxplot(
        all_lengths, vert=False, patch_artist=True,
        boxprops=dict(facecolor='#55A868', alpha=0.8),
        medianprops=dict(color='yellow', linewidth=2),
        flierprops=dict(marker='o', markersize=4)
    )
    
    # 添加统计标注
    stats = {
        '最短': min(all_lengths),
        '最长': max(all_lengths),
        '中位数': np.median(all_lengths),
        'Q1': np.percentile(all_lengths, 25),
        'Q3': np.percentile(all_lengths, 75)
    }
    
    plt.text(
        min(all_lengths)-50, 1.25, 
        f"最短: {stats['最短']} aa", fontsize=10, 
        fontproperties=font_prop
    )
    plt.text(
        max(all_lengths)+50, 1.25, 
        f"最长: {stats['最长']} aa", fontsize=10, ha='left', 
        fontproperties=font_prop
    )
    plt.text(
        stats['中位数'], 0.7, 
        f"中位数: {stats['中位数']:.0f} aa",
        fontsize=10, ha='center', 
        fontproperties=font_prop
    )
    
    plt.title('序列长度箱线图', fontsize=14, fontweight='bold', fontproperties=font_prop)
    plt.xlabel('序列长度 (aa)', fontsize=12, fontproperties=font_prop)
    
    # 保存图表
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_report(args, long_count, total_count, long_ids, all_lengths):
    """
    生成专业的分析报告
    参数:
        args: 命令行参数
        long_count: 长序列数量
        total_count: 总序列数量
        long_ids: 长序列ID列表
        all_lengths: 所有序列长度列表
    返回:
        报告保存路径的提示信息
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    percentage = long_count / total_count * 100 if total_count > 0 else 0
    
    # 创建输出目录
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    
    with open(args.output, "w", encoding='utf-8') as f:
        # ===== 报告头部 =====
        f.write(f"{'序列长度分析报告':^80}\n")
        f.write("=" * 80 + "\n")
        f.write(f"{'生成时间:':<15}{timestamp}\n")
        f.write(f"{'输入文件:':<15}{args.input}\n")
        f.write(f"{'输出报告:':<15}{args.output}\n")
        f.write(f"{'长度阈值:':<15}{args.min_length} aa\n")
        f.write("-" * 80 + "\n\n")
        
        # ===== 主要统计结果 =====
        f.write(f"{'总序列数:':<20}{total_count}\n")
        f.write(f"{'长序列数 (>阈值):':<20}{long_count}\n")
        f.write(f"{'长序列占比:':<20}{percentage:.2f}%\n")
        f.write("\n")
        
        # ===== 详细统计信息 =====
        f.write("[全局统计信息]\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'最短序列:':<15}{min(all_lengths)} aa\n")
        f.write(f"{'最长序列:':<15}{max(all_lengths)} aa\n")
        f.write(f"{'平均长度:':<15}{sum(all_lengths)/len(all_lengths):.1f} aa\n")
        f.write(f"{'长度中位数:':<15}{np.median(all_lengths):.1f} aa\n")
        
        # ===== 长度分布统计 =====
        f.write("\n[长度分布]\n")
        f.write("-" * 60 + "\n")
        bins = [0, 100, 300, args.min_length, 500, 1000, float('inf')]
        labels = [
            f"<100", "100-300", f"300-{args.min_length}", 
            f"{args.min_length}-500", "500-1000", ">1000"
        ]
        
        bin_counts = [0] * len(bins)
        for length in all_lengths:
            for i in range(len(bins)-1):
                if bins[i] <= length < bins[i+1]:
                    bin_counts[i] += 1
                    break
            else:
                if length >= bins[-1]:
                    bin_counts[-1] += 1
        
        for label, count in zip(labels, bin_counts):
            percent = count / total_count * 100
            f.write(f"{label+' aa:':<15}{count:>5} 序列 ({percent:.1f}%)\n")
        
        # ===== 长序列ID列表 =====
        f.write("\n[长序列ID列表]\n")
        f.write("-" * 60 + "\n")
        f.write(f"共找到 {len(long_ids)} 条长度超过阈值的序列:\n\n")
        for i, seq_id in enumerate(long_ids):
            f.write(f"{i+1}. {seq_id}\n")
            if (i+1) % 50 == 0:  # 每50个ID添加分隔符
                f.write("\n" + "-" * 60 + "\n\n")
    
    return f"报告已保存至: {args.output}"

def main():
    # ===================== 命令行参数解析 =====================[6,7](@ref)
    parser = argparse.ArgumentParser(
        description="FASTA序列长度分析工具：统计超过指定长度的序列并生成专业报告",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="输入FASTA文件路径（必需）"
    )
    parser.add_argument(
        "-o", "--output", 
        default="sequence_length_report.txt",
        help="分析报告输出路径（默认：sequence_length_report.txt）"
    )
    parser.add_argument(
        "-l", "--min-length", 
        type=int, 
        default=428,
        help="最小长度阈值（氨基酸数），示例：--min-length=500（默认：428）"
    )
    parser.add_argument(
        "-p", "--plot", 
        default="sequence_length_plot.png",
        help="长度分布图表输出路径（默认：sequence_length_plot.png）"
    )
    
    # 解析参数
    args = parser.parse_args()
    
    # 验证输入文件
    if not os.path.isfile(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在!")
        sys.exit(1)
    
    # 设置字体
    font_prop = set_chinese_font('./font/SimHei.ttf') 

    # 执行序列分析
    long_count, total_count, long_ids, all_lengths = count_long_sequences(
        args.input, args.min_length
    )
    
    # 生成统计图表
    plot_length_distribution(all_lengths, args.min_length, args.plot, font_prop)
    
    # 生成并保存报告
    report_msg = generate_report(args, long_count, total_count, long_ids, all_lengths)
    
    # ===================== 终端输出摘要 =====================
    print(f"\n{' 序列分析结果 ':=^60}")
    print(f"{'输入文件:':<15}{args.input}")
    print(f"{'序列总数:':<15}{total_count}")
    print(f"{'长度阈值:':<15}{args.min_length} aa")
    print(f"{'长序列数:':<15}{long_count} ({long_count/total_count*100:.2f}%)")
    print(f"{'最短序列:':<15}{min(all_lengths)} aa")
    print(f"{'最长序列:':<15}{max(all_lengths)} aa")
    print(f"{'平均长度:':<15}{sum(all_lengths)/len(all_lengths):.1f} aa")
    print(f"{'长度中位数:':<15}{np.median(all_lengths):.1f} aa")
    print(f"\n{report_msg}")
    print(f"{'分布图表:':<15}{args.plot}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
