#!/bin/bash

###############################################################################
# 软件名称：VCF SNP Select Tool
# 版本：1.1
# 著作权所有：(c) 2025 南京农业大学棉花研究所
# 许可证：GPLv3
# 
# 功能描述：
# 本软件用于对VCF文件进行SNP过滤和选择，包括两个步骤：
# 1. 基于窗口的过滤：在每个窗口内选择QUAL最高的SNP，并根据QUAL、AF阈值过滤
# 2. 基于LDblock的过滤：将LDblock划分为子窗口，在每个子窗口内选择QUAL最高的SNP
# 
# 使用方法：
# ./vcf_snp_filter.sh -i <input.vcf> -o <output.vcf> [选项]
# 
# 选项：
#   -i, --input         输入VCF文件 (必需)
#   -o, --output        输出VCF文件 (必需)
#   -w, --window        窗口大小 (默认: 10000)
#   -q, --qual          QUAL过滤阈值 (默认: 30)
#   -a, --af            AF过滤阈值 (默认: 0.05)
#   -l, --ldblock      LDblock文件 (可选)
#   -s, --ldwindow     LDblock内子窗口大小 (默认: 20000)
#   -e, --evaluate     输出简化SNP评估报告
#   -h, --help         显示帮助信息
###############################################################################

# 默认参数
WINDOW_SIZE=10000
QUAL_THRESHOLD=30
AF_THRESHOLD=0.05
LD_WINDOW_SIZE=20000
INPUT_VCF=""
OUTPUT_VCF=""
LD_BLOCK_FILE=""
DO_EVALUATE=0

# 显示帮助信息
show_help() {
    echo "VCF SNP Filter and Select Tool v1.1"
    echo "版权所有 (c) 2025 南京农业大学棉花研究所"
    echo ""
    echo "使用方法: $0 -i <input.vcf> -o <output.vcf> [选项]"
    echo ""
    echo "必需参数:"
    echo "  -i, --input     输入VCF文件"
    echo "  -o, --output    输出VCF文件"
    echo ""
    echo "可选参数:"
    echo "  -w, --window    窗口大小 (bp), 默认: 10000"
    echo "  -q, --qual      QUAL过滤阈值, 默认: 30"
    echo "  -a, --af        AF过滤阈值, 默认: 0.05"
    echo "  -l, --ldblock   LDblock文件 (如果提供，将执行LDblock过滤)"
    echo "  -s, --ldwindow  LDblock内子窗口大小 (bp), 默认: 20000"
    echo "  -e, --evaluate  输出简化SNP评估报告"
    echo "  -h, --help     显示此帮助信息"
    echo ""
    echo "示例:"
    echo "  $0 -i input.vcf -o output.vcf -w 10000 -q 30 -a 0.05"
    echo "  $0 -i input.vcf -o output.vcf -l ld_blocks.txt -s 20000 -e"
    echo ""
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_VCF="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_VCF="$2"
            shift 2
            ;;
        -w|--window)
            WINDOW_SIZE="$2"
            shift 2
            ;;
        -q|--qual)
            QUAL_THRESHOLD="$2"
            shift 2
            ;;
        -a|--af)
            AF_THRESHOLD="$2"
            shift 2
            ;;
        -l|--ldblock)
            LD_BLOCK_FILE="$2"
            shift 2
            ;;
        -s|--ldwindow)
            LD_WINDOW_SIZE="$2"
            shift 2
            ;;
        -e|--evaluate)
            DO_EVALUATE=1
            shift 1
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "错误: 未知参数 $1"
            show_help
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z "$INPUT_VCF" || -z "$OUTPUT_VCF" ]]; then
    echo "错误: 必须提供输入和输出文件"
    show_help
    exit 1
fi

# 检查输入文件是否存在
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "错误: 输入VCF文件不存在: $INPUT_VCF"
    exit 1
fi

if [[ -n "$LD_BLOCK_FILE" && ! -f "$LD_BLOCK_FILE" ]]; then
    echo "错误: LDblock文件不存在: $LD_BLOCK_FILE"
    exit 1
fi

# 创建临时目录
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT INT TERM

echo "=================================================="
echo "VCF SNP Filter and Select Tool v1.1"
echo "版权所有 (c) 2025 南京农业大学棉花研究所"
echo "开始处理..."
echo "=================================================="
echo "输入文件: $INPUT_VCF"
echo "输出文件: $OUTPUT_VCF"
echo "窗口大小: $WINDOW_SIZE bp"
echo "QUAL阈值: $QUAL_THRESHOLD"
echo "AF阈值: $AF_THRESHOLD"
if [[ -n "$LD_BLOCK_FILE" ]]; then
    echo "LDblock文件: $LD_BLOCK_FILE"
    echo "LDblock子窗口: $LD_WINDOW_SIZE bp"
fi
if [[ $DO_EVALUATE -eq 1 ]]; then
    echo "评估报告: 启用"
fi
echo "=================================================="

# 第一步：窗口过滤
echo "步骤1: 执行窗口过滤..."
window_filter() {
    local input=$1
    local output=$2
    local window_size=$3
    local qual_threshold=$4
    local af_threshold=$5
    
    # 提取表头
    grep '^#' "$input" > "$TMP_DIR/header.tmp"
    
    # 处理VCF主体 - 保持原脚本的核心逻辑不变，但移除DP和GQ过滤，添加AF过滤
    awk -v window_size="$window_size" \
        -v qual_threshold="$qual_threshold" \
        -v af_threshold="$af_threshold" '
    BEGIN {
        OFS = "\t"
    }

    # 只处理非注释行
    !/^#/ {
        # 检查基本过滤条件
        if ($6 <= qual_threshold) next
        if ($7 != "PASS" && $7 != ".") next
        
        # 检查INFO字段中的AF
        af_found = 0
        af_value = 0
        split($8, info_fields, ";")
        for (i in info_fields) {
            if (info_fields[i] ~ /^AF=/) {
                split(info_fields[i], af_parts, "=")
                af_value = af_parts[2]
                af_found = 1
                break
            }
        }
        
        # 如果没有找到AF信息，跳过该SNP
        if (!af_found) next
        
        # 检查AF是否大于阈值
        if (af_value <= af_threshold) next

        # 处理符合条件的SNP
        chr = $1
        pos = $2
        qual = $6
        
        # 计算窗口
        win_start = int((pos-1)/window_size) * window_size + 1
        win_end = win_start + window_size - 1
        win_key = chr "_" win_start
        
        # 处理新窗口
        if (win_key != current_win) {
            if (current_win != "" && best_qual > 0) {
                print best_line
            }
            current_win = win_key
            best_qual = 0
            best_line = ""
        }
        
        # 更新当前窗口的最佳SNP
        if (qual > best_qual) {
            best_qual = qual
            best_line = $0
        }
    }

    END {
        # 输出最后一个窗口的最佳SNP
        if (best_qual > 0) {
            print best_line
        }
    }
    ' "$input" > "$TMP_DIR/body.tmp"
    
    # 合并表头和结果
    cat "$TMP_DIR/header.tmp" "$TMP_DIR/body.tmp" > "$output"
}

# 执行窗口过滤
window_filter "$INPUT_VCF" "$TMP_DIR/window_filtered.vcf" "$WINDOW_SIZE" "$QUAL_THRESHOLD" "$AF_THRESHOLD"

window_filtered_count=$(grep -v '^#' "$TMP_DIR/window_filtered.vcf" | wc -l)
echo "窗口过滤完成，找到 $window_filtered_count 个SNP"

# 第二步：LDblock过滤（如果提供了LDblock文件）
if [[ -n "$LD_BLOCK_FILE" ]]; then
    echo "步骤2: 执行LDblock过滤..."
    
    # 预处理LDblock文件
    sort -k1,1 -k2,2n "$LD_BLOCK_FILE" > "$TMP_DIR/ld_blocks_sorted.txt"
    
    # LDblock过滤函数 - 保持原脚本的核心逻辑不变
    ld_filter() {
        local input=$1
        local output=$2
        local ld_window_size=$3
        
        awk -v ld_window_size="$ld_window_size" '
        BEGIN {
            # 读取LDblock区间
            while (getline < "'"$TMP_DIR/ld_blocks_sorted.txt"'" > 0) {
                chr = $1
                start = $2
                end = $3
                block_length = end - start
                # 存储LDblock信息
                ld_blocks[chr][++count[chr]] = start "|" end "|" block_length
            }
        }
        # 处理VCF文件头
        /^#/ {
            print
            next
        }

        {
            chrom = $1
            pos = $2
            qual = $6
            found_in_ld = 0

            # 检查SNP是否在LDblock区间内
            if (chrom in count) {
                for (i = 1; i <= count[chrom]; i++) {
                    split(ld_blocks[chrom][i], arr, "|")
                    start_ld = arr[1]
                    end_ld = arr[2]
                    if (pos >= start_ld && pos <= end_ld) {
                        found_in_ld = 1
                        ld_block_length = arr[3]
                        ld_index = i
                        break
                    }
                }
            }

            if (!found_in_ld) {
                # 不在LDblock区间内，直接输出
                print
            } else {
                # 在LDblock区间内，存储信息以备后续处理
                key = chrom "|" ld_index
                if (!(key in snp_count)) {
                    snp_count[key] = 0
                    ld_start[key] = start_ld
                    ld_end[key] = end_ld
                    ld_len[key] = ld_block_length
                }
                snp_count[key]++
                idx = snp_count[key]
                snp_qual[key][idx] = qual
                snp_line[key][idx] = $0
            }
        }
        END {
            # 处理每个LDblock内的SNP
            for (key in snp_count) {
                n = snp_count[key]
                split(key, arr, "|")
                chrom = arr[1]
                ld_idx = arr[2]
                start_ld = ld_start[key]
                end_ld = ld_end[key]
                block_length_ld = ld_len[key]

                if (block_length_ld <= ld_window_size) {
                    # LDblock长度≤子窗口大小
                    if (n > 1) {
                        # 保留QUAL最高的SNP
                        max_qual = -1
                        max_idx = -1
                        for (i = 1; i <= n; i++) {
                            if (snp_qual[key][i] > max_qual) {
                                max_qual = snp_qual[key][i]
                                max_idx = i
                            }
                        }
                        print snp_line[key][max_idx]
                    }
                    # 如果只有一个SNP，不输出（按原逻辑）
                } else {
                    # LDblock长度>子窗口大小，分成子窗口
                    num_bins = int((block_length_ld - 1) / ld_window_size) + 1
                    for (bin = 1; bin <= num_bins; bin++) {
                        bin_start = start_ld + (bin - 1) * ld_window_size
                        bin_end = (bin == num_bins) ? end_ld : bin_start + ld_window_size - 1
                        max_qual_bin = -1
                        max_line_bin = ""
                        for (i = 1; i <= n; i++) {
                            split(snp_line[key][i], fields, "\t")
                            pos_snp = fields[2]
                            if (pos_snp >= bin_start && pos_snp <= bin_end) {
                                if (snp_qual[key][i] > max_qual_bin) {
                                    max_qual_bin = snp_qual[key][i]
                                    max_line_bin = snp_line[key][i]
                                }
                            }
                        }
                        if (max_line_bin != "") {
                            print max_line_bin
                        }
                    }
                }
            }
        }
        ' "$input" > "$output"
    }
    
    # 执行LDblock过滤
    ld_filter "$TMP_DIR/window_filtered.vcf" "$OUTPUT_VCF" "$LD_WINDOW_SIZE"
    echo "LDblock过滤完成"
else
    # 如果没有LDblock文件，直接使用窗口过滤结果
    cp "$TMP_DIR/window_filtered.vcf" "$OUTPUT_VCF"
    echo "跳过LDblock过滤步骤"
fi

# 对输出文件进行排序
echo "对输出VCF文件进行排序..."

# 创建临时文件进行排序
TMP_SORTED="$TMP_DIR/sorted_output.vcf"

# 提取表头
grep '^#' "$OUTPUT_VCF" > "$TMP_SORTED"

# 对非表头行进行排序（按染色体和位置）
grep -v '^#' "$OUTPUT_VCF" | sort -k1,1 -k2,2n >> "$TMP_SORTED"

# 替换原文件
mv "$TMP_SORTED" "$OUTPUT_VCF"

# 输出最终统计信息
final_snp_count=$(grep -v '^#' "$OUTPUT_VCF" | wc -l)
echo "=================================================="
echo "处理完成!"
echo "最终输出包含 $final_snp_count 个SNP"
echo "输出文件: $OUTPUT_VCF"
echo "=================================================="

# 评估功能 - 这是新增和优化的部分
evaluate_representation() {
    local input_vcf=$1
    local output_vcf=$2
    local window_size=$3
    local ld_window_size=$4
    local ld_block_file=$5
    
    local report_file="${output_vcf%.vcf}_evaluation_report.txt"
    
    echo "=================================================="
    echo "简化基因组SNP代表性评估"
    echo "=================================================="
    
    # 创建报告文件
    {
        echo "简化基因组SNP代表性评估报告"
        echo "生成时间: $(date)"
        echo "输入文件: $input_vcf"
        echo "输出文件: $output_vcf"
        echo "窗口大小: $window_size bp"
        echo "QUAL阈值: $QUAL_THRESHOLD"
        echo "AF阈值: $AF_THRESHOLD"
        if [[ -n "$ld_block_file" ]]; then
            echo "LDblock文件: $ld_block_file"
            echo "LDblock子窗口: $ld_window_size bp"
        fi
        echo "=================================================="
        echo ""
    } > "$report_file"
    
    # 计算基本统计
    local total_snps=$(grep -v '^#' "$input_vcf" | wc -l)
    local selected_snps=$(grep -v '^#' "$output_vcf" | wc -l)
    local reduction_rate="0"
    if [[ $total_snps -gt 0 ]]; then
        reduction_rate=$(echo "scale=2; (1 - $selected_snps / $total_snps) * 100" | bc)
    fi
    
    {
        echo "基本统计:"
        echo "---------"
        echo "原始SNP数量: $total_snps"
        echo "选择后SNP数量: $selected_snps"
        echo "简化率: $reduction_rate%"
        echo ""
    } >> "$report_file"
    
    # 基因组覆盖评估
    genome_coverage_assessment "$input_vcf" "$output_vcf" "$window_size" "$report_file"
    
    # MAF分布比较
    maf_distribution_comparison "$input_vcf" "$output_vcf" "$report_file"
    
    # 染色体分布评估
    chromosome_distribution_assessment "$input_vcf" "$output_vcf" "$report_file"
    
    # 如果使用了LDblock，评估LDblock覆盖
    if [[ -n "$ld_block_file" ]]; then
        ld_block_coverage_assessment "$output_vcf" "$ld_block_file" "$ld_window_size" "$report_file"
    fi
    
    # 添加评估结论
    {
        echo "评估结论:"
        echo "========="
        echo "1. 理想的简化基因组应该保持:"
        echo "   - 基因组窗口覆盖率 > 90%"
        echo "   - MAF分布与原始数据集相似"
        echo "   - 各染色体保留比例均衡"
        echo "   - LDblock覆盖率 > 85% (如果适用)"
        echo ""
        echo "2. 建议:"
        echo "   - 如果覆盖率不足，可考虑减小窗口大小"
        echo "   - 如果MAF分布偏差较大，可调整过滤阈值"
        echo "   - 如果某些染色体保留比例异常，检查该染色体数据质量"
        echo ""
        echo "=================================================="
        echo "报告生成完成"
    } >> "$report_file"
    
    echo "评估报告已保存至: $report_file"
}

genome_coverage_assessment() {
    local input_vcf=$1
    local output_vcf=$2
    local window_size=$3
    local report_file=$4
    
    echo "1. 基因组覆盖度评估..."
    
    # 创建临时文件
    local tmp_input="$TMP_DIR/input_pos.txt"
    local tmp_output="$TMP_DIR/output_pos.txt"
    
    # 提取位置信息
    grep -v '^#' "$input_vcf" | cut -f1,2 > "$tmp_input"
    grep -v '^#' "$output_vcf" | cut -f1,2 > "$tmp_output"
    
    # 计算染色体长度和覆盖窗口，按指定顺序输出
    awk -v ws="$window_size" '
    BEGIN {
        # 定义染色体顺序
        chrom_order_index = 0
        chrom_order[chrom_order_index++] = "A01"
        chrom_order[chrom_order_index++] = "A02"
        chrom_order[chrom_order_index++] = "A03"
        chrom_order[chrom_order_index++] = "A04"
        chrom_order[chrom_order_index++] = "A05"
        chrom_order[chrom_order_index++] = "A06"
        chrom_order[chrom_order_index++] = "A07"
        chrom_order[chrom_order_index++] = "A08"
        chrom_order[chrom_order_index++] = "A09"
        chrom_order[chrom_order_index++] = "A10"
        chrom_order[chrom_order_index++] = "A11"
        chrom_order[chrom_order_index++] = "A12"
        chrom_order[chrom_order_index++] = "A13"
        chrom_order[chrom_order_index++] = "D01"
        chrom_order[chrom_order_index++] = "D02"
        chrom_order[chrom_order_index++] = "D03"
        chrom_order[chrom_order_index++] = "D04"
        chrom_order[chrom_order_index++] = "D05"
        chrom_order[chrom_order_index++] = "D06"
        chrom_order[chrom_order_index++] = "D07"
        chrom_order[chrom_order_index++] = "D08"
        chrom_order[chrom_order_index++] = "D09"
        chrom_order[chrom_order_index++] = "D10"
        chrom_order[chrom_order_index++] = "D11"
        chrom_order[chrom_order_index++] = "D12"
        chrom_order[chrom_order_index++] = "D13"
        
        # 创建染色体顺序映射
        for (i in chrom_order) {
            chrom_map[chrom_order[i]] = i
        }
        
        OFS = "\t"
    }
    {
        chr = $1
        pos = $2
        if (max_pos[chr] < pos) max_pos[chr] = pos
        # 计算窗口
        win = int((pos-1)/ws) * ws + 1
        windows[chr","win] = 1
    }
    END {
        total_windows = 0
        covered_windows = 0
        
        print "染色体\t长度(bp)\t总窗口数\t覆盖窗口数\t覆盖率"
        print "------\t--------\t--------\t----------\t------"
        
        # 首先按预定顺序输出已知染色体
        for (i = 0; i < chrom_order_index; i++) {
            chr = chrom_order[i]
            if (chr in max_pos) {
                chr_len = max_pos[chr]
                chr_windows = int((chr_len-1)/ws) + 1
                total_windows += chr_windows
                
                # 统计该染色体覆盖的窗口数
                chr_covered = 0
                for (j = 1; j <= chr_windows; j++) {
                    win_start = (j-1)*ws + 1
                    if (windows[chr","win_start] == 1) {
                        covered_windows++
                        chr_covered++
                    }
                }
                printf "%s\t%d\t%d\t%d\t%.2f%%\n", chr, chr_len, chr_windows, chr_covered, (chr_covered/chr_windows)*100
            }
        }
        
        # 然后输出其他染色体（如果有）
        for (chr in max_pos) {
            if (!(chr in chrom_map)) {
                chr_len = max_pos[chr]
                chr_windows = int((chr_len-1)/ws) + 1
                total_windows += chr_windows
                
                # 统计该染色体覆盖的窗口数
                chr_covered = 0
                for (j = 1; j <= chr_windows; j++) {
                    win_start = (j-1)*ws + 1
                    if (windows[chr","win_start] == 1) {
                        covered_windows++
                        chr_covered++
                    }
                }
                printf "%s\t%d\t%d\t%d\t%.2f%%\n", chr, chr_len, chr_windows, chr_covered, (chr_covered/chr_windows)*100
            }
        }
        
        print "------\t--------\t--------\t----------\t------"
        printf "总计\t-\t%d\t%d\t%.2f%%\n", total_windows, covered_windows, (covered_windows/total_windows)*100
    }
    ' "$tmp_output" > "$TMP_DIR/coverage_report.txt"
    
    cat "$TMP_DIR/coverage_report.txt" >> "$report_file"
    echo "" >> "$report_file"
}

maf_distribution_comparison() {
    local input_vcf=$1
    local output_vcf=$2
    local report_file=$3
    
    echo "2. MAF分布比较..."
    
    awk '
    BEGIN {
        print "MAF区间\t原始SNP数\t原始比例\t选择后SNP数\t选择后比例"
        print "------\t--------\t--------\t----------\t----------"
    }
    /^[^#]/ {
        # 尝试从INFO字段提取AF
        split($8, info, ";")
        af = 0.5  # 默认值
        for (i in info) {
            if (info[i] ~ /^AF=/) {
                split(info[i], af_arr, "=")
                af = af_arr[2]
                break
            }
        }
        # 计算MAF
        maf = (af > 0.5) ? 1 - af : af
        maf_bin = int(maf * 10) * 10  # 按10%分箱
        
        # 统计
        if (FILENAME == "'"$input_vcf"'") {
            input_maf[maf_bin]++
            input_total++
        } else {
            output_maf[maf_bin]++
            output_total++
        }
    }
    END {
        for (bin = 0; bin <= 50; bin += 10) {
            input_count = input_maf[bin] + 0
            output_count = output_maf[bin] + 0
            input_pct = (input_total > 0) ? (input_count / input_total) * 100 : 0
            output_pct = (output_total > 0) ? (output_count / output_total) * 100 : 0
            printf "%d-%d%%\t%d\t%.2f%%\t%d\t%.2f%%\n", 
                bin, bin+10, input_count, input_pct, output_count, output_pct
        }
        print "------\t--------\t--------\t----------\t----------"
        printf "总计\t%d\t100.00%%\t%d\t100.00%%\n", input_total, output_total
    }
    ' "$input_vcf" "$output_vcf" > "$TMP_DIR/maf_comparison.txt"
    
    cat "$TMP_DIR/maf_comparison.txt" >> "$report_file"
    echo "" >> "$report_file"
}

chromosome_distribution_assessment() {
    local input_vcf=$1
    local output_vcf=$2
    local report_file=$3
    
    echo "3. 染色体分布评估..."
    
    awk '
    BEGIN {
        # 定义染色体顺序
        chrom_order_index = 0
        chrom_order[chrom_order_index++] = "A01"
        chrom_order[chrom_order_index++] = "A02"
        chrom_order[chrom_order_index++] = "A03"
        chrom_order[chrom_order_index++] = "A04"
        chrom_order[chrom_order_index++] = "A05"
        chrom_order[chrom_order_index++] = "A06"
        chrom_order[chrom_order_index++] = "A07"
        chrom_order[chrom_order_index++] = "A08"
        chrom_order[chrom_order_index++] = "A09"
        chrom_order[chrom_order_index++] = "A10"
        chrom_order[chrom_order_index++] = "A11"
        chrom_order[chrom_order_index++] = "A12"
        chrom_order[chrom_order_index++] = "A13"
        chrom_order[chrom_order_index++] = "D01"
        chrom_order[chrom_order_index++] = "D02"
        chrom_order[chrom_order_index++] = "D03"
        chrom_order[chrom_order_index++] = "D04"
        chrom_order[chrom_order_index++] = "D05"
        chrom_order[chrom_order_index++] = "D06"
        chrom_order[chrom_order_index++] = "D07"
        chrom_order[chrom_order_index++] = "D08"
        chrom_order[chrom_order_index++] = "D09"
        chrom_order[chrom_order_index++] = "D10"
        chrom_order[chrom_order_index++] = "D11"
        chrom_order[chrom_order_index++] = "D12"
        chrom_order[chrom_order_index++] = "D13"
        
        # 创建染色体顺序映射
        for (i in chrom_order) {
            chrom_map[chrom_order[i]] = i
        }
        
        print "染色体\t原始SNP数\t选择后SNP数\t保留比例"
        print "------\t--------\t----------\t--------"
    }
    /^[^#]/ {
        chr = $1
        if (FILENAME == "'"$input_vcf"'") {
            input_count[chr]++
            input_total++
        } else {
            output_count[chr]++
            output_total++
        }
    }
    END {
        # 首先按预定顺序输出已知染色体
        for (i = 0; i < chrom_order_index; i++) {
            chr = chrom_order[i]
            if (chr in input_count) {
                input_cnt = input_count[chr]
                output_cnt = output_count[chr] + 0
                retention = (input_cnt > 0) ? (output_cnt / input_cnt) * 100 : 0
                printf "%s\t%d\t%d\t%.2f%%\n", chr, input_cnt, output_cnt, retention
            }
        }
        
        # 然后输出其他染色体（如果有）
        for (chr in input_count) {
            if (!(chr in chrom_map)) {
                input_cnt = input_count[chr]
                output_cnt = output_count[chr] + 0
                retention = (input_cnt > 0) ? (output_cnt / input_cnt) * 100 : 0
                printf "%s\t%d\t%d\t%.2f%%\n", chr, input_cnt, output_cnt, retention
            }
        }
        
        print "------\t--------\t----------\t--------"
        # 输出总计
        total_retention = (input_total > 0) ? (output_total / input_total) * 100 : 0
        printf "总计\t%d\t%d\t%.2f%%\n", input_total, output_total, total_retention
    }
    ' "$input_vcf" "$output_vcf" > "$TMP_DIR/chr_distribution.txt"
    
    cat "$TMP_DIR/chr_distribution.txt" >> "$report_file"
    echo "" >> "$report_file"
}

ld_block_coverage_assessment() {
    local output_vcf=$1
    local ld_block_file=$2
    local ld_window_size=$3
    local report_file=$4
    
    echo "4. LDblock覆盖评估..."
    
    awk -v window_size="$ld_window_size" '
    # 读取LDblock文件
    NR == FNR {
        chr = $1
        start = $2
        end = $3
        ld_blocks[++block_count] = chr "\t" start "\t" end
        next
    }
    
    # 处理VCF文件（非注释行）
    /^[^#]/ {
        chr = $1
        pos = $2
        for (i = 1; i <= block_count; i++) {
            split(ld_blocks[i], block, "\t")
            if (chr == block[1] && pos >= block[2] && pos <= block[3]) {
                block_snps[i]++
                # 标记该block有SNP
                block_covered[i] = 1
                break
            }
        }
    }
    END {
        total_blocks = block_count
        covered_blocks = 0
        total_snps_in_blocks = 0
        
        print "LDblock统计:"
        print "-------------"
        print "总LDblock数: " total_blocks
        
        for (i = 1; i <= total_blocks; i++) {
            if (block_covered[i]) {
                covered_blocks++
                total_snps_in_blocks += block_snps[i] + 0
            }
        }
        
        print "覆盖的LDblock数: " covered_blocks
        print "LDblock覆盖率: " sprintf("%.2f%%", (covered_blocks / total_blocks) * 100)
        if (covered_blocks > 0) {
            print "LDblock内平均SNP密度: " sprintf("%.3f", total_snps_in_blocks / covered_blocks) " SNP/block"
        } else {
            print "LDblock内平均SNP密度: 0 SNP/block"
        }
    }
    ' "$ld_block_file" "$output_vcf" > "$TMP_DIR/ld_coverage.txt"
    
    cat "$TMP_DIR/ld_coverage.txt" >> "$report_file"
    echo "" >> "$report_file"
}

# 在处理完成后添加评估调用
if [[ $DO_EVALUATE -eq 1 ]]; then
    evaluate_representation "$INPUT_VCF" "$OUTPUT_VCF" "$WINDOW_SIZE" "$LD_WINDOW_SIZE" "$LD_BLOCK_FILE"
fi
