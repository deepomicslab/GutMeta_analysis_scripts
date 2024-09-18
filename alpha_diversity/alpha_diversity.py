import csv, json

def read_tsv(file_path):
    result = []
    with open(file_path, 'r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            result.append(row)
    return result

def save_dict_as_json(data, file_path):
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=2)  # 将字典转换为 JSON 格式并写入文件

if __name__ == "__main__":
    print('start')
    # results = read_tsv('/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/Diversity table::Alpha_diversity.tsv')
    # group_info = read_tsv('/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/Group information::group_info.tsv')
    # p_values = read_tsv('/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/P-value result::Alpha_diversity_pvalue.tsv')
    results = read_tsv('/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/Result/output.alpha_diversity.tsv')
    group_info = read_tsv('/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/Result/merged_input.group_info.tsv')
    p_values = read_tsv('/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/Result/output.pvalue.tsv')

    results_json = {'data': results[1:], 'columns': results[0]}
    group_info_json = {'data': group_info[1:], 'columns': group_info[0]}
    p_values_json = {'data': p_values[1:], 'columns': p_values[0]}

    save_dict_as_json(results_json, '/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/results.json')
    save_dict_as_json(group_info_json, '/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/group_info.json')
    save_dict_as_json(p_values_json, '/Users/zhangyuanzheng/Downloads/tmp/Alpha diversity/p_values.json')
