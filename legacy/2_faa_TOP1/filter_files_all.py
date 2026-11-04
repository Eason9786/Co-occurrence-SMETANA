import os
import shutil

# 配置主文件夹路径
base_folder = "./2_faa_TOP1/UL_faa"  # 主文件夹路径
output_base_folder = "./2_faa_TOP1/Processed_Output"  # 输出主文件夹路径

# 创建输出主文件夹
if not os.path.exists(output_base_folder):
    os.makedirs(output_base_folder)

# 遍历主文件夹中的所有子文件夹
for sub_folder in os.listdir(base_folder):
    sub_folder_path = os.path.join(base_folder, sub_folder)
    
    # 确保当前路径是一个子文件夹
    if os.path.isdir(sub_folder_path):
        output_folder = os.path.join(output_base_folder, sub_folder + "_TOP1")  # 当前子文件夹的输出路径

        # 创建子文件夹的输出文件夹
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # 分组字典
        file_groups = {}

        # 遍历子文件夹中的文件
        for file_name in os.listdir(sub_folder_path):
            if file_name.endswith(".faa"):  # 只处理 faa 文件
                # 提取文件名中的分组标识符
                parts = file_name.split('s__')
                if len(parts) > 1 and '_' in parts[1]:
                    group_name = parts[1].split('_')[0]
                    file_path = os.path.join(sub_folder_path, file_name)
                    file_size = os.path.getsize(file_path)  # 获取文件大小
                    if group_name not in file_groups:
                        file_groups[group_name] = []
                    file_groups[group_name].append((file_name, file_size, file_path))

        # 按组处理文件
        for group_name, files in file_groups.items():
            # 按文件大小降序排序
            files_sorted = sorted(files, key=lambda x: x[1], reverse=True)
            # 保留前 1 个文件
            top_1_files = files_sorted[:1]
            for file_name, file_size, file_path in top_1_files:
                # 将文件复制到输出文件夹
                destination_path = os.path.join(output_folder, file_name)
                shutil.copy(file_path, destination_path)
                print(f"文件 {file_name}（大小：{file_size} 字节）已复制到 {output_folder}")

print("所有子文件夹处理完成！")
