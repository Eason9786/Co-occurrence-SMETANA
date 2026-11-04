import os
import csv

def read_ck1_tsv(ck1_path):
    species_dict = {}
    with open(ck1_path, 'r', encoding='utf-8') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            species = row['Species']
            seq_id = row['#Query']
            main_id = '_'.join(seq_id.split('_')[:-1]) if '_' in seq_id else seq_id
            if species and species != '-':
                if species not in species_dict:
                    species_dict[species] = set()
                species_dict[species].add(main_id)
    return species_dict

def read_fasta(fasta_path):
    fasta_dict = {}
    with open(fasta_path, 'r', encoding='utf-8') as fasta_file:
        seq_id = None
        full_id = None
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    fasta_dict[seq_id] = (full_id, ''.join(sequence))
                full_id = line[1:]
                parts = full_id.split('_')
                if len(parts) >= 2 and parts[-1].isdigit():
                    if parts[-2].isdigit():
                        main_id = '_'.join(parts[:-2])
                    else:
                        main_id = '_'.join(parts[:-1])
                else:
                    main_id = full_id
                seq_id = main_id
                sequence = []
            else:
                sequence.append(line)
        if seq_id:
            fasta_dict[seq_id] = (full_id, ''.join(sequence))
    return fasta_dict

def write_species_fasta(species_dict, fasta_dict, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for species, seq_ids in species_dict.items():
        safe_species = species.replace(' ', '_').replace('/', '_')
        output_path = os.path.join(output_dir, f"{safe_species}.faa")
        with open(output_path, 'w', encoding='utf-8') as out_fa:
            for seq_id in seq_ids:
                if seq_id in fasta_dict:
                    full_id, sequence = fasta_dict[seq_id]
                    out_fa.write(f">{full_id}\n")
                    out_fa.write(f"{sequence}\n")
                else:
                    print(f"警告: 序列ID {seq_id} 未在FASTA文件中找到。")
        print(f"已生成文件: {output_path}")

def main():
    ck1_path = './1_tsv_to_faa/UL_TS_3.tsv'
    fasta_path = './1_tsv_to_faa/gene.uniGeneset.faa'
    output_dir = 'UL_TS_3_species_faa'

    print("读取CK1.tsv文件...")
    species_dict = read_ck1_tsv(ck1_path)
    print("读取gene.uniGeneset.faa文件...")
    fasta_dict = read_fasta(fasta_path)

    total_seq_ids = sum(len(seq_ids) for seq_ids in species_dict.values())
    matched_seq_ids = sum(1 for species, seq_ids in species_dict.items() for seq_id in seq_ids if seq_id in fasta_dict)
    print(f"总共有 {total_seq_ids} 个序列ID，其中 {matched_seq_ids} 个匹配。")

    print("按物种写入FASTA文件...")
    write_species_fasta(species_dict, fasta_dict, output_dir)
    print("完成。")

if __name__ == "__main__":
    main()
