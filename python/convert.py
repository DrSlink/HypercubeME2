import csv

SEGMENT = 'S5'
aa_seq_idx = 0
mut_list_idx = 18
nogap_idx = 14
middle_idx = 17

if __name__ == '__main__':
    seq = []
    mut = []
    with open(SEGMENT + '_scaled_info_v2.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_file:
            row = row.split('\t')
            if row[nogap_idx] == '1' and row[middle_idx] == '1':
                seq.append(row[aa_seq_idx])
                mut.append(row[mut_list_idx])
    print(len(seq))
    with open(SEGMENT + '_cleared_seq.txt', 'w') as f:
        f.write('aa_seq\n')
        f.write('\n'.join(seq))
    with open(SEGMENT + '_cleared_mut.txt', 'w') as f:
        f.write('mut_list\n')
        f.write('\n'.join(mut))
