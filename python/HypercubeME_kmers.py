import argparse
import os
import sys
import time


def load_dataset(sequences_path: str) -> (dict, int):
    data_list = []
    with open(sequences_path) as fh:
        fh.readline()
        for line in fh:
            first = line.split()[0]
            data_list.append(first)
    data_list.sort()
    return {('',): data_list}, len(data_list[0])


def process_diagonal(diag: tuple, sequences: list, seq_len: int, result: dict):
    diag_end = int(diag[-1][1:-1]) + 1 if diag[-1] else 0
    for pos in range(diag_end, seq_len):
        same_kmers = {}
        for seq in sequences:
            kmer = seq[:pos] + seq[pos + 1:]
            if kmer in same_kmers:
                for prev_seq in same_kmers[kmer]:
                    new_diag = diag + (sys.intern(prev_seq[pos] + str(pos) +
                                                  seq[pos]),)
                    if new_diag in result:
                        result[new_diag].append(seq)
                    else:
                        result[new_diag] = [seq]
                same_kmers[kmer].append(seq)
            else:
                same_kmers[kmer] = [seq]


def process_diagonals(diagonals: dict, seq_len: int) -> dict:
    result = {}
    for diag, sequences in diagonals.items():
        process_diagonal(diag, sequences, seq_len, result)
    return result


def process_dataset(sequences_path: str, dst_path: str):
    dataset, seq_len = load_dataset(sequences_path)

    dimension = 1
    while True:
        dataset = process_diagonals(dataset, seq_len)
        number_diagonals = len(dataset)
        print('Number of %2d dimensional diagonals: %10d'
              % (dimension, number_diagonals))
        number_hypercubes = 0

        with open(os.path.join(dst_path, 'hypercubes_' + str(dimension) + '.txt'),
                  'w') as f:
            for diag, sequences in sorted(dataset.items(), key=lambda x: x[0]):
                number_hypercubes += len(sequences)
                diag_name = ':'.join(diag[1:])
                f.write(''.join(diag_name + '\t' + seq + '\n' for seq in sequences))

        print('Number of %2d dimensional hypercubes: %9d'
              % (dimension, number_hypercubes))
        dimension += 1
        if number_diagonals == number_hypercubes:
            break


if __name__ == '__main__':
    start_time = time.time()

    p = argparse.ArgumentParser()
    p.add_argument('-s', '--sequences',
                   help='the filename with the list of measured genotypes',
                   required=True)
    p.add_argument('-d', '--dst',
                   help='name of destination folder for hypercubes'
                        '"hypercubes" by default',
                   default='hypercubes')
    args = p.parse_args()

    if not os.path.exists(args.dst):
        os.makedirs(args.dst)

    process_dataset(args.sequences, args.dst)

    print(time.time() - start_time)
