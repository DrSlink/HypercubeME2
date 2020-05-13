import argparse
import os
import sys
import time


def process_next_dimension(omics_dict):
    new_dict = {}
    for key, value in omics_dict.items():
        number_sequences = len(value)
        for i in range(number_sequences - 1):
            first_len = len(value[i])
            for j in range(i + 1, number_sequences):
                if abs(first_len - len(value[j])) > 1:
                    continue
                diff_set = value[i].symmetric_difference(value[j])
                if len(diff_set) == 1:
                    mutation = diff_set.pop()
                    new_diag_elem = mutation[-1] + mutation[:-1] + 'Z'
                    if new_diag_elem > key[-1]:
                        new_diag = key + (sys.intern(new_diag_elem),)
                        new_set = value[i].copy()
                        new_set.discard(mutation)
                        if new_diag in new_dict:
                            new_dict[new_diag].append(new_set)
                        else:
                            new_dict[new_diag] = [new_set]
                elif len(diff_set) == 2:
                    first_mutation = diff_set.pop()
                    second_mutation = diff_set.pop()
                    if first_mutation[:-1] == second_mutation[:-1]:
                        if first_mutation[-1] < second_mutation[-1]:
                            new_diag_elem = first_mutation[-1] + \
                                            first_mutation[:-1] + \
                                            second_mutation[-1]
                        else:
                            new_diag_elem = second_mutation[-1] + \
                                            first_mutation[:-1] + \
                                            first_mutation[-1]
                        if new_diag_elem > key[-1]:
                            new_diag = key + (sys.intern(new_diag_elem),)
                            new_set = value[i].copy()
                            new_set.discard(first_mutation)
                            new_set.discard(second_mutation)
                            if new_diag in new_dict:
                                new_dict[new_diag].append(new_set)
                            else:
                                new_dict[new_diag] = [new_set]
    return new_dict


if __name__ == '__main__':
    start_time = time.time()

    p = argparse.ArgumentParser()
    p.add_argument('-m', '--mutations',
                   help='the filename with the list of measured mutations',
                   required=True)
    p.add_argument('-d', '--dst',
                   help='name of destination folder for hypercubes'
                        '"hypercubes_mut" by default',
                   default='hypercubes_mut')
    args = p.parse_args()

    with open(args.mutations) as f:
        omics_dict = {('',): []}
        f.readline()
        for line in f:
            first = line.split('\t')[0].replace('\n', '')
            if first == '' or first == 'wt':
                omics_dict[('',)].append(set())
            else:
                omics_dict[('',)].append(
                    set(sys.intern(mut) for mut in first.split(':')))
    if not os.path.exists(args.dst):
        os.makedirs(args.dst)
    dimension = 1
    while True:
        omics_dict = process_next_dimension(omics_dict)
        number_diagonals = len(omics_dict)
        print('Number of %2d dimensional diagonals: %10d'
              % (dimension, number_diagonals))
        number_hypercubes = 0
        with open(os.path.join(args.dst, 'hypercubes_' + str(dimension) + '.txt'),
                  'w') as f:
            for diagonal, offsets in sorted(omics_dict.items(), key=lambda x: x[0]):
                number_hypercubes += len(offsets)
                diagonal_name = ':'.join(diagonal)[1:]
                f.write(''.join(diagonal_name + '\t' +
                                ':'.join(val for val in offset) + '\n'
                                for offset in offsets))
        print('Number of %2d dimensional hypercubes: %9d'
              % (dimension, number_hypercubes))
        dimension += 1
        if number_diagonals == number_hypercubes:
            break

    print(time.time() - start_time)
