import argparse
import os
import sys
import time


def load_dataset(path: str) -> dict:
    with open(path) as fh:
        omics_list = []
        fh.readline()
        for line in fh:
            first = line.split()[0]
            omics_list.append(first)
    omics_list.sort()  # Just to make sure it is already sorted. It is fast.
    return {('',): omics_list}


def apply_mutation(genotype: str, mutation: tuple) -> str:
    return genotype[:mutation[0]] + mutation[1] + genotype[mutation[0] + 1:]


def process_genotype(key: tuple, genotype: str, genes: set, mutations: set,
                     result: dict):
    for mutation in mutations:
        if genotype[mutation[0]] != mutation[1]:
            new_diagonal_point = mutation[1] + \
                                 str(mutation[0]) + \
                                 genotype[mutation[0]]  # only if genes were sorted
            if new_diagonal_point > key[-1]:
                mutated = apply_mutation(genotype, mutation)
                if mutated in genes:
                    new_diagonal = key + (sys.intern(new_diagonal_point),)
                    offset = sys.intern(genotype[:mutation[0]] + '*' +
                                        genotype[mutation[0] + 1:])
                    if new_diagonal in result:
                        result[new_diagonal].append(offset)
                    else:
                        result[new_diagonal] = [offset]
    genes.add(genotype)
    for i in range(len(genotype)):
        if genotype[i] != '*':
            mutations.add((i, genotype[i]))


def process_diagonals(diagonals: dict) -> dict:
    result = dict()
    for key, value in diagonals.items():
        genes = set()
        mutations = set()
        for genotype in value:
            process_genotype(key, genotype, genes, mutations, result)
    return result


if __name__ == '__main__':
    start_time = time.time()

    p = argparse.ArgumentParser()
    p.add_argument('-g', '--genotypes',
                   help='the filename with the list of measured genotypes',
                   required=True)
    p.add_argument('-d', '--dst',
                   help='name of destination folder for hypercubes'
                        '"hypercubes" by default',
                   default='hypercubes')
    args = p.parse_args()

    dataset = load_dataset(args.genotypes)
    if not os.path.exists(args.dst):
        os.makedirs(args.dst)

    dimension = 1
    while True:
        dataset = process_diagonals(dataset)
        number_diagonals = len(dataset)
        print('Number of %2d dimensional diagonals: %10d'
              % (dimension, number_diagonals))
        number_hypercubes = 0

        with open(os.path.join(args.dst, 'hypercubes_' + str(dimension) + '.txt'),
                  'w') as f:
            for diagonal, offsets in sorted(dataset.items(), key=lambda x: x[0]):
                number_hypercubes += len(offsets)
                diagonal_name = ':'.join(diagonal)[1:]
                f.write(''.join(diagonal_name + '\t' + offset + '\n'
                                for offset in offsets))

        print('Number of %2d dimensional hypercubes: %9d'
              % (dimension, number_hypercubes))
        dimension += 1
        if number_diagonals == number_hypercubes:
            break
    print(time.time() - start_time)
