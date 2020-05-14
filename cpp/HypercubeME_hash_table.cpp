#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <algorithm>


int get_last_diag_start(const std::string &key) {
    for (int i = (int) key.size() - 1; i >= 0; i--) {
        if (key[i] == ':') {
            return i + 1;
        }
    }
    return 0;
}


void process_genotype(const std::string &key, const int last_diag_start, const std::string &genotype,
                      std::unordered_set<std::string> &genes, std::unordered_set<std::string> &mutations,
                      std::unordered_map<std::string, std::vector<std::string>> &result) {
    for (const std::string &mutation : mutations) {
        int mut_idx = mutation[0];
        char mutation_name = mutation.back();
        if (genotype[mut_idx] != mutation_name) {
            std::string new_diagonal_point{mutation_name, (char) mut_idx, genotype[mut_idx]};
            if (key.compare(last_diag_start, 3, new_diagonal_point) < 0) {
                std::string mutated = genotype;
                mutated[mut_idx] = mutation_name;
                if (genes.find(mutated) != genes.end()) {
                    std::string diagonal = key + ':' + new_diagonal_point;
                    mutated[mut_idx] = '*';
                    if (result.find(diagonal) != result.end()) {
                        result[diagonal].push_back(mutated);
                    } else {
                        result[diagonal] = {mutated};
                    }
                }
            }
        }
    }
    genes.insert(genotype);
    for (int i = 0; i < genotype.size(); i++) {
        if (genotype[i] != '*') {
            mutations.insert({(char) i, genotype[i]});
        }
    }
}


void process_diagonals(const std::unordered_map<std::string, std::vector<std::string>> &dataset,
                       std::unordered_map<std::string, std::vector<std::string>> &result) {
    for (const auto &diagonal : dataset) {
        std::unordered_set<std::string> genes;
        std::unordered_set<std::string> mutations;
        int last_diag_start = get_last_diag_start(diagonal.first);
        for (const std::string &genotype : diagonal.second) {
            process_genotype(diagonal.first, last_diag_start, genotype, genes, mutations, result);
        }
    }
}


int main(int argc, char *argv[]) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string input = argv[1];
    std::ifstream my_file(input);
    std::string line;
    std::vector<std::string> genomes;
    if (my_file.is_open()) {
        std::getline(my_file, line);
        while (std::getline(my_file, line)) {
            genomes.push_back(line);
        }
        std::unordered_map<std::string, std::vector<std::string>> dataset[2];
        int dimension = 1;
        dataset[dimension][""] = genomes;
        while (true) {
            process_diagonals(dataset[dimension % 2], dataset[(dimension + 1) % 2]);
            std::cout << dataset[(dimension + 1) % 2].size() << std::endl;
            int number_diagonals = dataset[(dimension + 1) % 2].size();
            int number_hypercubes = 0;

            std::vector<std::pair<std::string, std::vector<std::string>>> res_data(dataset[(dimension + 1) % 2].begin(),
                                                                                   dataset[(dimension + 1) % 2].end());
            std::sort(res_data.begin(), res_data.end());

            std::string buff_string;
            for (const auto &diagonal : res_data) {
                number_hypercubes += diagonal.second.size();
                std::string new_diagonal;
                for (int i = 4; i <= diagonal.first.size(); i += 4) {
                    new_diagonal.append({diagonal.first[i - 3]});
                    new_diagonal.append(std::to_string(diagonal.first[i - 2]));
                    new_diagonal.append(diagonal.first, i - 1, 2);
                }
                new_diagonal.append({'\t'});
                for (const std::string &offset : diagonal.second) {
                    buff_string.append(new_diagonal + offset);
                    buff_string.append({'\n'});
                }
            }
            std::ofstream res_file;
            res_file.open("hypercubes_" + std::to_string(dimension) + ".txt");
            res_file << buff_string;
            res_file.close();

            std::cout << number_hypercubes << std::endl;
            dataset[dimension % 2].clear();
            dimension += 1;
            if (number_diagonals == number_hypercubes) {
                break;
            }
        }
        my_file.close();
    }
    end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Time: " << elapsed_seconds << "ms\n";
    return 0;
}