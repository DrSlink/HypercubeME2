#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <algorithm>


namespace std {
    template<>
    struct hash<pair<unsigned long, char>> {
        std::size_t operator()(const pair<unsigned long, char> &mut) const {
            static std::hash<unsigned long> hashFn1;
            static std::hash<char> hashFn2;
            std::size_t temp = hashFn1(mut.first);
            std::size_t result = 31 + (temp ^ (temp >> 32));
            temp = hashFn2(mut.second);
            return 31 * result + (temp ^ (temp >> 32));
        }
    };
}


unsigned long get_diag_end(const std::string &diag) {
    for (long i = (long) diag.size() - 4; i >= 0; i--) {
        if (diag[i] == ':') {
            return std::stoul(diag.substr(i + 2, diag.size() - i - 3)) + 1;
        }
    }
    return 0;
}


void process_diagonals(std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &diagonals,
                       std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &result,
                       const unsigned long seq_len) {
    auto diagonal_it = diagonals.begin();
    while (diagonal_it != diagonals.end()){
        std::size_t num_seq = diagonal_it->second.size();
        if (num_seq > 1) {
            unsigned long diag_end = get_diag_end(diagonal_it->first);
            std::unordered_set<std::string> seq_set;
            seq_set.reserve(num_seq);
            std::unordered_set<std::pair<unsigned long, char>> mutations;
            for (const auto &seq : diagonal_it->second) {
                for (const auto &mut : mutations) {
                    if ((*seq)[mut.first] != mut.second) {
                        std::string mutated = *seq;
                        mutated[mut.first] = mut.second;
                        if (seq_set.find(mutated) != seq_set.end()) {
                            result[diagonal_it->first + ':' + mut.second + std::to_string(mut.first) +
                                   (*seq)[mut.first]].push_back(seq);
                        }
                    }
                }
                seq_set.insert(*seq);
                for (unsigned long i = diag_end; i < seq_len; i++) {
                    mutations.insert({i, (*seq)[i]});
                }
            }
        }
        diagonal_it = diagonals.erase(diagonal_it);
    }
}

int comp_str_shared_ptr(const std::shared_ptr<const std::string> &first, const std::shared_ptr<const std::string> &second) {
    return *first < *second;
}

int comp_str_ptr(const std::string* first, const std::string* second) {
    return *first < *second;
}

int main(int argc, char *argv[]) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string input = argv[1];
    std::ifstream my_file(input);
    if (my_file.is_open()) {
        std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> dataset[2];
        int dimension = 1;
        std::string line;
        std::getline(my_file, line);
        while (std::getline(my_file, line)) {
            dataset[dimension][""].push_back(std::make_shared<const std::string>(line));
        }
        my_file.close();
        std::sort(dataset[dimension][""].begin(), dataset[dimension][""].end(), comp_str_shared_ptr);
        unsigned long seq_len = (*dataset[dimension][""][0]).size();
        while (true) {
            process_diagonals(dataset[dimension % 2], dataset[(dimension + 1) % 2], seq_len);

            std::size_t number_diagonals = dataset[(dimension + 1) % 2].size();
            std::cout << "Number of " << dimension << " dimensional diagonals: " << number_diagonals << std::endl;

            std::vector<const std::string*> res_data;
            res_data.reserve(number_diagonals);
            for (const auto &elem : dataset[(dimension + 1) % 2]) {
                res_data.push_back(&elem.first);
            }
            std::sort(res_data.begin(), res_data.end(), comp_str_ptr);

            std::size_t number_hypercubes = 0;
            std::ofstream res_file;
            res_file.open("hypercubes_" + std::to_string(dimension) + ".txt");
            for (const auto &diagonal : res_data) {
                number_hypercubes += dataset[(dimension + 1) % 2][*diagonal].size();
                std::string new_diagonal = (*diagonal).substr(1) + '\t';
                for (const auto &offset : dataset[(dimension + 1) % 2][*diagonal]) {
                    res_file << new_diagonal << *offset << '\n';
                }
            }
            res_file.close();

            std::cout << "Number of " << dimension << " dimensional hypercubes: " << number_hypercubes << std::endl;
            dimension += 1;
            if (number_diagonals == number_hypercubes) {
                break;
            }
        }
    }
    end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Time: " << elapsed_seconds << " ms\n";
    return 0;
}
