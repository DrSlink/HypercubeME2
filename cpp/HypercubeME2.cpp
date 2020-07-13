#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <algorithm>


//=====================================================================================================================
// COMMON
//=====================================================================================================================

unsigned long get_diag_end(const std::string &diag) {
    for (long i = (long) diag.size() - 4; i >= 0; i--) {
        if (diag[i] == ':') {
            return std::stoul(diag.substr(i + 2, diag.size() - i - 3)) + 1;
        }
    }
    return 0;
}

int comp_str_shared_ptr(const std::shared_ptr<const std::string> &first,
                        const std::shared_ptr<const std::string> &second) {
    return *first < *second;
}

int comp_str_ptr(const std::string *first, const std::string *second) {
    return *first < *second;
}

//=====================================================================================================================
// ORIGINAL ALGORITHM O(N^2)
//=====================================================================================================================

long long compare_seqs(const std::string &first, const std::string &second, const unsigned long seq_len) {
    long long pos = -1;
    for (long long i = 0; i < seq_len; i++) {
        if (first[i] != second[i]) {
            if (pos == -1) {
                pos = i;
            } else {
                return -1;
            }
        }
    }
    return pos;
}

void process_diagonals_orig(std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &diagonals,
                            std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &result,
                            const unsigned long seq_len) {
    auto diagonal_it = diagonals.begin();
    while (diagonal_it != diagonals.end()) {
        std::size_t num_seq = diagonal_it->second.size();
        if (num_seq > 1) {
            unsigned long diag_end = get_diag_end(diagonal_it->first);
            for (unsigned long i = 1; i < num_seq; i++) {
                for (unsigned long j = 0; j < i; j++) {
                    long long pos = compare_seqs(*diagonal_it->second[i], *diagonal_it->second[j], seq_len);
                    if (pos != -1 && pos >= diag_end) {
                        result[diagonal_it->first + ':' + (*diagonal_it->second[j])[pos] + std::to_string(pos) +
                               (*diagonal_it->second[i])[pos]].push_back(diagonal_it->second[i]);
                    }
                }
            }
        }
        diagonal_it = diagonals.erase(diagonal_it);
    }
}

//=====================================================================================================================
// K-MERS ALGORITHM O(N*L)
//=====================================================================================================================

bool kmer_equal(const std::shared_ptr<const std::string> &first, const std::shared_ptr<const std::string> &second,
                const unsigned long ex_pos) {
    std::string_view str_view_f(*first);
    std::string_view str_view_s(*second);
    return str_view_f.substr(0, ex_pos) == str_view_s.substr(0, ex_pos) &&
           str_view_f.substr(ex_pos + 1) == str_view_s.substr(ex_pos + 1);
}

std::size_t kmer_hash(const std::shared_ptr<const std::string> &first, const unsigned long ex_pos) {
    static std::hash<std::string_view> hashFn;
    std::string_view str_view(*first);

    std::size_t temp = hashFn(str_view.substr(0, ex_pos));
    temp ^= hashFn(str_view.substr(ex_pos + 1)) + 0x9e3779b9 + (temp << 6) + (temp >> 2);
    return temp;
}

void process_kmer(const std::string *diagonal, const std::vector<std::shared_ptr<const std::string>> *group,
                  std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> *result,
                  const unsigned long ex_pos) {
    auto hash = [ex_pos](const std::shared_ptr<const std::string> &n) { return kmer_hash(n, ex_pos); };
    auto equal = [ex_pos](const std::shared_ptr<const std::string> &l, const std::shared_ptr<const std::string> &r) {
        return kmer_equal(l, r, ex_pos);
    };
    std::unordered_map<std::shared_ptr<const std::string>, std::vector<std::shared_ptr<const std::string>>,
            decltype(hash), decltype(equal)> m(8, hash, equal);

    for (const auto &seq : *group) {
        if (m.find(seq) != m.end()) {
            auto group_it = m.find(seq);
            for (const auto &elem : group_it->second) {
                (*result)[*diagonal + ':' + (*elem)[ex_pos] + std::to_string(ex_pos) +
                          (*seq)[ex_pos]].push_back(seq);
            }
        }
        m[seq].push_back(seq);
    }
}

void
process_diagonals_kmers(std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &diagonals,
                        std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &result,
                        const unsigned long seq_len) {
    for (const auto &diagonal : diagonals) {
        std::size_t num_seq = diagonal.second.size();
        if (num_seq > 1) {
            unsigned long diag_end = get_diag_end(diagonal.first);
            for (unsigned long ex_pos = diag_end; ex_pos < seq_len; ex_pos++) {
                process_kmer(&diagonal.first, &diagonal.second, &result, ex_pos);
            }
        }
    }
}

//=====================================================================================================================
// HASH-TABLE ALGORITHM O(N*A*L)
//=====================================================================================================================

namespace std {
    template<>
    struct hash<pair<unsigned long, char>> {
        std::size_t operator()(const pair<unsigned long, char> &mut) const {
            static std::hash<unsigned long> hashFn1;
            static std::hash<char> hashFn2;
            std::size_t temp = hashFn1(mut.first);
            temp ^= hashFn2(mut.second) + 0x9e3779b9 + (temp << 6) + (temp >> 2);
            return temp;
        }
    };
}

void process_diagonals_ht(std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &diagonals,
                          std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &result,
                          const unsigned long seq_len) {
    auto diagonal_it = diagonals.begin();
    while (diagonal_it != diagonals.end()) {
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

//=====================================================================================================================
// RECURSIVE HASH-TABLE ALGORITHM O(N*Log(L))
//=====================================================================================================================

bool substr_equal(const std::shared_ptr<const std::string> &first, const std::shared_ptr<const std::string> &second,
                  const unsigned long start, const unsigned long len) {
    return (*first).compare(start, len, (*second), start, len) == 0;
}

std::size_t substr_hash(const std::shared_ptr<const std::string> &first,
                        const unsigned long start, const unsigned long len) {
    static std::hash<std::string_view> hashFn;
    return hashFn(std::string_view((*first).c_str() + start, len));
}

void add_combinations(const std::string &diagonal, const std::vector<std::shared_ptr<const std::string>> &group,
                      std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &result,
                      const unsigned long rest_start, const unsigned long group_size) {
    std::string res_start_str = std::to_string(rest_start);
    for (int i = 1; i < group_size; i++) {
        for (int j = 0; j < i; j++) {
            std::string new_diagonal = diagonal;
            new_diagonal += ':';
            new_diagonal += (*group[j])[rest_start];
            new_diagonal += res_start_str;
            new_diagonal += (*group[i])[rest_start];
            result[new_diagonal].push_back(group[i]);
        }
    }
}

void proc_rec(const std::string *diagonal, const std::vector<std::shared_ptr<const std::string>> *group,
              std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> *result,
              const unsigned long start, const unsigned long end, const unsigned long rest_start,
              const unsigned long rest_end) {
    unsigned long len = end - start;
    auto hash = [start, len](const std::shared_ptr<const std::string> &n) {
        return substr_hash(n, start, len);
    };
    auto equal = [start, len](const std::shared_ptr<const std::string> &l,
                              const std::shared_ptr<const std::string> &r) {
        return substr_equal(l, r, start, len);
    };
    std::unordered_map<std::shared_ptr<const std::string>, std::vector<std::shared_ptr<const std::string>>,
            decltype(hash), decltype(equal)> subgroups(8, hash, equal);
    for (const auto &seq : (*group)) {
        subgroups[seq].push_back(seq);
    }

    unsigned long half = (rest_start + rest_end) / 2;

    auto elem_it = subgroups.begin();
    while (elem_it != subgroups.end()) {
        unsigned long group_size = elem_it->second.size();
        if (group_size > 1) {
            if (rest_start == half && rest_end - half == 1) {
                add_combinations(*diagonal, elem_it->second, *result, half, group_size);
            } else if (half != rest_end) {
                proc_rec(diagonal, &elem_it->second, result, rest_start, half, half, rest_end);
            }
            if (rest_end == half && half - rest_start == 1) {
                add_combinations(*diagonal, elem_it->second, *result, rest_start, group_size);
            } else if (half != rest_start) {
                proc_rec(diagonal, &elem_it->second, result, half, rest_end, rest_start, half);
            }
        }
        elem_it = subgroups.erase(elem_it);
    }
}

void process_diagonals_rht(std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &diagonals,
                           std::unordered_map<std::string, std::vector<std::shared_ptr<const std::string>>> &result,
                           const unsigned long seq_len) {
    auto diagonal_it = diagonals.begin();
    while (diagonal_it != diagonals.end()) {
        std::size_t num_seq = diagonal_it->second.size();
        if (num_seq > 1) {
            unsigned long diag_end = get_diag_end(diagonal_it->first);
            proc_rec(&diagonal_it->first, &diagonal_it->second, &result, 0, diag_end, diag_end, seq_len);
        }
        diagonal_it = diagonals.erase(diagonal_it);
    }
}

//=====================================================================================================================
// PROCESSING
//=====================================================================================================================

int main(int argc, char *argv[]) {
    if (argc <= 1) {
        std::cout << "First argument should be algorithm: cmp, k-mer, hash or recursive\n";
        std::cout << "Second argument should be name of existing file with sequences";
        return 0;
    }
    int algorithm;
    if (std::strcmp(argv[1], "cmp") == 0) {
        algorithm = 0;
        std::cout << "ORIGINAL ALGORITHM O(N^2)" << std::endl;
    } else if (std::strcmp(argv[1], "k-mer") == 0) {
        algorithm = 1;
        std::cout << "K-MERS ALGORITHM O(N*L)" << std::endl;
    } else if (std::strcmp(argv[1], "hash") == 0) {
        algorithm = 2;
        std::cout << "HASH-TABLE ALGORITHM O(N*A*L)" << std::endl;
    } else if (std::strcmp(argv[1], "recursive") == 0) {
        algorithm = 3;
        std::cout << "RECURSIVE HASH-TABLE ALGORITHM O(N*Log(L))" << std::endl;
    } else {
        std::cout << "First argument should be algorithm: cmp, k-mer, hash or recursive";
        return 0;
    }
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    if (argc <= 2) {
        std::cout << "Second argument should be name of existing file with sequences";
        return 0;
    }
    std::string input = argv[2];
    std::cout << input << std::endl;
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
            if (algorithm == 0) {
                process_diagonals_orig(dataset[dimension % 2], dataset[(dimension + 1) % 2], seq_len);
            } else if (algorithm == 1) {
                process_diagonals_kmers(dataset[dimension % 2], dataset[(dimension + 1) % 2], seq_len);
                dataset[dimension % 2].clear();
            } else if (algorithm == 2) {
                process_diagonals_ht(dataset[dimension % 2], dataset[(dimension + 1) % 2], seq_len);
            } else {
                process_diagonals_rht(dataset[dimension % 2], dataset[(dimension + 1) % 2], seq_len);
            }
            std::size_t number_diagonals = dataset[(dimension + 1) % 2].size();
            std::cout << "Number of " << dimension << " dimensional diagonals: " << number_diagonals << std::endl;

            std::vector<const std::string *> res_data;
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
                std::sort(dataset[(dimension + 1) % 2][*diagonal].begin(),
                          dataset[(dimension + 1) % 2][*diagonal].end(),
                          comp_str_shared_ptr);
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
    } else {
        std::cout << "File " + input + "not found";
        return 0;
    }
    end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Time: " << elapsed_seconds << " ms\n";
    return 0;
}
