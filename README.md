# HypercubeME2
Reimplementation of the algorithm for generating all combinatorially complete datasets of sequences, proposed in the article: ["HypercubeME: two hundred million combinatorially complete datasets from a single experiment"](https://doi.org/10.1093/bioinformatics/btz841)  
A new algorithm is proposed based on hash tables, which is asymptotically faster.

| Script | Algorithm | S7 time | S5 time |
| ------------- | :-------------: | :-------------: | :-------------: |
| [`HypercubeME.py (orig)`](https://github.com/ivankovlab/HypercubeME/blob/master/HypercubeME.py) | proposed in paper |   1h58min54s (≈2h)* | not mesured (≈10d)* |
| `HypercubeME_re.py` | proposed in paper | 7min01s (17x)** | not mesured | 
| `HypercubeME_hash_table.py` | hash table | 0min26s (274x)** | 1min16s (11368x)** |
| `HypercubeME_hash_table.cpp` |  hash table | 0min19s (375x)** | 1min00s (14400x)** |

(*) results from paper  
(**) speed increase, compared to original script

All data was taken from [this](https://github.com/Lcarey/HIS3InterspeciesEpistasis/tree/master/Data) repository and cleared as suggested by the paper.

`HypercubeME_re.py` works on the same input format as the original script.  
`python3 HypercubeME_re.py -m mutations_file.txt -d destination_folder`

Where `mutations_file.txt` should look like:  
mut_list  
0C:7A:20F:21V:26S:27S  
0C:7G:20F:26S:27L  

The script creates `destination_folder` if it doesn't exist. You may not specify this parameter, then the script will create the folder with name `hypercubes_mut` 

`HypercubeME_hash_table.py` and `HypercubeME_hash_table.cpp` works on full sequenses but not mutations.  
`python3 HypercubeME_hash_table.py -s sequences_file.txt -d destination_folder`

Where `genotypes_file.txt` should look like:  
aa_seq  
CHALAKHAGWSLIVECIGDLFVDDHHSSED  
CHALAKHGGWSLIVECIGDLFIDDHHSLED


##### Time for generation of all 1-dimentional hypercubes in milliseconds

| Script | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 | S9 | S10 | S11 | S12 |
| ------ | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| `HypercubeME.py (orig)`| - | - | - | - | - | - | 1116662 | - | - | - | 5600298 | 17012930 |
| `HypercubeME_re.py` |935761 | 1606538 | 1153554 | 1012247 | 1121166 | 948298 | 63736 | 803821 | 1359025 | 993878 | 312357 | 812896 |
| `HypercubeME_hash_table.py` | 11809 | 13447 | 8424 | 9510 | 11920 | 8500 | 1917 | 7189 | 6200 | 8803 | 2515 | 5082 |
| `HypercubeME_hash_table.cpp` | 2165 | 2677 | 2016 | 1662 | 2425 | 1882 | 376 | 1521 | 1348 | 1370 | 521 | 1040 |

| original vs reimplemented | reimplemented vs python hash table | python hash table vs c++ hash table |
| :-: | :-: | :-: |
| ![orig_re](/plots/orig_vs_re.png) | ![all](/plots/re_vs_pyhash.png) | ![hash_py_cpp](/plots/pyhash_vs_cpphash.png) |
