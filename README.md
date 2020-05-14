# HypercubeME2
Reimplementation of the algorithm for generating all combinatorially complete datasets of sequences, proposed in the article: ["HypercubeME: two hundred million combinatorially complete datasets from a single experiment"](https://doi.org/10.1093/bioinformatics/btz841)  
Here proposed a new algorithm based on hash tables, which is asymptotically faster.

| Script | Algorithm | S7 time | S5 time |
| ------------- | :-------------: | :-------------: | :-------------: |
| [`HypercubeME.py (orig)`](https://github.com/ivankovlab/HypercubeME/blob/master/HypercubeME.py) | proposed in paper |   1h58min54s (≈2h)* | not mesured (≈10d)* |
| `HypercubeME_re.py` | proposed in paper | 7min01s (17x)** | not mesured | 
| `HypercubeME_hash_table.py` | hash table | 1min39s (72x)** | 6min49s (2112x)** |
| `HypercubeME_hash_table.cpp` |  hash table | 0min19s (375x)** | 1min00s (14400x)** |

(*) results from paper  
(**) speed increase, compared to original script

All data was taken from [this](https://github.com/Lcarey/HIS3InterspeciesEpistasis/tree/master/Data) repository and cleared as it proposed in the paper.

`HypercubeME_re.py` works on the same input format as the original script.  
`python3 HypercubeME_re.py -m mutations_file.txt -d destination_folder`

Where `mutations_file.txt` should look like:  
mut_list  
0C:7A:20F:21V:26S:27S  
0C:7G:20F:26S:27L  

The script would create `destination_folder` if it not existed before. You may not specify this parameter, then the script will create the folder with name hypercubes_mut 

`HypercubeME_hash_table.py` and `HypercubeME_hash_table.cpp` works on full sequenses but not mutations.  
`python3 HypercubeME_hash_table.py -g genotypes_file.txt -d destination_folder`

Where `genotypes_file.txt` should look like:  
aa_seq  
CHALAKHAGWSLIVECIGDLFVDDHHSSED  
CHALAKHGGWSLIVECIGDLFIDDHHSLED
