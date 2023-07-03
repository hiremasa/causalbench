import os
from io import BytesIO
from typing import List, Set, Tuple
import pkg_resources
import gdown

data_directory = "storage/"

def download_if_not_exist(url, output_directory, filename):
    path = os.path.join(output_directory, filename)
    if not os.path.exists(path):
        gdown.download(url, path, quiet=False)
    return path

def download_corum(output_directory):
    # NOTE: Temporarily providing a mirror for the CORUM file due to the below page being impacted by downtime.
    # URL = "https://mips.helmholtz-muenchen.de/corum/download/releases/current/humanComplexes.txt.zip"
    URL = "https://github.com/causalbench/causalbench-mirror/raw/main/corum_complexes.txt.zip"
    filename = "corum_complexes.txt.zip"
    path = download_if_not_exist(URL, output_directory, filename)
    return path


def download_ligand_receptor_pairs(output_directory):
    URL = "https://raw.githubusercontent.com/LewisLabUCSD/Ligand-Receptor-Pairs/ba44c3c4b4a3e501667309dd9ce7208501aeb961/Human/Human-2020-Shao-LR-pairs.txt"
    filename = "human_lr_pair.txt"
    path = download_if_not_exist(URL, output_directory, filename)
    return path


def download_string_network(output_directory):
    URL = "https://stringdb-static.org/download/protein.links.detailed.v11.5/9606.protein.links.detailed.v11.5.txt.gz"
    filename = "protein.links.txt.gz"
    path = download_if_not_exist(URL, output_directory, filename)
    return path


def download_string_physical(output_directory):
    URL = "https://stringdb-static.org/download/protein.physical.links.detailed.v11.5/9606.protein.physical.links.detailed.v11.5.txt.gz"
    filename = "protein.physical.links.txt.gz"
    path = download_if_not_exist(URL, output_directory, filename)
    return path


def download_string_protein_info(output_directory):
    URL = "https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz"
    filename = "protein.info.txt.gz"
    path = download_if_not_exist(URL, output_directory, filename)
    return path


# dataset_corum
map_name_to_ensembl = GeneNameMapLoader("storage/").load()
path_corum = download_corum("storage/")
df = pd.read_csv(
    path_corum,
    index_col="subunits(Gene name)",
    sep="\t",
    compression="zip",
)
row_names = df.index.values.tolist()
gene_from, gene_to = [], []
for row_name in row_names:
    complex_members = row_name.split(";")
    complex_members = [
        map_name_to_ensembl.get(gene)
        for gene in complex_members
        if gene in map_name_to_ensembl
    ]
    for i in range(len(complex_members)):
        for other in complex_members[:i] + complex_members[i + 1 :]:
            gene_to.append(other)
            gene_from.append(complex_members[i])
dataset_corum = set(zip(gene_from, gene_to))


# dataset_lr
path_lr_pairs = download_ligand_receptor_pairs(
    "storage/"
)
df = pd.read_csv(path_lr_pairs, index_col="ligand_gene_symbol", sep="\t")
dataset_lr = set()
for row in df[["ligand_ensembl_gene_id", "receptor_ensembl_gene_id"]].values:
    if not pd.isna(row[0]) and not pd.isna(row[1]):
        dataset_lr.add((row[0], row[1]))
        dataset_lr.add((row[1], row[0]))


# string 
path_string_network_pairs = download_string_network(
    data_directory
)
path_string_physical_pairs = download_string_physical(
    data_directory
)
filename_info = download_string_protein_info(
    data_directory
)
info = pd.read_csv(filename_info, sep="\t", compression="gzip")
protein_id_to_gene_name = dict()
for gene_name, string_id in info[
    ["preferred_name", "#string_protein_id"]
].to_numpy():
    protein_id_to_gene_name[string_id] = gene_name
df1 = pd.read_csv(path_string_network_pairs, sep=" ", compression="gzip")
df2 = pd.read_csv(path_string_physical_pairs, sep=" ", compression="gzip")
map_name_to_ensembl = GeneNameMapLoader(data_directory).load()
def create_dataset(df):
    dataset_string = set()
    for row in df[["protein1", "protein2"]].values:
        prot1, prot2 = row
        prot1, prot2 = (
            protein_id_to_gene_name[prot1],
            protein_id_to_gene_name[prot2],
        )
        if prot1 in map_name_to_ensembl and prot2 in map_name_to_ensembl:
            prot1, prot2 = (
                map_name_to_ensembl[prot1],
                map_name_to_ensembl[prot2],
            )
            dataset_string.add((prot1, prot2))
            dataset_string.add((prot2, prot1))
    return dataset_string
dataset_string_network, dataset_string_physical = create_dataset(df1), create_dataset(df2)


# chipseq
map_name_to_ensembl = GeneNameMapLoader(data_directory).load()
dataset_name = "weissmann_k562"
if dataset_name == "weissmann_k562":
    byte_data = pkg_resources.resource_string(__name__, "data/K562_ChipSeq.csv")
else:
    byte_data = pkg_resources.resource_string(
        __name__, "data/Hep_G2_ChipSeq.csv"
    )
df = pd.read_csv(BytesIO(byte_data))
dataset_chip_seq = set()
for row in df[["source", "target"]].values:
    if not pd.isna(row[0]) and not pd.isna(row[1]):
        if row[0] in map_name_to_ensembl and row[1] in map_name_to_ensembl:
            gene1, gene2 = (
                map_name_to_ensembl[row[0]],
                map_name_to_ensembl[row[1]],
            )
            dataset_chip_seq.add((gene1, gene2))

class UnionFind:
    def __init__(self, elements):
        self.parent = {element: element for element in elements}
        self.rank = {element: 0 for element in elements}

    def find(self, element):
        if self.parent[element] != element:
            self.parent[element] = self.find(self.parent[element])
        return self.parent[element]

    def union(self, element1, element2):
        root1 = self.find(element1)
        root2 = self.find(element2)

        if root1 != root2:
            if self.rank[root1] > self.rank[root2]:
                self.parent[root2] = root1
            else:
                self.parent[root1] = root2
                if self.rank[root1] == self.rank[root2]:
                    self.rank[root2] += 1

# Initialize UnionFind with all genes
data_set = dataset_corum 

def check_data_statistics(data_set):
  all_genes = set()
  for gene_pair in data_set:
      all_genes.add(gene_pair[0])
      all_genes.add(gene_pair[1])
  uf = UnionFind(all_genes)
  
  # Merge genes that are in the same complex
  for gene1, gene2 in data_set:
      uf.union(gene1, gene2)
  
  # Find the number of unique complexes
  unique_complexes = set(uf.find(gene) for gene in all_genes)
  num_complexes = len(unique_complexes)
  print("# complexes: ", num_complexes)
  
  undirected_edges = set()
  for gene_pair in data_set:
      # sort the gene pair so that the edge (A, B) becomes the same as (B, A)
      undirected_edge = tuple(sorted(gene_pair))
      undirected_edges.add(undirected_edge)
  num_undirected_edges = len(undirected_edges)
  print("# undirected edges: ", num_undirected_edges)
  
  for gene_pair in data_set:
      reversed_pair = (gene_pair[1], gene_pair[0])
      if reversed_pair not in data_set:
          print(f"Unidirectional edge detected: {gene_pair}")
  
  for gene_pair in data_set:
      if gene_pair[0] == gene_pair[1]:
          print(f"Loop edge detected: {gene_pair}")

  all_genes = set()
  for gene_pair in data_set:
      # add each gene in the pair to the set of all genes
      all_genes.add(gene_pair[0])
      all_genes.add(gene_pair[1])
  num_genes = len(all_genes)
  print("# genes: ", num_genes)

