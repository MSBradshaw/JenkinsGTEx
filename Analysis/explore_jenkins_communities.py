import pandas as pd
import pickle
from os import path
import random
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import obonet

"""
Returns a pandas DataFrame of the GTEx gene expression data
"""


def load_gtex():
    gtex_pickle_path = 'GTEx/GTEx.pickle'
    gtex = None
    if path.exists(gtex_pickle_path):
        print('loading GTEx')
        gtex = pickle.load(open(gtex_pickle_path, 'rb'))
    else:
        print('Reading GTEx: This could take a while...')
        gtex = pd.read_csv('GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', skiprows=2)
        # protocol 4 is to allow for large file handling (over 4gb)
        pickle.dump(gtex, open(gtex_pickle_path, 'wb'), protocol=4)
    return gtex


def get_jenkins_communities():
    becketts_pickle_path = 'CachedFiles/jenkins_communities.pickle'
    if path.exists(becketts_pickle_path):
        print('Loading Jenkins Beckett Communities')
        return pickle.load(open(becketts_pickle_path, 'rb'))
    else:
        print('Generating Jenkins Becketts Communities')
        base_name = 'Jenkins-raw/outfile'
        best_q = 0.00
        best_file = None

        # find the permutation with the best modularity score
        for i in range(1, 100):
            for j in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']:
                file = base_name + str(i) + j + '.txt'
                if not path.exists(file):
                    continue
                with open(file) as f:
                    first_line = f.readline()
                    if float(first_line.strip()) > best_q:
                        best_q = float(first_line.strip())
                        best_file = file

        # read in the communities associated with that best modularity score.
        # read it in as a 2D list
        first = True
        coms = []
        file = open(best_file, "r")
        for line in file:
            if first:
                first = False
            else:
                coms.append(line.strip().split(','))
        file.close()

        # fix a formatting error introduced from R
        for i in range(len(coms)):
            for j in range(len(coms[i])):
                if coms[i][j][0:3] == 'HP.':
                    coms[i][j] = coms[i][j].replace('.', ':')
        pickle.dump(coms, open(becketts_pickle_path, 'wb'))
        return coms


# how many tissues are there that are specific to every gene in the community
def get_tissues_specific_to_all_genes(tissue_df, com):
    # get just the genes in the community
    genes = [x for x in com if x[:3] != 'HP:']
    if len(genes) <= 1:
        # print('Warning, only one gene in the community')
        return []
    tissues = []
    # find the tissues that gene is specific for
    tissue_specific_genes_count = 0
    for g in genes:
        temp_tissues = list(tissue_df[tissue_df['Description'] == g]['tissue'])
        if len(temp_tissues) > 0:
            tissue_specific_genes_count += 1
        tissues += temp_tissues

    if tissue_specific_genes_count == 1:
        # print('Only one tissues specific gene')
        return []

    common_tissues = []
    for t in list(set(tissues)):
        if tissues.count(t) >= tissue_specific_genes_count:
            common_tissues.append(t)
    return common_tissues


# how many tissues are there that are specific to every gene in the community
def get_tissues_specific_to_half_genes(tissue_df, com):
    # get just the genes in the community
    genes = [x for x in com if x[:3] != 'HP:']
    if len(genes) <= 1:
        # print('Warning, only one gene in the community')
        return []
    tissues = []
    # find the tissues that gene is specific for
    tissue_specific_genes_count = 0
    for g in genes:
        temp_tissues = list(tissue_df[tissue_df['Description'] == g]['tissue'])
        if len(temp_tissues) > 0:
            tissue_specific_genes_count += 1
        tissues += temp_tissues

    if int(tissue_specific_genes_count * .5) <= 1:
        # print('Only one tissues specific gene')
        return []

    common_tissues = []
    for t in list(set(tissues)):
        if tissues.count(t) >= int(tissue_specific_genes_count * .5):
            common_tissues.append(t)
    return common_tissues


"""
Given
tsg: tissue specific gene dataframe
communites: list of lists of communities
strict: True: find communities where all genes are specific for a common tissue
strict: False: find communities where half the genes are specific for a common tissue
Returns a diction of of the community ID and a list of tissues it is specific for
"""


def count_tissues_in_communities(tsg, communities, strict):
    tissue_specific_communites = {}
    community_index = 0
    for comm in communities:
        tis = None
        if strict:
            tis = get_tissues_specific_to_all_genes(tsg, comm)
        else:
            tis = get_tissues_specific_to_half_genes(tsg, comm)
        if len(tis) > 0:
            tissue_specific_communites[community_index] = tis
        community_index += 1
    return tissue_specific_communites


"""
Given:
communities: a list of lists of communities
size: how many communities ot generate
Return a list of list of randomly sampled and created communities containing only genes
"""


def make_monte_carlo_communities(communities, size):
    genes = []
    for com in communities:
        genes += com
    genes = [x for x in genes if x[:3] != 'HP:']

    com_sizes = [len(com) for com in communities]
    monte_carlo_sizes = random.choices(com_sizes, k=size)

    # make 1000 fake communities
    monte_carlo_communities = []
    for i in range(len(monte_carlo_sizes)):
        monte_carlo_communities.append(random.sample(genes, monte_carlo_sizes[i]))
    return monte_carlo_communities


"""
Given:
tsg: tissue specific gene dataframe
communites: list of lists of communities
repititions: the number of repititions to run
Returns: 2 lists containing the number of communities in each repetition where genes: 
    List 1) all shared a common tissue specificity 
    List 2) half shared a common tissue specific gene   
"""


def run_monte_carlo_simulations(tsg, communities, repetitions):
    full_sizes = []
    half_sizes = []
    for i in range(repetitions):
        print(i)
        com = make_monte_carlo_communities(communities, len(communities))
        full = count_tissues_in_communities(tsg, com, True)
        half = count_tissues_in_communities(tsg, com, False)
        full_sizes.append(full)
        half_sizes.append(half)
    return full_sizes, half_sizes


def get_mean(data):
    return sum(data) / len(data)


random.seed(0)
jenkins_coms = get_jenkins_communities()

# get the number of genes per com in jenkins
jen_com_sizes = [[x for x in com if x[:3] == 'HP:' ] for com in jenkins_coms]
filtered_jen_com_sizes = [com for com in jen_com_sizes if len(com) > 1]
print((len(jen_com_sizes)))
print((len(filtered_jen_com_sizes)))
print(filtered_jen_com_sizes)

tsg = pd.read_csv('TissueSpecificGenes/tissue_specific_genes.csv')

jen_full = count_tissues_in_communities(tsg, jenkins_coms, True)
print('Jenkins full match count: ' + str(len(jen_full)))

jen_half = count_tissues_in_communities(tsg, jenkins_coms, False)
print('Jenkins half match count: ' + str(len(jen_half)))

jf = pickle.load(open('CachedFiles/tissue_expression_jenkins_full.pickle', 'rb'))
jh = pickle.load(open('CachedFiles/tissue_expression_jenkins_half.pickle', 'rb'))

monte_jenkins_full_lengths = [len(x) for x in jf]
monte_jenkins_half_lengths = [len(x) for x in jh]

jenkins_X_full = np.random.normal(2, 0.2, len(monte_jenkins_full_lengths))
jenkins_X_half = np.random.normal(4, 0.2, len(monte_jenkins_half_lengths))

jenkins_Y_full = monte_jenkins_full_lengths
jenkins_Y_half = monte_jenkins_half_lengths

plt.scatter(jenkins_X_full, jenkins_Y_full, color='red')
plt.scatter(jenkins_X_half, jenkins_Y_half, color='green')
plt.scatter([2], [len(jen_full)], color='blue')
plt.scatter([4], [len(jen_half)], color='blue')
plt.ylabel('Number of communities sharing tissue specificity')
plt.xticks([2, 4], ['Jenkins Full','Jenkins Half'])
plt.plot((1.5,2.5),(get_mean(jenkins_Y_full),get_mean(jenkins_Y_full)),'black')
plt.plot((3.5,4.5),(get_mean(jenkins_Y_half),get_mean(jenkins_Y_half)),'black')
plt.savefig('jenkins_tissue_exprs_mc_bee_hive.png')

# make a list of communities that are reconstructed and ones that are not
url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
graph = obonet.read_obo(url)

# just the communities that are reconstructed
recon_com_indexes = list(jen_full.keys())

non_recon_com_indexes = [x for x in list(range(len(jenkins_coms))) if x not in recon_com_indexes]


def write_com_stats(indexes,filename):
    # for each community
    recon_out = open(filename,'w')
    good_count = 0
    total_count = 0
    for com_i in indexes:
        print(com_i)
        com = jenkins_coms[com_i]
        try:
            recon_out.write(str(jen_full[com_i]) + '\n')
        except KeyError:
            recon_out.write(str(com_i) + '\n')
            pass
        for hpo in [x for x in com if x[:3] == 'HP:']:
            total_count+=1
            try:
                graph.nodes[hpo]
            except KeyError:
                print('Key Error: ' + hpo)
                continue

            recon_out.write(hpo + '\t')
            recon_out.write(graph.nodes[hpo]['name'] + '\t')
            if 'def' in graph.nodes[hpo].keys():
                recon_out.write(graph.nodes[hpo]['def'])
                good_count+=1
            elif 'comment' in graph.nodes[hpo].keys():
                print('Comment: ' + hpo)
                recon_out.write(graph.nodes[hpo]['comment'])
            else:
                print('No info: ' + hpo)
            recon_out.write('\n')
        for gene in [x for x in com if x[:3] != 'HP:']:
            recon_out.write(gene + '\n')
        recon_out.write('\n')

write_com_stats(recon_com_indexes,'reconstructed_communities.txt')
write_com_stats(non_recon_com_indexes, 'not_reconstructed_communities.txt')


