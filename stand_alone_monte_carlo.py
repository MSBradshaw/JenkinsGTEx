import pandas as pd
import random
import pickle
import sys


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


def run_monte_carlo_simulations(tsg, communities, repetitions):
    full_sizes = []
    half_sizes = []
    coms = []
    for i in range(repetitions):
        print(i)
        com = make_monte_carlo_communities(communities, len(communities))
        coms.append(com)
        full = count_tissues_in_communities(tsg, com, True)
        half = count_tissues_in_communities(tsg, com, False)
        full_sizes.append(full)
        half_sizes.append(half)
    return full_sizes, half_sizes, coms


tsg = pickle.load(open('tsg.pickle', 'rb'))
jenkins_coms = pickle.load(open('jenkins_communities.pickle', 'rb'))
# beckett_communities = pickle.load(open('beckett_communities.pickle', 'rb'))

index = sys.argv[1]

j_full, j_half, j_coms = run_monte_carlo_simulations(tsg, jenkins_coms, 10)
# m_full, m_half = run_monte_carlo_simulations(tsg, beckett_communities, 10)

pickle.dump(j_full, open('monte_carlo_com_results_jenkins_full_' + index + '.pickle','wb'))
pickle.dump(j_half, open('monte_carlo_com_results_jenkins_half_' + index + '.pickle','wb'))
pickle.dump(j_coms, open('monte_carlo_com_results_jenkins_communities_' + index + '.pickle','wb'))
# pickle.dump(m_full, open('monte_carlo_com_results_mg2_full_' + index + '.pickle','wb'))
# pickle.dump(m_half, open('monte_carlo_com_results_mg2_half_' + index + '.pickle','wb'))
