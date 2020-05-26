import pickle

jf = []
jh = []
mf = []
mh = []
for i in range(1, 101):
    jf += pickle.load(open('MonteCarloSims/monte_carlo_com_results_jenkins_full_' + str(i) + '.pickle', 'rb'))
    jh += pickle.load(open('MonteCarloSims/monte_carlo_com_results_jenkins_half_' + str(i) + '.pickle', 'rb'))
    mf += pickle.load(open('MonteCarloSims/monte_carlo_com_results_mg2_full_' + str(i) + '.pickle', 'rb'))
    mh += pickle.load(open('MonteCarloSims/monte_carlo_com_results_mg2_half_' + str(i) + '.pickle', 'rb'))

pickle.dump(jf, open('CachedFiles/tissue_expression_jenkins_full.pickle', 'wb'))
pickle.dump(jh, open('CachedFiles/tissue_expression_jenkins_half.pickle', 'wb'))
pickle.dump(mf, open('CachedFiles/tissue_expression_mygene2_full.pickle', 'wb'))
pickle.dump(mh, open('CachedFiles/tissue_expression_mygene2_half.pickle', 'wb'))
