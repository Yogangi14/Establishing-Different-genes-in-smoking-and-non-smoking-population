#!/usr/bin/env python
# coding: utf-8

# In[662]:


import os
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import f
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.stats.multitest as smm


# In[664]:


#loading the data
data_path = "../data/"
file_name = "C:/Users/yogan/Downloads/Raw Data_GeneSpring.txt"
gene_profile = pd.read_csv(( file_name), delimiter = '\t')
gene_profile = pd.DataFrame(gene_profile)


# In[666]:


print (gene_profile)


# In[668]:


gene_profile.iloc[:, 1:49] = gene_profile.iloc[:, 1:49] ** 2

print(gene_profile)


# In[670]:


N = np.array([[1, 0, 1, 0], [1, 0, 0, 1], [0, 1, 1, 0], [0,1,0,1]])
N = np.repeat(N, repeats=12, axis=0)


# In[672]:


print(len(N))


# In[674]:


D = np.array([[1, 0, 0, 0], [0, 1, 0,0 ], [0, 0, 1, 0], [0,0,0,1]])

D = np.repeat(D, repeats=12, axis=0)


# In[676]:


print(len(D))


# In[678]:


D = np.matrix(D) # returns the list of lists into a matrix form
N = np.matrix(N)
rank_N = np.linalg.matrix_rank(N)
rank_D = np.linalg.matrix_rank(D) 


# In[680]:


# degrees of freedom
degree_of_freedoms = (48-rank_D) / (rank_D-rank_N)


# In[682]:


#computing F-statistic
def compute_F_statistics(N, D,gene_profile, degree_of_freedoms):
    I = np.identity(48)
    F_statistics=[]

    numerator = I - (np.matmul(np.matmul(N,np.linalg.pinv(np.matmul(N.T,N))),N.T))
    denominator = I - (np.matmul(np.matmul(D,np.linalg.pinv(np.matmul(D.T,D))),D.T))

    #estimating m, f, s, ns and m_S, m_ns, f_s, f_ns
    for index, row in gene_profile.iterrows():
        #loop through each datset
        X = np.array(row.iloc[1:49].to_numpy().tolist())
        ## calclate f statistic for each row
        x1 = np.matmul(X.T, numerator)
        x2 = np.matmul(X.T, denominator)
        f = ((np.matmul(x1, X) / np.matmul(x2, X)) - 1) * degree_of_freedoms
        F_statistics.append(f)

    return F_statistics


# In[684]:


F_statistics = compute_F_statistics(N, D, gene_profile, degree_of_freedoms)
#print(F_statistics)


# In[686]:


F_statistics = np.array(F_statistics)


# In[702]:


#print(F_statistics)


# In[704]:


len(gene_profile)


# In[706]:


# calculation of p_values for each row
dfn = rank_D-rank_N
dfd = 48 - rank_D
cdf_value = stats.f.cdf(F_statistics, dfn, dfd)
p_values = 1 - cdf_value
print(p_values)


# In[708]:


p_values = p_values.flatten()


# In[710]:


## histogram plot
plt.hist(p_values, bins=1000)
plt.savefig("histogram of p_values.png", format="png")
plt.show()
print("Histogram of plotted of p_values")


# In[712]:


threshold = 4.08  # from the f table with df 1 and 44
significant_genes = [i for i, f in enumerate(F_statistics) if f > threshold]


# In[714]:


#genes of concern extra
significant_genes = [i for i, p in enumerate(p_values) if p < 0.05]


# In[716]:


print (significant_genes)


# In[718]:


print (len(significant_genes))


# In[720]:


plt.scatter(range(len(F_statistics)), F_statistics)
plt.xlabel('Gene Index')
plt.ylabel('F-statistic')
plt.title('Scatter Plot of F-statistics')
plt.savefig("F_stats freqs.png", format="png")
plt.show()






