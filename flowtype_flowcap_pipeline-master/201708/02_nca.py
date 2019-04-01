import metric_learn
import numpy as np
import pandas as pd
from sklearn.datasets import load_iris

# visualisation imports
import metric_learn
import numpy as np
import pandas as pd
from sklearn.datasets import load_iris

# visualisation imports
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



# loading our dataset

iris_data = load_iris()
# this is our data
X = iris_data['data']
X = pd.read_csv('./result_py/matrixCountAdj_cortrim.csv',sep=',').as_matrix()
# these are our constraints
Y = iris_data['target']
Y = np.squeeze(np.asarray( pd.read_csv('./result_py/aml.csv',sep=',').as_matrix() ))



# function to plot the results
def plot(X, Y):
    x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
    y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
    plt.figure(2, figsize=(8, 6))

    # clean the figure
    plt.clf()

    plt.scatter(X[:, 0], X[:, 1], c=Y, cmap=plt.cm.Paired)
    plt.xlabel('Sepal length')
    plt.ylabel('Sepal width')

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xticks(())
    plt.yticks(())

    plt.show()


# plotting the dataset as is.

plot(X, Y)


# ---------------------------------------

# setting up LMNN
lmnn = metric_learn.LMNN(k=5, learn_rate=1e-6)
# fit the data!
lmnn.fit(X, Y)
# transform our input space
X_lmnn = lmnn.transform()
plot(X_lmnn, Y)

# ITML; http://www.cs.utexas.edu/users/pjain/pubs/metriclearning_icml.pdf
itml = metric_learn.ITML_Supervised(num_constraints=200)
X_itml = itml.fit_transform(X, Y)
plot(X_itml, Y)

# SDML; http://lms.comp.nus.edu.sg/sites/default/files/publication-attachments/icml09-guojun.pdf
sdml = metric_learn.SDML_Supervised(num_constraints=200)
X_sdml = sdml.fit_transform(X, Y, random_state = np.random.RandomState(1234))
plot(X_sdml, Y)

# LSML; http://web.cs.ucla.edu/~weiwang/paper/ICDM12.pdf
lsml = metric_learn.LSML_Supervised(num_constraints=200)
X_lsml = lsml.fit_transform(X, Y)
plot(X_lsml, Y)

# NCA; https://papers.nips.cc/paper/2566-neighbourhood-components-analysis.pdf

nca = metric_learn.NCA(max_iter=1000, learning_rate=0.01)
nca.fit(X, Y)
X_nca = nca.transform()
plot(X_nca, Y)

# LFDA; http://www.machinelearning.org/proceedings/icml2006/114_Local_Fisher_Discrim.pdf
lfda = metric_learn.LFDA(k=2, dim=2)
X_lfda = lfda.fit_transform(X, Y)
plot(X_lfda, Y)

# RCA; https://www.aaai.org/Papers/ICML/2003/ICML03-005.pdf
rca = metric_learn.RCA_Supervised(num_chunks=30, chunk_size=2)
X_rca = rca.fit_transform(X, Y)
plot(X_rca, Y)





