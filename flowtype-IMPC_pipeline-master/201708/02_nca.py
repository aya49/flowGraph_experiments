
import metric_learn
import numpy as np
import pandas as pd
from sklearn.datasets import load_iris
import os
import sys
import random as rd

# visualisation imports
import metric_learn
import numpy as np
import pandas as pd
from sklearn.datasets import load_iris

# visualisation imports
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import TSNE
import bokeh.palettes as bp



min_samples = 2 #there must be >min_samples per class
plot_dml_dir = "./result/P1/Sanger_SPLEEN/plots/dml"
#plot_dml_dir = "/home/ayue/projects/IMPC/result/P1/Sanger_SPLEEN/plots/dml"
if not os.path.exists(plot_dml_dir):
    os.makedirs(plot_dml_dir)
#matrix_dir = '/home/ayue/projects/IMPC/result/P1/Sanger_SPLEEN/matrix'
matrix_dir = './result/P1/Sanger_SPLEEN/matrix'
matrix_type = ['PvalTRIM_CountAdj', 'LogFoldTRIM_CountAdj']



# function to plot the results
def plot(X, Y, fname, Y_original=None, title="Metric_Learning", classperplot=12, ext="png"):
    TSNEfig = ""
    cmap = plt.cm.get_cmap("hsv", classperplot + 1)
    if (X.shape[1] > 2):
        tsne = TSNE(n_components=2, random_state=0)
        X = tsne.fit_transform(X)
        TSNEfig = " TSNE"
    if (Y_original is None):
        Y_original = Y
    x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
    y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5

    # plot.

    plt.figure(2, figsize=(8, 6))
    plt.clf()  # clean the figure

    plt.scatter(X[:, 0], X[:, 1], c=Y, cmap=plt.cm.Paired, label=Y_original)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xticks(())
    plt.yticks(())

    plt.title(title + TSNEfig)
    #plt.show()
    plt.savefig(fname + " " + "." + ext)


    # plot more if there are too many classes to fit in one plot
    Y_classes, Y_copy = np.unique(Y_original, return_inverse=True)

    if len(Y_classes)>classperplot:
        plot_no = 1
        colour = bp.Category20_12
        while len(Y_classes)>0:
            Y_class = Y_classes[:min(len(Y_classes),classperplot)]
            Y_classes = Y_classes[min(len(Y_classes),classperplot):]

            plt.figure(2, figsize=(8, 6))
            plt.clf()  # clean the figure
            ax = plt.subplot(111)

            # rows = [i for i, e in enumerate(Y_original) if e in set(Y_class)]
            # Y_class_classes, Y_class_factors = np.unique(Y_original[rows], return_inverse=True)
            # plt.scatter(X[rows, 0], X[rows, 1], c=Y_class_factors, label=Y_class)

            col_ind = 0
            for i in range(len(Y_class)):
                Y_classi = Y_class[i]
                rows = [i for i, e in enumerate(Y_original) if e in Y_classi]
                ax.scatter(X[rows, 0], X[rows, 1], c=[colour[col_ind]]*len(rows), label=Y_classi)
                col_ind += 1
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])

            plt.title(title + TSNEfig + " " + str(plot_no))
            ax.legend(loc='center left', fontsize=10, bbox_to_anchor=(1, 0.5), fancybox=True).get_frame().set_alpha(0.5)

            plt.xlim(x_min, x_max)
            plt.ylim(y_min, y_max)
            plt.xticks(())
            plt.yticks(())

            #plt.show()
            plt.savefig(fname + " " + str(plot_no) + "." + ext)

            plot_no += 1



for mcp in matrix_type:
    # loading our dataset

    #iris_data = load_iris()
    # this is our data
    #X = iris_data['data']
    X = pd.read_csv(matrix_dir + mcp+'.csv',sep=',').as_matrix()


    # these are our constraints
    #Y = iris_data['target']
    Y_original = X[:,0]
    Y_classes, Y = np.unique(Y_original, return_inverse=True)
    Y_good = np.array([y if len(np.flatnonzero(Y==y))>min_samples else -1 for y in Y])
    Y_goodind = np.flatnonzero(Y_good!=-1)
    Y_original = Y_original[Y_goodind]
    Y_classes, Y = np.unique(Y_original, return_inverse=True)
    # delete classes with only 1 sample

    X = X[Y_goodind,1:]
    #X = X.reshape(X.shape[1:]) #squeeze out extra dimensions


    # plotting the dataset as is.
    plot(X, Y, Y_original=Y_original, title=mcp+"Original", fname=plot_dml_dir+"/"+mcp+"_Orig")




    # ---------------------------------------

    # setting up LMNN
    try:
        lmnn = metric_learn.LMNN(k=2, learn_rate=1e-6)
        # fit the data!
        lmnn.fit(X, Y)
        # transform our input space
        X_lmnn = lmnn.transform()
        plot(X_lmnn, Y, Y_original=Y_original, title=mcp+" LMNN", fname=plot_dml_dir+"/"+mcp+"_LMNN")
    except:
        e = sys.exc_info()[0]
        print e

    # ITML; http://www.cs.utexas.edu/users/pjain/pubs/metriclearning_icml.pdf
    try:
        itml = metric_learn.ITML_Supervised(num_constraints=200)
        itml.fit(X, Y)
        X_itml = itml.transform()
        plot(X_itml, Y, Y_original=Y_original, title=mcp+" ITML", fname=plot_dml_dir+"/"+mcp+"_ITML")
    except:
        e = sys.exc_info()[0]
        print e


    # SDML; http://lms.comp.nus.edu.sg/sites/default/files/publication-attachments/icml09-guojun.pdf
    try:
        sdml = metric_learn.SDML_Supervised(num_constraints=200)
        sdml.fit(X, Y)#, num_constraints=200)
        X_sdml = sdml.transform()
        plot(X_sdml, Y, Y_original=Y_original, title=mcp+" SDML", fname=plot_dml_dir+"/"+mcp+"_SDML")
    except:
        e = sys.exc_info()[0]
        print e

    # LSML; http://web.cs.ucla.edu/~weiwang/paper/ICDM12.pdf
    try:
        lsml = metric_learn.LSML_Supervised(num_constraints=200)
        lsml.fit(X, Y)
        X_lsml = lsml.transform()
        plot(X_lsml, Y, Y_original=Y_original, title=mcp+" LSML", fname=plot_dml_dir+"/"+mcp+"_LSML")
    except:
        e = sys.exc_info()[0]
        print e

    # NCA; https://papers.nips.cc/paper/2566-neighbourhood-components-analysis.pdf
    try:
        nca = metric_learn.NCA(max_iter=1000, learning_rate=0.01)
        nca.fit(X, Y)
        X_nca = nca.transform()
        plot(X_nca, Y, Y_original=Y_original, title=mcp+" NCA", fname=plot_dml_dir+"/"+mcp+"_NCA")
    except:
        e = sys.exc_info()[0]
        print e

    # LFDA; http://www.machinelearning.org/proceedings/icml2006/114_Local_Fisher_Discrim.pdf
    try:
        lfda = metric_learn.LFDA(k=2, dim=2)
        lfda.fit(X, Y)
        X_lfda = lfda.transform()
        plot(X_lfda, Y, Y_original=Y_original, title=mcp+" LFDA", fname=plot_dml_dir+"/"+mcp+"_LFDA")
    except:
        e = sys.exc_info()[0]
        print e

    # RCA; https://www.aaai.org/Papers/ICML/2003/ICML03-005.pdf
    try:
        rca = metric_learn.RCA_Supervised(num_chunks=30, chunk_size=2)
        rca.fit(X, Y)
        X_rca = rca.transform()
        plot(X_rca, Y, Y_original=Y_original, title=mcp+" RCA", fname=plot_dml_dir+"/"+mcp+"_RCA")
    except:
        e = sys.exc_info()[0]
        print e










