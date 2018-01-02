#!/usr/bin/python
import sys
import numpy
from sklearn import svm
from sklearn.linear_model import Lasso
import os

train_size = int(sys.argv[1])

#the following function fits a linear regression on stability scores data, x are the features and y are the binary predictions of binding
def scikit_programs_regression(x_train, x_test, y_train, y_test, reg_weight):
    linear_model = Lasso(alpha=reg_weight).fit(x_train, y_train)
    corr_coef = numpy.corrcoef(y_test, linear_model.predict(x_test))
    intercept = []
    intercept.append(linear_model.intercept_)
    corr_coef_var = []
    corr_coef_var.append(corr_coef[0, 1])

    return linear_model.coef_.ravel(), corr_coef_var, intercept

# the following function reads the features for regression and fits a linear model and saves the model parameters to file
# train to test ratio specifies the percentage of data used for training, M is the size of the dataset, model is a variable that keeps track of classification/regression globally, num_f is the size of the feature vector, C is the regularization parameter and x is the target variable (the stability scores)) 
#the function was originally meant for train: test ratio <1 and 10-fold cross validation for parameter tuning. 
def train_regression(train_to_test_ratio, M, model, num_f, C, x):
    features = open("regression_features_%d.txt" %train_size)
    f = features.readlines()
    K = int(train_to_test_ratio * M)
    overall_cum_coef = []
    for param in range(0, 1):
        cum_coef = 0
        for cv in range(0, 1):
            arr = numpy.arange(M)
            numpy.random.shuffle(arr)
            c = numpy.zeros((M, num_f + 1))
            for j in range(0, M):
                index = arr[j]
                line = list(f[index])
                one_line = []
                for jj in range(0, len(line) - 1):
                    one_line.append(int(line[jj]))
                one_line = numpy.asarray(one_line)
                current = numpy.hstack((one_line, x[index]))
                c[j, :] = current
            ff = c[:, 0:num_f]
            y = c[:, num_f:num_f + 1]
            y = numpy.reshape(y, (M,))
            x_train = ff[0:K]
            x_test = ff[K:M]
            y_train = y[0:K]
            y_test = y[K:M]
            coef1, coef, intercept = scikit_programs_regression(x_train, x_train, y_train,y_train, C)

        overall_cum_coef.append(float(cum_coef / 10))
    numpy.savetxt('coefficients_%d_%d.txt' % (train_size,model), coef1)
    numpy.savetxt('intercept_%d_%d.txt' % (train_size,model), numpy.array(intercept))
    numpy.savetxt('corr_coef_%d_%d.txt' % (train_size,model), coef)



# the following function reads the features for classifcation, fits a SVM and saves the model parameters to file
# train to test ratio specifies the percentage of data used for training, M is the size of the dataset, model is a variable that keeps track of classification/regression globally, num_f is the size of the feature vector, C is the regularization parameter and x is the target variable (the stability scores)) 
#the function was originally meant for train: test ratio <1 and 10-fold cross validation for parameter tuning. 
def train_classification(train_to_test_ratio, M, model, num_f, weight, x):
    features = open("classification_features_%d.txt" %train_size)
    f = features.readlines()

    K = int(train_to_test_ratio * M)
    overall_train_error = []
    overall_test_error = []
    overall_test_error_ones = []
    overall_test_error_minus_ones = []
    reg_weights = [weight]
    for param in range(0, len(reg_weights)):
        C = reg_weights[param]
        train_error = 0
        test_error = 0
        test_error_ones = 0
        test_error_minus_ones = 0
        for cv in range(0, 1):
            arr = numpy.arange(L)
            numpy.random.shuffle(arr)
            c = numpy.zeros((M, num_f + 1))
            for j in range(0, M):
                index = arr[j]
                line = list(f[index])
                one_line = []
                for jj in range(0, len(line) - 1):
                    one_line.append(int(line[jj]))
                one_line = numpy.asarray(one_line)
                current = numpy.hstack((one_line, x[index]))
                c[j, :] = current
            ff = c[:, 0:num_f]
            y = c[:, num_f:num_f + 1]
            y = numpy.reshape(y, (M,))
            x_train = ff[0:K]
            x_test = ff[K:M]
            y_train = y[0:K]
            y_train = numpy.array(y_train).reshape((-1, 1))
            y_test = y[K:M]
            y_test = numpy.array(y_test).reshape((-1, 1))

            e1, e2, e3, e4, coef, intercept = scikit_programs_classification(x_train, x_train, y_train,y_train, C)

            train_error = train_error + e1
            test_error = test_error + e2
            test_error_ones = test_error_ones + e3
            test_error_minus_ones = test_error_minus_ones + e4
        overall_train_error.append(float(train_error))
        overall_test_error.append(float(test_error))
        overall_test_error_ones.append(float(test_error_ones))
        overall_test_error_minus_ones.append(float(test_error_minus_ones))
        numpy.savetxt('coefficients_%d_%d.txt' % (train_size,model), coef)
        numpy.savetxt('intercept_%d_%d.txt' % (train_size,model), intercept)
        numpy.savetxt('train_error_%d_%s_%d.txt' % ( train_size, C, model),overall_train_error)
        numpy.savetxt('test_error_%d_%s_%d.txt' % (  train_size, C, model),overall_test_error)
        numpy.savetxt('test_error_ones_%d_%s_%d.txt' % (  train_size, C, model), overall_test_error_ones)
        numpy.savetxt('test_error_minus_ones_%d_%s_%d.txt' % (train_size, C, model), overall_test_error_minus_ones)

# this function learns a SVM classification model on binding data, x are the features and y are the binary predictions of binding
# it returns the overall and per class prediction error
def scikit_programs_classification(x_train, x_test, y_train, y_test, reg_weight):
    classifier = svm.LinearSVC(tol=0.0001, multi_class='ovr', fit_intercept=True, intercept_scaling=1,
                             class_weight="auto", C=reg_weight, verbose=0, random_state=None).fit(x_train, y_train)
    
    #uncomment the following command to learn a non-linear SVM with rbf kernel
    #logistic=svm.SVC(C=reg_weight, kernel="rbf",gamma=0.001,   class_weight="auto", verbose=0, random_state=None).fit(x_train,y_train)
    train_error = 0
    for i in range(0, len(x_train)):
        if (y_train[i] - classifier.predict(x_train[i,])) != 0:
            train_error = train_error + 1
    train_error = float(train_error) / len(x_train)

    test_error = 0
    for i in range(0, len(x_test)):
        if (y_test[i] - classifier.predict(x_test[i,])) != 0:
            test_error = test_error + 1
    test_error = float(test_error) / len(x_test)

    test_error_ones = 0
    count = 0
    for i in range(0, len(x_test)):
        if y_test[i] == 1:
            count = count + 1
            if (y_test[i] - classifier.predict(x_test[i,])) != 0:
                test_error_ones = test_error_ones + 1

    test_error_ones = float(test_error_ones) / count
    test_error_minus_one = 0
    count = 0
    for i in range(0, len(x_test)):
        if y_test[i] == -1:
            count = count + 1
            if (y_test[i] - classifier.predict(x_test[i,])) != 0:
                test_error_minus_one = test_error_minus_one + 1

    test_error_minus_one = float(test_error_minus_one) / count


    return train_error, test_error, test_error_ones, test_error_minus_one, classifier.coef_.ravel(), classifier.intercept_

#the following function maps the 20 amino acids to number indices
def obtain_index(acid):
    acid_vector = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    for i in range(0, len(acid_vector)):
        if acid == acid_vector[i]:
            return i

###########################################################################################################################
# the following program constructs the feature vectors from data
###########################################################################################################################
# open the train antibody and virus sequence files
a_file = open("train_set_ab_%d.txt" %train_size)
v_file = open("train_set_v_%d.txt" %train_size)
lines_a = a_file.readlines()
lines_v = v_file.readlines()
L = len(lines_a)

N_a = 27
N_v = 32

# open a new file to to write the feature vectors for classification (based on the ab and the v side and the pairwise interactions)
feature_file = open("classification_features_%d.txt" %train_size,"w")


#for each data point in the training set create a feature vector and write to the file
for data_point in range(0, len(lines_a)):
    features_a = lines_a[data_point]
    features_v = lines_v[data_point]
    count = 0

    for i in range(0, N_a):
        # the antibody side
        acid = features_a[i]
        acid_index = obtain_index(acid)
        for j in range(0, 20):
            if acid_index == j:
                feature_file.write("%d" % 1)
            else:
                feature_file.write("%d" % 0)
            count = count + 1

    for i in range(0, N_v):
        # the virus side
        acid = features_v[i]
        acid_index = obtain_index(acid)
        for j in range(0, 20):
            if acid_index == j:
                feature_file.write("%d" % 1)
            else:
                feature_file.write("%d" % 0)
            count = count + 1
    #the interaction pairs
    for i in range(0, N_a):
        acid_a = features_a[i]
        acid_index_a = obtain_index(acid_a)
        for m in range(0, N_v):
            acid_v = features_v[m]
            acid_index_v = obtain_index(acid_v)
            for j in range(0, 20):
                for k in range(0, 20):
                    if acid_index_a == j and acid_index_v == k:
                        feature_file.write("%d" % 1)
                    else:
                        feature_file.write("%d" % 0)
                    count = count + 1

    feature_file.write("\n")
feature_file.close()

#open a new feature file for regression (this feature vector is a function of the antibody sequence only)

feature_file = open("regression_features_%d.txt" %train_size,"w")

for data_point in range(0, len(lines_a)):
    features_a = lines_a[data_point]
    features_v = lines_v[data_point]
    count = 0

    for i in range(0, N_a):
        # the antibody side
        acid = features_a[i]
        acid_index = obtain_index(acid)
        for j in range(0, 20):
            if acid_index == j:
                feature_file.write("%d" % 1)
            else:
                feature_file.write("%d" % 0)
            count = count + 1

    feature_file.write("\n")
feature_file.close()


###########################################################################################################
# the following part trains the models on data using the functions defined above, we have set the regularization parameters to values mentioned in the paper
###########################################################################################################
z = numpy.loadtxt("train_set_scores_binding_%d.txt" %train_size)
y = numpy.loadtxt("train_set_scores_stability_%d.txt" %train_size)
train_classification(1, L, 1, N_a * 20 + N_v * 20 + N_a * N_v * 20 * 20, 0.001, z)
train_regression(1, L, 2, N_a * 20, 0.01, y)
