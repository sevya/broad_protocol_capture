#!/usr/bin/python
import numpy
import sys


train_size = int(sys.argv[1])

acid_vector = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X"]
M = 20


def obtain_index(acid):
    acid_vector = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    for i in range(0, len(acid_vector)):
        if acid == acid_vector[i]:
            return i


N_ab = 27
N_v = 32


#the following loads the saved parameters of the model trained on all the data (in this case 30 virus sequences)
virus_file = open("v_30.txt")
virus = virus_file.readlines()
weights = numpy.loadtxt("coefficients_%d_1.txt" %(30))
intercept = numpy.loadtxt("intercept_%d_1.txt" % (30))
weights_s = numpy.loadtxt("coefficients_%d_2.txt" % (30))
intercept_s = numpy.loadtxt("intercept_%d_2.txt" % (30))


# open the optimized (breadth maximized) antibody and test it against the saved model
ab_file = open("BM_ab_%d.txt" % (train_size))
ab = ab_file.readlines()
native_antibody = ab[0]
native_antibody = native_antibody[0:N_ab]


# map the optimized antibody to feature space and evaluate the prediction of binding by a) multiplying the feature with saved parameter vector
#and b) adding the intercept 
bind_count = 0
for v in range(0, 30):
    current_v = virus[v]
    current_v = current_v[0:32]
    count = 0
    feature = numpy.zeros((346780, 1))
    for ab_pos in range(0, 27):
        acid = native_antibody[ab_pos]
        acid_index = obtain_index(acid)
        for x in range(0, 20):
            if acid_index == x:
                feature[count] = 1
            else:
                feature[count] = 0
            count = count + 1
    for w in range(0, 32):
        acid = current_v[w]
        acid_index = obtain_index(acid)
        for x in range(0, 20):
            if acid_index == x:
                feature[count] = 1
            else:
                feature[count] = 0
            count = count + 1
    for u in range(0, 27):
        acid_ab = native_antibody[u]
        acid_index_ab = obtain_index(acid_ab)
        for m in range(0, 32):
            acid_v = current_v[m]
            acid_index_v = obtain_index(acid_v)
            for w in range(0, 20):
                for x in range(0, 20):
                    if acid_index_ab == w and acid_index_v == x:
                        feature[count] = 1
                    else:
                        feature[count] = 0
                    count = count + 1

    product = 0
    for c in range(0, 346780):
        product = product + weights[c] * feature[c]

    if product + intercept < 0:
        bind_count = bind_count + 1



# compute breadth
breadth = (float(bind_count)) / 30

# write breadth to file 
b_file = open("predicted_breadth_%d.txt" % (train_size), "w")
b_file.write("%f" % breadth)
b_file.close()


# map the optimized antibody to the feature space of the regression model  

count = 0
feature = numpy.zeros((20 * 27, 1))
for ab_pos in range(0, 27):
    acid = native_antibody[ab_pos]
    acid_index = obtain_index(acid)
    for x in range(0, 20):
        if acid_index == x:
            feature[count] = 1
        else:
            feature[count] = 0
        count = count + 1


# evaluate predicted stability score
product = 0
for c in range(0, 20 * 27):
    product = product + weights_s[c] * feature[c]

stability = product + intercept_s

# write the predicted score to file
b_file = open("predicted_stability_%d.txt" % (train_size), "w")
b_file.write("%d" % stability)
b_file.close()
