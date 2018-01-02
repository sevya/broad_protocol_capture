#!/usr/bin/python
import numpy
import sys
import cplex
from cplex.exceptions import CplexError

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

# following is a constant used in the ILP
Z = 50

# read the data and the saved model parameters to construct the ILP
virus_file = open("v_%d.txt" % (train_size))
virus = virus_file.readlines()
weights = numpy.loadtxt("coefficients_%d_1.txt" % (train_size))
intercept = numpy.loadtxt("intercept_%d_1.txt" % (train_size))
weights_s = numpy.loadtxt("coefficients_%d_2.txt" % (train_size))
intercept_s = numpy.loadtxt("intercept_%d_2.txt" % (train_size))
ab_weights = weights[0:N_ab * 20]
v_weights = weights[N_ab * 20:N_v * 20 + N_ab * 20]
train_size = len(virus)

# write the ILP in cplex format
variables = []
objective = []
upper_bounds = []
v_types = []

for i in range(0, N_ab):
    for j in range(0, M):
        variables.append("a_" + str(i) + "_" + str(j))
        objective.append(weights_s[i * 20 + j])
        upper_bounds.append(1.0)
        v_types.append("I")

prob = cplex.Cplex()
prob.objective.set_sense(prob.objective.sense.minimize)
prob.variables.add(obj=objective, ub=upper_bounds, names=variables, types=v_types)

var = []
coef = []
right_hand_side = intercept_s

right_hand_side = -1 * right_hand_side - Z

for t in range(0, train_size):

    native_virus = virus[t]
    native_virus = native_virus[0:32]
    
    native_feature_vector = []
    for i in range(0, N_v):
        acid = native_virus[i]
        acid_index = obtain_index(acid)
        for j in range(0, 20):
            if acid_index == j:
                native_feature_vector.append(1)
            else:
                native_feature_vector.append(0)


    var = []
    coef = []
    right_hand_side = intercept

    for i in range(0, N_v):
        for j in range(0, M):
            index = i * M + j
            right_hand_side = right_hand_side + native_feature_vector[index] * v_weights[index]
    right_hand_side = -0.0000001 - right_hand_side
    for k in range(0, N_ab):

        for u in range(0, M):
            var.append('a_' + str(k) + '_' + str(u))
            index = k * 20 + u
            antibody_coef_sum = ab_weights[index]
            for i in range(0, N_v):
                for j in range(0, M):
                    matrix = weights[((N_ab * 20 + N_v * 20) + k * 20 * 20 * N_v + 20 * 20 * i):(
                    (N_ab * 20 + N_v * 20) + k * 20 * 20 * N_v + 20 * 20 * (i + 1))]
                    Q_ki = numpy.reshape(matrix, (20, 20))
                    antibody_coef_sum = antibody_coef_sum + native_feature_vector[20 * i + j] * Q_ki[u][j]
            coef.append(antibody_coef_sum)
    prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=var, val=coef)], senses=["L"], rhs=[right_hand_side])


for i in range(0, N_ab):
    var = []
    coef = []
    for j in range(0, M):
        var.append('a_' + str(i) + '_' + str(j))
        coef.append(1)
    prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=var, val=coef)], senses=["E"], rhs=[1])


try:
    prob.write("BM_%d.lp" %train_size)
    prob.solve()
except CplexError, exc:
    print exc

print ()
# solution.get_status() returns an integer code
print "Solution status = ", prob.solution.get_status(), ":",
# the following line prints the corresponding string
print(prob.solution.status[prob.solution.get_status()])
print("Solution value  = ", prob.solution.get_objective_value())

x = prob.solution.get_values()

numcols = prob.variables.get_num()
#for j in range(numcols):
#    print("Column %d:  Value = %10f %s " % (j, x[j], variables[j]))



# the following writes the optimized ab sequence to file
write_file = open("BM_ab_%d.txt" % (train_size), "w")
best_ab = []
for i in range(0, N_ab):
    for j in range(0, M):
       if int(x[i * M + j]) == 1:
            acid = acid_vector[j]
            best_ab.append(acid)
            write_file.write("%s" % acid)

write_file.close()

print best_ab
