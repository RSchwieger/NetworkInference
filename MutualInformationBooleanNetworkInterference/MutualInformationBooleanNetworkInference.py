"""
Implementation of the algorithm described in

A novel mutual information-based Boolean network inference method from time-series gene expression data

"""

from copy import deepcopy
from itertools import product, repeat
from math import log2 as log

from MutualInformationBooleanNetworkInterference.test_case import v,K


def time_series_to_random_variable(time_series):
    random_variable = {}
    for elem in time_series:
        if elem in random_variable.keys():
            random_variable[elem] += 1/len(time_series)
        else:
            random_variable[elem] = 1 / len(time_series)
    return random_variable

X = [time_series_to_random_variable(elem) for elem in v]

def gene_wise_dynamics_consistency(observed_variable, predicted_variable):
    return sum([(x == y) for (x,y) in zip(observed_variable, predicted_variable)])/(len(v)-1) # ToDo: -1 richtig?

def entropy(X):
    return sum([X[x]*log(X[x]) for x in X.keys()])*-(1)

def joined_entropy(X,Y):
    Z = {(keyX, keyY): X[keyX] * Y[keyY] for (keyX, keyY) in product(X.keys(), Y.keys())}
    return entropy(Z)

def mutual_information(X,Y):
    return entropy(X)+entropy(Y)-joined_entropy(X,Y)

def argmax(dict):
    max_key, max_value = None, None
    for key in dict.keys():
        if max_key is None:
            max_key = key
            max_value = dict[key]
        elif dict[key] > max_value:
            max_key = key
            max_value = dict[key]
    return max_key


def MIFS(v_0, W, k):
    """

    :param v_0: index
    :param W: list of indices
    :param k: natural number
    :return: set of regulators
    """
    S = []
    for i in range(k):
        dictio = {w: mutual_information(X[v_0], X[w])-sum([mutual_information(X[w], X[s]) for s in S]) for w in W} # ToDo Besser ordnen
        v = argmax(dictio)
        S += [v]
        W.remove(v)
    return S

def boolean_conjunction_generator(S):
    """
    Generator for boolean conjunctions of variables in S.
    """
    for elem in product( *list(repeat([1,0],len(S)+1)) ):
        def boolean_function(x):
            atoms = [x[boolean_function.S[i]] if boolean_function.elem[i]
                     else (not x[boolean_function.S[i]]) for i in range(len(boolean_function.S))]
            result = True
            for temp in atoms:
                result = result and temp
            #print("S = "+str(boolean_function.S))
            #print("x "+str(x))
            #print("elem = "+str(boolean_function.elem))
            #print("atoms "+str(atoms)+" disjunction? "+str(elem[-1]))
            #print("result "+str(result))
            # If the last digit is False negate the function to turn it into a disjunction
            if elem[-1]:
                return result
            else:
                return not result



        boolean_function.elem = elem
        boolean_function.S = S
        yield boolean_function

def search_update_rule(v_0, S):
    """

    :param v_0: index
    :param S: set of regulators
    :return:
    """
    # ToDo: Disjunktionen S+1 blabla
    E_MAX = -1 # any value <0
    boolean_function_max = None
    for boolean_function in boolean_conjunction_generator(S):
        #print(S)
        #print("v_0 "+str(v_0)+" update")
        time_series = [boolean_function(x) for x in zip(*v)]
        E_temp = gene_wise_dynamics_consistency(v[v_0], time_series)
        #print("sim_time_series["+str(v_0)+"] = "+str(time_series))
        #print("time_series[" + str(v_0) + "] = " + str(v[v_0]))
        #print("boolean_function = " + str(boolean_function.elem))
        #print("E_temp = " + str(E_temp))
        if E_temp > E_MAX:
            E_MAX = E_temp
            boolean_function_max = boolean_function
    #print(10*"#####")
    #print("E_MAX = "+str(E_MAX))
    #print("boolean_function_max = " + str(boolean_function_max.elem))
    #print(10 * "#####")
    return E_MAX, boolean_function_max


def SWAP(v_0, S, W):
    """

    :param v_0: index of target variable
    :param S:  set of selected Boolean variables
    :param W: set of unselected Boolean variables
    :return: S, E_MAX, boolean_function_max
    """
    E_MAX, boolean_function_max = search_update_rule(v_0, S)
    print(10*"-------")
    print("E_MAX = "+str(E_MAX))
    print("boolean_function_max = " + str(boolean_function_max.elem))
    print(10 * "-------")
    for i in range(len(W)):
        for j in range(len(S)):
            print("S = "+str(S[0:j]+[W[i]]+S[j+1:]))
            E_temp, boolean_function_temp = search_update_rule(v_0, S[0:j]+[W[i]]+S[j+1:])
            if E_temp>E_MAX:
                temp_S_j = S[j]
                S = S[0:j]+[W[i]]+S[j+1:]
                W = W[0:i]+[temp_S_j]+W[i+1:]
                E_MAX, boolean_function_max = E_temp, boolean_function_temp
    print(10*"#####")
    print("E_MAX = "+str(E_MAX))
    print("boolean_function_max = " + str(boolean_function_max.elem))
    print(10 * "#####")
    return S, E_MAX, boolean_function_max

def get_i_th_constant(i):
    def constant_function(x):
        return 1  # ToDo: Kein Regulator => constant?
    return constant_function

def MIBNI():
    bool_func = [None for _ in v]
    S = list(range(len(v)))
    W = list(range(len(v)))
    for v_0 in range(len(v)):
        if entropy(X[v_0]) == 0:
            bool_func[v_0] = get_i_th_constant(v_0)
            continue
        else:
            for k in range(1,K+1):
                #print(str(v_0)+" nr. reg. "+" -> "+str(k))
                #print("S before "+str(S))
                S = MIFS(v_0, deepcopy(W), k) # ToDo: deepcopy correct?
                #print("S after MIFS "+str(S))
                #print("S before " + str(S))
                S, E_MAX, bool_func[v_0] = SWAP(v_0, S, [w for w in W if w not in S]) # ToDo: Muss ich hier W dauerhaft ändern?
                #print("S after SWAP " + str(S))
                #print(str(v_0)+" E_MAX " + " -> " + str(E_MAX))
                if E_MAX == 1:
                    break
    return bool_func

if __name__ == "__main__":
    bool_func = MIBNI()
    # --------------------------
    resulting_time_series = [[li[0] for li in v]]
    for i in range(len(v)):
        resulting_time_series += [[bool_func[i](resulting_time_series[-1]) for i in range(len(v))]]
    print(resulting_time_series)



#print(SWAP(0, [1,2], [0,3]))
#S, E_MAX, boolean_function_max = SWAP(0, [1,2], [0,3])
#print(E_MAX)
#print([boolean_function_max(x) for x in zip(*v)])
#print(MIFS(0, [1,2, 3], 2))
