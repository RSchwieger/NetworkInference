"""
Implementation of the algorithm described in

A novel mutual information-based Boolean network inference method from time-series gene expression data

"""

from copy import deepcopy
from itertools import product, repeat
from math import log2 as log

from MutualInformationBooleanNetworkInterference.test_case import observed_state_sequence,K


def time_series_to_random_variable(time_series):
    """
    Converts the time series data into a random variable
    :param time_series: list of states
    :return: random variable saved as dict of the form event->probability
    """
    random_variable = {}
    for elem in time_series:
        if elem in random_variable.keys():
            random_variable[elem] += 1/len(time_series)
        else:
            random_variable[elem] = 1 / len(time_series)
    print(random_variable)
    assert (max(random_variable.values()) <= 1 and min(random_variable.values())>=0)
    return random_variable

# Convert the observed sequence of states into a list of random variables (saved as list of dictionaries)
X = [time_series_to_random_variable(elem) for elem in observed_state_sequence]

def gene_wise_dynamics_consistency(observed_variable, predicted_variable):
    """
    This corresponds to E(v,v') in the paper
    :param observed_variable: list of states in {0,1}
    :param predicted_variable: list of states in {0,1}
    :return: gene_wise_dynamics_consistency
    """
    return sum([(x == y) for (x,y) in zip(observed_variable, predicted_variable)])\
           / (len(observed_state_sequence) - 1) # ToDo: -1 richtig?

def entropy(X):
    """
    This corresponds to the definition of H(X)=-sum p(x) * log p(x) in the paper
    :param X: random variable (dictionary of the form event->probability)
    :return: entropy of a random variable
    """
    return sum([X[x]*log(X[x]) for x in X.keys()])*-(1)

def joined_entropy(X,Y):
    """
    Corresponds to H(X,Y) in the paper
    :param X: random variable
    :param Y: random variable
    :return: joined entropy
    """
    # ToDo: Check
    # construct a new two dimensional random variable and use the entropy function
    Z = {(keyX, keyY): X[keyX] * Y[keyY] for (keyX, keyY) in product(X.keys(), Y.keys())}
    return entropy(Z)

def mutual_information(X,Y):
    """
    Corresponds to I(X;Y) in the paper
    :param X: random variable
    :param Y: random variable
    :return: entropy
    """
    return entropy(X)+entropy(Y)-joined_entropy(X,Y)

def argmax(dict):
    """
    Finds any key corresponding to one of the biggest values of the dictionary
    :param dict: dictionary
    :return: the key corresponding to one of the biggest values
    """
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
    Corresponds to the Pseudo-code of the MIFS function in the paper
    :param v_0: index
    :param W: list of indices
    :param k: natural number
    :return: set of regulators
    """
    S = []
    for i in range(k):
        # Create a dictionary corresponding to the expression I(v_0, w) - sum I(w,s) in the paper
        dictio = {w: mutual_information(X[v_0], X[w])-sum([mutual_information(X[w], X[s]) for s in S]) for w in W} # ToDo Besser ordnen
        # Find argmax
        v = argmax(dictio)
        S += [v] # S union {v}
        W.remove(v) # W\{v}
    return S

class BooleanFunction(object):
    """
    Functor that is used in the boolean function generator
    """
    def __init__(self, elem, S):
        """
        :param elem: Specifies which elements are negated in the function
        :param S: regulators
        """
        self.elem = elem
        self.S = S
    def __call__(self, x):
        """
        :param x: any boolean state {0,1}^n
        :return: 0 or 1
        """
        # Make a list. For each regulator check if elem[i] is zero or one. If zero negate the entry
        atoms = [x[self.S[i]] if self.elem[i]
                 else (not x[self.S[i]]) for i in range(len(self.S))]

        # Now compute the disjunction of atoms, which looks somehow like this: atoms = [True, False, ....]
        result = True
        for temp in atoms:
            result = result and temp
        if self.elem[-1]:
            return result
        else:
            return not result

def boolean_conjunction_generator(S):
    """
    Generator for boolean conjunctions of variables in S.
    """
    # The last digit [0,1], says us if we negate the complete function or not
    for elem in product( *list(repeat([1,0],len(S))), [0,1]):
        yield BooleanFunction(elem=elem, S=S)

def search_update_rule(v_0, S):
    """

    :param v_0: index
    :param S: set of regulators
    :return: (result of the fitting, best boolean function)
    """
    # ToDo: Disjunktionen S+1 blabla
    E_MAX = -1 # any value <0
    boolean_function_max = None
    for boolean_function in boolean_conjunction_generator(S):
        # Compute the time series that the current boolean function would produce
        # observed_state_sequence has the form [states of 1. component, states of 2. component, ...]
        # Each state is a list of zeros and ones - i.e. smth like this [0,1,1,0,1]
        # saying which value the component has in which time step

        # zip(*observed_state_sequence) rotates the list into [state1, state2, ...]
        # The states have the form state_i = [first component, second component,...]
        time_series = [boolean_function(x) for x in zip(*observed_state_sequence)]
        # Compare it with observed time series
        E_temp = gene_wise_dynamics_consistency(observed_state_sequence[v_0], time_series)
        if E_temp > E_MAX:
            E_MAX = E_temp
            boolean_function_max = boolean_function
    return E_MAX, boolean_function_max


def SWAP(v_0, S, W):
    # ToDo: Remove self regulators
    """
    Corresponds to the SWAP method in the pseudo-code of the paper
    :param v_0: index of target variable
    :param S:  set of selected Boolean variables
    :param W: set of unselected Boolean variables
    :return: S, E_MAX, boolean_function_max
    """
    E_MAX, boolean_function_max = search_update_rule(v_0, S)
    for i in range(len(W)):
        for j in range(len(S)):
            # search the "best" boolean function wtíth the regulators S\W[i]
            E_temp, boolean_function_temp = search_update_rule(v_0, S[0:j]+[W[i]]+S[j+1:])
            if E_temp>E_MAX:
                temp_S_j = S[j]
                S = S[0:j]+[W[i]]+S[j+1:]
                W = W[0:i]+[temp_S_j]+W[i+1:]
                E_MAX, boolean_function_max = E_temp, boolean_function_temp
    return S, E_MAX, boolean_function_max

def get_i_th_constant(i):
    def constant_function(x):
        return 1  # ToDo: Kein Regulator => constant?
    return constant_function

def MIBNI():
    """
    The main function of the algorithm
    :return: boolean function
    """
    bool_func = [None for _ in observed_state_sequence]
    W = list(range(len(observed_state_sequence)))     # indices of all Boolean variables
    for v_0 in range(len(observed_state_sequence)):   # iterate over all Boolean variables
        if entropy(X[v_0]) == 0: # if the entropy is zero => no regulatory gene for v_0 ( = constant function)
            bool_func[v_0] = get_i_th_constant(v_0)
            continue
        else:
            for k in range(1,K+1): # Iterator over the potential number of regulators
                # Get candidate regulators (actually indices of the regulators) with MIFS
                S = MIFS(v_0, deepcopy(W), k) # ToDo: deepcopy correct?
                # ToDo: Muss ich hier W dauerhaft ändern?
                S, E_MAX, bool_func[v_0] = SWAP(v_0, S, [w for w in W if w not in S]) # SWAP(v0, S, W\S)
                if E_MAX == 1:
                    break
    return bool_func

if __name__ == "__main__":
    bool_func = MIBNI()
    # --------------------------
    resulting_time_series = [[li[0] for li in observed_state_sequence]]
    for i in range(len(observed_state_sequence)):
        resulting_time_series += [[bool_func[i](resulting_time_series[-1]) for i in range(len(observed_state_sequence))]]
    print(resulting_time_series)