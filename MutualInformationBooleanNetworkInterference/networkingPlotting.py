from PyBoolNet import InteractionGraphs as IGs
from PyBoolNet import QuineMcCluskey as QMC
from PyBoolNet import StateTransitionGraphs as STGs
from MutualInformationBooleanNetworkInterference.MutualInformationBooleanNetworkInference import MIBNI

bool_func = MIBNI()
new_bool_func = [None for func in bool_func]

def get_f_i(bool_func, i):
    # x0, x1, x2, x3
    # x0, x1, x2, x3, x4, x5, x6, x7, x8, x9
    def f_i(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9):
        return bool_func[i]([x0, x1, x2, x3, x4, x5, x6, x7, x8, x9])
    return f_i

for i in range(len(new_bool_func)):
    new_bool_func[i] = get_f_i(bool_func,i)

variable_to_boolean_function = {"x"+str(i): new_bool_func[i] for i in range(len(bool_func))}
print(5*"###")
prime_implicants = QMC.functions2primes(variable_to_boolean_function)
# Create the interaction graph
igraph = IGs.primes2igraph(prime_implicants)

# Construct the state transition graph
state_transition_graph = STGs.primes2stg(prime_implicants, "synchronous")

STGs.stg2image(state_transition_graph, "stg.pdf", LayoutEngine="dot")
IGs.create_image(prime_implicants, "igraph.pdf")
