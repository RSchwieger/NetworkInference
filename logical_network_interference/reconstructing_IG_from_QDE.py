from constraint import *
from itertools import product
from PyBoolNet import InteractionGraphs as IGs
import networkx as nx
from logical_network_interference.example2 import edges, ig_size

"""
Given a subset of edges of the qde-graph edges of the interaction graph are computed via a constraint solver.
Each edge leads to one constraints.
Actually also the non-existence leads to a constraint. This is not
taken into account, since it is unlikely to obtain empirical proofs for this.
Instead all solutions are searched for consistent parts of the solutions (edges which occur with the same label in
all solutions).
"""


problem = Problem()


problem.addVariables([(i,j) for (i,j) in product(range(ig_size), range(ig_size)) if i != j], [-1,1])

def add_constraint(edge, problem):
    v = edge[0]
    w = edge[1]
    diff = [i for i in range(len(w)) if v[i] != w[i]]
    equal = [i for i in range(len(v)) if v[i] == w[i]]
    # print(str(v)+", "+str(w)+" equal: "+str(equal)+" diff: "+str(diff))
    if len(equal) == 0:
        print("EXCEPTION: At least one entry must be equal on the edge "+str(edge))

    for index_diff in range(len(diff)):
        # print(str(diff)+", index: "+str(index_diff))
        def condition(*row_equal):
            # The index gives us also the corresponding object in "equal"
            # print("v = "+str(v)+", w = "+str(w))
            for index_equal in range(len(row_equal)):
                # print("index = "+str((diff[condition.index_diff],equal[index_equal]))+", equal = "+str(equal)+", diff = "+str(diff))
                # print("arg = " + str(condition.arg) +", row_i = " + str(row_equal))
                if row_equal[index_equal] == w[diff[condition.index_diff]]*v[equal[index_equal]]:
                    return True
            return False
        condition.index_diff = index_diff
        condition.arg = tuple([(diff[condition.index_diff],j) for j in equal])

        # print("Add constraint with "+str(condition.arg)+" for "+str(v)+", "+str(w))
        problem.addConstraint(condition,
                              condition.arg)

for edge in edges:
    add_constraint(edge, problem)
#print("constraints added")
solutions = problem.getSolutions()
print(len(solutions))
print(solutions)

def test_if_all_solutions_have_same_value(solutions, key):
    val = solutions[0][key]
    for solution in solutions:
        if solution[key] != val:
            return False
    return True


recontructed_edges = []
if len(solutions) > 0:
    keys = solutions[0].keys()
    for key in keys:
        if test_if_all_solutions_have_same_value(solutions, key):
            recontructed_edges += [key]
    print([[edge, solutions[0][edge]] for edge in recontructed_edges])

    iG = nx.DiGraph()
    iG.add_nodes_from([str(i) for i in range(ig_size)])
    for (v,w) in recontructed_edges:
        iG.add_edge(str(w), str(v))
        iG.edge[str(w)][str(v)]["sign"] = {solutions[0][(v,w)]}
        if solutions[0][(v,w)] == -1:
            iG.edge[str(w)][str(v)]["color"] = "red"
    IGs.igraph2image(iG, "./reconstructed_interaction_graph.pdf")
