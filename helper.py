import os
import re

# =============================================================================
# USEFUL FUNCTIONS
# =============================================================================
def setwd():
    os.chdir(os.path.dirname(__file__))

def flash():
    from IPython import get_ipython
    get_ipython().magic('%reset -sf')
    get_ipython().magic('%clear')

def reader():
    file = open('input.txt', 'r')
    read = file.readlines()
    data = []
    for line in read:
        data.append(line.strip())
    file.close()
    return data

def writer(data):
    with open('out.txt', 'w') as file:
        strings = '\n'.join(data)
        strings = data
        file.write(str(strings))
        file.close()
        
def list_writer(lst):
    with open('out.txt', 'w') as file:
        out = ''
        for ele in lst:
            out += str(ele) + '\n'
        file.write(out)
        file.close()

def dict_values(d):
    out = []
    for key in d:
        for value in d[key]:
            out.append(value)
    return out

def pickle_dict(mydict, filename):
    import pickle
    output = open(filename, 'wb')
    pickle.dump(mydict, output)
    output.close()

def graph_from_adj(filename):
    import networkx as nx
    '''
    creates graph from text file of format:
        0->4:11
        
    where each line shows, an edge in the form out_node->in_node:weight
        
    '''
    file = open(filename, 'r')
    read = file.readlines()
    data = []
    for line in read:
        temp = []
        for ele in re.split('->|:', line):
            temp.append(ele.strip())
        data.append(temp)
    file.close()
    G = nx.DiGraph()
    for edge in data:
        G.add_edge(edge[0], edge[1], weight=int(edge[2]))
    return G
