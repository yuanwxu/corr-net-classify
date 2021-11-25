import pickle

def load_imp_res(file, edge=None):
  with open(file, 'rb') as pickle_in:
    node_res = pickle.load(pickle_in)
    if edge:
        edge_res = pickle.load(pickle_in)
        return node_res, edge_res
  return node_res

