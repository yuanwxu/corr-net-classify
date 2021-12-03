# DGCNN classification of microbial correlation networks and
# node importance calculation

import stellargraph as sg

try:
    sg.utils.validate_notebook_version("1.2.1")
except AttributeError:
    raise ValueError(
        f"This notebook requires StellarGraph version 1.2.1, but a different version {sg.__version__} is installed.  Please see <https://github.com/stellargraph/stellargraph/issues/1172>."
    ) from None

import pandas as pd
import numpy as np
from stellargraph.mapper import PaddedGraphGenerator
from stellargraph.layer import DeepGraphCNN
from stellargraph import StellarGraph

from sklearn import model_selection

from tensorflow.keras import Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Conv1D, MaxPool1D, Dropout, Flatten
from tensorflow.keras.losses import binary_crossentropy, categorical_crossentropy
from tensorflow.keras.metrics import Precision, Recall, AUC
import tensorflow as tf
import pickle


def read_graphs2(W, node_features=None):
  """Read graphs into list of StellarGraph instances
     Args:
          W: dataframe of graphs with 4 columns: graph_id, source, target, weight
          node_features: node feature matrix, if None will use identity
  """
  out = list()
  gid0 = W.iloc[0].at['graph_id']
  g0 = W[W.graph_id == gid0]
  nodes = list(set(g0.source).union(set(g0.target)))
  if not node_features:
    node_features = sg.IndexedArray(np.identity(len(nodes)), index=list(nodes))
  
  for _, g in W.groupby('graph_id'):
    out.append(StellarGraph(nodes=node_features,
                            edges=g.drop(columns='graph_id'),
                            node_type_default='microbe',
                            edge_type_default='correlation'))
  
  # Check all graphs have the same number of nodes
  nn = [g.number_of_nodes() for g in out]
  if not all(nn[0] == x for x in nn):
    raise ValueError(
        "Not all graphs have same number of nodes."
    )
  return out

def train_dgcnn(graphs, graph_labels, n_epochs=50):
  """ Build and train DGCNN model """
  generator = PaddedGraphGenerator(graphs=graphs)
  k = graphs[0].number_of_nodes()  # the number of rows for the output tensor, no truncation
                                 # done here because all graphs have same number of nodes
  layer_sizes = [32, 32, 32, 1]

  dgcnn_model = DeepGraphCNN(
      layer_sizes=layer_sizes,
      activations=["tanh", "tanh", "tanh", "tanh"],
      k=k,
      bias=False,
      generator=generator,
  )
  x_inp, x_out = dgcnn_model.in_out_tensors()
  x_out = Conv1D(filters=16, kernel_size=sum(layer_sizes), strides=sum(layer_sizes))(x_out)
  x_out = MaxPool1D(pool_size=2)(x_out)
  x_out = Conv1D(filters=32, kernel_size=5, strides=1)(x_out)
  x_out = Flatten()(x_out)
  x_out = Dense(units=128, activation="relu")(x_out)
  x_out = Dropout(rate=0.5)(x_out)
  predictions = Dense(units=1, activation="sigmoid")(x_out)
  model = Model(inputs=x_inp, outputs=predictions)
  model.compile(
      optimizer=Adam(learning_rate=0.0001), loss=binary_crossentropy, 
      metrics=["accuracy", Precision(name='precision'), Recall(name='recall'), AUC(name='auc')],
  )

  train_graphs, test_graphs = model_selection.train_test_split(
    graph_labels, train_size=0.8, test_size=None, stratify=graph_labels
  )

  gen = PaddedGraphGenerator(graphs=graphs)

  # if use symmetric normalization, problem arise in negative degree values (because of negative correlations), 
  # and so can't take square root of those.
  train_gen = gen.flow(
      list(train_graphs.index),
      targets=train_graphs.values,
      batch_size=50,
      symmetric_normalization=False, 
      weighted=True,
  )

  test_gen = gen.flow(
      list(test_graphs.index),
      targets=test_graphs.values,
      batch_size=1,
      symmetric_normalization=False,
      weighted=True,
  )

  history = model.fit(
    train_gen, epochs=n_epochs, verbose=0, validation_data=test_gen, shuffle=True
  )
  # Print test set metrics
  test_metrics = model.evaluate(test_gen, verbose=0)
  print(f'Test Set Metrics: ')
  for name, val in zip(model.metrics_names, test_metrics):
      print("\t{}: {:0.4f}".format(name, val))
  
  return model, test_metrics, history

class ImportanceDGCNN:  

  def __init__(self, W, model, node_features=None):
    """Initialize object for computing importance in the DGCNN graph classification model
       Args:
            W: dataframe of graphs with 4 columns: graph_id, source, target, weight
            node_features: used to build StellarGraph graph instance, same as used in read_graphs
            model: the trained keras model of DGCNN
    """
    # Take any graph from W to find its nodes and edges
    gid0 = W.iloc[0].at['graph_id']
    g0 = W[W['graph_id'] == gid0]
    self.nodes = list(set(g0.source).union(set(g0.target)))
    self.edges = list(zip(g0.source, g0.target))
    self.ngraphs = W.groupby('graph_id').ngroups
    self.model = model

    if not node_features:
        node_features = sg.IndexedArray(np.identity(len(self.nodes)), index=list(self.nodes))
    self._W = W
    self._node_features = node_features

    # Check if all graphs have same set of edges
    for _, g in W.groupby('graph_id'):
      if set(zip(g.source, g.target)) != set(self.edges):
        raise ValueError("Not all graphs have the same set of edges. This case is not implemented.") 

  def _null_edge_graphs(self, val=0):
    """ Generator of StellarGraph graphs with exactly one edge set to 'val' (default 0)
    """
    for src, tar in self.edges:
      cond = (self._W['source'] == src) & (self._W['target'] == tar)
      W2 = self._W.copy()
      W2['weight'].mask(cond, val, inplace=True) # set weights corresonding to edge to 0
      for _, g in W2.groupby('graph_id'):
        yield StellarGraph(nodes=self._node_features, edges=g.drop(columns='graph_id'),
                           node_type_default='microbe', edge_type_default='correlation')


  def _null_node_graphs(self):
    """ Generator of StellarGraph graphs with all edges incident to a node set to 0
    """
    for n in self.nodes:
      cond = (self._W['source'] == n) | (self._W['target'] == n)
      W2 = self._W.copy()
      W2['weight'].mask(cond, 0, inplace=True)
      for _, g in W2.groupby('graph_id'):
        yield StellarGraph(nodes=self._node_features, edges=g.drop(columns='graph_id'),
                           node_type_default='microbe', edge_type_default='correlation')
        
        
  def _null_2nodes_graphs(self):
    """Generator of StellarGraph graphs with all edges incident to two nodes set to 0
    """
    for n1, n2 in self.edges:
      cond1 = (self._W['source'] == n1) | (self._W['target'] == n1)
      cond2 = (self._W['source'] == n2) | (self._W['target'] == n2)
      W2 = self._W.copy()
      W2['weight'].mask(cond1 | cond2, 0, inplace=True)
      for _, g in W2.groupby('graph_id'):
        yield StellarGraph(nodes=self._node_features, edges=g.drop(columns='graph_id'),
                           node_type_default='microbe', edge_type_default='correlation')
        
  def _null_nnodes_graphs(self, nlist):
    """ Generator of StellarGraph graphs with all edges incident to n nodes set to 0,
        Assume the first n-1 nodes are given as nlist, the generator then generates 
        graphs where each of the remaining len(self.nodes) - len(nlist) is added to 
        the given n-1 nodes, and edges linked to the resulting n nodes are set to 0.
    """
    from functools import reduce
    import operator

    if not set(nlist).issubset(self.nodes):
      raise ValueError("Not all provided nodes are found in the graph")

    conds = [(self._W['source'] == nd) | (self._W['target'] == nd) for nd in nlist]
    for n in self.nodes:
      if n in nlist:
        continue
      combined_cond = conds + [(self._W['source'] == n) | (self._W['target'] == n)]
      reduced_cond = reduce(operator.or_, combined_cond)
      W2 = self._W.copy()
      W2['weight'].mask(reduced_cond, 0, inplace=True)
      for _, g in W2.groupby('graph_id'):
        yield StellarGraph(nodes=self._node_features, edges=g.drop(columns='graph_id'),
                           node_type_default='microbe', edge_type_default='correlation')
    
        

  @staticmethod
  def _batch(iterable, n):
    """ Generate prediction batch of size n, using the grouper idiom """
    iters = [iter(iterable)] * n
    return zip(*iters)

  @staticmethod
  def compute_lor(pred, P_new):
    """ Compute log-odds ratio between new and original predicted probs 
        Args:
            pred: prediction on the original graphs, output of model.predict(),
                  shape N-by-1, where N number of graph instances
            P_new: predicition on new graphs, shape N-by-K, where K = number of
                   edges/nodes depending on edge or node importance 
        Returns:
            numpy array same shape as P_new
    """
    eps = 1e-6
    lo1 = np.log(P_new+eps) - np.log(1-P_new+eps)
    lo2 = np.log(pred+eps) - np.log(1-pred+eps)
    return lo1 - lo2

  def read_sg(self):
    """ Read graphs into list of StellarGraph instances """
    out = list()
    for _,g in self._W.groupby('graph_id'):
      out.append(StellarGraph(nodes=self._node_features,     
                              edges=g.drop(columns='graph_id'),
                              node_type_default='microbe',
                              edge_type_default='correlation')) 
    return out 

  def predict_graph(self, graphs):
    """Use the model to predict the probability of positive class
       Args:
          graphs: list of StellarGraph graph instances
    """
    fl = PaddedGraphGenerator(graphs=graphs).flow(range(len(graphs)), 
                                                  batch_size=len(graphs), 
                                                  symmetric_normalization=False, 
                                                  weighted=True)
    return self.model.predict(fl)

  def edge_imp(self, set_wt=0):
    """Calclate edge importance by change in log-odds between the original graph and one 
       where an edge weight is set to 'set_wt' (default 0)
    """
    sg = self.read_sg()
    pred = self.predict_graph(sg)
    P_new = np.empty((self.ngraphs, len(self.edges)))
    
    gen = self._null_edge_graphs(set_wt)
    for i, bch in enumerate(ImportanceDGCNN._batch(gen, self.ngraphs)):
      pred_new = self.predict_graph(list(bch)).reshape(-1)
      P_new[:,i] = pred_new
      print(f'{i}: EDGE {self.edges[i]} DONE.')

    LR = ImportanceDGCNN.compute_lor(pred, P_new)
    stats = self.summary_stats(LR, 'edge')
    self.LR_edge, self.LR_edge_stats = LR, stats
    return stats, LR

  def node_imp(self):
    """Calclate node importance by change in log-odds between the original graph and one 
       where all edges linked to a node is set to 0.
    """
    sg = self.read_sg()
    pred = self.predict_graph(sg)
    P_new = np.empty((self.ngraphs, len(self.nodes)))

    gen = self._null_node_graphs()
    for i, bch in enumerate(ImportanceDGCNN._batch(gen, self.ngraphs)):
      pred_new = self.predict_graph(list(bch)).reshape(-1)
      P_new[:,i] = pred_new
      print(f'{i}: NODE {self.nodes[i]} DONE.')

    LR = ImportanceDGCNN.compute_lor(pred, P_new)
    stats = self.summary_stats(LR, 'node')
    self.LR_node, self.LR_node_stats = LR, stats
    return stats, LR

  def node_pair_imp(self):
    """Calculate node pair importance by knocking out each pair of nodes
    """
    sg = self.read_sg()
    pred = self.predict_graph(sg)
    P_new = np.empty((self.ngraphs, len(self.edges)))

    gen = self._null_2nodes_graphs()
    for i, bch in enumerate(ImportanceDGCNN._batch(gen, self.ngraphs)):
      pred_new = self.predict_graph(list(bch)).reshape(-1)
      P_new[:,i] = pred_new
      print(f'{i}: NODES {self.edges[i][0], self.edges[i][1]} DONE.')

    LR = ImportanceDGCNN.compute_lor(pred, P_new)
    stats = self.summary_stats(LR, 'node2')
    self.LR_node2, self.LR_node2_stats = LR, stats
    return stats, LR

  def nnode_imp(self, n):
    """Calculate n-node importance by knocking out n nodes, using a greedy search
       strategy where after the first node resulting in maximum change in log-odds 
       is found, the node from the remaining node set resulting in maximun change 
       in log-odds is added to form 2-node, and it continues until n nodes are added

       Returns:
            List of tuples. The first component is the names of knocked-out nodes,
            the second component is LR of shape (n_graphs, n_nodes-k+1), 
            with n_nodes-k+1 equal to the length of the first component. 
            Return all k-node importance from k = 1 to n
    """
    sg = self.read_sg()
    pred = self.predict_graph(sg)
    n_full = list(self.nodes)
    nlist = []

    out = [] 
    for k in range(1, n+1):
      P_new = np.empty((self.ngraphs, len(self.nodes)-k+1))
      if k == 1:
        gen = self._null_node_graphs()
      else:
        gen = self._null_nnodes_graphs(nlist)

      for i, bch in enumerate(ImportanceDGCNN._batch(gen, self.ngraphs)):
        pred_new = self.predict_graph(list(bch)).reshape(-1)
        P_new[:,i] = pred_new

      # Find which node to add to nlist
      LR = ImportanceDGCNN.compute_lor(pred, P_new)
      maxi = np.argmax(np.median(np.abs(LR), axis=0)) # index of node with max median absolute LR 
      n_remain = [nn for nn in n_full if nn not in nlist]
      out.append(([tuple(nlist + [x]) for x in n_remain], LR))
      nlist = nlist + [n_remain[maxi]]
    
    return out


  def summary_stats(self, LR, which):
    """ Get mean, median and std err of log-odds ratio """
    lor_mean, lor_med = np.mean(LR, axis=0), np.median(LR, axis=0)
    lor_std = np.std(LR, axis=0)
    df = pd.DataFrame({'lor_mean': lor_mean,
                       'lor_med': lor_med, 
                       'std_err': lor_std/np.sqrt(LR.shape[0])})
    if which == 'edge':
      df['source'] = [e[0] for e in self.edges]
      df['target'] = [e[1] for e in self.edges]

    if which == 'node':
      df['node'] = self.nodes
    
    if which == 'node2':
      df['node1'] = [e[0] for e in self.edges]
      df['node2'] = [e[1] for e in self.edges]
    
    return df

def save_imp_res(node_res, edge_res=None, runid='run'):
    with open('imp_' + runid + '.pkl', 'wb') as pickle_out:
        pickle.dump(node_res, pickle_out)  # serialize node importance result
        if edge_res:
            pickle.dump(edge_res, pickle_out)


# Main script
def main():
    K = 10 # number of DGCNN runs
    N_EPOCH = 100 # number of epochs per run
    IMP_SUBSAMPLE = 0.3 # fraction of graphs sampled per DGCNN run used to calculate importance,
                        # decrease this if out-of-memory
    N = 21 # up to N-node importance 

    W_ctrl = pd.read_csv('W_mbqc_ctrl.csv', dtype={'graph_id': 'int'})
    W_case = pd.read_csv('W_mbqc_case.csv', dtype={'graph_id': 'int'})

    g_ctrl = read_graphs2(W_ctrl)
    g_case = read_graphs2(W_case)
    graphs = g_ctrl + g_case
    graph_labels = pd.Series(len(g_ctrl) * ['0'] + len(g_case) * ['1'])
    graph_labels = pd.get_dummies(graph_labels, drop_first=True)

    nn_imp_result = []
    history = []
    test_metrics = []
    for k in range(K):
        print(f"DGCNN Run {k+1} ...")
        model, tm, hist = train_dgcnn(graphs, graph_labels, n_epochs=N_EPOCH)
        history.append(hist)
        test_metrics.append(tm)

        # Subsample a fraction of graphs to compute n-node importance
        rind = 1 + np.random.choice(len(g_case), size=round(len(g_case)*IMP_SUBSAMPLE), replace=False)
        W_case_sub = W_case[W_case['graph_id'].isin(rind)]
        imp_calculator = ImportanceDGCNN(W_case_sub, model)
        nn_imp = imp_calculator.nnode_imp(N) # N-node imp

        # Each model can result in different selected n-nodes (n=1,...,N), 
        # so will not attempt to combine the LR arrays
        nn_imp_result.append(nn_imp) 
    
    # Save n-node importance result as pickle
    save_imp_res(nn_imp_result, runid="mbqc")

    # Save loss and accuracy trajectory 
    fig = sg.utils.plot_history(history, return_figure=True)
    fig.savefig('hist_mbqc.png')

    # Save test set metrics across all DGCNN runs
    test_metrics = pd.DataFrame(np.array(test_metrics), columns=['loss', 'accuracy', 'precision', 'recall', 'auc'])
    test_metrics.agg(['mean', 'std']).to_csv("test_metrics_mbqc.csv")


if __name__ == "__main__":
    main()
