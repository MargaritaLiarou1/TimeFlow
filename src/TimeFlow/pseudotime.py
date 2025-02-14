from __future__ import division
import pandas as pd
import numpy as np
from sklearn.neighbors import kneighbors_graph
import igraph as ig


np.random.seed(1)


class Pseudotime:

    """
    Pseudotime class represents the data with a connected k-Nearest Neighbour Graph. The graph is
    constructed based on the cell pairwise Euclidean distances. The graph is initially unweighted, but weights
    are later added based on the absolute difference of the log probability density at the endpoints of each edge.
    The user should specify the index of the root cell.
    """

    def __init__(self, markers, markers_array, density_values, root_id, neighbours_number):

        # pandas Dataframe with markers
        self.markers = markers
        # Array with CD markers
        self.markers_arr = markers_array
        # DataFrame with density values at each cell
        self.pdf = density_values
        # Index id for the root cell
        self.root_id = root_id
        # Number of nearest neighbours
        self.k = neighbours_number

        # Variables for operations on the k-NN graph
        self.weighted_euc_knnG = None
        self.adjacency_matrix = None
        self.weighted_knnG = None
        self.weighted_graph_info = None
        self.graph_edges = None
        self.graph_info = None
        self.paths_edges_ids = None
        self.paths_nodes_ids = None
        self.kNN_Graph = None
        self.kNN_Graph_summary = None
        self.all_paths_edges = []
        self.pairwise_dist_mat = None
        self.pairs_euc = None

        # Variables for euclidean distance calculations
        self.euc_node_pairs_list = None
        self.euc_edges_weights = None
        self.euc_edges_weights_df = None
        self.euc_weighted_euc_knnG = None
        self.edges_indices = []
        self.edges_weights_df = None
        self.all_nodes_traversed = []
        self.all_edges_traversed = None

        # Variables for pseudotime computation
        self.pseudotime = []
        self.costs_df = None

    def graph_construction(self):

        """
        Construction of the k-Nearest Neighbour Graph and adjacency matrix

        Args:
            self.markers: markers to construct the graph from
            self.k: number of nearest neighbors to consider for each marker

        Returns:
            self.adjacency_matrix: adjacency matrix of the k-nearest neighbors graph
                                   in COOrdinate (COO) sparse format
        """

        knn_graph = kneighbors_graph(self.markers, self.k, mode='connectivity', include_self=False)

        # Save adjacency matrix in COOrdinate sparse format
        self.adjacency_matrix = knn_graph.tocoo()

        return self.adjacency_matrix

    def graph_density_weighting_scheme(self):

        """
        Implementation of graph weighting scheme: calculation of edge weights for a k-nearest neighbors graph based on
        the absolute log-probability difference of node pairs.

        Instead of representing a pair of nodes connected by an edge as (x_i, x_j) (manuscript version),
        here we represent the pair as (node_1, node_2).
        1. Convert the DataFrame 'self.pdf' to a NumPy array 'pdf_arr'
        2. Find all pairs of nodes connected by an edge using 'self.adjacency_matrix'
        3. Initialize zero weights for the edges
        4. Iterates through each edge and computes the absolute log-probability difference
           between its endpoints using 'pdf_arr'
        5. Construct a DataFrame 'self.weighted_knnG' containing 'node_1', 'node_2', and
           'edges_weights' columns, where 'edges_weights' represent the calculated weights

        Returns:
        self.weighted_knnG': output weighted k-NN graph

        """

        pdf_arr = self.pdf.to_numpy()

        # Find all pairs of nodes connected by an edge
        paired_nodes_df = pd.DataFrame({"node_1": self.adjacency_matrix.row, "node_2": self.adjacency_matrix.col})
        node_pairs_list = list(zip(paired_nodes_df.iloc[:, 0], paired_nodes_df.iloc[:, 1]))

        # Initialize zero weights for the edges
        edges_weights = np.zeros(paired_nodes_df.shape[0])

        # Calculate edge weights based on the absolute log-probability difference of node pairs
        for i in range(edges_weights.shape[0]):
            x = node_pairs_list[i][0]
            y = node_pairs_list[i][1]
            edges_weights[i] = np.abs(np.take(pdf_arr, x) - np.take(pdf_arr, y))
        edges_weights_df = pd.DataFrame(edges_weights, columns=["edges_weights"])

        # Return the weighted k-NN graph
        self.weighted_knnG = pd.concat([paired_nodes_df, edges_weights_df], axis=1)

        return self.weighted_knnG

    def cell_graph_info(self):

        """
        Graph information DataFrame with edge list for the weighted k-NN graph.

        Returns:
        self.graph_edges: a list of tuples representing edges of the weighted k-NN graph
        """

        self.graph_info = pd.concat([self.markers, self.pdf], axis=1)
        self.graph_info['node_1'] = self.graph_info.index
        self.weighted_graph_info = pd.merge(self.graph_info, self.weighted_knnG, on='node_1', how='outer')
        # List all edges of the graph
        self.graph_edges = list(zip(self.weighted_graph_info.node_1, self.weighted_graph_info.node_2))

        return self.graph_edges

    def graph_object(self):

        """
        Edge weight assingment on the k-NN graph.

        Args:
            self.graph_edges: list tuples representing edges between nodes
            self.weighted_graph_info: edge weights for the graph

        Returns:
            self.kNN_Graph: k-NN graph object
            self.kNN_Graph_summary: summary dataframe containing edge information of the graph
        """

        # Create as many nodes as data rows
        nodes = self.markers.shape[0]
        self.kNN_Graph = ig.Graph(nodes, self.graph_edges)
        g_edge_weights = self.weighted_graph_info.edges_weights.tolist()
        # Assign edge weights
        self.kNN_Graph.es["weights"] = g_edge_weights
        self.kNN_Graph_summary = self.kNN_Graph.get_edge_dataframe()
        # Return k-NN graph and dataframe with edge weights
        return self.kNN_Graph, self.kNN_Graph_summary

    def shortest_path_edges_traversed(self):

        """
           Set of edges traversed in each shortest path starting from a root node.

            Args:
                self.root_id: id of the root node from which shortest paths are computed
                self.kNN_Graph: k-NN graph used to compute shortest paths

            Returns:
                self.all_edges_traversed: a list of lists, where each inner list contains edges found in the shortest path
                                        from the root node to each target node (cell)
        """

        # Find the set of edges involved in each shortest path
        self.all_edges_traversed = self.kNN_Graph.get_shortest_paths(self.root_id, to=None, weights="weights",
                                                                     output='epath')

        return self.all_edges_traversed

    def shortest_path_nodes_traversed(self):

        """
        Set of nodes traversed in each shortest path starting from a root node.

        Args:
            self.root_id: id of the root node from which shortest paths are computed
            self.kNN_Graph: k-NN graph used to compute shortest paths

        Returns:
            self.all_nodes_traversed: a list of lists, where each set contains nodes found in
                                      the shortest path from the root node to each target node (cell)
        """

        # Find the set of nodes involved in each shortest path
        self.all_nodes_traversed = self.kNN_Graph.get_shortest_paths(self.root_id, to=None, weights="weights",
                                                                     output='vpath')
        return self.all_nodes_traversed

    def graph_construction_euclidean(self):

        """
        Construction of the k-Nearest Neighbour Graph and Euclidean distance pairwise matrix.

        Args:
            self.markers: markers to construct the graph from
            self.k: number of nearest neighbors to consider for each marker

        Returns:
            self.pairwise_dist_mat: the pairwise Euclidean distance matrix of the k-nearest neighbors graph
                                    in COOrdinate (COO) sparse format
        """

        knn_graph = kneighbors_graph(self.markers, self.k, mode='distance', include_self=False)

        # Save pairwise euclidean distances in COOrdinate sparse format
        self.pairwise_dist_mat = knn_graph.tocoo()

        return self.pairwise_dist_mat

    def euclidean_weights(self):

        """
        Calculation of edge weights for the k-NN graph based on Euclidean distances between node pairs.

        Instead of representing a pair of nodes connected by an edge as (x_i, x_j) (manuscript version),
        here we represent the pair as (node_a, node_b).
        1. Construct a DataFrame 'self.pairs_euc' with columns 'node_a' and 'node_b' from 'self.pairwise_dist_mat'
        2. Create a list 'self.euc_node_pairs_list' containing tuples of node pairs
        3. Initialize zero weights for the edges in 'self.euc_edges_weights'
        4. Iterates through each edge and computes the Euclidean distance between its endpoints using 'markersArr'
        5. Construct a DataFrame 'self.euc_weighted_euc_knnG' containing 'node_a', 'node_b',
           and 'euc_edges_weights' columns, where 'euc_edges_weights' represent the calculated weights

        Returns:
        self.euc_weighted_euc_knnG: weighted k-NN graph based on Euclidean distances
        """

        # self.pairwise_dist_mat_arr = self.pairwsise_dist_mat.toarray()
        self.pairs_euc = pd.DataFrame({"node_a": self.pairwise_dist_mat.row, "node_b": self.pairwise_dist_mat.col})
        self.euc_node_pairs_list = list(zip(self.pairs_euc.iloc[:, 0], self.pairs_euc.iloc[:, 1]))
        self.euc_edges_weights = np.zeros(self.pairs_euc.shape[0])
        markersArr = self.markers.to_numpy()

        for i in range(self.euc_edges_weights.shape[0]):
            x = self.euc_node_pairs_list[i][0]
            y = self.euc_node_pairs_list[i][1]
            # self.euc_edges_weights[i] = np.abs(euclidean_distances(markersArr[x].reshape(1, -1), markersArr[y].reshape(1,       -1))).reshape(1,)
            self.euc_edges_weights[i] = np.abs(
                np.linalg.norm(markersArr[x].reshape(1, -1) - markersArr[y].reshape(1, -1))).reshape(1, )

        self.euc_edges_weights_df = pd.DataFrame(self.euc_edges_weights, columns=["euc_edges_weights"])
        self.euc_weighted_euc_knnG = pd.concat([self.pairs_euc, self.euc_edges_weights_df], axis=1)

        return self.euc_weighted_euc_knnG

    def pseudotime_(self):

        """
        Pseudotime computation as evolution of markers from the root of the trajectory.

        Returns:
            pseudotime_costs: pseudotime value for each cell (not scaled in [0,1]
        """

        # edges_indices: list of edges found by iterating over 'all_edges_traversed'.
        # for each index 'i' in the range of the length of 'all_edges_traversed',
        # add the element at 'all_edges_traversed[i]' to 'edges_indices'.
        self.edges_indices = [self.all_edges_traversed[i] for i in range(len(self.all_edges_traversed))]

        # Sum the corresponding Euclidean costs of all edges across a path (edges of each path pre-defined based on
        # density weighting scheme)
        self.pseudotime_costs = [sum(self.euc_weighted_euc_knnG.iloc[j, 2]) for j in self.edges_indices]

        return self.pseudotime_costs



