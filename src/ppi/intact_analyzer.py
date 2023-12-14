"""This is a module to graphically represent protein-protein interactions"""


import networkx as nx
import matplotlib.pyplot as plt


class IntActAnalyzer:
    """This class represents a protein-protein interaction graph analyzer."""

    def __init__(self, graph: nx.MultiGraph):
        """Initializes an instance of IntActAnalyzer class.

        Args:
            graph (nx.MultiGraph): A networkx MultiGraph.
        """
        self.graph: nx.MultiGraph = graph

    def get_protein_with_highest_bc(self):
        """Finds the protein with the highest betweenness score.

        Returns:
            node: The node with the maximum betweeness score.
        """
        # Initialize betweenness centrality value as 0
        max_bc = 0

        # Initialize the id of the node with the maximum betweenness score
        max_bc_node = None

        # Calculate the betweeness score of all the nodes
        bc_dict = nx.betweenness_centrality(self.graph)

        for node, bc in bc_dict.items():
            self.graph.nodes[node]["node_id"] = node
            self.graph.nodes[node]["bc_value"] = bc

            # If the betweennes score is greater than the maximum score,
            # assign the new score to max_bc
            if bc > max_bc:
                max_bc: float = bc
                max_bc_node = node

        node = self.graph.nodes[max_bc_node]
        return node

    def get_neighbors_name(self, name: str) -> list[str] | bool:
        """Get names of the neighbors of the node with the specified name.

        Args:
            name (str): Name of the node

        Returns:
            list[str] | bool: List of names of the neighboring nodes if
            a node with the specified name is found. False otherwise.
        """

        node_name: str = name

        # Initialize node_id as None
        node_id = None

        # Find the nodes with the specified name
        for node, attrs in self.graph.nodes(data=True):
            if attrs.get("name") == node_name:
                # Assign node_id when a match is found
                node_id = node
                break

        # Check if a node with the specified name was found
        if node_id is not None:
            neighbor_ids = self.graph.neighbors(node_id)
            neighbor_names: list[str] = [
                self.graph.nodes[n]["name"]
                for n in neighbor_ids  # Access node attributes directly
            ]
            return neighbor_names

        else:
            print("No node with the specified name found")
            return False

    def draw_graph(self, edge_label="id", node_label="id", figsize=(10, 5)):
        """Shows the graph.

        Arguments `edge_label` and `node_label` allows to change the labels in the graph.

        node_label: Any of the keys in node data
        edge_label: Any of the keys in edge data

        Args:
            edge_label (str, optional): Label to be shown on edges. Defaults to "id".
            node_label (str, optional): Label to be shown on nodes. Defaults to "id".
            figsize (tuple, optional): Size of the graph. Defaults to (10, 5).
        """
        plt.figure(figsize=figsize)
        if node_label == "id":
            node_labels = {x: x for x in self.graph.nodes}
        else:
            node_labels = {x: self.graph.nodes[x][node_label] for x in self.graph.nodes}
        edge_labels = {}
        for edge in self.graph.edges:
            eid = str(self.graph.edges[edge][edge_label])
            if edge[:2] not in edge_labels:
                edge_labels[edge[:2]] = eid
            else:
                edge_labels[edge[:2]] += f",{eid}"
        pos = nx.spring_layout(self.graph)  # Layout for the graph
        nx.draw_networkx_nodes(self.graph, pos)
        nx.draw_networkx_edges(self.graph, pos)
        nx.draw_networkx_labels(self.graph, pos, labels=node_labels)
        nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_labels)

        plt.show()

    def number_of_nodes(self) -> int:
        """Number of nodes in the graph.

        Returns:
            number_of_nodes (int): Number of nodes in the graph.
        """

        number_of_nodes: int = self.graph.number_of_nodes()
        return number_of_nodes
