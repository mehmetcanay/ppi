import networkx as nx


class IntActAnalyzer:
    def __init__(self, graph: nx.MultiGraph):
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
