import networkx as nx
from ppi.database import Database
from ppi.intact_analyzer import IntActAnalyzer

PATH_TO_DATA = "./tests/data/test_ppi.tsv"


class TestIntActAnalyzer:
    def setup_method(self, method):
        """Setup method to provide always a GraphAnalyzer instance as self.ga."""
        db = Database()
        db.set_path_to_data_file(path=PATH_TO_DATA)
        db.import_data()
        graph = db.get_graph()
        self.iaa: IntActAnalyzer = IntActAnalyzer(graph)

    def test_multi_graph(self):
        """Test the type."""
        assert isinstance(self.iaa.graph, nx.MultiGraph)

    def test_get_protein_with_highest_bc(self):
        """Test to get the protein with the highest betweenness centrality in th graph.

        Note: that additionally to the node data following keys are in the result:
        1. bc_value: calculated betweenness centrality value
        2. node_id: ID of node
        """
        assert self.iaa.get_protein_with_highest_bc() == {
            "accession": "node_id2",
            "name": "name_2",
            "taxid": 1,
            "node_id": 2,
            "bc_value": 0.8,
        }

    def test_get_neighbors_name(self):
        """Test getting names of neighbors."""
        assert self.iaa.get_neighbors_name("name_5") == ["name_2", "name_6"]
