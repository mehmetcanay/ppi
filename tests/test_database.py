import os

from pandas import DataFrame
from ppi.database import Database
import pytest

HOME: str = os.path.expanduser("~")
PROJECT_FOLDER: str = os.path.join(HOME, ".ppi")
PATH_TO_DB: str = os.path.join(PROJECT_FOLDER, "ppi.sqlite")
PATH_TO_DATA = "./tests/data/test_ppi.tsv"
PATH_TO_EMPTY_DATA = "./tests/data/empty_ppi.tsv"

expected_columns: list[str] = [
    "confidence_value",
    "detection_method",
    "a_uniprot_id",
    "b_uniprot_id",
    "interaction_type",
    "pmid",
    "a_name",
    "a_taxid",
    "b_name",
    "b_taxid",
]

expected_proteins: dict = {
    "index": [1, 2, 3, 4, 5, 6, 7],
    "columns": ["accession", "name", "taxid"],
    "data": [
        ["node_id1", "name_1", 1],
        ["node_id2", "name_2", 1],
        ["node_id3", "name_3", 1],
        ["node_id4", "name_4", 1],
        ["node_id5", "name_5", 1],
        ["node_id6", "name_6", 2],
        ["node_id7", "name_7", 1],
    ],
    "index_names": ["id"],
    "column_names": [None],
}

expected_interactions: dict = {
    "index": [1, 2, 3, 4, 5, 6, 7, 8],
    "columns": [
        "confidence_value",
        "detection_method",
        "interaction_type",
        "pmid",
        "protein_a_id",
        "protein_b_id",
    ],
    "data": [
        [0.1, "dm1", "it1", "pmid1", 1, 2],
        [0.2, "dm2", "it2", "pmid1", 2, 3],
        [0.3, "dm3", "it2", "pmid1", 2, 4],
        [0.4, "dm1", "it2", "pmid1", 2, 5],
        [0.5, "dm4", "it2", "pmid1", 5, 6],
        [0.6, "dm1", "it3", "pmid1", 2, 6],
        [0.7, "dm5", "it3", "pmid1", 6, 7],
        [0.8, "dm2", "it3", "pmid1", 2, 3],
    ],
    "index_names": ["id"],
    "column_names": [None],
}

expected_graph_edges = {
    (1, 2, 0),
    (2, 3, 0),
    (2, 3, 1),
    (2, 4, 0),
    (2, 5, 0),
    (2, 6, 0),
    (5, 6, 0),
    (6, 7, 0),
}


class TestDatabase:
    def setup_method(self, method):
        """Setup method to provide always a Database instance as self.db."""
        self.db = Database()
        self.db.set_path_to_data_file(path=PATH_TO_DATA)

    def test_drop_database(self):
        """Test if database is dropped (SQLite file not exists)"""
        self.db.drop_database()
        assert os.path.isfile(PATH_TO_DB) == False

    def test_exists(self):
        """Test if database exists (SQLite file is created)."""
        self.db.import_data()
        assert self.db.exists == True

    def test_empty(self):
        """Test if interaction table has data."""
        self.db.import_data()
        assert self.db.has_data == True
        self.db.set_path_to_data_file(PATH_TO_EMPTY_DATA)
        self.db.import_data()
        assert self.db.has_data == False

    def test_folder_create(self):
        """Test if project folder (where database is stored) exists.

        The creation of the project folder should be implemented in Database.__init__
        """
        assert os.path.isdir(PROJECT_FOLDER)

    def test_set_not_existing_path(self):
        """Test if an no existing path is used for the data file."""
        with pytest.raises(FileNotFoundError):
            assert self.db.set_path_to_data_file(path="stupid_path")

    def test_read_data(self):
        """Test read the raw data and if return DataFrame has the expected columns."""
        df_all: DataFrame = self.db.read_data()
        assert df_all.shape[0] == 8
        assert list(df_all.columns) == expected_columns

    def test_get_proteins_df(self):
        """Test if a protein DataFrame is created."""
        df: DataFrame = self.db.get_proteins()
        assert isinstance(df, DataFrame)
        assert df.index.name == "id"
        assert df.to_dict(orient="tight") == expected_proteins

    def test_get_interaction(self):
        """Test if a interaction DataFrame is created."""
        df_interactions: DataFrame = self.db.get_interactions()
        calculated_interactions = df_interactions.to_dict(orient="tight")
        assert calculated_interactions == expected_interactions

    def test_tables_exists(self):
        """Test if all tables after import of data exists in database exists."""
        self.db.import_data()
        tables: list[str] = self.db.get_table_names()
        assert set(tables) == set(["interaction", "protein"])

    def test_protein_table_column_names(self):
        """Test if protein table has all needed columns."""
        self.db.import_data()
        columns: list[str] = self.db.get_columns(table="protein")
        assert set(columns) == set(["id", "accession", "name", "taxid"])

    def test_interaction_table_column_names(self):
        """Test if interaction table has all needed columns."""
        self.db.import_data()
        columns: list[str] = self.db.get_columns(table="interaction")
        expected_columns = set(
            [
                "id",
                "confidence_value",
                "detection_method",
                "interaction_type",
                "pmid",
                "protein_a_id",
                "protein_b_id",
            ]
        )
        assert set(columns) == expected_columns

    def test_detection_method_statistics(self):
        """Test distribution of detection_method in interactions."""
        df = self.db.get_detection_method_statistics()
        assert df.to_dict() == {
            "number": {"dm1": 3, "dm2": 2, "dm5": 1, "dm4": 1, "dm3": 1}
        }
        assert df.index.name == "detection_method"

    def test_pmid_statistics(self):
        """Test distribution of pmid in interactions."""
        df = self.db.get_pmid_statistics()
        assert df.to_dict() == {"number": {"pmid1": 8}}
        assert df.index.name == "pmid"

    def test_interaction_type_statistics(self):
        """Test distribution of interaction_type in interactions."""
        df = self.db.get_interaction_type_statistics()
        assert df.to_dict() == {"number": {"it2": 4, "it3": 3, "it1": 1}}
        assert df.index.name == "interaction_type"

    def test_confidence_value_statistics(self):
        """Test distribution of confidence_value in interactions."""
        df = self.db.get_confidence_value_statistics()
        assert df.to_dict() == {
            "number": {0.8: 1, 0.7: 1, 0.6: 1, 0.5: 1, 0.4: 1, 0.3: 1, 0.2: 1, 0.1: 1}
        }
        assert df.index.name == "confidence_value"

    def test_get_where_pmid(self):
        """Test if the WHERE SQL part for column pmid is correctly created."""
        where = " WHERE pmid = 'pmid1'"
        assert self.db.get_where(pmid="pmid1") == where

    def test_get_where_confidence_value(self):
        """Test if the WHERE SQL part for column confidence_value is correctly created."""
        where = " WHERE confidence_value >= 0.5"
        assert self.db.get_where(confidence_value_gte=0.5) == where

    def test_get_where_detection_method(self):
        """Test if the WHERE SQL part for column detection_method is correctly created."""
        where = " WHERE detection_method = 'my_method'"
        assert self.db.get_where(detection_method="my_method") == where

    def test_get_where_interaction_type(self):
        """Test if the WHERE SQL part for column interaction_type is correctly created."""
        where = " WHERE interaction_type = 'my_method'"
        assert self.db.get_where(interaction_type="my_method") == where

    def test_get_where_disallow_self_interaction(self):
        """Test if the WHERE SQL part for disallow_self_interaction is correctly created."""
        where = " WHERE protein_a_id != protein_b_id"
        assert self.db.get_where(disallow_self_interaction=True) == where

    def test_get_where_all(self):
        """Test if the WHERE SQL part correctly created."""
        where = " WHERE pmid = 'pmid' AND detection_method = 'dt' AND interaction_type = 'it' AND confidence_value >= 0.8 AND protein_a_id != protein_b_id"
        result = self.db.get_where(
            pmid="pmid",
            detection_method="dt",
            interaction_type="it",
            confidence_value_gte=0.8,
            disallow_self_interaction=True,
        )
        assert result == where

    def test_get_graph(self):
        """Tests if the graph created from a database query has the correct values."""
        g = self.db.get_graph()
        assert g.number_of_edges() == 8
        assert set(g.edges) == expected_graph_edges
        assert set(g.nodes) == {1, 2, 3, 4, 5, 6, 7}

    def test_data_in_node(self):
        """Test if data from the protein table are loaded into the nodes."""
        graph = self.db.get_graph()
        assert graph.nodes[1] == {"accession": "node_id1", "name": "name_1", "taxid": 1}

    def test_data_in_graphs(self):
        """Test if data from the interaction table are loaded into the edges."""
        graph = self.db.get_graph()
        # id in the result is equivalent to the id in the interaction table
        assert graph.edges[(1, 2, 0)] == {
            "id": 1,
            "confidence_value": 0.1,
            "pmid": "pmid1",
            "interaction_type": "it1",
            "detection_method": "dm1",
        }

    def test_filter_confidence_value(self):
        """Test the confidence_value_gte filter for the get_graph method.

        This filter gets all entries with a confidence_value >= confidence_value_gte
        """
        g = self.db.get_graph(confidence_value_gte=0.7)
        assert set(g.nodes) == {2, 3, 6, 7}
        assert set(g.edges) == {(2, 3, 0), (6, 7, 0)}

    def test_filter_detection_method(self):
        """Test the detection_method filter for the get_graph method."""
        g = self.db.get_graph(detection_method="dm3")
        assert set(g.nodes) == {2, 4}
        assert set(g.edges) == {(2, 4, 0)}

    def test_filter_interaction_type(self):
        """Test the interaction_type filter for the get_graph method."""
        g = self.db.get_graph(interaction_type="it1")
        assert set(g.nodes) == {1, 2}
        assert set(g.edges) == {(1, 2, 0)}

    def test_filter_pmid(self):
        """Test the pmid filter for the get_graph method."""
        g = self.db.get_graph(pmid="pmid1")
        assert g.number_of_nodes() == 7
        assert g.number_of_edges() == 8

    def test_filter_disallow_self_interaction(self):
        """Test the filter disallow self interaction.

        This can not been tested, because no self interaction are in the test dataset.
        """
        g = self.db.get_graph(disallow_self_interaction=True)
        assert g.number_of_edges() == 8
