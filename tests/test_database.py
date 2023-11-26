import os

from pandas import DataFrame
from ppi.database import Database
import pytest

HOME: str = os.path.expanduser("~")
PROJECT_FOLDER: str = os.path.join(HOME, ".ppi")
PATH_TO_DB: str = os.path.join(PROJECT_FOLDER, "ppi.sqlite")
PATH_TO_DATA = "./tests/data/test_ppi.tsv"

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


class TestDatabase:
    def setup_method(self, method):
        """Setup method to provide always a Database instance as self.db."""
        self.db = Database()
        self.db.set_path_to_data_file(path=PATH_TO_DATA)

    def test_folder_create(self):
        """Test if project folder (where database is stored) exists.

        The creation of the project folder should be implemented in Database.__init__
        """
        assert os.path.isdir(PROJECT_FOLDER)

    def test_set_not_existing_path(self):
        """Test if an non existing path is used for the data file."""
        with pytest.raises(FileNotFoundError):
            assert self.db.set_path_to_data_file(path="stupid_path")

    def test_read_data(self):
        """Test read the raw data and return DataFrame."""
        df_all: DataFrame = self.db.read_data()
        assert df_all.shape[0] == 8
        assert list(df_all.columns) == expected_columns

    def test_get_proteins_df(self):
        """Test id a protein dataframe is created."""
        df: DataFrame = self.db.get_proteins()
        assert isinstance(df, DataFrame)
        assert df.index.name == "id"
        assert df.to_dict(orient="tight") == expected_proteins

    def test_get_interaction(self):
        """Test if a interaction dataframe is created."""
        df_interactions: DataFrame = self.db.get_interactions()
        calculated_interactions = df_interactions.to_dict(orient="tight")
        assert calculated_interactions == expected_interactions

    def test_tables_exists(self):
        """Test if all tables in database exists."""
        self.db.import_data()
        tables: list[str] = self.db.get_table_names()
        assert set(tables) == set(["interaction", "protein"])

    def test_protein_table_column_names(self):
        """Test if protein table have all needed columns."""
        self.db.import_data()
        columns: list[str] = self.db.get_columns(table="protein")
        assert set(columns) == set(["id", "accession", "name", "taxid"])

    def test_interaction_table_column_names(self):
        """Test if interaction table have all needed columns."""
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
