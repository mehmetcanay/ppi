import os
import errno
import pandas as pd
from sqlalchemy import create_engine


class Database:
    def __init__(
        self,
        path: str = "",
        tables={"protein": pd.DataFrame(), "interaction": pd.DataFrame()},
    ) -> None:
        """Initializes an instance of Database class with default
        values if no values are passed. Creates a PROJECT_FOLDER in
        HOME named ".ppi" if not already exists. Creates an SQL database
        named ppi.sqlite in the PROJECT_FOLDER. If such a database already
        exists, it will be replaced.

        Args:
            path (str, optional): Path to the .tsv file containing data. Defaults to "".
            tables (dict, optional): A dictionary of dataframes. 
            Defaults to {"protein": pd.DataFrame(), "interaction": pd.DataFrame()}.
        """

        self.path: str = path
        self.tables: dict[str, pd.DataFrame] = tables

        HOME: str = os.path.expanduser("~")
        PROJECT_FOLDER: str = os.path.join(HOME, ".ppi")

        # If the folder does not exists, create one
        if not os.path.isdir(PROJECT_FOLDER):
            os.mkdir(PROJECT_FOLDER)

        DB_PATH: str = os.path.join(PROJECT_FOLDER, "ppi.sqlite")
        self.engine = create_engine("sqlite:///" + DB_PATH)

    def set_path_to_data_file(self, path: str) -> None:
        """Changes self.path value if the given path to the data file is valid.

        Args:
            path (str): A path to a data file

        Raises:
            FileNotFoundError: Given path is not valid
        """

        # Set self.path to path if the given path is valid,
        # raise an error otherwise
        if os.path.isfile(path):
            self.path = path
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)

    def read_data(self) -> pd.DataFrame:
        """Reads the data from the .tsv file in the given path.

        Returns:
            df (pd.DataFrame): Pandas DataFrame containing information from the file in the given path
        """

        df: pd.DataFrame = pd.read_csv(self.path, sep="\t")
        return df

    def get_proteins(self) -> pd.DataFrame:
        """Creates a proteins dataframe from the data file in the given path.

        Returns:
            proteins (pd.DataFrame): A dataframe containing information
            of unique proteins in the data file in the given path.
        """

        # Read the data in the given path
        df: pd.DataFrame = self.read_data()

        # Extract informations on protein a and rename columns
        protein_a: pd.DataFrame = df.loc[:, ["a_uniprot_id", "a_name", "a_taxid"]]
        protein_a.rename(
            columns={"a_uniprot_id": "accession", "a_name": "name", "a_taxid": "taxid"},
            inplace=True,
        )

        # Extract informations on protein b and rename columns
        protein_b: pd.DataFrame = df.loc[:, ["b_uniprot_id", "b_name", "b_taxid"]]
        protein_b.rename(
            columns={"b_uniprot_id": "accession", "b_name": "name", "b_taxid": "taxid"},
            inplace=True,
        )

        # Concatane informations of protein a and protein b
        proteins: pd.DataFrame = pd.concat([protein_a, protein_b])

        # Drop duplicated informations
        proteins.drop_duplicates(inplace=True)

        # Sort the data by the accession number
        proteins.sort_values(by=["accession"], inplace=True)

        # Reset the index and drop the previous index column
        proteins.reset_index(inplace=True, drop=True)

        # Increminate indices by 1 for consistency in SQL format
        proteins.index += 1

        # Rename index to id for consistency in SQL format
        proteins.index.name = "id"

        proteins.to_sql(name="protein", con=self.engine, if_exists="replace")

        return proteins

    def get_interactions(self) -> pd.DataFrame:
        """Generates an interactions table from the data in the given path.

        Returns:
            interactions (pd.DataFrame): A dataframe containing information on
            unique interaction entries in the data file in the given path.
        """
        # Read the data file in the given path
        df: pd.DataFrame = self.read_data()

        # Extract informations about unique proteins
        proteins: pd.DataFrame = self.get_proteins()

        # Columns to be used in filtration
        columns: list[str] = [
            "confidence_value",
            "detection_method",
            "interaction_type",
            "pmid",
            "a_uniprot_id",
            "b_uniprot_id",
        ]

        # Filter the data file
        interactions: pd.DataFrame = df.loc[:, columns]

        # Drop the duplicated entries
        interactions.drop_duplicates(inplace=True)

        # Change the uniprot_ids with ids from proteins dataframe
        interactions["protein_a_id"] = interactions.apply(
            lambda x: proteins.loc[proteins["accession"] == x["a_uniprot_id"]].index[0],
            axis=1,
        )
        interactions["protein_b_id"] = interactions.apply(
            lambda x: proteins.loc[proteins["accession"] == x["b_uniprot_id"]].index[0],
            axis=1,
        )

        # Drop the columns that are no longer needed
        interactions.drop(["a_uniprot_id", "b_uniprot_id"], axis=1, inplace=True)

        # Reset the index
        interactions.reset_index(inplace=True, drop=True)

        # Increment the indices by 1 and rename it as id for consistency in SQL format
        interactions.index += 1
        interactions.index.name = "id"

        interactions.to_sql(name="interaction", con=self.engine, if_exists="replace")

        return interactions

    def import_data(self) -> None:
        """Imports SQL database as Pandas DataFrame and stores it in
        self.tables dictionary. The SQL database must contain two
        tables named "protein" and "interaction".
        """

        # Reading tables from the SQL database
        protein: pd.DataFrame = pd.read_sql("protein", con=self.engine)
        interaction: pd.DataFrame = pd.read_sql("interaction", con=self.engine)

        # Storing them in self.tables dictionary of dataframes.
        self.tables["protein"] = protein
        self.tables["interaction"] = interaction

    def get_table_names(self) -> list[str]:
        """Retrieves table names in the SQL database.

        Returns:
            table_names (list[str]): A list containing the name of
            the tables in the SQL database.
        """

        table_names: list[str] = list(self.tables.keys())
        return table_names

    def get_columns(self, table: str) -> list[str]:
        """Returns the list of column names in Pandas
        DataFrame generated from an SQL database.

        Args:
            table (str): Name of the table

        Returns:
            columns (list[str]): A list containing names of the columns in a dataframe
        """
        df: pd.DataFrame = self.tables[table]
        columns: list[str] = list(df.columns)

        return columns
