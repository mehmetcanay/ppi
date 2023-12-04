import os
import errno
from typing import Any
import pandas as pd
from sqlalchemy import create_engine, inspect

HOME: str = os.path.expanduser("~")
PROJECT_FOLDER: str = os.path.join(HOME, ".ppi")
PATH_TO_DB: str = os.path.join(PROJECT_FOLDER, "ppi.sqlite")


class Database:
    def __init__(self, path: str = ""):
        """Initializes an instance of Database class with default
        values if no values are passed. Creates a PROJECT_FOLDER in
        HOME named ".ppi" if not already exists.

        Args:
            path (str, optional): Path to the .tsv file containing
            data. Defaults to "".
        """

        self.path: str = path

        # If the folder does not exists, create one
        if not os.path.isdir(PROJECT_FOLDER):
            os.mkdir(PROJECT_FOLDER)

        self.engine = None
        self.exists: bool = False
        self.has_data: bool = False

    def drop_database(self) -> bool:
        """Drops the existing database.

        Raises:
            FileNotFoundError: There is no existing database.

        Returns:
            bool: True if database exists and dropped.
        """

        # If the database exists, drop it.
        if os.path.isfile(PATH_TO_DB):
            os.remove(PATH_TO_DB)
            self.exists = False
            return True

        # Raise an error if there is no database
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), PATH_TO_DB)

    def set_path_to_data_file(self, path: str) -> bool:
        """Changes self.path value if the given path to the data file is valid.

        Args:
            path (str): A path to a data file.

        Raises:
            FileNotFoundError: Given path is not valid.

        Returns:
            bool: True, if the path is valid.
        """

        # Set self.path to path if the given path is valid,
        # raise an error otherwise
        if os.path.isfile(path):
            self.path = path
            return True
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)

    def read_data(self) -> pd.DataFrame:
        """Reads the data from the .tsv file in the given path.

        Returns:
            df (pd.DataFrame): Pandas DataFrame containing information
            from the file in the given path.
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

        # Column rename map of protein a
        column_names_a: dict[str, str] = {
            "a_uniprot_id": "accession",
            "a_name": "name",
            "a_taxid": "taxid",
        }

        # Column rename map of protein b
        column_names_b: dict[str, str] = {
            "b_uniprot_id": "accession",
            "b_name": "name",
            "b_taxid": "taxid",
        }

        # Extract informations on protein a and rename columns
        columns_a: list[str] = ["a_uniprot_id", "a_name", "a_taxid"]
        protein_a: pd.DataFrame = df.loc[:, columns_a]
        protein_a.rename(
            columns=column_names_a,
            inplace=True,
        )

        # Extract informations on protein b and rename columns
        columns_b: list[str] = ["b_uniprot_id", "b_name", "b_taxid"]
        protein_b: pd.DataFrame = df.loc[:, columns_b]
        protein_b.rename(
            columns=column_names_b,
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

        return proteins

    def get_interactions(self) -> pd.DataFrame:
        """Generates an interactions table from the data in the given path.

        Returns:
            interactions (pd.DataFrame): A dataframe containing information on
            unique interaction entries in the data file in the given path.
        """
        # Read the data file in the given path
        df: pd.DataFrame = self.read_data()

        # Check if the dataframe is empty
        if df.empty:
            # Columns to be dropped
            columns: list[str] = ["a_name", "a_taxid", "b_name", "b_taxid"]

            # Column names to be renamed
            rename_dict: dict[str, str] = {
                "a_uniprot_id": "protein_a_id",
                "b_uniprot_id": "protein_b_id",
            }

            # Drop the unnecessary columns
            interactions = df.drop(columns, axis=1)

            # Rename the columns
            interactions.rename(rename_dict, inplace=True)

            # Increment the indices by 1 and rename it as
            # id for consistency in SQL format
            interactions.index += 1
            interactions.index.name = "id"

            return interactions

        else:
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
            accession = proteins["accession"]
            interactions["protein_a_id"] = interactions.apply(
                lambda x: proteins.loc[accession == x["a_uniprot_id"]].index[0],
                axis=1,
            )
            interactions["protein_b_id"] = interactions.apply(
                lambda x: proteins.loc[accession == x["b_uniprot_id"]].index[0],
                axis=1,
            )

            # Drop the columns that are no longer needed
            interactions.drop(["a_uniprot_id", "b_uniprot_id"], axis=1, inplace=True)

            # Reset the index
            interactions.reset_index(inplace=True, drop=True)

            # Increment the indices by 1 and rename it as
            # id for consistency in SQL format
            interactions.index += 1
            interactions.index.name = "id"

            return interactions

    def import_data(self) -> None:
        """Exports protein and interaction tables created by get_proteins() and
        get_interactions() functions into the SQL database.
        """

        # Creating an engine and assigning it to engine attribute
        self.engine = create_engine("sqlite:///" + PATH_TO_DB)

        # Changing exists attribute to True
        self.exists = True

        # Getting the tables
        protein: pd.DataFrame = self.get_proteins()
        interaction: pd.DataFrame = self.get_interactions()

        # If the interaction table is empty, change has_data
        # attribute to False, else to True.
        if interaction.empty:
            self.has_data = False

        else:
            self.has_data = True

        # Exporting the tables to an SQL database
        protein.to_sql(name="protein", con=self.engine, if_exists="replace")
        interaction.to_sql(name="interaction", con=self.engine, if_exists="replace")

    def get_table_names(self) -> list[str]:
        """Retrieves table names in the SQL database.

        Returns:
            table_names (list[str]): A list containing the name of
            the tables in the SQL database.
        """

        # Check if the database is created
        assert self.engine is not None, "Database is not created"

        # Inspect the connection of the connected engine
        # to retrieve the table names
        with self.engine.connect() as connection:
            inspector = inspect(connection)
            table_names: list[str] = inspector.get_table_names()

        return table_names

    def get_columns(self, table: str) -> list[str]:
        """Returns the list of column names in a given table
        in the SQL database.

        Args:
            table (str): Name of the table.

        Returns:
            columns (list[str]): A list containing names
            of the columns in a table.
        """

        # Check if the database is created
        assert self.engine is not None, "Database is not created"

        # Inspect the connection of the connected engine
        # to retrieve the column names
        with self.engine.connect() as connection:
            inspector = inspect(connection)

            # Get columns from a specified table
            columns = inspector.get_columns(table)

            # Retrieve name information from ReflectedColumn object
            columns = [column["name"] for column in columns]

        return columns

    def get_detection_method_statistics(self) -> pd.DataFrame:
        """Generates a table containing number of used detection methods to
        detect interactions.

        Returns:
            statistics_df (pd.DataFrame): A dataframe containing detection
            method statistics.
        """

        # Get the interactions table
        df: pd.DataFrame = self.get_interactions()

        # Retrieve the unique detection methods
        methods: list[str] = list(df.detection_method.unique())

        # Initialize a dictionary to hold statistic values
        statistics: dict[str, list[Any]] = {"detection_method": methods, "number": []}

        # Iterate over the detection methods and retrieve the occurences of them. Store
        # the number of occurences in statistics dictionary
        for method in methods:
            number: int = df.detection_method.value_counts()[method]
            statistics["number"].append(number)

        # Create a dataframe and set detection_method column as the index
        statistics_df = pd.DataFrame(statistics)
        statistics_df.set_index("detection_method", inplace=True)

        return statistics_df

    def get_pmid_statistics(self) -> pd.DataFrame:
        """Generates a table of pmid statistics.

        Returns:
            statistics_df (pd.DataFrame): A dataframe containing pmid
            statistics.
        """

        # Get the interactions table
        df: pd.DataFrame = self.get_interactions()

        # Retrieve the unique PMIDs
        pmids: list[str] = list(df.pmid.unique())

        # Initialize a dictionary to hold statistic values
        statistics: dict[str, list[Any]] = {"pmid": pmids, "number": []}

        # Iterate over the pmids and retrieve the occurences of them. Store the
        # number of occurences in statistics dictionary
        for pmid in pmids:
            number: int = df.pmid.value_counts()[pmid]
            statistics["number"].append(number)

        # Create a dataframe and set pmid column as the index
        statistics_df = pd.DataFrame(statistics)
        statistics_df.set_index("pmid", inplace=True)

        return statistics_df

    def get_interaction_type_statistics(self) -> pd.DataFrame:
        """Generates a table of interaction type statistics.

        Returns:
            statistics_df (pd.DataFrame): A dataframe containing interaction
            type statistics.
        """

        # Get the interactions table
        df: pd.DataFrame = self.get_interactions()

        # Retrieve the unique interaction types
        interaction_types: list[str] = list(df.interaction_type.unique())

        # Initialize a dictionary to hold statistic values
        statistics: dict[str, list[Any]] = {
            "interaction_type": interaction_types,
            "number": [],
        }

        # Iterate over the interaction types and retrieve the occurences of
        # them. Store the number of occurences in statistics dictionary
        for interaction_type in interaction_types:
            number: int = df.interaction_type.value_counts()[interaction_type]
            statistics["number"].append(number)

        # Create a dataframe and set interaction_type column as the index
        statistics_df = pd.DataFrame(statistics)
        statistics_df.set_index("interaction_type", inplace=True)

        return statistics_df
    
    def get_confidence_value_statistics(self) -> pd.DataFrame:
        """Generates a table of confidence value statistics.

        Returns:
            statistics_df (pd.DataFrame): A dataframe containing confidence
            value statistics.
        """

        # Get the interactions table
        df: pd.DataFrame = self.get_interactions()

        # Retrieve the unique interaction types
        confidence_values: list[str] = list(df.confidence_value.unique())

        # Initialize a dictionary to hold statistic values
        statistics: dict[str, list[Any]] = {
            "confidence_value": confidence_values,
            "number": [],
        }

        # Iterate over the confidence values and retrieve the occurences of
        # them. Store the number of occurences in statistics dictionary
        for confidence_value in confidence_values:
            number: int = df.confidence_value.value_counts()[confidence_value]
            statistics["number"].append(number)

        # Create a dataframe and set confidence_value column as the index
        statistics_df = pd.DataFrame(statistics)
        statistics_df.set_index("confidence_value", inplace=True)

        return statistics_df
