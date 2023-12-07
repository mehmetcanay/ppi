import os
import errno
from typing import Any
import pandas as pd
from sqlalchemy import create_engine, inspect
import networkx as nx

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

        # Extract information on protein a and rename columns
        columns_a: list[str] = ["a_uniprot_id", "a_name", "a_taxid"]
        protein_a: pd.DataFrame = df.loc[:, columns_a]
        protein_a.rename(
            columns=column_names_a,
            inplace=True,
        )

        # Extract information on protein b and rename columns
        columns_b: list[str] = ["b_uniprot_id", "b_name", "b_taxid"]
        protein_b: pd.DataFrame = df.loc[:, columns_b]
        protein_b.rename(
            columns=column_names_b,
            inplace=True,
        )

        # Concatenate information of protein a and protein b
        proteins: pd.DataFrame = pd.concat([protein_a, protein_b])

        # Drop duplicated information
        proteins.drop_duplicates(inplace=True)

        # Sort the data by the accession number
        proteins.sort_values(by=["accession"], inplace=True)

        # Reset the index and drop the previous index column
        proteins.reset_index(inplace=True, drop=True)

        # Incriminate indices by 1 for consistency in SQL format
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
            # Extract information about unique proteins
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

        # Iterate over the detection methods and retrieve the occurrences of them. Store
        # the number of occurrences in statistics dictionary
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

        # Iterate over the pmids and retrieve the occurrences of them. Store the
        # number of occurrences in statistics dictionary
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

        # Iterate over the interaction types and retrieve the occurrences of
        # them. Store the number of occurrences in statistics dictionary
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

        # Iterate over the confidence values and retrieve the occurrences of
        # them. Store the number of occurrences in statistics dictionary
        for confidence_value in confidence_values:
            number: int = df.confidence_value.value_counts()[confidence_value]
            statistics["number"].append(number)

        # Create a dataframe and set confidence_value column as the index
        statistics_df = pd.DataFrame(statistics)
        statistics_df.set_index("confidence_value", inplace=True)

        return statistics_df

    def get_where(
        self,
        pmid: str = "",
        detection_method: str = "",
        interaction_type: str = "",
        confidence_value_gte: float = -1,
        disallow_self_interaction: bool = False,
    ) -> str:
        """Generates a WHERE SQL query with given input.

        Args:
            pmid (str, optional): Specified PMID. Defaults to "".

            confidence_value_gte (float, optional): Specified minimum
            confidence value. Defaults to -1.

            detection_method (str, optional): Specified detection method.
            Defaults to "".

            interaction_type (str, optional): Specified interaction type.
            Defaults to "".

            disallow_self_interaction (bool | None, optional): True if self
            interaction is allowed, False otherwise. Defaults to False.

        Returns:
            sql_query (str): WHERE SQL query with specified parameters.
        """

        sql_query: str = ""

        if pmid:
            if sql_query:
                sql_query += f" AND pmid = '{pmid}'"
            else:
                sql_query = f" WHERE pmid = '{pmid}'"

        if detection_method:
            if sql_query:
                sql_query += f" AND detection_method = '{detection_method}'"
            else:
                sql_query = f" WHERE detection_method = '{detection_method}'"

        if interaction_type:
            if sql_query:
                sql_query += f" AND interaction_type = '{interaction_type}'"
            else:
                sql_query = f" WHERE interaction_type = '{interaction_type}'"

        if confidence_value_gte >= 0:
            if sql_query:
                sql_query += f" AND confidence_value >= {confidence_value_gte}"
            else:
                sql_query = f" WHERE confidence_value >= {confidence_value_gte}"

        if disallow_self_interaction is True:
            if sql_query:
                sql_query += f" AND protein_a_id != protein_b_id"
            else:
                sql_query: str = f" WHERE protein_a_id != protein_b_id"

        return sql_query

    def get_graph(
        self,
        pmid: str = "",
        detection_method: str = "",
        interaction_type: str = "",
        confidence_value_gte: float = -1,
        disallow_self_interaction: bool = False,
    ) -> nx.MultiGraph | bool:
        if self.engine:
            # Initialize the MultiGraph instance
            G = nx.MultiGraph()

            query: str = self.get_where(
                pmid=pmid,
                detection_method=detection_method,
                interaction_type=interaction_type,
                confidence_value_gte=confidence_value_gte,
                disallow_self_interaction=disallow_self_interaction,
            )

            # Read the dataframes
            protein: pd.DataFrame = pd.read_sql("protein", self.engine, index_col="id")
            interaction: pd.DataFrame = pd.read_sql(
                f"SELECT * from interaction{query}", self.engine, index_col="id"
            )

            # Extracting protein ids existing in the interaction database
            proteins_a = set(interaction.protein_a_id.unique())
            proteins_b = set(interaction.protein_b_id.unique())
            proteins = proteins_a.union(proteins_b)
            proteins = list(sorted(proteins))

            # Filter the protein db based on the existing protein ids
            # in the interaction database
            protein = protein.loc[proteins]

            # Iterate over the indices of filtered protein dataframe
            # and add them as nodes with their corresponding attributes
            for node in protein.index:
                G.add_node(
                    node,
                    accession=protein.loc[node, "accession"],
                    name=protein.loc[node, "name"],
                    taxid=protein.loc[node, "taxid"],
                )

            # Iterate over the indices of interaction dataframe
            # and add them as edges with their corresponding attributes
            for i in interaction.index:
                G.add_edge(
                    u_for_edge=interaction.loc[i, "protein_a_id"],
                    v_for_edge=interaction.loc[i, "protein_b_id"],
                    id=i,
                    confidence_value=interaction.loc[i, "confidence_value"],
                    pmid=interaction.loc[i, "pmid"],
                    interaction_type=interaction.loc[i, "interaction_type"],
                    detection_method=interaction.loc[i, "detection_method"],
                )

            return G

        else:
            print("Please create a database first by using import_data() function")
            return False
