from typing import Any
import click
from ppi import database, intact_analyzer
from networkx import MultiGraph


@click.group()
def main():
    pass


@main.command()
@click.option(
    "--path",
    "--always_required",
    required=True,
    type=str,
    help="A path to a datafile is required.",
)
def bcentrality(path: str):
    """Finds the protein with the highest betweenness score.

    Args:
        path (str): Path to the data file.

    Returns:
        dict[str, Any]: Dictionary of attributes of the protein with highest betweenness score.
    """
    db = database.Database()
    db.set_path_to_data_file(path=path)
    db.import_data()
    graph: MultiGraph | bool = db.get_graph()
    iaa = intact_analyzer.IntActAnalyzer(graph=graph)
    the_highest_bc: dict[str, Any] = iaa.get_protein_with_highest_bc()
    click.echo(message=the_highest_bc)


@main.command()
@click.option(
    "--path",
    "--always_required",
    required=True,
    help="A path to a datafile is required.",
)
def number_of_nodes(path: str):
    """Number of nodes in the graph.

    Args:
        path (str): Path to the data file.

    Returns:
        num_of_nodes (int): Number of nodes in the graph.
    """
    db = database.Database()
    db.set_path_to_data_file(path=path)
    db.import_data()
    graph: MultiGraph | bool = db.get_graph()
    iaa = intact_analyzer.IntActAnalyzer(graph=graph)
    num_of_nodes: int = iaa.number_of_nodes()
    click.echo(message=num_of_nodes)


if __name__ == "__main__":
    main()
