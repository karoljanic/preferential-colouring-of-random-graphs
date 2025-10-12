from typing import Union, Callable
import sympy
import networkx as nx

import expected_number_of_subgraphs

""" Contains functions to compute the expected number of subdivisions of specific subgraphs
    within a larger Barabási---Albert graph. Results are symbolic expressions in terms of
    the number of vertices (n) and the mean degree (m) of the graph. The expected number
    of subdivisions is in the notation O(*).
"""


def expected_number_of_subdivisions(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    subgraph_vertices_number: int,
    subgraph_edges_number: int,
    subgraph_validator: Callable[[nx.Graph], bool],
    subdivision_size: int
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of subdivisions of a specific subgraph within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        subgraph_vertices_number (int): Number of vertices in the subgraph.
        subgraph_edges_number (int): Number of edges in the subgraph.
        subgraph_validator (Callable[[nx.Graph], bool]): A function that takes a NetworkX graph and returns True if it matches the desired subgraph structure.

    Returns:
        sympy.Expr: The expected number of occurrences of the subgraph subdivisions - O(*).

    Note:
        All degrees in the subgraph must be at least 3.
    """


# Usage:
if __name__ == "__main__":
    expected_number_of_k5_subdivisions = expected_number_of_subdivisions(
        sympy.symbols("n"),
        sympy.symbols("m"),
        3,
        3,
        lambda g: len(g.nodes) == 3 and all(d == 2 for _, d in g.degree()),
    )
    print("Expected number of K5 subdivisions:", expected_number_of_k5_subdivisions)

    expected_number_of_k33_subdivisions = expected_number_of_subdivisions(
        sympy.symbols("n"),
        sympy.symbols("m"),
        6,
        9,
        lambda g: len(g.nodes) == 6
        and all(d == 3 for _, d in g.degree())
        and nx.is_bipartite(g),
    )
    print("Expected number of K3,3 subdivisions:", expected_number_of_k33_subdivisions)
