from typing import Union, Callable
import sympy
import networkx as nx

from expected_number_of_subgraphs import expected_number_of_subgraphs

""" Contains functions to compute the expected number of subdivisions of specific subgraphs
    within a larger Barabási---Albert graph. Results are symbolic expressions in terms of
    the number of vertices (n) and the mean degree (m) of the graph. The expected number
    of subdivisions is in the Theta asymptotic.
"""


def expected_number_of_subdivisions_with_size(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    subgraph_vertices_number: int,
    subgraph_edges_number: int,
    subgraph_validator: Callable[[nx.Graph], bool],
    subdivision_size: Union[sympy.Integer, sympy.Symbol],
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of subdivisions of a specific subgraph within a larger graph and given the size of the subdivisions.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        subgraph_vertices_number (int): Number of vertices in the subgraph.
        subgraph_edges_number (int): Number of edges in the subgraph.
        subgraph_validator (Callable[[nx.Graph], bool]): A function that takes a NetworkX graph and returns True if it matches the desired subgraph structure.
        subdivision_size (int): The size of the subdivisions to consider.

    Returns:
        sympy.Expr: The expected number of occurrences of the subgraph subdivisions - Theta asymptotic.

    Note:
        All degrees in the subgraph must be at least 3.
    """

    subdivisions_number = (
        sympy.binomial(subdivision_size, subgraph_vertices_number)
        * sympy.binomial(
            subdivision_size - subgraph_vertices_number + subgraph_edges_number - 1,
            subgraph_edges_number - 1,
        )
        * sympy.factorial(subdivision_size - subgraph_vertices_number)
    )

    degree_distributions_number = (
        graph_average_degree ** (subdivision_size - subgraph_vertices_number)
        + (
            graph_average_degree
            ** (
                sympy.Integer(2)
                / sympy.Integer(3)
                * (subdivision_size - subgraph_vertices_number)
            )
        )
        * (
            sympy.Integer(2)
            ** ((subdivision_size - subgraph_vertices_number) / sympy.Integer(3))
        )
    ) / (sympy.Integer(2) ** (subdivision_size - subgraph_vertices_number))

    crossings_number = sympy.E ** (
        subdivision_size**3
        * sympy.log(graph_vertices_number * graph_average_degree)
        / graph_vertices_number
        / graph_average_degree
    )

    index_choices = sympy.harmonic(graph_vertices_number) ** (
        subdivision_size - subgraph_vertices_number
    )

    return (
        subdivisions_number
        * degree_distributions_number
        * crossings_number
        * expected_number_of_subgraphs(
            graph_vertices_number,
            graph_average_degree,
            subgraph_vertices_number,
            subgraph_edges_number,
            subgraph_validator,
        )
        * index_choices
    )


def expected_number_of_subdivisions(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    subgraph_vertices_number: int,
    subgraph_edges_number: int,
    subgraph_validator: Callable[[nx.Graph], bool],
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
        sympy.Expr: The expected number of occurrences of the subgraph subdivisions - Theta asymptotic.

    Note:
        All degrees in the subgraph must be at least 3.
    """

    subdivision_size = sympy.symbols("k", integer=True, positive=True)

    expected_value = expected_number_of_subdivisions_with_size(
        graph_vertices_number,
        graph_average_degree,
        subgraph_vertices_number,
        subgraph_edges_number,
        subgraph_validator,
        subdivision_size,
    )

    return sympy.summation(
        expected_value,
        (subdivision_size, subgraph_vertices_number, graph_vertices_number),
    )


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
