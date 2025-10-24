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
        / 1000
    )

    index_choices = sympy.harmonic(graph_vertices_number) ** (
        subdivision_size - subgraph_vertices_number
    )

    print(
        float(subdivisions_number),
        float(degree_distributions_number),
        float(crossings_number),
        float(index_choices),
    )

    return (
        subdivisions_number
        / 1000
        * degree_distributions_number
        / 1000
        * crossings_number
        / 100
        * expected_number_of_subgraphs(
            graph_vertices_number,
            graph_average_degree,
            subgraph_vertices_number,
            subgraph_edges_number,
            subgraph_validator,
        )
        / 1000
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


def expected_number_of_clique_subdivisions_with_size(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    clique_size: int,
    subdivision_size: Union[sympy.Integer, sympy.Symbol],
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of subdivisions of a complete graph within a larger graph and given the size of the subdivisions.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        clique_size (int): Number of vertices in the complete graph.
        subdivision_size (Union[sympy.Integer, sympy.Symbol]): The size of the subdivisions to consider.

    Returns:
        sympy.Expr: The expected number of occurrences of the complete graph subdivisions - Theta asymptotic.
    """

    def clique_validator(g: nx.Graph) -> bool:
        return len(g.nodes) == clique_size and all(
            d == clique_size - 1 for _, d in g.degree()
        )

    subgraph_vertices_number = clique_size
    subgraph_edges_number = clique_size * (clique_size - 1) // 2

    return expected_number_of_subdivisions_with_size(
        graph_vertices_number,
        graph_average_degree,
        subgraph_vertices_number,
        subgraph_edges_number,
        clique_validator,
        subdivision_size,
    )


def expected_number_of_complete_bipartite_subdivisions_with_size(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    part_size_a: int,
    part_size_b: int,
    subdivision_size: Union[sympy.Integer, sympy.Symbol],
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of subdivisions of a complete bipartite graph within a larger graph and given the size of the subdivisions.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        part_size_a (int): Number of vertices in the first part of the bipartite graph.
        part_size_b (int): Number of vertices in the second part of the bipartite graph.
        subdivision_size (Union[sympy.Integer, sympy.Symbol]): The size of the subdivisions to consider.

    Returns:
        sympy.Expr: The expected number of occurrences of the complete bipartite graph subdivisions - Theta asymptotic.
    """

    def bipartite_validator(g: nx.Graph) -> bool:
        return (
            len(g.nodes) == part_size_a + part_size_b
            and all(d == part_size_b for _, d in g.degree() if d >= part_size_b)
            and all(d == part_size_a for _, d in g.degree() if d >= part_size_a)
            and nx.is_bipartite(g)
        )

    subgraph_vertices_number = part_size_a + part_size_b
    subgraph_edges_number = part_size_a * part_size_b

    return expected_number_of_subdivisions_with_size(
        graph_vertices_number,
        graph_average_degree,
        subgraph_vertices_number,
        subgraph_edges_number,
        bipartite_validator,
        subdivision_size,
    )


def expected_number_of_clique_subdivisions(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    clique_size: int,
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of subdivisions of a complete graph within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        clique_size (int): Number of vertices in the complete graph.

    Returns:
        sympy.Expr: The expected number of occurrences of the complete graph subdivisions - Theta asymptotic.
    """

    def clique_validator(g: nx.Graph) -> bool:
        return len(g.nodes) == clique_size and all(
            d == clique_size - 1 for _, d in g.degree()
        )

    subgraph_vertices_number = clique_size
    subgraph_edges_number = clique_size * (clique_size - 1) // 2

    return expected_number_of_subdivisions(
        graph_vertices_number,
        graph_average_degree,
        subgraph_vertices_number,
        subgraph_edges_number,
        clique_validator,
    )


def expected_number_of_complete_bipartite_subdivisions(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    part_size_a: int,
    part_size_b: int,
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of subdivisions of a complete bipartite graph within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        part_size_a (int): Number of vertices in the first part of the bipartite graph.
        part_size_b (int): Number of vertices in the second part of the bipartite graph.

    Returns:
        sympy.Expr: The expected number of occurrences of the complete bipartite graph subdivisions - Theta asymptotic.
    """

    def bipartite_validator(g: nx.Graph) -> bool:
        return (
            len(g.nodes) == part_size_a + part_size_b
            and all(d == part_size_b for _, d in g.degree() if d >= part_size_b)
            and all(d == part_size_a for _, d in g.degree() if d >= part_size_a)
            and nx.is_bipartite(g)
        )

    subgraph_vertices_number = part_size_a + part_size_b
    subgraph_edges_number = part_size_a * part_size_b

    return expected_number_of_subdivisions(
        graph_vertices_number,
        graph_average_degree,
        subgraph_vertices_number,
        subgraph_edges_number,
        bipartite_validator,
    )


# Usage:
if __name__ == "__main__":
    expected_number_of_k5_subdivisions_with_k_size = (
        expected_number_of_clique_subdivisions_with_size(
            sympy.symbols("n"),
            sympy.symbols("m"),
            5,
            sympy.symbols("k"),
        )
    )
    print(
        "Expected number of K5 subdivisions with k size:",
        expected_number_of_k5_subdivisions_with_k_size,
    )

    expected_number_of_k33_subdivisions_with_k_size = (
        expected_number_of_complete_bipartite_subdivisions_with_size(
            sympy.symbols("n"),
            sympy.symbols("m"),
            3,
            3,
            sympy.symbols("k"),
        )
    )
    print(
        "Expected number of K3,3 subdivisions with k size:",
        expected_number_of_k33_subdivisions_with_k_size,
    )

    # expected_number_of_k5_subdivisions = expected_number_of_clique_subdivisions(
    #     sympy.symbols("n"),
    #     sympy.symbols("m"),
    #     5,
    # )
    # print("Expected number of K5 subdivisions:", expected_number_of_k5_subdivisions)

    # expected_number_of_k33_subdivisions = (
    #     expected_number_of_complete_bipartite_subdivisions(
    #         sympy.symbols("n"),
    #         sympy.symbols("m"),
    #         3,
    #         3,
    #     )
    # )
    # print("Expected number of K3,3 subdivisions:", expected_number_of_k33_subdivisions)
