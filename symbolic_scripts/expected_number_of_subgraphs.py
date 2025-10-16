from typing import Union, List, Callable, Iterator
import sympy
import networkx as nx
import itertools

""" Contains functions to calculate the expected number of occurrences of specific subgraphs 
    within a larger Barabási---Albert graph. Results are symbolic expressions in terms of
    the number of vertices (n) and the mean degree (m) of the graph and with assymptotic accuracy (1 + o(1)).
    
    Based on the method described in the paper:
    https://www.stat.berkeley.edu/~aldous/Networks/boll1.pdf """


def generate_labeled_graphs(
    vertices_number: int, edges_number: int
) -> Iterator[nx.Graph]:
    """
    Generates all labeled graphs with a given number of vertices and edges.

    Parameters:
        vertices_number (int): Number of vertices in the graph.
        edges_number (int): Number of edges in the graph.

    Yields:
        Iterator[nx.Graph]: An iterator over all labeled graphs with the specified vertices and edges.
    """

    vertices = list(range(1, vertices_number + 1))
    possible_edges = list(itertools.combinations(vertices, 2))

    for edges in itertools.combinations(possible_edges, edges_number):
        G = nx.Graph()
        G.add_nodes_from(vertices)
        G.add_edges_from(edges)
        yield G


def integer_partitions(n: int) -> Iterator[List[int]]:
    """
    Generate all integer partitions.

    Parameters:
        n (int): The integer to partition.

    Yields:
        Iterator[List[int]]: An iterator over all partitions of n.
    """

    dp = [[] for _ in range(n + 1)]
    dp[0] = [[]]
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            for p in dp[i - j]:
                if not p or j <= p[-1]:
                    dp[i].append(p + [j])

    for partition in dp[n]:
        yield partition


def expected_number_of_subgraphs(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    subgraph_vertices_number: int,
    subgraph_edges_number: int,
    subgraph_validator: Callable[[nx.Graph], bool],
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of a specific subgraph within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        subgraph_vertices_number (int): Number of vertices in the subgraph.
        subgraph_edges_number (int): Number of edges in the subgraph.
        subgraph_validator (Callable[[nx.Graph], bool]): A function that takes a NetworkX graph and returns True if it matches the desired subgraph structure.

    Returns:
        sympy.Expr: The expected number of occurrences of the subgraph - (1 + o(1)).
    """

    total_probability = sympy.Integer(0)

    for subgraph in generate_labeled_graphs(
        subgraph_vertices_number, subgraph_edges_number
    ):
        if not subgraph_validator(subgraph):
            continue

        out_degrees = {
            node: sum(1 for u, v in subgraph.edges if max(u, v) == node)
            for node in subgraph.nodes
        }

        in_degrees = {
            node: sum(1 for u, v in subgraph.edges if min(u, v) == node)
            for node in subgraph.nodes
        }

        out_edges_distribution = {
            node: sympy.prod(
                [graph_average_degree - i for i in range(out_degrees[node])],
                start=sympy.Integer(1),
            )
            for node in subgraph.nodes
        }

        in_edges_distribution = {node: {} for node in subgraph.nodes}
        for node in subgraph.nodes:
            for degree_partition in integer_partitions(in_degrees[node]):
                degree_counts = {
                    degree: degree_partition.count(degree)
                    for degree in set(degree_partition)
                }

                in_edges_distribution[node][tuple(degree_partition)] = (
                    sympy.prod(
                        [graph_average_degree - i for i in range(len(degree_partition))]
                    )
                    * sympy.factorial(in_degrees[node])
                    / sympy.prod([sympy.factorial(r) for r in degree_partition])
                    / sympy.prod([sympy.factorial(c) for _, c in degree_counts.items()])
                )

        in_edges_distribution_list = [
            list(v.items()) for k, v in in_edges_distribution.items()
        ]
        out_edges_distribution_list = list(out_edges_distribution.values())

        for combination in itertools.product(*in_edges_distribution_list):
            single_probability = sympy.prod(out_edges_distribution_list)
            for term in combination:
                degs, count = term
                single_probability *= count * sympy.prod(
                    [sympy.factorial(d) for d in degs]
                )
            single_probability /= 2**subgraph_edges_number
            single_probability /= graph_average_degree**subgraph_edges_number

            total_probability += single_probability

    return (
        (sympy.Integer(1) / sympy.factorial(subgraph_vertices_number))
        ** (
            sympy.Integer(subgraph_edges_number)
            / sympy.Integer(subgraph_vertices_number)
        )
        * total_probability.factor()
        * sympy.harmonic(
            graph_vertices_number,
            sympy.Integer(subgraph_edges_number)
            / sympy.Integer(subgraph_vertices_number),
        )
        ** sympy.Integer(subgraph_vertices_number)
    )


def expected_number_of_cycle_subgraphs(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    cycle_length: int,
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of cycle subgraphs of a given length within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        cycle_length (int): Length of the cycle subgraph.

    Returns:
        sympy.Expr: The expected number of occurrences of the cycle subgraph - (1 + o(1)).
    """

    return expected_number_of_subgraphs(
        graph_vertices_number,
        graph_average_degree,
        cycle_length,
        cycle_length,
        lambda g: len(g.nodes) == cycle_length and all(d == 2 for _, d in g.degree()),
    )


def expected_number_of_clique_subgraphs(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    clique_size: int,
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of clique subgraphs of a given size within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        clique_size (int): Size of the clique subgraph.

    Returns:
        sympy.Expr: The expected number of occurrences of the clique subgraph - (1 + o(1)).
    """

    return expected_number_of_subgraphs(
        graph_vertices_number,
        graph_average_degree,
        clique_size,
        clique_size * (clique_size - 1) // 2,
        lambda g: len(g.nodes) == clique_size
        and all(d == clique_size - 1 for _, d in g.degree()),
    )


def expected_number_of_bipartite_complete_subgraphs(
    graph_vertices_number: Union[sympy.Integer, sympy.Symbol],
    graph_average_degree: Union[sympy.Integer, sympy.Symbol],
    part_size: int,
) -> sympy.Expr:
    """
    Calculates the expected number of occurrences of complete bipartite subgraphs K_{part_size, part_size} within a larger graph.

    Parameters:
        graph_vertices_number (Union[sympy.Integer, sympy.Symbol]): Total number of vertices in the larger graph.
        graph_average_degree (Union[sympy.Integer, sympy.Symbol]): Mean degree of the larger graph.
        part_size (int): Size of each part in the bipartite subgraph.

    Returns:
        sympy.Expr: The expected number of occurrences of the complete bipartite subgraph - (1 + o(1)).
    """

    return expected_number_of_subgraphs(
        graph_vertices_number,
        graph_average_degree,
        2 * part_size,
        part_size * part_size,
        lambda g: nx.algorithms.bipartite.is_bipartite(g)
        and len(g.nodes) == 2 * part_size
        and all(d == part_size for _, d in g.degree()),
    )


# Usage:
if __name__ == "__main__":
    expected_number_of_k3_subgraphs = expected_number_of_cycle_subgraphs(
        sympy.symbols("n"),
        sympy.symbols("m"),
        3,
    )
    print(f"Expected number of K3 subgraphs: {expected_number_of_k3_subgraphs}")

    expected_number_of_k5_subgraphs = expected_number_of_clique_subgraphs(
        sympy.symbols("n"),
        sympy.symbols("m"),
        5,
    )
    print(f"Expected number of K5 subgraphs: {expected_number_of_k5_subgraphs}")

    expected_number_of_k33_subgraphs = expected_number_of_bipartite_complete_subgraphs(
        sympy.symbols("n"),
        sympy.symbols("m"),
        3,
    )
    print(f"Expected number of K3,3 subgraphs: {expected_number_of_k33_subgraphs}")
