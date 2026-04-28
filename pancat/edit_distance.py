"Computes edition between pangenome graphs"
from json import dump
from gfagraphs import Graph
from logging import info
from os import getpid
from psutil import Process
from sys import stderr
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import cpu_count
from typing import Callable
from collections.abc import Iterable


def futures_collector(
    func: Callable,
        argslist: list,
        kwargslist: list[dict] | None = None,
        num_processes: int = cpu_count(),
) -> list:
    """
    Spawns len(arglist) instances of func and executes them at num_processes instances at time.

    * func : a function
    * argslist (list): a list of tuples, arguments of each func
    * kwargslist (list[dict]) a list of dicts, kwargs for each func
    * num_processes (int) : max number of concurrent instances.
        Default : number of available logic cores
    * memory (float|None) : ratio of memory to be used, ranging from .05 to .95. Will not work if *resource* is incompatible.
    """
    if kwargslist is None or len(kwargslist) == len(argslist):
        with ThreadPoolExecutor(max_workers=num_processes) as executor:
            futures = [
                executor.submit(
                    func,
                    *args if isinstance(args, Iterable) else args
                ) if kwargslist is None else
                executor.submit(
                    func,
                    *args if isinstance(args, Iterable) else args,
                    **kwargslist[i]
                ) for i, args in enumerate(argslist)
            ]
        return [f.result() for f in futures]
    else:
        raise ValueError(
            f"""Positionnal argument list length ({len(argslist)})
            does not match keywords argument list length ({len(kwargslist)}).""")


def perform_edition(
        gfa_A: str,
        gfa_B: str,
        output_path: str,
        graph_level: bool,
        selection: list[str] | bool = True,
        cores: int = 1,
        regular_expression: str = ".*",
        trace_memory: bool = False,
) -> tuple:
    """
    In this function, we do calculate the distance between G1 and G2, by trying to modify G2 into G1.
    Note the two graphs can be freely swapped, we just need to invert scores for reciprocal events and operations

    Args:
        gfa_A (str): a path to a gfa file
        gfa_B (str): a path to another gfa file
        output_path (str): path where to store editions
        graph_level (bool): asks to perform graph-level edition
        selection (list[str] | bool, optional):
            - if true, compute edition on intersection of paths at path-level.
            - if a list[str], compute edition on specified paths
            Defaults to True.

    Returns:
        tuple: results from the edition
    """
    base_interpeter_memory: int = Process(
        getpid()).memory_info().rss / 1024 ** 2
    graph_A: Graph = Graph(
        gfa_file=gfa_A,
        with_sequence=False,
        low_memory=True,
        regexp=regular_expression
    )
    info(f"Loaded graph {gfa_A} in memory")
    print('Paths of Graph_A', ', '.join(graph_A.paths.keys()))
    graph_B: Graph = Graph(
        gfa_file=gfa_B,
        with_sequence=False,
        low_memory=True,
        regexp=regular_expression
    )
    info(f"Loaded graph {gfa_B} in memory")
    print('Paths of Graph_B', ', '.join(graph_B.paths.keys()))
    if trace_memory:
        info(
            f"Graph loading Python interpreter memory consumption: {(Process(getpid()).memory_info().rss / 1024 ** 2)-base_interpeter_memory} Mb used")

    # Prints out names of paths (for debugging purposes)

    results: dict = dict()
    if graph_level:
        if cores > 1:
            info(f"Computing edition for graphs with {cores} threads")
            results = graph_level_edition_multiprocessing(graph_A, graph_B)
        else:
            info("Computing edition for graphs on single thread")
            results = graph_level_edition(graph_A, graph_B)
    else:
        # We compute the intersection of paths in both graphs
        path_intersect: set[str] = set(
            graph_A.paths.keys()).intersection(set(graph_B.paths.keys()))
        if isinstance(selection, list):
            # We perform edition on selected paths, if all paths are in both graph
            if not all([x in path_intersect for x in selection]):
                raise ValueError()
            if cores > 1:
                info(f"Computing edition for graphs with {cores} threads")
                results = path_level_edition_multiprocessing(
                    graph_A, graph_B, set(selection))
            else:
                info(f"Computing edition for graphs on single thread")
                results = path_level_edition(graph_A, graph_B, set(selection))
        else:
            # We perform edition on shared paths, hoping the best for non-common paths \o/
            # (Best practice is to validate before if all paths are shared)
            if cores > 1:
                info(f"Computing edition for graphs with {cores} threads")
                results = path_level_edition_multiprocessing(
                    graph_A, graph_B, path_intersect)
            else:
                info(f"Computing edition for graphs on single thread")
                results = path_level_edition(graph_A, graph_B, path_intersect)
    info(f"Saving results in {output_path}")
    if output_path.endswith('.json'):
        dump(results, open(output_path, 'w', encoding='utf-8'))
    else:
        with open(output_path, 'w', encoding='utf-8') as logwriter:
            logwriter.write(
                f"# Editions results for {gfa_A} against {gfa_B} performed at {('path level','graph level')[graph_level]}\n"
            )
            logwriter.write(
                f'# Selected paths {(path_intersect,selection)[isinstance(selection, list)]}\n'
            )
            if graph_level:
                for category in ['merges', 'splits']:
                    logwriter.write(f'## {category}\n')
                    for edition in results[category]:
                        logwriter.write(
                            '\t'.join([str(__) for __ in edition])+'\n')
            else:
                for path_name, edit_category in results.items():
                    logwriter.write(f'## {path_name}\n')
                    for category in ['merges', 'splits']:
                        logwriter.write(f'## {category}\n')
                        for edition in edit_category[category]:
                            logwriter.write(
                                '\t'.join([str(__) for __ in edition])+'\n')
    info(f"Job terminated sucessfully!")


def path_level_edition(graph_A: Graph, graph_B: Graph, selected_paths: set[str]) -> dict:
    """Compute edition, path by path, between the two graphs.
    The graph_A will be used as reference

    Args:
        graph_A (Graph): pangenome graph
        graph_B (Graph): pangenome graph
        selected_paths (set[str]): the paths where the edition needs to be computed

    Returns:
        dict: results of edition
    """
    edition_results: dict = dict()

    print(f"PATH SELECTION: {selected_paths}")
    # Iterating on each pair of paths
    for path_name in selected_paths:
        print(f"Processing path: {path_name}", file=stderr)

        i: int = 0  # counter of segmentations on graph_A
        j: int = 0  # counter of segmentations on graph_B

        equivalences_count: int = 0  # counter of equivalences
        splits_count: int = 0  # counter of splits
        merges_count: int = 0  # counter of merges

        merges: set[tuple[int, tuple[str]]] = set()  # set for merges
        splits: set[tuple[int, tuple[str]]] = set()  # set for merges

        pos_A: int = 0  # Absolute position in BP on A
        pos_B: int = 0  # Absolute position in BP on B

        global_pos: int = 0  # Position across both genomes
        global_list: list[int] = list()

        # Iterating until we did not go through both segmentations
        while i < len(graph_A.paths[path_name]['path']) and j < len(graph_B.paths[path_name]['path']):
            # Currently evaluated nodes
            current_node_A: str = graph_A.paths[path_name]['path'][i][0]
            current_node_B: str = graph_B.paths[path_name]['path'][j][0]

            # We compute the next closest breakpoint
            global_pos = min(
                global_pos +
                (
                    graph_A.segments[current_node_A]
                    ['length']-(global_pos-pos_A)
                ),
                global_pos +
                (
                    graph_B.segments[current_node_B]
                    ['length']-(global_pos-pos_B)
                )
            )

            # We added the interval to current positions
            match (global_pos-pos_A == graph_A.segments[current_node_A]['length'], global_pos-pos_B == graph_B.segments[current_node_B]['length']):
                case (True, True):
                    equivalences_count += 1
                    # Iterating on both, no edition needed
                    pos_A += graph_A.segments[current_node_A]['length']
                    pos_B += graph_B.segments[current_node_B]['length']
                    i += 1
                    j += 1
                case (True, False):
                    # Iterating on top, split required
                    print(
                        f"{path_name}\t{global_pos}\tS\t{current_node_A}\t{current_node_B}", file=stderr)
                    splits_count += 1
                    splits.add(
                        (
                            global_pos,
                            (
                                current_node_A,
                                current_node_B,
                            )
                        )
                    )
                    pos_A += graph_A.segments[current_node_A]['length']
                    i += 1
                case (False, True):
                    # Iterating on bottom, merge required
                    print(
                        f"{path_name}\t{global_pos}\tM\t{current_node_A}\t{current_node_B}", file=stderr)
                    merges_count += 1
                    merges.add(
                        (
                            global_pos,
                            (
                                current_node_A,
                                current_node_B,
                            )
                        )
                    )
                    pos_B += graph_B.segments[current_node_B]['length']
                    j += 1
                case (False, False):
                    raise ValueError()

        edition_results[path_name] = {
            'merges': sorted(list(merges)),
            'splits': sorted(list(splits))
        }
        print(
            f"# Path {path_name} has {equivalences_count} equivalences, {merges_count} merges and {splits_count} splits")
    return edition_results


def graph_level_edition(graph_A: Graph, graph_B: Graph) -> set:
    """Compute edition, at graph level, between the two graphs.
    The graph_A will be used as reference.
    Key idea is to store a series of path-pos where to merge/split, and filter those editions

    Args:
        graph_A (Graph): pangenome graph
        graph_B (Graph): pangenome graph

    Returns:
        set: results of edition: contains sets of tuples path_name <-> position where the edition is
    """
    # Computing offsets on both graphs
    graph_A.sequence_offsets()
    graph_B.sequence_offsets()

    edition_results: dict[list] = dict()

    merges: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for merges
    splits: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for splits
    equiv: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for splits

    intersect: set = set(graph_A.paths.keys()).intersection(
        set(graph_B.paths.keys()))
    print(f"PATH INTERSECTION: {intersect}")
    # Iterating on each pair of paths
    for path_name in intersect:

        i: int = 0  # counter of segmentations on graph_A
        j: int = 0  # counter of segmentations on graph_B

        pos_A: int = 0  # Absolute position in BP on A
        pos_B: int = 0  # Absolute position in BP on B

        global_pos: int = 0  # Position across both genomes

        # Iterating until we did not go through both segmentations
        while i < len(graph_A.paths[path_name]['path']) and j < len(graph_B.paths[path_name]['path']):
            # Currently evaluated nodes
            current_node_A: str = graph_A.paths[path_name]['path'][i][0]
            current_node_B: str = graph_B.paths[path_name]['path'][j][0]

            # We compute the next closest breakpoint
            global_pos = min(
                global_pos +
                (
                    graph_A.segments[current_node_A]['length'] -
                    (global_pos-pos_A)
                ),
                global_pos +
                (
                    graph_B.segments[current_node_B]['length'] -
                    (global_pos-pos_B)
                )
            )

            # We added the interval to current positions
            match (global_pos-pos_A == graph_A.segments[current_node_A]['length'], global_pos-pos_B == graph_B.segments[current_node_B]['length']):
                case (True, True):
                    # Iterating on both, no edition needed
                    equiv.add(
                        frozenset(
                            (
                                (
                                    path_name,
                                    frozenset(
                                        x for x, _, _ in list_of_positions
                                    ),
                                    (
                                        current_node_A,
                                        current_node_B,
                                    )
                                )
                                for path_name, list_of_positions in graph_B.segments[current_node_B]['PO'].items()
                            )
                        )
                    )
                    pos_A += graph_A.segments[current_node_A]['length']
                    pos_B += graph_B.segments[current_node_B]['length']
                    i += 1
                    j += 1
                case (True, False):
                    # Iterating on top, split required
                    splits.add(
                        frozenset(
                            (
                                (
                                    path_name,
                                    frozenset(
                                        x+global_pos-pos_A for x, _, _ in list_of_positions
                                    ),
                                    (
                                        current_node_A,
                                        current_node_B,
                                    )
                                )
                                for path_name, list_of_positions in graph_B.segments[current_node_B]['PO'].items()
                            )
                        )
                    )
                    pos_A += graph_A.segments[current_node_A]['length']
                    i += 1
                case (False, True):
                    # Iterating on bottom, merge required
                    merges.add(
                        frozenset(
                            (
                                (
                                    path_name,
                                    frozenset(
                                        x for _, x, _ in list_of_positions
                                    ),
                                    (
                                        current_node_A,
                                        current_node_B,
                                    )
                                )
                                for path_name, list_of_positions in graph_B.segments[current_node_B]['PO'].items()
                            ),
                        )
                    )
                    pos_B += graph_B.segments[current_node_B]['length']
                    j += 1
                case (False, False):
                    raise ValueError()

    edition_results['merges'] = [
        [
            (path_name, [x for x in pos], nodes) for path_name, pos, nodes in ext_fset
        ] for ext_fset in merges
    ]
    edition_results['splits'] = [
        [
            (path_name, [x for x in pos], nodes) for path_name, pos, nodes in ext_fset
        ] for ext_fset in splits
    ]
    edition_results['equiv'] = [
        [
            (path_name, [x for x in pos], nodes) for path_name, pos, nodes in ext_fset
        ] for ext_fset in equiv
    ]

    return edition_results


def graph_level_edition_multiprocessing(graph_A: Graph, graph_B: Graph) -> set:
    """Compute edition, at graph level, between the two graphs.
    The graph_A will be used as reference.
    Key idea is to store a series of path-pos where to merge/split, and filter those editions

    Args:
        graph_A (Graph): pangenome graph
        graph_B (Graph): pangenome graph

    Returns:
        set: results of edition: contains sets of tuples path_name <-> position where the edition is
    """
    # Computing offsets on both graphs
    graph_A.sequence_offsets()
    graph_B.sequence_offsets()

    edition_results: dict[list] = dict()

    merges: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for merges
    splits: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for splits

    intersect: set = set(graph_A.paths.keys()).intersection(
        set(graph_B.paths.keys()))
    print(f"PATH INTERSECTION: {intersect}")
    editions: list[tuple[set, set]] = futures_collector(
        func=edit_single_path_graph_level,
        argslist=[
            (
                path_name,
                graph_A,
                graph_B
            ) for path_name in intersect
        ]
    )

    for ml, sl in editions:
        merges = merges.union(ml)
        splits = splits.union(sl)

    edition_results['merges'] = [
        [
            (path_name, [x for x in pos], nodes) for path_name, pos, nodes in ext_fset
        ] for ext_fset in merges
    ]
    edition_results['splits'] = [
        [
            (path_name, [x for x in pos], nodes) for path_name, pos, nodes in ext_fset
        ] for ext_fset in splits
    ]

    return edition_results


def edit_single_path_graph_level(path_name: str, graph_A: Graph, graph_B: Graph) -> tuple[set]:
    """Perform the edition calculation on a single path, and returns positions of merges and splits

    Args:
        path_name (str): the path being investigated
        graph_A (Graph): first graph
        graph_B (Graph): second graph

    Raises:
        ValueError: if no edition can be determined

    Returns:
        tuple[set]: (merges, splits)
    """
    merges: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for merges
    splits: set[frozenset[tuple[str, frozenset, tuple[str]]]
                ] = set()  # set for splits

    i: int = 0  # counter of segmentations on graph_A
    j: int = 0  # counter of segmentations on graph_B

    pos_A: int = 0  # Absolute position in BP on A
    pos_B: int = 0  # Absolute position in BP on B

    global_pos: int = 0  # Position across both genomes

    # Iterating until we did not go through both segmentations
    while i < len(graph_A.paths[path_name]['path']) and j < len(graph_B.paths[path_name]['path']):
        # Currently evaluated nodes
        current_node_A: str = graph_A.paths[path_name]['path'][i][0]
        current_node_B: str = graph_B.paths[path_name]['path'][j][0]

        # We compute the next closest breakpoint
        global_pos = min(
            global_pos +
            (
                graph_A.segments[current_node_A]['length'] -
                (global_pos-pos_A)
            ),
            global_pos +
            (
                graph_B.segments[current_node_B]['length'] -
                (global_pos-pos_B)
            )
        )

        # We added the interval to current positions
        match (global_pos-pos_A == graph_A.segments[current_node_A]['length'], global_pos-pos_B == graph_B.segments[current_node_B]['length']):
            case (True, True):
                # Iterating on both, no edition needed
                pos_A += graph_A.segments[current_node_A]['length']
                pos_B += graph_B.segments[current_node_B]['length']
                i += 1
                j += 1
            case (True, False):
                # Iterating on top, split required
                splits.add(
                    frozenset(
                        (
                            (
                                path_name,
                                frozenset(
                                    x+global_pos-pos_A for x, _, _ in list_of_positions
                                ),
                                (
                                    current_node_A,
                                    current_node_B,
                                )
                            )
                            for path_name, list_of_positions in graph_B.segments[current_node_B]['PO'].items()
                        )
                    )
                )
                pos_A += graph_A.segments[current_node_A]['length']
                i += 1
            case (False, True):
                # Iterating on bottom, merge required
                merges.add(
                    (
                        frozenset(
                            (
                                path_name,
                                frozenset(
                                    x for _, x, _ in list_of_positions
                                ),
                                (
                                    current_node_A,
                                    current_node_B,
                                )
                            )
                            for path_name, list_of_positions in graph_B.segments[current_node_B]['PO'].items()
                        )
                    )
                )
                pos_B += graph_B.segments[current_node_B]['length']
                j += 1
            case (False, False):
                raise ValueError()
    return (merges, splits)


def path_level_edition_multiprocessing(graph_A: Graph, graph_B: Graph, selected_paths: set[str]) -> dict:
    """Compute edition, path by path, between the two graphs.
    The graph_A will be used as reference

    Args:
        graph_A (Graph): pangenome graph
        graph_B (Graph): pangenome graph
        selected_paths (set[str]): the paths where the edition needs to be computed

    Returns:
        dict: results of edition
    """

    paths_of_graph: list[str] = list(selected_paths)
    # Iterating on each pair of paths
    print(f"PATH SELECTION: {paths_of_graph}")

    editions: list[dict] = futures_collector(
        func=edit_single_path_path_level,
        argslist=[
            (path_name, graph_A, graph_B) for path_name in paths_of_graph
        ]
    )

    return {
        paths_of_graph[i]: edits for i, edits in enumerate(editions)
    }


def edit_single_path_path_level(path_name: str, graph_A: Graph, graph_B: Graph) -> dict:
    """Performs the edition over a single path

    Args:
        path_name (str): name of the path
        graph_A (Graph): first graph
        graph_B (Graph): second graph

    Raises:
        ValueError: if no edition is forseeable

    Returns:
        dict: editions
    """
    i: int = 0  # counter of segmentations on graph_A
    j: int = 0  # counter of segmentations on graph_B

    merges: set[tuple[int, tuple[str]]] = set()  # set for merges
    splits: set[tuple[int, tuple[str]]] = set()  # set for merges

    pos_A: int = 0  # Absolute position in BP on A
    pos_B: int = 0  # Absolute position in BP on B

    global_pos: int = 0  # Position across both genomes

    # Iterating until we did not go through both segmentations
    while i < len(graph_A.paths[path_name]['path']) and j < len(graph_B.paths[path_name]['path']):
        # Currently evaluated nodes
        current_node_A: str = graph_A.paths[path_name]['path'][i][0]
        current_node_B: str = graph_B.paths[path_name]['path'][j][0]

        # We compute the next closest breakpoint
        global_pos = min(
            global_pos +
            (
                graph_A.segments[current_node_A]
                ['length']-(global_pos-pos_A)
            ),
            global_pos +
            (
                graph_B.segments[current_node_B]
                ['length']-(global_pos-pos_B)
            )
        )

        # We added the interval to current positions
        match (global_pos-pos_A == graph_A.segments[current_node_A]['length'], global_pos-pos_B == graph_B.segments[current_node_B]['length']):
            case (True, True):
                # Iterating on both, no edition needed
                pos_A += graph_A.segments[current_node_A]['length']
                pos_B += graph_B.segments[current_node_B]['length']
                i += 1
                j += 1
            case (True, False):
                # Iterating on top, split required
                splits.add(
                    (
                        global_pos,
                        (
                            current_node_A,
                            current_node_B,
                        )
                    )
                )
                pos_A += graph_A.segments[current_node_A]['length']
                i += 1
            case (False, True):
                # Iterating on bottom, merge required
                merges.add(
                    (
                        global_pos,
                        (
                            current_node_A,
                            current_node_B,
                        )
                    )
                )
                pos_B += graph_B.segments[current_node_B]['length']
                j += 1
            case (False, False):
                raise ValueError()
    return {
        'merges': sorted(list(merges)),
        'splits': sorted(list(splits))
    }
