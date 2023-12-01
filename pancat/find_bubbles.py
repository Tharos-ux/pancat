from typing import Generator
from pgGraphs import Graph
from tharospytools.bio_tools import revcomp
from tharospytools.list_tools import flatten
from tharospytools.path_tools import path_allocator


def grouper(iterable, n=2, m=1):
    """Collect data into overlapping fixed-length chunks or blocks"""
    return [iterable[i:i+n] for i in range(0, len(iterable)-1, n-m)]


def common_members(elements: list[set]):
    first_path, other_paths = elements[0], elements[1:]
    return sorted(list(first_path.intersection(*other_paths)))


def bubble_caller(gfa_graph: Graph) -> list[dict]:
    """Calls out the stritly disjoint (super)bubbles in the graph.
    A bubble can be defined as having a starting and an ending node
    with a in and out node with degree equal to the number of paths
    for superbubble level we don't have to watch the order, as 

    Args:
        gfa_file (str): path to a gfa-like file

    Returns:
        list[dict]: a list of mappings between paths names and the subchain in the bubble
                    one element per bubble
    """
    all_sets = {
        path_name:
            [
                node_name for node_name, _ in path_datas['path']
            ]
        for path_name, path_datas in gfa_graph.paths.items()
    }

    all_maps = {
        path_name: path_datas['path'] for path_name, path_datas in gfa_graph.paths.items()
    }

    bubbles_endpoints: list = sorted(common_members(
        list(
            set(x) for x in all_sets.values()
        )
    ), key=int)

    # loop until convergence
    convergence: int = 0
    while (len(bubbles_endpoints) != convergence):
        convergence: int = len(bubbles_endpoints)
        bubbles: list[dict] = [
            {} for _ in range(len(bubbles_endpoints)-1)
        ]
        for path_name in gfa_graph.paths.keys():
            # Computing endpoint positions in list for each path
            endpoints_indexes: list = grouper(
                [
                    all_sets[
                        path_name
                    ].index(
                        endpoint
                    ) for endpoint in bubbles_endpoints
                ],
                2
            )
            # Getting bubble chains
            for i, (start, end) in enumerate(endpoints_indexes):
                bubbles[i][path_name] = all_sets[path_name][start:end]

            # Decompressing all paths
            embed_nodes: set = set(flatten(
                [chain_content[1:-1] for bubble in bubbles for _, chain_content in bubble.items()]))
            # Search for endpoints that are in the set
            # if in the middle of the chain we notice a endpoint, THIS IS NO ENDPOINT and we need to clear it
            bubbles_endpoints = [
                endpoint for endpoint in bubbles_endpoints if endpoint not in embed_nodes]
            # If we need to suppress endpoints, we will have different length, so we will loop
    #  Extracting reading way
    oriented_bubbles: list[dict] = [{}
                                    for _ in range(len(bubbles_endpoints)-1)]
    for path_name in gfa_graph.paths.keys():
        endpoints_indexes: list = grouper(
            [
                all_sets[
                    path_name
                ].index(
                    endpoint
                ) for endpoint in bubbles_endpoints
            ],
            2
        )
        for i, (start, end) in enumerate(endpoints_indexes):
            oriented_bubbles[i][path_name
                                ] = all_maps[path_name][start:end+1]
    return oriented_bubbles


def call_variants(gfa_file: str, gfa_type: str) -> Generator:
    """Given a GFA file and a path name, calls all rank 1 variants against it

    Args:
        gfa_file (str): path to a gfa file
        gfa_type (str): subformat
    """
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
    for bubble in bubbles:
        yield {
            path_name: ''.join(
                [
                    gfa_graph.segments[node]['seq'] if orientation.value == '+' else revcomp(
                        gfa_graph.segments[node]['seq']
                    )
                    for node, orientation in path_chain
                ]
            ) for path_name, path_chain in bubble.items()
        }


def flattenable_bubbles(bubbles: list[dict]) -> Generator:
    """Only returns superbubbles given a list of bubbles

    Args:
        bubbles (list[dict]): a list of mixed bubbles

    Yields:
        list[dict]: a list of superbubbles
    """
    yield from [bubble for bubble in bubbles if any([len(chain) > 3 for chain in bubble.values()])]


def linearize_bubbles(gfa_file: str, gfa_type: str, output: str) -> Generator:
    """Given a GFA file, flattens the bubbles

    Args:
        gfa_file (str): path to a gfa file
        gfa_type (str): subformat
        output (str): output file path
    """
    output_path: str = path_allocator(
        output, particle=".gfa", default_name="graph")
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
    if any([not len(x) for bubble in bubbles for x in bubble.values()]):
        raise ValueError(
            "Unexepected event detected. A zero-sized bubble was detected, indicating a potential loop."
        )
    # Init return graph
    output_graph: Graph = Graph(
        gfa_file=None,
        with_sequence=True
    )
    output_graph.headers = gfa_graph.headers
    contained_nodes: set = set()
    path_builder: dict = {
        path_name: []
        for path_name in gfa_graph.paths.keys()
    }
    # For each bubble, we compute new nodes
    for bubble in bubbles:
        for path_name, path_chain in bubble.items():
            # min bubble is of size 2, if two nodes are one next to another
            (source, ori_source), chained, (sink, ori_sink) = path_chain[0], path_chain[1:len(
                path_chain)-1], path_chain[-1]
            for node, ori in [(source, ori_source), (sink, ori_sink)]:
                if node not in contained_nodes:
                    output_graph.add_node(node, gfa_graph.segments[node]['seq'] if ori.value == '+' else revcomp(
                        gfa_graph.segments[node]['seq']))
                    contained_nodes.add(node)

            # if theres a central chain
            if len(chained) > 0:
                # node + edge between source + new node and new node + sink
                target: str = chained[0][0]
                output_graph.add_node(target, ''.join([gfa_graph.segments[node]['seq'] if orientation.value == '+' else revcomp(
                    gfa_graph.segments[node]['seq'])for node, orientation in chained]))
                output_graph.add_edge(source, ori_source.value, target, '+')
                output_graph.add_edge(target, '+', sink, ori_sink.value)
                if len(path_builder[path_name]) == 0:
                    path_builder[path_name] = [
                        (source, ori_source), (target, '+'), (sink, ori_sink)]
                else:
                    path_builder[path_name] = path_builder[path_name] + \
                        [(target, '+'), (sink, ori_sink)]
            else:
                # edge between source and sink
                output_graph.add_edge(
                    source, ori_source.value, sink, ori_sink.value)
                if len(path_builder[path_name]) == 0:
                    path_builder[path_name] = [
                        (source, ori_source), (sink, ori_sink)]
                else:
                    path_builder[path_name] = path_builder[path_name] + \
                        [(sink, ori_sink)]
    # Adding paths
    for path_name, chain in path_builder.items():
        output_graph.add_path(path_name, chain)
    output_graph.save_graph(output_path)
