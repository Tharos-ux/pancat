from typing import Generator
from gfagraphs import Graph
from tharospytools import path_allocator, flatten, revcomp


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
    gfa_paths: list = gfa_graph.get_path_list()

    all_sets = {
        path.datas['name']:
            [
                node_name for node_name, _ in path.datas['path']
        ]
        for path in gfa_paths
    }

    all_maps = {
        path.datas['name']: path.datas['path'] for path in gfa_paths
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
        bubbles: list[dict] = [{}
                               for _ in range(len(bubbles_endpoints)-1)]
        for path in gfa_paths:
            # Computing endpoint positions in list for each path
            endpoints_indexes: list = grouper(
                [
                    all_sets[
                        path.datas['name']
                    ].index(
                        endpoint
                    ) for endpoint in bubbles_endpoints
                ],
                2
            )
            # Getting bubble chains
            for i, (start, end) in enumerate(endpoints_indexes):
                bubbles[i][path.datas['name']
                           ] = all_sets[path.datas['name']][start:end]

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
    for path in gfa_paths:
        endpoints_indexes: list = grouper(
            [
                all_sets[
                    path.datas['name']
                ].index(
                    endpoint
                ) for endpoint in bubbles_endpoints
            ],
            2
        )
        for i, (start, end) in enumerate(endpoints_indexes):
            oriented_bubbles[i][path.datas['name']
                                ] = all_maps[path.datas['name']][start:end+1]
    return oriented_bubbles


def call_sequences(gfa_file: str, gfa_type: str) -> Generator:
    """Given a GFA file and a path name, extracts all chains from it

    Args:
        gfa_file (str): path to a gfa file
        gfa_type (str): subformat
    """
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_type,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
    for bubble in bubbles:
        yield {
            path_name: ''.join(
                [gfa_graph.get_segment(node=node).datas['seq'] if orientation.value == '+' else revcomp(gfa_graph.get_segment(node=node).datas['seq'])
                 for node, orientation in path_chain]
            ) for path_name, path_chain in bubble.items()
        }


if __name__ == "__main__":
    gfa_file: str = "/home/sidubois/Workspace/Notes/graph_cactus_sandra.gfa"
    gfa_ver: str = "GFA1.1"

    all_sequences: Generator = call_sequences(
        gfa_file=gfa_file,
        gfa_type=gfa_ver
    )

    bubble_count: int = 0
    for sequence in all_sequences:
        bubble_count += 1
    print(bubble_count)
