"Creates a graph we can navigate in."
from os import path, remove
from json import load
from networkx import MultiDiGraph, compose
from pyvis.network import Network
from gfagraphs import Graph, GFANetwork, GFAParser, GFAFormat
from os.path import exists, isdir
from pathlib import Path
from typing import Iterable


def find_interval_on_path(
    graph: Graph,
    path_name: str,
    start: int,
    end: int,
) -> tuple[int, int]:
    """Finds the interval on a path (in node positions in this path) that corresponds to start and end positions.
    Ignore PO offset but compute it on the fly for each path, and do not store it in the graph
    It will help to compute the path offset for each path in the graph without loops getting in the way

    Args:
        graph (Graph): a gfagraphs Graph object
        path_name (str): the path name
        start (int): start position
        end (int): end position

    Returns:
        tuple[int,int]: a tuple of two integers
    """
    # Handling the case where the start and end are inverted
    if start > end:
        start, end = end, start
    # Handling the case where the start is negative
    if start < 0:
        start = 0

    # Searching for the start and end node indexes in the reference path
    start_node_index: int = -1
    end_node_index: int = -1
    current_position: int = graph.paths[path_name]['start_offset'] or 0
    for i, (node, _) in enumerate(graph.paths[path_name]['path']):
        # We search for the first and last node indexes in path corresponding to the start-stop interval
        if current_position < start:
            start_node_index = i
        if current_position >= end:
            end_node_index = i
            break
        current_position += graph.segments[node]['length']

    # Handling the case where the start or end of the interval is not reached
    if start_node_index == -1:
        start_node_index = 0
    if end_node_index == -1:
        end_node_index = len(graph.paths[path_name]['path'])-1

    print(f"Found start and endpoints〈{graph.paths[path_name]['path'][start_node_index][0]},{graph.paths[path_name]['path'][end_node_index][0]}〉")

    return (start_node_index, end_node_index)


def extract_subgraph(gfa_path: str, x: int, y: int, sequence: str, output: str, save_to_file: bool = True) -> None:
    """Gathers the subgraph from two points on a reference

    Parameters
    ----------
    gfa_path : str
        path to file
    x : int
        start point
    y : int
        end point
    sequence : str
        reference genome for positions
    output : str
        path to output file
    """
    # Load the graph
    gfa_graph: Graph = load_graph(gfa_path)

    # If no start or end is given, we return the whole graph
    if x is None or y is None or sequence is None:
        return (gfa_graph, set(gfa_graph.segments.keys()))

    # We use the find_interval_on_path function to get the start and end nodes
    source_index, sink_index = find_interval_on_path(
        graph=gfa_graph,
        path_name=sequence,
        start=x,
        end=y,
    )

    source:str = gfa_graph.paths[sequence]['path'][source_index][0]
    sink:str = gfa_graph.paths[sequence]['path'][sink_index][0]
    
    path_dict: dict = dict()
    # We use the intersection of the paths crossing the nodes to get the subgraph
    path_intersect: set[str] = set(gfa_graph.segments[source]['PO'].keys()).intersection(
        set(gfa_graph.segments[sink]['PO'].keys())
    )
    # We iterate on the intersection to seek source and sink nodes
    for path_name in path_intersect:
        # We will only keep non spurious node pairs
        if len(gfa_graph.segments[sink]['PO'][path_name]) == len(gfa_graph.segments[source]['PO'][path_name]) == 1:
            source_index, sink_index = -1, -1
            for i, (node, _) in enumerate(gfa_graph.paths[path_name]['path']):
                if node == source:
                    source_index = i
                if node == sink:
                    sink_index = i
            if source_index > sink_index:
                source_index, sink_index = sink_index, source_index

            # We write the new path
            new_path: list[tuple] = gfa_graph.paths[path_name]['path'][source_index:sink_index+1]

    for path_name, path_data in gfa_graph.paths.items():
        loops_over_path: int = 0
        source_index, sink_index = None, None
        for i, (seg_name, _) in enumerate(path_data['path']):
            try:
                if seg_name == source:
                    source_index: int = i
                    source_offsets: tuple = gfa_graph.segments[seg_name]['PO'][sequence][loops_over_path]

                if seg_name == sink:
                    sink_index: int = i
                    sink_offsets: tuple = gfa_graph.segments[seg_name]['PO'][sequence][loops_over_path]

                if source_index is not None and sink_index is not None:
                    # Write new path
                    if source_index > sink_index:
                        new_path: list[tuple] = path_data['path'][sink_index:source_index+1]
                    else:
                        new_path: list[tuple] = path_data['path'][source_index:sink_index+1]
                    start_position: int = min(
                        source_offsets[0], sink_offsets[0], source_offsets[1], sink_offsets[1])
                    end_position: int = max(
                        source_offsets[0], sink_offsets[0], source_offsets[1], sink_offsets[1])
                    path_dict[f"{path_name}\0{loops_over_path}"] = {
                        "name": f"{path_name}\0{loops_over_path}",
                        "start_offset": start_position,
                        "stop_offset": end_position,
                        "path": new_path,
                    }
                    # Reset indexes
                    source_index, sink_index = None, None
                    loops_over_path += 1
            except IndexError:
                pass

    gfa_graph.paths = path_dict

    selected_nodes_set: set[str] = {
        x for z in path_dict.values() for x, _ in z['path']
    }

    if not save_to_file:
        return (gfa_graph, selected_nodes_set)
    GFAParser.save_subgraph(
        graph=gfa_graph,
        output_path=output,
        nodes=selected_nodes_set,
        force_format=False,
        minimal_graph=True
    )


def grouper(iterable: Iterable, n: int = 2, m: int = 1) -> list:
    """Collect data into possibly overlapping fixed-length chunks or blocks

    Args:
        iterable (Iterable): the iterable you want to pack
        n (int, optional): size of the chunks. Defaults to 2.
        m (int, optional): overlaps of the chunks. Defaults to 1.

    Returns:
        list: a list of n-sized chunks with an overlap of m
    """
    return [iterable[i:i+n] for i in range(0, len(iterable)-1, n-m)]


def path_allocator(
    path_to_validate: str,
    particle: str | None = None,
    default_name: str = 'file',
    always_yes: bool = True
) -> str:
    """Checks if a file exists in this place, and arborescence exists.
    If not, creates the arborescence

    Args:
        path_to_validate (str): a string path to the file
        particle (str | None, optional): file extension. Defaults to None.
        default_name (str): a name if name is empty
        always_yes (bool, optional): if file shall be erased by default. Defaults to True.

    Returns:
        str: the path to the file, with extension
    """
    if ('/') in path_to_validate:
        if isdir(path_to_validate) and not path_to_validate.endswith('/'):
            path_to_validate = path_to_validate + '/'
        folder_path, sep, file_name = path_to_validate.rpartition('/')
    else:
        folder_path = ""
        sep = ""
        file_name = path_to_validate
    if file_name == "":
        file_name = default_name
    if particle and not file_name.endswith(particle):
        file_name = file_name+particle
    full_path: str = folder_path+sep+file_name
    if not always_yes and exists(full_path):
        if not input('File already exists. Do you want to write over it? (y/n): ').lower().strip() == 'y':
            raise OSError("File already exists. Aborting.")
    if folder_path != "":
        Path(folder_path).mkdir(parents=True, exist_ok=True)
    return full_path


def display_graph(graph: MultiDiGraph, colors_paths: dict[str, str], annotations: dict, output_path: str) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        graph (MultiDiGraph): a graph combining multiple pangenomes to highlight thier similarities
        name (str): output name for graph render
        colors_paths (dict[str, str]): a set of colors to keep path colors consistent
    """
    output_path: str = path_allocator(
        output_path, particle='.html', default_name='graph')
    output_path_temp: str = path_allocator(
        output_path, particle='.tmp.html', default_name='graph')
    graph_visualizer: Network = Network(
        height='1000px', width='100%', directed=True, select_menu=False, filter_menu=False, bgcolor='#ffffff')
    graph_visualizer.set_template_dir(path.dirname(__file__), 'template.html')
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    html: str = graph_visualizer.generate_html()
    legend: str = '\n'.join(
        [f"<li><span class='{key}'></span> <a href='#'>{key}</a></li>" for key in colors_paths.keys()])
    with open(output_path_temp, "w+", encoding='utf-8') as out:
        out.write(html)
    with open(output_path, "w", encoding="utf-8") as html_writer:
        with open(output_path_temp, "r", encoding="utf-8") as html_file:
            for line in html_file:
                if "<div class='sidenav'>" in line:
                    html_writer.write(
                        f"""{line}{''.join(
                            [
                                "<a href='#' title=''>"+str(key)+" : <b>"+str(value)+"</b></a>" for key,value in annotations.items()
                                ]
                                )}\n<ul class='legend'>{legend}</ul>"""
                    )
                elif "/* your colors */" in line:
                    html_writer.write(
                        ''.join(
                            [".legend ."+key+" { background-color: "+val +
                                "; }" for key, val in colors_paths.items()]
                        )
                    )
                else:
                    html_writer.write(line)
    if path.exists(output_path_temp):
        remove(output_path_temp)


def compute_stats(
    graph: Graph,
    length_classes: tuple[list] = (
        [0, 1], [2, 10], [11, 50], [51, 200], [201, 500], [
            501, 1000], [1001, 10000], [10001, float('inf')]
    )
) -> dict:
    """Computes some basic metrics for the graph

    Args:
        graph (Graph): a gfagraphs Graph object

    Returns:
        dict: a container for metrics
    """

    stats: dict = {}
    stats["Number of segments"] = len(graph.segments)
    stats["Number of edges"] = len(graph.lines)

    # segment_sizes[size_class] = [number,cum_length]
    segment_sizes: dict = {x: [0, 0] for x in range(len(length_classes))}
    total_size: int = 0
    for seg_datas in graph.segments.values():
        for x in segment_sizes.keys():
            low_bound, high_bound = length_classes[x]
            if seg_datas['length'] >= low_bound and seg_datas['length'] <= high_bound:
                segment_sizes[x] = [segment_sizes[x]
                                    [0]+1, segment_sizes[x][1]+seg_datas['length']]
                total_size += seg_datas['length']
    for x, (number, cum_size) in segment_sizes.items():
        low_bound, high_bound = length_classes[x]
        stats[f'Number of {low_bound} bp - {high_bound} bp'] = f"{number} ({round((number/len(graph.segments))*100,ndigits=2)}%)"
        stats[f'Size of {low_bound} bp - {high_bound} bp'] = f"{cum_size} ({round((cum_size/total_size)*100,ndigits=2)}%)"
    stats["Total size of segments"] = total_size
    try:
        stats["Is graph acyclic"] = all([len(set([x for x, _ in path_datas["path"]])) == len(
            path_datas["path"]) for path_datas in graph.paths.values])
    except:
        pass
    return stats


def multigraph_viewer(
    file_A: str,
    file_B: str,
    file_editions: str,
    boundaries: list,
    output: str,
    reference: str | None = None,
    start: int | None = None,
    end: int | None = None,
) -> None:
    """Renders a digraph alignment and explicits their connexions inbetween them.

    Args:
        file_A (str): first gfa file
        file_B (str): second gfa file
        file_editions (str): json file containing editions, must be graph-level
        boundaries (list): list of node class sizes
        output (str): output for html file
    """
    bounds: list[tuple[int, int]] = grouper(
        [0] + [bound+x for bound in boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
    )[::2]

    # Creating gfagraphs object
    gfa_graph_A: Graph = clean_graph(
        *extract_subgraph(
            gfa_path=file_A,
            x=start,
            y=end,
            sequence=reference,
            output='temp_A.gfa',
            save_to_file=False,
        )
    )
    gfa_graph_B: Graph = clean_graph(
        *extract_subgraph(
            gfa_path=file_B,
            x=start,
            y=end,
            sequence=reference,
            output='temp_B.gfa',
            save_to_file=False,
        )
    )

    # Computing NetworkX structure
    pangenome_graph_A: MultiDiGraph = GFANetwork.compute_networkx(
        graph=gfa_graph_A,
        node_prefix='A',
        node_size_classes=bounds
    )
    pangenome_graph_B: MultiDiGraph = GFANetwork.compute_networkx(
        graph=gfa_graph_B,
        node_prefix='B',
        node_size_classes=bounds
    )
    nodes_a:list[str] = [node[2:] for node in pangenome_graph_A.nodes.keys()]
    nodes_b:list[str] = [node[2:] for node in pangenome_graph_B.nodes.keys()]

    editions: dict = load(open(file_editions, 'r', encoding='utf-8'))

    edge_bank: dict[tuple, str] = dict()
    for basename in editions.keys():
        for style in ['merges', 'splits']:
            for _, [node_A, node_B] in editions[basename][style]:
                if node_A in nodes_a and node_B in nodes_b:
                    v_1,v_2 = edge_bank.get((node_A,node_B),(0,0))
                    if style == 'merges':
                        v_1 += 1
                    elif style == 'splits':
                        v_2 += 1
                    edge_bank[(node_A,node_B)] = (v_1,v_2)


    full_graph: MultiDiGraph = compose(
        pangenome_graph_A,
        pangenome_graph_B
    )


    for (a, b), (v_1,v_2) in edge_bank.items():
        full_graph.add_edge(
            f"A_{a}",
            f"B_{b}",
            arrows='',
            alpha=.5,
            color='blue',
            label=f"{v_1}:{v_2}",
        )

    stats_A = compute_stats(
        graph=gfa_graph_A,
        length_classes=tuple(bounds)
    )
    stats_B = compute_stats(
        graph=gfa_graph_B,
        length_classes=tuple(bounds)
    )

    # Computing stats and displaying
    display_graph(
        graph=full_graph,
        colors_paths={
            **gfa_graph_A.metadata['colors'], **gfa_graph_B.metadata['colors']},
        annotations={
            **{f'A_{key}': val for key, val in stats_A.items()},
            **{f'B_{key}': val for key, val in stats_B.items()}
        },
        output_path=output
    )


def graph_viewer(
    file: str,
    output: str,
    boundaries: list,
    force_rgfa: bool,
    start: int | None = None,
    end: int | None = None,
    reference: str | None = None,
) -> None:
    bounds: list[tuple[int, int]] = grouper(
        [0] + [bound+x for bound in boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
    )[::2]

    # Creating gfagraphs object
    if start is None and end is None:
        gfa_graph: Graph = Graph(
            gfa_file=file,
            with_sequence=True
        )
    else:
        gfa_graph: Graph = clean_graph(
            *extract_subgraph(
                gfa_path=file,
                x=start,
                y=end,
                sequence=reference,
                output='temp.gfa',
                save_to_file=False,
            )
        )

    # Computing NetworkX structure
    pangenome_graph: MultiDiGraph = GFANetwork.compute_networkx(
        graph=gfa_graph,
        enforce_format=GFAFormat.RGFA if force_rgfa else GFAFormat.GFA1,
        node_prefix=None,
        node_size_classes=bounds,
    )

    # Computing stats and displaying
    stats = compute_stats(
        graph=gfa_graph,
        length_classes=tuple(bounds)
    )
    display_graph(
        graph=pangenome_graph,
        colors_paths=gfa_graph.metadata['colors'],
        annotations=stats,
        output_path=output,
    )


def graph_stats(file: str, boundaries: list) -> None:
    bounds = grouper(
        [0] + [bound+x for bound in boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
    )

    # Creating gfagraphs object
    gfa_graph: Graph = Graph(
        gfa_file=file,
        with_sequence=True
    )

    # Computing stats and displaying
    graph_stats = compute_stats(
        graph=gfa_graph,
        length_classes=tuple(bounds)
    )
    for key, value in graph_stats.items():
        print(f"{key}: {value}")


def load_graph(graph_path: str) -> Graph:
    """Load a graph in memory and computes its offsets

    Parameters
    ----------
    graph_path : str
        path to a valid gfa file

    Returns
    -------
    Graph
        annotated gfagraphs object
    """
    gfa: Graph = Graph(
        gfa_file=graph_path,
        with_sequence=True,
        with_reverse_edges=False,
        low_memory=False,
        regexp='.*',
    )
    gfa.compute_neighbors()
    gfa.sequence_offsets()
    return gfa


def find_anchors(gfa: Graph, interval: tuple[int, int], reference: str) -> tuple[str, str]:
    """Finds the pair of nodes that corresponds to start and end positions

    Parameters
    ----------
    gfa : Graph
        loaded and anotated gfa graph
    positions : tuple[int,int]
        a couple of (start,stop) positions
    path_name : str
        name of the path for positions

    Returns
    -------
    tuple[str,str]
        a tuple of node names
    """
    if interval[0] > interval[1]:
        interval = (interval[1], interval[0])
    if interval[0] < 0:
        interval = (0, interval[1])
    if interval[1] > (ln := sum(gfa.segments[x]['length'] for x, _ in gfa.paths[reference]['path'])) or interval[0] > interval[1]:
        interval = (interval[0], ln)

    encountered_nodes: dict[str, int] = {}
    start: bool | str = False
    stop: bool | str = False
    for node, _ in gfa.paths[reference]['path']:
        if gfa.segments[node]['PO'][reference][encountered_nodes.get(node, 0)][1] >= interval[0] and not start:
            start = node
            stac: int = encountered_nodes.get(node, 0)
        if gfa.segments[node]['PO'][reference][encountered_nodes.get(node, 0)][1] >= interval[1] and start and not stop:
            stop = node
            stoc: int = encountered_nodes.get(node, 0)
        encountered_nodes[node] = encountered_nodes.get(node, 0) + 1
    return ((start, stac), (stop, stoc))


def clean_graph(gfa_graph: Graph, node_between_set: set[str] | None) -> Graph:
    """Cleans the graph by removing any node that is not in the set.
    Edits paths to reflect the changes.
    Edits edges to reflect the changes.

    Parameters
    ----------
    gfa_graph : Graph
        a gfagraphs objet
    node_between_set : set[str]
        nodes to keep

    Returns
    -------
    Graph
        a graph that has been cleaned
    """
    if node_between_set:
        gfa_graph.segments = {
            node_name: node_data for node_name, node_data in gfa_graph.segments.items() if node_name in node_between_set}
        for path_data in gfa_graph.paths.values():
            path_data['path'] = [
                (node, vect) for node, vect in path_data['path'] if node in node_between_set]
        gfa_graph.lines = {
            (source, sink): line_data for (source, sink), line_data in gfa_graph.lines.items() if source in node_between_set and sink in node_between_set
        }
    return gfa_graph
