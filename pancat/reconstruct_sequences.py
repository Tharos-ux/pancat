'Extract sequences given a graph'
from Bio import SeqIO
from tharospytools.bio_tools import revcomp
from pgGraphs import Graph


def reconstruct_paths(gfa_file: str, output: str, selected_paths: list | bool = False, split: bool = True) -> dict:
    """Given a GFA file with paths, follows those paths to reconstruct the linear sequences

    Args:
        gfa_file (str): gfa graph file, with paths
        gfa_version (str): the version of the file

    Raises:
        ValueError: if no path in graph

    Returns:
        dict: mapping between path name and sequence
    """
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        with_sequence=True
    )
    if len(gfa_graph.paths) > 0:
        genomes: dict = gfa_graph.reconstruct_sequences()
        selected_paths: list = (selected_paths, sorted(list(gfa_graph.paths.keys())))[
            isinstance(selected_paths, bool)]
        for i, label in enumerate(selected_paths):
            with open(f'{output}_{i}.fasta' if split else f'{output}.fasta', 'w' if split else 'a', encoding="utf-8") as writer:
                writer.write(
                    f">{label}\n{''.join(genomes[label])}\n"
                )
                del genomes[label]
    else:
        raise ValueError(
            "Graph has no embed paths. Could'nt determine what to reconstruct!"
        )


def graph_against_fasta(gfa_graph: str, pipeline_txt: str) -> bool:
    with open(pipeline_txt, 'r', encoding='utf-8') as pipeline:
        lines = pipeline.readlines()
        reconstruction: dict = reconstruct_paths(gfa_graph)
        complete_pangenome_graph: list[bool] = [
            False for _ in range(len(lines))]
        for y, line in enumerate(lines):
            # We assume the pipeline file is in the cactus format
            seq_name, seq_path = line.strip().split('\t')
            if seq_name in reconstruction:
                # Do comparison
                # We load the file
                with open(seq_path, 'r', encoding='utf-8') as reader:
                    for record in SeqIO.parse(reader, "fasta"):
                        seq1: str = record.seq.lower()
                        seq2: str = ''.join(reconstruction[seq_name]).lower()
                        if len(seq1) != len(seq2):
                            print(
                                f"[{seq_name}] Length of sequence with header {record.id} does not match the path {seq_name} from a size-perspective (resp. {len(seq1)} vs {len(seq2)}).")
                        else:
                            densities: list = [0 for _ in range(101)]
                            for i in range(len(seq1)):
                                if seq1[i] != seq2[i]:
                                    densities[round((i/len(seq1))*100)] += 1
                            if sum(densities) != 0:
                                print(
                                    f"[{seq_name}] Sequence {record.id} and path have the same length ({len(seq1)}). {round((sum(densities)/len(seq1))*100,2)}% of base pairs are not shared.")
                            else:
                                print(
                                    f"[{seq_name}] Sequence {record.id} and path represents the same sequence.")
                                complete_pangenome_graph[y] = True
            else:
                print(
                    f"[{seq_name}] Can't find {seq_name} in the paths of the graph.")
        if all(complete_pangenome_graph):
            print(f"{gfa_graph} is a complete pangenome graph")
        else:
            print(f"{gfa_graph} is not a complete pangenome graph")
    return all(complete_pangenome_graph)
