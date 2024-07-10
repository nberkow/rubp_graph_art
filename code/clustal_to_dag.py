import sys
import networkx as nx
import matplotlib.pyplot as plt

class clustal_dag:

    def __init__(self):

        self.node_subset = set()
        self.alignments = {}
        self.total_length = -1

        self.min_aln_len = 10
        self.DG = nx.MultiDiGraph()
        self.DG.add_node("root", **{"member_sequences" : set()})
        
    def load_node_subset(self, subset_file):

        """
        optionally provide a file with one node name per line

        this subset of nodes is used to build the graph. Other aln inputs are ignored

        if ommitted the subset is the full set
        """

        with open(subset_file) as sfile:
            for x in sfile:
                self.node_subset.add(x.rstrip())

    def parse_alignment_file(self, clustal_file):

        unfiltered_aln = {}

        with open(clustal_file) as aln_file:
            for aln in aln_file:
                elements = aln.split("      ")
                if len(elements) == 2:
                    seq_name = elements[0]
                    if seq_name in self.node_subset or len(self.node_subset) == 0:
                        align_str = elements[1].rstrip()
                        if not seq_name in unfiltered_aln:
                            unfiltered_aln[seq_name] = ""
                            self.DG.nodes['root']["member_sequences"].add(seq_name)
                        unfiltered_aln[seq_name] = f"{unfiltered_aln[seq_name]}{align_str}"

        aln_by_len = {}
        for s in unfiltered_aln:
            aln_len = len(unfiltered_aln[s])
            if aln_len >= self.min_aln_len:
                if not aln_len in aln_by_len:
                    aln_by_len[aln_len] = []
                aln_by_len[aln_len].append(s)

        filtered_names = []
        for n in aln_by_len:
            
            aln_count = len(aln_by_len[n])
            if aln_count > len(filtered_names):
                filtered_names = aln_by_len[n]
                self.total_length = n

        for n in filtered_names:
            self.alignments[n] = unfiltered_aln[n]        

    def make_graph(self):
    
        seq_names = list(self.alignments.keys())
        current_nodes_by_member_tuple = {
            tuple(seq_names) : self.DG['root']
        }

        for n in range(self.total_length):
            aln_by_residue = {}
            for s in seq_names:
                aa_residue = self.alignments[s][n]

                # don't store gaps as nodes
                if aa_residue != "-":
                    if not aa_residue in aln_by_residue:
                        aln_by_residue[aa_residue] = []
                    aln_by_residue[aa_residue].append(s)

            if self.check_member_set_identical(aln_by_residue, current_nodes_by_member_tuple):
                pass
            else:
                current_nodes_by_member_tuple = self.update_graph(aln_by_residue, n)

    def update_graph(self, aln_by_residue, n):

        updated_nodes_by_member_tuple = {}

        for a in aln_by_residue:

            seq_names = aln_by_residue[a]
            seq_names.sort()

            new_node_name = f"{n}{a}"
            attrib = {}
            attrib['member_sequences'] = set(seq_names)
            self.DG.add_node(new_node_name, **attrib)
            updated_nodes_by_member_tuple[tuple(seq_names)] = self.DG[new_node_name]

            # link the newly made node to upstream nodes that share sequences
            for graph_node in self.DG.nodes:

                # skip self edges
                if graph_node != new_node_name:

                    # find the common sequences between the new node and an existing graph node
                    shared_sequences = self.DG.nodes[graph_node]["member_sequences"].intersection(set(seq_names))

                    # if the sequence is alread represented by an out edge, nothing to do.
                    # otherwise create an edge from the old node to the new node with the sequence name as an attrib
                    graph_node_out_sequences = set() 
                    for out_edge in self.DG.out_edges(graph_node, data=True):
                        seq_name = out_edge[2]['seq_name']
                        graph_node_out_sequences.add(seq_name)

                    shared_sequences_to_add = shared_sequences - graph_node_out_sequences
                    for seq_name in shared_sequences_to_add:
                        if seq_name in seq_names:
                            self.DG.add_edge(graph_node, new_node_name, seq_name=seq_name)

        return updated_nodes_by_member_tuple

    def check_member_set_identical(self, aln_by_residue, current_nodes_by_member_tuple):

        """
        for a given position in the alignment, check if the set of sequences included is
        the same as the previous position

        input: 
            aln_by_residue - dict with input seq names indexed by the aa residue at a single position
            e.g. {
                'E' : ('seq1', 'seq3'),
                'M' : ('seq2', 'seq4', 'seq5)
            }
            current_node_by_member_tuple - dict with seq names tuples as keys, graph_nodes as values

        output:
            boolean - does the value set in aln_by_residue perfectly match the key 
            set in current_aln_by_member_tuple?
        
        """

        # sequence grouping by current residue
        abr_set = set()
        for m in aln_by_residue.values():
            m.sort()
            abr_set.add(tuple(m))

        # sequence grouping by previous residue
        cbr_set = set()
        for m in current_nodes_by_member_tuple:
            cbr_set.add(tuple(m))

        return(abr_set == cbr_set)

if __name__ == "__main__":
    cdag = clustal_dag()
    #cdag.load_node_subset("/Users/NBerkowitz/Projects/rubp_graph_art/test_subset.txt")

    cdag.parse_alignment_file(sys.argv[1])
    cdag.make_graph()

    for layer, nodes in enumerate(nx.topological_generations(cdag.DG)):
        # `multipartite_layout` expects the layer as a node attribute, so add the
        # numeric layer value as a node attribute
        for node in nodes:
            cdag.DG.nodes[node]["layer"] = layer


    pos = nx.multipartite_layout(cdag.DG, subset_key="layer")
    last_x = pos['root'][0]
    shift_y = 0.02
    
    for p in pos:
        if pos[p][0] != last_x:
            shift_y = 0 - shift_y
            last_x = pos[p][0]
        pos[p][1] += shift_y

    fig, ax = plt.subplots()
    nx.draw_networkx(cdag.DG, pos=pos, ax=ax)
    fig.tight_layout()
    plt.savefig("/Users/NBerkowitz/Projects/rubp_graph_art/test_dag_multip.png")

    """
    nx.draw(cdag.DG)
    kwargs = {}
    plt.savefig("/Users/NBerkowitz/Projects/rubp_graph_art/test_dag.png")
    plt.clf()

    nx.draw_circular(cdag.DG, with_labels=True, **kwargs)
    plt.savefig("/Users/NBerkowitz/Projects/rubp_graph_art/test_dag_circ.png")
    plt.clf()

    nx.draw_kamada_kawai(cdag.DG, with_labels=True, **kwargs)
    plt.savefig("/Users/NBerkowitz/Projects/rubp_graph_art/test_dag_kk.png")
    plt.clf()

    nx.draw_planar(cdag.DG, with_labels=True, **kwargs)
    plt.savefig("/Users/NBerkowitz/Projects/rubp_graph_art/test_dag_plan.png")
    plt.clf()

    nx.draw_spring(cdag.DG, with_labels=True, **kwargs)
    plt.savefig("/Users/NBerkowitz/Projects/rubp_graph_art/test_dag_spring.png")
    plt.clf()"""





            
