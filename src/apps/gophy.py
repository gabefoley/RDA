import godata
import phylo
import sys
import annotations
import os

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Usage: gophy <obo-file> <gaf-file> <nwk-file> <p-value-threshold> where', \
        "\n\t<obo-file> is the Gene Ontology OBO-file (or other OBO-file to define ontology of terms)",\
        "\n\t<gaf-file> is the Gene Annotation File (linking gene names to terms)",\
        "\n\t<nwk-file> is a Newick file defining the phylogenetic tree (with leaves labelled with gene names)",\
        "\n\t<p-value-threshold> is statistical significance level required for enrichment to be reported"
              "\n\t -o <outfile> is an optional filename to save the annotated tree to")
        sys.exit(1)
    pvalue = float(sys.argv[4])
    print('P-value is ', pvalue)
    go = godata.GO(sys.argv[2], sys.argv[1])

    # If it's a Newick file then read it in and create a fresh NexusAnnotations file
    if (sys.argv[3].split(".")[1] == "nwk"):
        tree = phylo.readNewick(sys.argv[3])
        nexus_annotations = annotations.NexusAnnotations(tree=tree)
        tree.putAnnotations(nexus_annotations)

    # Else if it's a Nexus file read it in and check if it has an existing Uniprot annotation
    elif (sys.argv[3].split(".")[1] == "nxs" or sys.argv[3].split(".")[1] == "nex"):
        tree = phylo.read_nexus(sys.argv[3])
        tree.swap_annotations("UNIPROT")
    else:
        raise RuntimeError("Tree file should end in .nwk, .nex, or .nxs")

    # Set the outpath to write the annotated tree to
    if "-o" in sys.argv:
        outfile = sys.argv[(sys.argv.index("-o") +1)]
        out_path = (os.getcwd() + "/" + outfile)

    # for every node
    nodes = tree.getNodes()
    background = []
    for node in nodes:
        if node.isLeaf():
            background.append(node.label)
    print('Number of nodes:', len(nodes))


    for node in nodes:
        members = node.getDescendants(transitive=True)
        members.insert(0, node)
        if len(members) > 1:
            cnt = 1
            foreground = []
            for member in members:
                if member.isLeaf():
                    foreground.append(member.label)
#                    print('\t\t', cnt, member)
                    cnt += 1
            report = go.getEnrichmentReport(foreground, background, threshold = pvalue)

            if len(report) > 0:
                text_annotations = []
                print('============================\n', node.label, '\t', str(node), '\n\thas', len(members),
                      'leaf members (incl itself)')
                for row in report:
                    print('\t\t\t', row)
                    text_annotations.append(row[4])

                tree.nexus_annotations.add_annotation(node=node, key="GO_Terms", annotation=text_annotations)
                tree.nexus_annotations.add_colour_annotation(node=node)


    # If we were using UNIPROT annotations from an annotated Nexus file, switch the tree back to using the original
    # annotations
    tree.swap_annotations("Original")


    # If user defined a file to write the tree to, save it there
    if out_path:
        tree.write_to_nexus(out_path, use_symbols=True)
        print(tree.nexus_annotations.annotation_symbols)
        with open(out_path.split(".")[0] +"_annotations.txt", "w+") as legends_file:
            for k, v in tree.nexus_annotations.annotation_symbols.items():
                legends_file.write(k + " : " + v + "\n")



    # tree.write_to_nexus(out_path_filtered, exclude_annotations=["PDB"], use_symbols=True)
