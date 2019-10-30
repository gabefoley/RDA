"""
evodiv.py Module to explore diversity in sequences in an evolutionary branch
By Mikael Boden
"""
import sequence
import phylo
import prob
import sys

def saveConsensus(aln, theta2 = 0.01, countgaps = False, filename = None):
    """ Display a table with rows for each alignment column, showing
        symbols in order of decreasing probability.
        theta2 is the percent threshold (0.01 is 1 percent) for inclusion (symbols below are ignored).
        countgaps, if true, count gaps (default false).
        filename is name of file to save the output to (default stdout)."""
    if filename == None:
        f = sys.stdout
    else: # assume filename is a textstring, which is first cleaned from strange characters
        filename = ''.join(e for e in filename if e.isalnum() or e == '_' or e == '.')
        f = open(filename, 'w')
    for col in range(aln.alignlen):
        # collect probabilities for column, with or without gap
        myalpha = aln.alphabet
        if countgaps:
            alist = list(aln.alphabet)
            alist.append('-')
            myalpha = sequence.Alphabet(alist)
        d = prob.Distrib(myalpha)
        for seq in aln.seqs:
            if seq[col] in myalpha:
                d.observe(seq[col])
        symprobs = d.getProbsort() # the symbols sorted by probability
        ninclusions = 0
        for (s, p) in symprobs:
            if p >= theta2:
                ninclusions += 1
            else:
                break
        if ninclusions > 1:
            f.write("%d:" % (col + 1))
            for (s, p) in symprobs:
                if p >= theta2:
                    f.write("%c%02d" % (s, int((p*100))))
            f.write(" ")
    f.write('\n')
    f.close()

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: evodiv <tree> <alignment> <nodes> where', \
        "\n\t<tree> is a completely labelled Newick file of a phylogenetic tree", \
        "\n\t<alignment> is a FASTA or Clustal file with a sequence for each label in tree", \
        "\n\t<nodes> is a FASTA file with a sequence entry for each label for which a variability is determined", \
        "\n\tVariability is saved to file with the sequence entry's name .txt")
        sys.exit(1)
    tree = phylo.readNewick(sys.argv[1])
    try:
        seqs = sequence.readFastaFile(sys.argv[2], alphabet=sequence.Protein_wX, gappy=True)
        aln = sequence.Alignment(seqs)
    except:
        aln = sequence.readClustalFile(sys.argv[2], sequence.Protein_Alphabet_wX)
    tree.putAlignment(aln)
    select = sequence.readFastaFile(sys.argv[3], alphabet=sequence.Protein_wX, gappy=True)
    nodes = []
    for selected in select:
        nodename = selected.name
        nodes.append(nodename)
    for nodename in nodes:
        node = tree.findLabel(nodename)
        if node == None:
            print('Warning: ', nodename, 'is not found in tree')
            break
        if node.sequence == None:
            print('Warning: ', nodename, 'has not got a sequence')
        all_children = tree.getDescendantsOf(node, transitive=True)
        if all_children == None:
            print('Warning: ', nodename, 'has no children')
        direct_ancestors = tree.getAncestorsOf(node, transitive=True)
        if direct_ancestors == None:
            print('Warning: ', nodename, 'has no ancestors')
            direct_ancestors = []
        relevant = [node.sequence]
        for child in all_children:
            if not child.sequence:
                pass
                #print 'Warning: ', nodename, 'has a child', child.label, 'with a sequence which is None'
            else:
                allgaps = True
                for pos in child.sequence:
                    if pos != '-':
                        allgaps = False
                        break
                if not allgaps:
                    relevant.append(child.sequence)
        for parent in direct_ancestors:
            if not parent.sequence:
                print('Warning: ', nodename, 'has an ancestor', parent.label, 'with a sequence which is None')
            else:
                relevant.append(parent.sequence)
        relevant_aln = sequence.Alignment(relevant)
        saveConsensus(relevant_aln, countgaps = True, filename = nodename + ".txt")
        sequence.writeFastaFile(nodename + ".fa", relevant_aln.seqs)