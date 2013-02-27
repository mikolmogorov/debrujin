
import sys

import build as build1
import build2

if __name__ == "__main__":
    try:
        seq_path, k = sys.argv[1:]
        k = int(k)
    except ValueError:
        print("Usage: check.py SEQ_PATH DOT_PATH")
    else:
        seq = build1.read_fasta(seq_path)
        g1 = build1.build_compressed_graph(seq, k)
        g2 = build2.build_compressed_graph(seq, k)
        print set(g1.values())
        print set(g2.values())
        assert g1 == g2
