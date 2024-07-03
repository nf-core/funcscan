import sys

f = sys.argv[1]

longest_contig = 0
shortest_contig = 0

with open(f) as f_handle:
    for line in f_handle:
        if not line.startswith(">"):
            c_len = len(line)
            print(c_len)
            if c_len > longest_contig:
                longest_contig = c_len
            if shortest_contig == 0:
                shortest_contig = c_len
            elif c_len < shortest_contig:
                shortest_contig = c_len

print("Shortest contig: " + str(shortest_contig))
print("Longest contig: " + str(longest_contig))