#!/usr/bin/env python3
"""Deterministic HiFi-like read simulator for checkpoint testing.
Generates a ~1.2Mb pseudo-genome with segmental duplications, then
samples ~25x coverage of 8-14kb reads with ~0.8% substitution errors.
"""
import random, sys

random.seed(42)
BASES = "ACGT"

def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))

# usage: gen_reads.py out.fa [genome_bases=1200000] [coverage=25]
G   = int(sys.argv[2]) if len(sys.argv) > 2 else 1_200_000
COV = int(sys.argv[3]) if len(sys.argv) > 3 else 25
genome = rand_seq(G)
# two duplicated segments (forces non-trivial overlaps)
dup1 = genome[G//12 : G//12 + G//40]
dup2 = genome[G//2  : G//2  + G//25]
genome = genome + dup1 + dup2

out = open(sys.argv[1] if len(sys.argv) > 1 else "reads.fa", "w")
RL_MIN, RL_MAX, ERR = 8000, 14000, 0.008
n_bases = int(len(genome) * COV)
rid, emitted = 0, 0
while emitted < n_bases:
    rl = random.randint(RL_MIN, RL_MAX)
    st = random.randint(0, len(genome) - rl)
    seq = list(genome[st:st+rl])
    for k in range(rl):
        if random.random() < ERR:
            seq[k] = random.choice([b for b in BASES if b != seq[k]])
    out.write(f">r{rid}\n{''.join(seq)}\n")
    emitted += rl
    rid += 1
out.close()
print(f"{rid} reads, {emitted} bases, genome {len(genome)}", file=sys.stderr)
