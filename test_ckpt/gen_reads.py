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

# genome: 1.2Mb unique + two duplicated segments (forces non-trivial overlaps)
G = 1_200_000
genome = rand_seq(G)
dup1 = genome[100_000:130_000]
dup2 = genome[500_000:560_000]
genome = genome + dup1 + dup2          # ~1.29Mb total

out = open(sys.argv[1] if len(sys.argv) > 1 else "reads.fa", "w")
COV, RL_MIN, RL_MAX, ERR = 25, 8000, 14000, 0.008
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
