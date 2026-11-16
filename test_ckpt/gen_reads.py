#!/usr/bin/env python3
"""Deterministic read simulator for checkpoint testing.
mode=hifi: ~0.8% substitution errors, 8-14kb reads -> .fa
mode=ont : ~10% errors (indel-biased), 6-20kb reads   -> .fq
Generates a pseudo-genome with segmental duplications so that
overlaps are non-trivial.
"""
import random, sys

random.seed(42)
BASES = "ACGT"

def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))

# usage: gen_reads.py out.fa [genome_bases=1200000] [coverage=25] [mode=hifi]
OUT  = sys.argv[1] if len(sys.argv) > 1 else "reads.fa"
G    = int(sys.argv[2]) if len(sys.argv) > 2 else 1_200_000
COV  = int(sys.argv[3]) if len(sys.argv) > 3 else 25
MODE = sys.argv[4] if len(sys.argv) > 4 else "hifi"

genome = rand_seq(G)
# two duplicated segments (forces non-trivial overlaps)
dup1 = genome[G//12 : G//12 + G//40]
dup2 = genome[G//2  : G//2  + G//25]
genome = genome + dup1 + dup2

if MODE == "ont":
    RL_MIN, RL_MAX = 6000, 20000
    P_SUB, P_INS, P_DEL = 0.04, 0.035, 0.035   # ~11% total, indel-heavy
else:
    RL_MIN, RL_MAX = 8000, 14000
    P_SUB, P_INS, P_DEL = 0.008, 0.0, 0.0

out = open(OUT, "w")
n_bases = int(len(genome) * COV)
rid, emitted = 0, 0
while emitted < n_bases:
    rl = random.randint(RL_MIN, RL_MAX)
    st = random.randint(0, len(genome) - rl)
    seq = list(genome[st:st+rl])
    k = 0
    while k < len(seq):
        r = random.random()
        if r < P_SUB:
            seq[k] = random.choice([b for b in BASES if b != seq[k]])
        elif r < P_SUB + P_INS:
            seq[k] = seq[k] + random.choice(BASES)
        elif r < P_SUB + P_INS + P_DEL:
            seq[k] = ""
        k += 1
    s = "".join(seq)
    if MODE == "ont":
        out.write(f"@r{rid}\n{s}\n+\n{'I'*len(s)}\n")
    else:
        out.write(f">r{rid}\n{s}\n")
    emitted += rl
    rid += 1
out.close()
print(f"{rid} reads ({MODE}), {emitted} bases, genome {len(genome)}", file=sys.stderr)
