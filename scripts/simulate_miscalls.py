"""Randomly mutate k% of read positions in a paired FASTQ (uniformly to a different
base) to simulate basecalling/miscall errors. Vectorized with numpy for speed; N and
other non-ACGT symbols are left untouched. Output is plain (uncompressed) FASTQ."""
import sys
import gzip
import numpy as np

R1_IN, R2_IN, K_PERCENT, R1_OUT, R2_OUT, SEED = sys.argv[1:7]
RATE = float(K_PERCENT) / 100.0
SEED = int(SEED)
BATCH = 50000

BASE_TO_IDX = np.full(256, -1, dtype=np.int8)
for i, b in enumerate(b"ACGT"):
    BASE_TO_IDX[b] = i
IDX_TO_BASE = np.frombuffer(b"ACGT", dtype=np.uint8)


def _opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def mutate_batch(seqs, rng):
    lens = np.array([len(s) for s in seqs])
    maxlen = int(lens.max())
    n = len(seqs)
    arr = np.zeros((n, maxlen), dtype=np.uint8)
    for i, s in enumerate(seqs):
        arr[i, :len(s)] = np.frombuffer(s.encode(), dtype=np.uint8)
    idx = BASE_TO_IDX[arr]
    valid = idx >= 0
    pos_mask = np.arange(maxlen)[None, :] < lens[:, None]
    rand = rng.random((n, maxlen))
    mutate_mask = (rand < RATE) & valid & pos_mask
    shift = rng.integers(1, 4, size=(n, maxlen))
    new_idx = (idx + shift) % 4
    new_base = IDX_TO_BASE[new_idx]
    arr_out = np.where(mutate_mask, new_base, arr)
    return [arr_out[i, :lens[i]].tobytes().decode() for i in range(n)]


def process(in_path, out_path, rng):
    n_mutated = 0
    n_bases = 0
    with _opener(in_path) as fin, open(out_path, "wt") as fout:
        headers, seqs, pluses, quals = [], [], [], []

        def flush():
            nonlocal n_mutated, n_bases
            if not seqs:
                return
            mutated = mutate_batch(seqs, rng)
            for h, orig, m, p, q in zip(headers, seqs, mutated, pluses, quals):
                fout.write(h)
                fout.write(m + "\n")
                fout.write(p)
                fout.write(q)
                n_mutated += sum(a != b for a, b in zip(orig, m))
                n_bases += len(orig)
            headers.clear(); seqs.clear(); pluses.clear(); quals.clear()

        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline().rstrip("\n")
            p = fin.readline()
            q = fin.readline()
            headers.append(h); seqs.append(s); pluses.append(p); quals.append(q)
            if len(seqs) >= BATCH:
                flush()
        flush()
    return n_mutated, n_bases


rng1 = np.random.default_rng(SEED)
rng2 = np.random.default_rng(SEED + 1)
m1, b1 = process(R1_IN, R1_OUT, rng1)
m2, b2 = process(R2_IN, R2_OUT, rng2)
print(f"k={K_PERCENT}% rate={RATE} seed={SEED} "
      f"R1: {m1}/{b1} bases mutated ({m1 / b1:.5%}) "
      f"R2: {m2}/{b2} bases mutated ({m2 / b2:.5%})", file=sys.stderr)
