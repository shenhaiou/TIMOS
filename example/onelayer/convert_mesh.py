#!/usr/bin/env python3
"""
Convert a LightPropagation .mesh file to TIMOS format.

The two formats share the same structure:
  Line 1:          numNodes
  Line 2:          numElems
  Next numNodes:   x y z
  Next numElems:   n1 n2 n3 n4 medIdx

The only difference is that LightPropagation medium indices are arbitrary
integers (e.g. 3), while TIMOS expects 1-based sequential indices.
This script remaps them in order of first appearance: first seen -> 1,
second seen -> 2, etc.

Usage:
  python3 convert_mesh.py <input.mesh> <output.mesh>
"""

import sys


def convert(src: str, dst: str) -> None:
    with open(src) as f:
        lines = f.read().split()

    idx = 0

    def read_int():
        nonlocal idx
        v = int(lines[idx]); idx += 1
        return v

    def read_float():
        nonlocal idx
        v = float(lines[idx]); idx += 1
        return v

    num_nodes = read_int()
    num_elems = read_int()

    nodes = []
    for _ in range(num_nodes):
        x, y, z = read_float(), read_float(), read_float()
        nodes.append((x, y, z))

    elems = []
    for _ in range(num_elems):
        n1, n2, n3, n4 = read_int(), read_int(), read_int(), read_int()
        med = read_int()
        elems.append((n1, n2, n3, n4, med))

    # Remap medium indices to 1-based sequential in order of first appearance.
    med_map: dict[int, int] = {}
    for _, _, _, _, med in elems:
        if med not in med_map:
            med_map[med] = len(med_map) + 1

    with open(dst, 'w') as f:
        f.write(f"{num_nodes}\n")
        f.write(f"{num_elems}\n")
        for x, y, z in nodes:
            f.write(f"{x} {y} {z}\n")
        for n1, n2, n3, n4, med in elems:
            f.write(f"{n1} {n2} {n3} {n4} {med_map[med]}\n")

    unique = sorted(med_map.items(), key=lambda kv: kv[1])
    print(f"Converted {src} -> {dst}")
    print(f"  Nodes: {num_nodes}, Elements: {num_elems}")
    print(f"  Medium index mapping: { {k: v for k, v in unique} }")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.mesh> <output.mesh>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
