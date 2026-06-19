# Graph Construction

## Seeding

`PoaGraph::new(seed, config)` builds an initial linear chain: one node per base in the seed
read, with each node connected to the next by an edge of weight 1 and coverage 1. This chain
is the starting backbone; every subsequent read is aligned into it.

The choice of seed matters. A seed that is much longer or shorter than the read set creates
unnecessary early branches and inflates the band width the aligner needs to track them. A
median-length spanning read is the right default (see [Seed Selection](../library/seed-selection.md)).

## Adding a read

`PoaGraph::add_read(read)` runs three steps:

### 1. Topological sort (Kahn's algorithm)

The graph is sorted into a linear ordering where every node appears before all nodes it has
an outgoing edge to. This is required by the DP: the aligner processes nodes in topological
order so that all predecessors of a node have been scored before the node itself is scored.

The sort runs in O(V + E) and is cached between reads via the **stale spine** mechanism
described in [Banded DP Alignment](banded-dp.md).

### 2. DP alignment

The read is aligned against the sorted graph using banded affine-gap DP. The alignment
returns an `AlignOp` sequence: `Match(node)`, `Insert(base)`, or `Delete(node)`.

### 3. Graph update (`add_to_graph`)

The traceback is replayed to update the graph:

| Op | Action |
|---|---|
| `Match(node)` | Increment `node.coverage` and the incoming edge weight |
| `Insert(base)` | Create a new node; connect it to the previous and next nodes in the traceback |
| `Delete(node)` | Traverse `node` without incrementing its coverage; increment only the traversal edge weight |

**Insert nodes** are allocated with coverage 1 (the current read is the first to traverse
them). Subsequent reads that match the same insert base will merge into the existing node
rather than creating a new one, because the DP naturally finds the highest-scoring path and
a match to an existing node scores +1 while a new insert scores 0 (open + extend penalty).

**Delete traversals** increment the outgoing edge weight at each skipped node, preserving
the read count on the traversal path for the heaviest-path DP, but do not count as evidence
that the node's base is correct. This is the key design choice that makes boundary trim work:
nodes that are skipped by the majority of reads have low `coverage` even if many reads
traversed the graph position.

## Graph growth over a read set

A typical graph on a clean 20-read set of 150 bp STR reads:

- Seed: 150 nodes, 149 edges
- After 20 reads (5% substitution error): ~165 nodes (15 error insert nodes), ~185 edges
- Heaviest path: 150 nodes (the true consensus)

Error insert nodes sit off the spine with edge weight 1 (or 2 if two reads share the same
error). The `(weight - 1)` normalisation in the heaviest path DP gives them a negative
contribution, keeping the heaviest path on the true spine.
