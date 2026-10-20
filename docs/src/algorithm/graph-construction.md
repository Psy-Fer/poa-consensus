# Graph Construction

## Seeding

The engine builds an initial linear chain from a seed read: one node per base, with each
node connected to the next by an edge of weight 1 and coverage 1. This chain is the starting
backbone; every subsequent read is aligned into it.

The choice of seed matters, so `consensus()` self-seeds on the median-length read internally
rather than trusting the caller's index. A seed that is much longer or shorter than the read
set creates unnecessary early branches and biases periodic repeats toward its length; a
median-length read avoids that (see [Seed Selection](../library/seed-selection.md)).

## Adding a read

Each subsequent read is folded into the graph in three steps:

### 1. Topological sort (Kahn's algorithm)

The graph is sorted into a linear ordering where every node appears before all nodes it has
an outgoing edge to. This is required by the DP: the aligner processes nodes in topological
order so that all predecessors of a node have been scored before the node itself is scored.

The sort runs in O(V + E) and is recomputed for each read that is added.

### 2. DP alignment

The read is aligned against the sorted graph using banded affine-gap DP (see
[Banded DP Alignment](banded-dp.md)). The alignment returns a sequence of operations:
match (the read base is aligned to a node, whether or not the bases agree), insert (a read
base with no graph node), or delete (a graph node the read skips).

### 3. Graph update

The traceback is replayed to update the graph:

| Op | Action |
|---|---|
| Match, same base | Increment the node's `cov`; increment the edge from the read's previous node |
| Match, different base | Node fusion: reuse the sibling in that column with the read's base, or create a new substitution node linked as an `aligned` sibling; increment its `cov` and the incoming edge |
| Insert | Create a new node with `cov` 1; connect it from the read's previous node |
| Delete | Increment the skipped node's `del` count only |

**Node fusion** (as in SPOA and abPOA) keeps all alternative bases at a position as
siblings in one alignment column. Later reads carrying the same alternative base reuse the
sibling instead of fragmenting into fresh nodes, so allele counts accumulate correctly and
the majority-frequency consensus can read each column directly.

**Insert nodes** start with coverage 1. A later read with the same inserted base usually
merges into the existing node, because matching an existing node scores positively while
opening a fresh insert pays the gap penalty.

**Deletes** touch the skipped node only through its `del` count. When the read rejoins the
graph after a delete run, the edge it takes gets its weight incremented like any other edge,
but that adjacency is left out of the read's *matched* adjacencies (`read_matched_edges`).
Bubble reporting uses the matched view, so a length variant's short-allele bypass does not
look like a competing arm; multi-allele phasing uses the full `read_paths`, where the bypass
is exactly the shorter allele's arm. A deleted node gains no `cov`, so nodes skipped by most
reads keep low coverage and boundary trim can remove them.

## Example: two-allele SNV graph

After building a graph from 20 reads split equally between two alleles with a SNV at
position 11 (A vs G), the graph has a two-node bubble. Both arms are supported by 10
reads each (plus one noisy read creating a third minor arm). The spine goes through the
first allele's arm; the second allele's arm is shown in grey.

![Two-allele SNV bubble: spine in blue, second allele arm in grey](../diagrams/poa_network.svg)

The same graph with a specific allele-B read overlaid (orange nodes are the path that read
takes):

![Same graph with allele-B read path highlighted in orange](../diagrams/poa_network_with_read.svg)

## Graph growth over a read set

Measured on 20 reads of a random 150 bp sequence with 5% substitution error (default
configuration):

- Seed: 150 nodes, 149 edges
- After 20 reads: 280 nodes, 403 edges. Almost all of the extra nodes are substitution
  siblings created by node fusion (about 150 substitution errors in total, spread over 150
  columns with three possible alternative bases each).
- Consensus: 150 bp, identical to the truth

Error siblings sit off the spine with edge weight 1 (or 2 if two reads share the same
error). The heaviest path always follows the largest out-edge weight, so it stays on the
majority base at every column.
