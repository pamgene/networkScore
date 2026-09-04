# networkScore

Computes golden score (paired kinase+sensitivity) and kinase-only score for
a condition, via permutation testing against networks built by `networkGen`.
See the suite-wide [Context Map](../CONTEXT-MAP.md) for shared vocabulary
(network-result object, condition, PPI network) and how this package relates
to the others.

## Language

**Golden score**:
The permutation-test significance score for a condition using paired
kinase-activity + sensitivity data. Not the same "golden" as any per-node
value — it describes the paired-data scoring path specifically, as opposed
to kinase-only.
_Avoid_: network score (ambiguous with the topology statistic it's derived
from)

**Kinase-only score**:
The same permutation-test significance scoring, using only kinase-activity
data (no paired sensitivity data).

**Significance score**:
The p-value-like output of a permutation test: the fraction of permuted
networks whose statistic was at least as extreme as the observed network's.
Low = the observed signal is unlikely to arise by chance. This is what
"score" means everywhere in this package — the per-node percentile value
from `networkGen` is called `percentile_score` specifically to avoid
colliding with this term.
_Avoid_: p-value (close but not statistically a true p-value), score (too
ambiguous outside this package)

**Observed network / observed build**:
The network built from the real, unshuffled input data for a condition —
compared against that condition's permuted networks to compute its
significance score.

**Permutation / permuted build**:
A network built from a shuffled version of the input data (kinase-activity
values randomly reassigned across kinases via `sample()`), used to build a
null distribution for the significance score. A condition's score is
computed from its observed build plus (typically 50) permuted builds.
_Avoid_: random network, null network

**Task list**:
The full flat list of build requests — one observed + all permutations, for
every condition being scored in one run — submitted to
`networkGen::generate_networks_batch()` in a single call. Each task is
tagged with its condition, role (observed/permutation), and permutation
index, so results can be regrouped by condition afterward.
_Avoid_: batch (that's networkGen's term for the mechanism; this is the data
networkScore builds to feed it)
