# 010 · Where NR5A2 acts — an assumption, tested rather than chosen

*Date:* 2026-08-15 · *Status:* accepted

## Question

NR5A2 is in the intended payload but had **no representation anywhere in the
model** (gap B5). Stage 3 must rank it against E47 and MIST1 for the u₃ slot,
and a candidate with no representation scores zero by construction. Where does
it act?

## Positions considered

**Position A — enhancer co-activator.** NR5A2 modulates `alpha_auto`, boosting
the 2.3-kb autoregulatory enhancer. This is the natural reading of "NR5A2 helps
rebuild acinar identity", and it puts NR5A2 directly into the loop that the
whole model is about.

Against: **the evidence does not show this.** Holmstrom 2011 (PMID 21852532,
GSE34295) shows LRH-1/NR5A2 and PTF1-L co-regulating an **exocrine
transcriptional network** — co-occupancy on digestive-enzyme genes. It does
**not** demonstrate binding at the *Ptf1a* autoregulatory enhancer. Position A
is an extrapolation from "co-regulates acinar genes" to "co-activates the master
regulator's own enhancer", and those are different claims.

**Position B — acinar-output co-activator.** NR5A2 acts on `A`, the digestive-
enzyme output, which is what the data actually shows. Conservative and directly
supported.

Against: if NR5A2 only boosts output downstream of the loop, it cannot help
*re-close* the loop, and it would be a poor payload component — which may be
true, but should be a finding rather than an assumption.

*Where they actually disagree:* on whether to encode the mechanism we hope is
true or the one the data supports. Both are guesses about an unmeasured link.

## Decision

**Neither. Both, as competing topologies — T6a and T6b.**

```
T6a  nr5a2_mode="enhancer"   enhancer *= (1 + kappa·NR5A2)     Position A
T6b  nr5a2_mode="output"     output_gain += kappa·NR5A2        Position B
```

Composable topologies (decision 003) make an extra variant nearly free, so
there is no reason to pick. **The Q-values decide.** Under a five-copied-files
architecture this would have cost a sixth file and the assumption would simply
have gone untested — which is a concrete example of decision 003 paying for
itself immediately.

**The assumption is labelled as an assumption wherever it appears** — in the
master plan's topology table, in `model.py` at the point of use, and here. A
placement that is presented as a finding when it is an extrapolation is exactly
the kind of thing a hostile reader finds and a careful one flags first.

This also converts the u₃ slot into a genuine three-way mechanistic comparison
rather than a relabelling: **E47 helps by relieving titration, NR5A2 by boosting
transcription, MIST1 by raising secretory capacity.** Three different mechanisms
competing for one payload slot is a real ranking question.

**Partial empirical check, qualitative only.** GSE34295 can answer *"is there
LRH-1 signal at the Ptf1a locus, yes or no"*. It is 2 samples, single-replicate,
2011-era ChIP-seq with no biological replication, so it can support presence or
absence and **nothing more** — no number may be reported from it, and it cannot
set `kappa`. Added as a task in the optional sequence module. A "no" would not
kill T6a outright (absence of evidence in one weak dataset), but it would mean
T6a must win on Q-value alone with its assumption stated plainly.

## What would reverse this

1. **If T6a and T6b give indistinguishable Q-values**, the placement does not
   matter at the resolution the data supports — report that, and stop treating
   it as an open question.
2. **If GSE34295 shows clear LRH-1 signal at the *Ptf1a* locus**, T6a graduates
   from assumption to supported, and the caveat can be dropped.
3. **If NR5A2 ranks last in Stage 3 under *both* placements**, it should come
   out of the payload — and that is a real result about a component Luqmaan
   currently intends to test at the bench, which is precisely the kind of
   finding the two-arm structure (decision 001) exists to make possible.
