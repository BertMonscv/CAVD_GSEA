# AMENDMENT: CAVD GSEA Pre-Analysis Plan v1 → v2

**Amendment date:** 2026-05-26
**Author:** Zongyue Li (LZY)
**Supersedes:** v1.0 of the plan (`preregistration/CAVD_GSEA_PreAnalysisPlan.md`),
frozen 2026-05-23, git commit `ed2d35d`, tag `prereg-frozen`.
**Replaces with:** v2 of the plan (`pre/CAVD_GSEA_PreAnalysisPlan_v2.md`),
frozen 2026-05-26, git commit `c76f64489e94180e77741d8cbf0aca9a37134e9c`.

---

## 1. Status declaration

**No GSEA / fgsea call has been executed against any of the three cohorts
(GSE51472, GSE51472+GSE12644, GSE83453) between the v1 freeze on 2026-05-23
and this amendment on 2026-05-26.** No differential expression ranking has
been generated; no fgsea output exists locally or in any committed file.
Every change documented below is therefore a *pre-running correction*,
committed before any analysis execution.

## 2. Scope of v2 (and why it is not "v1 rewritten from scratch")

v2 is a minor revision of v1, not a wholesale replacement. The
hypothesis, the five gene sets, the three cohorts, the ranking metric
(`logFC`), the significance threshold, `set.seed(42)`, `minSize=8 /
maxSize=500`, the commitment to report all five results regardless of
significance, and the directional-consistency requirement are unchanged.
What changes is a small set of specifications that fall into three
honest categories:

- **methodologically necessary** (without them v1 either could not run
  or could produce misleading results) — changes #2, #3, #4a, #4c,
- **necessary for reproducibility and audit** but not affecting the
  scientific conclusion — changes #4b, #5,
- **clarifying / weakly necessary** — changes #1, #6.

A self-review of each change's honest necessity is given in Section 4.

v1 is retained verbatim at its original location
(`preregistration/CAVD_GSEA_PreAnalysisPlan.md`) and is not modified.

## 3. Summary of changes

| # | Change | Honest necessity |
|---|---|---|
| 1 | Ranking metric stated as "`logFC` column from limma `topTable()`" (was "log2FC") | Wording precision; near-zero real gain |
| 2 | Probe-to-gene collapse rule specified (retain probe with highest IQR per gene) | Methodologically necessary |
| 3 | HGNC symbol normalization added before gene set matching, with alias-resolution log committed | Methodologically necessary (especially for this study's gene set composition) |
| 4a | fgsea function changed from `fgsea(... nperm=10000)` to `fgseaMultilevel()` | Methodologically necessary (v1's parameter combination was not executable) |
| 4b | `set.seed(42)` called immediately before each fgsea invocation rather than once | Reproducibility, not scientific impact (3rd-decimal NES differences only) |
| 4c | BH multiple-testing pool split: 5 hypothesis-driven sets and 50 Hallmark sets corrected separately | Methodologically necessary (joint pool dilutes hypothesis-driven adjusted P from /5 to /55) |
| 5 | Analysis script `run_gsea_v2.R` committed at freeze; no edits after | Compliance and audit; not scientific impact |
| 6 | Gene set overlap check committed before freeze with explicit Jaccard > 0.30 threshold | Useful curation QC; threshold made explicit (see §4) |

## 4. Per-change rationale (with honest self-review)

The rationales below are reproduced from the author's self-review of
each change, retaining its honest framing — including, where applicable,
explicit acknowledgement of limited necessity.

### Change 1: ranking metric stated as `logFC` column

V1 wrote "log2FC", which is informally unambiguous because limma's
`topTable()` output has exactly one column named `logFC`. V2 names
that column explicitly. The real informational gain is close to zero;
v1 was already executable as written. This change is preserved because
precise wording does no harm, but it is *not* claimed as a substantive
correction. The v2 header describes this as "clarifies the ranking
metric," using "clarifies" (not "adds") for that reason.

### Change 2: probe-to-gene collapse rule (highest IQR)

GPL570 commonly assigns 2–10 probes to a single gene. Different
collapse rules can yield different `logFC` values for the same gene
(e.g. the same gene varying from 1.2 to 2.8 across plausible probe
selections; Miller et al. 2011 *PLoS One* analyzed this). V1 did not
specify a collapse rule, leaving this open. V2 specifies "retain
the probe with the highest interquartile range per gene," applied
identically across all three cohorts.

Honest caveat: "highest IQR" is not the uniquely correct answer.
Alternatives such as `max(mean(expression))` or jetset selection are
equally defensible. The choice of IQR follows the WGCNA tutorial
recommendation (Langfelder & Horvath). The point that matters for
pre-registration is *that a rule is committed in advance*, not which
specific rule it is.

### Change 3: HGNC symbol normalization

Several genes in the custom gene sets have widely used aliases:
ADRP / ADFP → PLIN2, LC3B → MAP1LC3B, p62 → SQSTM1, S6K1 → RPS6KB1,
OPG → TNFRSF11B. If GEO annotation in a cohort uses the alias and the
custom set uses the current symbol (or vice versa), that gene silently
drops out of the leading edge — no error is raised. V2 commits to
normalizing all symbols to current HGNC official names before
matching, and committing the alias-resolution log.

Honest framing: this is treated as "background preprocessing" in many
pre-registrations and is often not written down. For this study's
specific gene set composition (which contains several aliased genes
in the leptin / autophagy / lipophagy axes), making the rule explicit
has more real value than it would for a generic GSEA pre-registration.

### Change 4a: `fgseaMultilevel()` in place of `fgsea(... nperm=10000)`

V1 wrote `nperm=10000`, which is a parameter of the classical
`fgsea()` function. `fgseaMultilevel()` does not accept `nperm` and
the two functions are not interchangeable. V1's parameter combination
is therefore not executable as written; v2 picks one of the two
options explicitly.

The choice of `fgseaMultilevel()` (rather than fixing v1 to use
classical `fgsea()` with `nperm`) is driven by the precision of small
adjusted P values: pilot exploratory work for this study produced
adjusted P values on the order of 1e-4, which classical permutation
with `nperm=10000` cannot resolve (its floor is approximately 1e-4).

### Change 4b: per-call `set.seed(42)`

R's RNG carries state across function calls. With four fgsea
invocations under v2 (three cohorts × hypothesis-driven pool +
one Hallmark descriptive call), calls 2/3/4 would inherit
unpredictable RNG state if the seed is set only once at the script's
start.

Honest framing: for `fgseaMultilevel()`, which performs heavy internal
sampling, the practical effect of state inheritance is small —
typically affecting only the 3rd decimal of NES values. This change
is better understood as a *formal requirement for bit-identical
reproducibility* than as a correction with scientific impact. It is a
0/1 question for an external reproducer ("did your run match mine
exactly?"), not a question about NES interpretation. Distinct from
change 4a (which addresses non-executability), this change addresses
small numerical drift across reproductions.

### Change 4c: BH pool separation

Pooling the 5 hypothesis-driven gene sets together with the 50
Hallmark sets under a single Benjamini–Hochberg correction would
dilute the hypothesis-driven adjusted P values by changing the
denominator from /5 to /55. Hallmark P values, which serve a
descriptive purpose only, would compete in rank with the
hypothesis-driven set P values. This is a structural multiple-testing
error: tests answering different questions should not be pooled.

V2 corrects this by declaring two independent pools — 5
hypothesis-driven, 50 Hallmark descriptive — corrected separately,
never merged. Without this change, several of the five
hypothesis-driven sets might become non-significant solely because of
pool dilution rather than because of biological signal.

### Change 5: analysis script committed at freeze

V1 stated that gene sets and parameters were frozen, but did not
require the script that runs them to be committed before execution.
"Frozen analysis" without a frozen script is incomplete: a future
reader cannot verify what was actually run.

V2 requires `run_gsea_v2.R` to be committed at freeze, with no edits
afterwards except via a further dated amendment.

Honest framing: this is a compliance and audit requirement. If the
script were not committed, the analysis would still be scientifically
valid as long as it was run once without modification and the numbers
match. What committing the script protects is the *ability of an
external reviewer to verify* "what was run is what was planned."

### Change 6: gene set overlap check with explicit Jaccard > 0.30 threshold

V1 did not require a gene set overlap check. V2 adds the requirement
that pairwise Jaccard overlap among the four custom gene sets is
computed before the freeze commit (via `verify_gene_set_overlap.R`,
logged in `gene_set_overlap_log.tsv`), and that any pair with
Jaccard > 0.30 is addressed by reassigning shared genes to their
most specific set.

Honest framing — on the threshold: an earlier draft of this clause
used the word "substantial" rather than naming a specific threshold.
Fuzzy wording in pre-registration plans is the exact phrasing that
becomes ambiguous under post-hoc pressure — if results turn out
borderline, "substantial" could be read either as meeting or not
meeting the threshold depending on which reading is convenient. The
threshold was made explicit (0.30) in the v2 freeze for this reason.

For this study specifically, the concern is moot: the overlap check
was executed before the v2 freeze, with all six primary pairs
producing Jaccard ≤ 0.0244 (well below the 0.30 threshold). The
result is logged in `gene_set_overlap_log.tsv`, verified by
independent Python and R implementations. The R script
(`verify_gene_set_overlap.R`) is committed at freeze.

## 5. What is unchanged from v1

For the record, these are identical between v1 and v2:

- The biological hypothesis (leptin–mTORC1–autophagy/lipophagy–
  osteogenesis axis engagement in human CAVD bulk transcriptome).
- The three cohorts and their roles (discovery / sensitivity /
  cross-platform validation).
- The five pre-specified gene sets, by name and N. (A supplementary
  sixth gene set, CUSTOM_S100A9_RAGE_CALCIFICATION, exists in
  `Supp_Table_S1_v2.tsv` but is not counted among the five
  pre-specified sets for purposes of the BH correction pool.)
- Primary ranking metric: limma `logFC`, decreasing order.
- `minSize = 8`, `maxSize = 500`, `set.seed(42)`.
- Significance threshold: adjusted P < 0.05.
- Commitment to report all five pre-specified results regardless
  of significance.
- The directional-consistency requirement (same NES sign across the
  three cohorts) as a secondary criterion.

## 6. Manuscript cross-references

Future manuscript versions (v19 and onward) should cite all three
artifacts in Methods §2.4 and the Data Availability statement:

- v1 plan: git commit `ed2d35d` (tag `prereg-frozen`),
  path `preregistration/CAVD_GSEA_PreAnalysisPlan.md`
- v2 plan: git commit `c76f64489e94180e77741d8cbf0aca9a37134e9c`,
  path `pre/CAVD_GSEA_PreAnalysisPlan_v2.md`
- This amendment: git commit `c76f64489e94180e77741d8cbf0aca9a37134e9c`,
  path `pre/AMENDMENT_2026-05-26_plan_v1_to_v2.md`

Suggested Methods §2.4 wording for v19+:

> "The GSEA analysis was pre-registered. The original plan (v1) was
> frozen on 2026-05-23 (git commit `ed2d35d`, tag `prereg-frozen`).
> Prior to any GSEA execution, the plan was amended on 2026-05-26
> to v2 (git commit `c76f64489e94180e77741d8cbf0aca9a37134e9c`) for the following reasons:
> (i) v1's stated `fgsea(nperm=10000)` parameter combination was not
> executable as written and was replaced by `fgseaMultilevel()`;
> (ii) v1 did not specify a probe-to-gene collapse rule;
> (iii) v1 did not specify HGNC alias normalization, which is
> material for several genes in the custom sets;
> (iv) v1 pooled all gene sets under a single BH correction, which
> would have diluted hypothesis-driven adjusted P values by mixing
> them with descriptive Hallmark tests. The full list of changes
> and their rationales is documented in
> `AMENDMENT_2026-05-26_plan_v1_to_v2.md` (git commit `c76f64489e94180e77741d8cbf0aca9a37134e9c`).
> No GSEA call was executed against any cohort between the v1 freeze
> and the v2 amendment."

---

*End of amendment.*
