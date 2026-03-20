GEO provides access to NCBI Gene Expression Omnibus expression studies and platform annotations, with tasks to search datasets, inspect metadata, retrieve expression matrices, translate probes to genes, compute differential expression, derive ranked signatures, and score those signatures against user gene sets.

The workflow downloads and parses GEO SOFT resources on demand, then reuses the cached dataset and platform files in later runs. Most analysis tasks accept standard GEO dataset identifiers such as `GDS4513` or `GSE5846`. Platform annotation tables are used to infer the identifier format of each matrix and, when `to_gene` is enabled, to translate or collapse probe-level data into gene-level data. GEO annotations are not perfectly uniform across platforms, so agents that require a specific identifier namespace should inspect the key field of the returned TSV and not assume that every translated output uses exactly the same gene label.

A common usage pattern is to search with `query`, inspect the result with `dataset_info` and `sample_info`, fetch the matrix with `matrix`, run `differential`, and then either extract `up_genes` and `down_genes` or build `signatures` for downstream `rank_query` and `rank_query_batch` analyses. The workflow also includes `barcode` for discretized expression activity calls.

The `main` and `contrast` inputs used by differential-analysis tasks are intentionally flexible. In practice they are commonly passed as comma-separated sample identifiers such as `GSM136321,GSM136322,GSM136323`, as factor selections such as `disease state=relapse`, or, for series datasets with sample titles, as regular expressions written between slashes such as `/melanoma/`. The workflow passes these values to the matrix layer and additionally resolves slash-delimited regular expressions against `sample_info` titles for GSE series.

A minimal Ruby session looks like this:

```ruby
require 'rbbt/workflow'
Workflow.require_workflow "GEO"

datasets = GEO.job(:query, nil, query: "NCI60[Title]").run
dataset  = datasets.first

info    = GEO.job(:dataset_info, nil, dataset: dataset).run
samples = GEO.job(:sample_info, nil, dataset: dataset).run
matrix  = GEO.job(:matrix, nil, dataset: dataset, to_gene: true).run
```

A typical differential-expression run looks like this:

```ruby
diff = GEO.job(:differential, "relapse",
  dataset: "GDS4513",
  main: "disease state=relapse",
  to_gene: true
).run

up = GEO.job(:up_genes, "relapse",
  dataset: "GDS4513",
  main: "disease state=relapse",
  to_gene: true,
  threshold: 0.005
).run
```

The same workflow can be called from the command line:

```bash
scout workflow task GEO matrix --dataset GDS4513 --to_gene
scout workflow task GEO differential --dataset GDS4513 --main "disease state=relapse" --to_gene
```

The repository also includes a small worked example under `examples/differential/NCI60_SK-MEL`, where the dataset is `GSE5846`, the main samples are `GSM136321,GSM136322,GSM136323`, and `to_gene` is enabled.

# Tasks

## dataset_info
Return the parsed metadata for a GEO dataset

This task loads the YAML metadata associated with a dataset code. The returned structure typically includes the GEO platform, the reported value type, subset definitions, and, for GSE series, per-sample title information. It is the best entry point when an agent needs to discover how a dataset is organized before choosing groups for downstream analysis.

Because many other tasks depend on the same metadata, `dataset_info` is also useful for debugging. If an analysis behaves unexpectedly, inspect this task first to see which platform, subset names, and sample annotations were inferred from GEO.

## platform
Return the platform accession used by a dataset

This task is a small convenience wrapper around dataset metadata. It extracts the platform accession, such as a `GPL` code, from a dataset so that downstream logic can retrieve platform annotation tables or inspect the identifier system used by the expression matrix.

Use it when an agent needs to branch from a dataset-level analysis to a platform-level one without parsing the whole metadata hash manually.

## platform_info
Return the parsed metadata for a GEO platform

This task loads the YAML metadata generated for a GEO platform. Platform metadata is useful when you need to know the organism, the size of the annotation table, or the annotation file that will be used to translate probes into gene identifiers.

In this workflow, platform information is especially relevant for understanding the behavior of `matrix`, `differential`, `up_genes`, `down_genes`, and `signatures` when `to_gene` is enabled.

## query
Search GEO DataSets and return matching dataset identifiers

This task sends the provided query to the NCBI GEO DataSets endpoint and returns the matching dataset accessions as an array. The implementation targets the `gds` database, so this task returns `GDS` identifiers. If an agent already knows a `GSE` series accession, it can pass that accession directly to the other tasks without using `query` first.

This is the discovery task of the workflow. It is useful for building dataset lists that can then be inspected individually with `dataset_info`, `sample_info`, and `matrix`.

## sample_info
Return a sample annotation table for a dataset

This task produces a TSV keyed by sample identifier. For GSE series it uses the parsed sample metadata and returns at least a `Title` field. For GDS datasets, where GEO often represents sample groupings as subsets, the task reconstructs a factor table by assigning each sample the corresponding subset values.

Use this task to understand how samples can be grouped before calling `differential`. It is also the task to inspect when building regular-expression selectors against sample titles.

## matrix
Return the expression matrix for a dataset

This task retrieves the expression matrix for a GEO dataset and returns it as a TSV with features as rows and samples as fields. On first use the underlying GEO files are downloaded and parsed; later calls reuse the cached matrix.

When `to_gene` is false, the matrix remains in the original platform identifier space. When `to_gene` is true, the workflow uses the platform annotation table to translate probes and collapse them into gene-level rows. This option is often the most convenient choice for downstream gene-set analysis, but agents that need the original probe space should leave it disabled.

## differential
Compute differential expression between sample groups

This task runs differential analysis on a dataset matrix. The `main` input selects the target group and `contrast` optionally selects the reference group. Common input styles are comma-separated sample identifiers, factor selections such as `disease state=relapse`, and slash-delimited regular expressions that are matched against GSE sample titles.

If `to_gene` is enabled, the task first translates the matrix to gene-level identifiers and then runs the differential calculation. The result is a TSV produced by the matrix differential implementation and typically contains signed statistics, including p-values and t-values, that can be consumed by `up_genes`, `down_genes`, and `signatures`.

## up_genes
Return significant hits in the positive direction of a differential result

This task depends on `differential` and reuses its inputs. It filters the differential result with the provided `threshold` and returns the feature identifiers that are significant in the positive direction according to the signed adjusted p-value field generated by the differential analysis.

When the originating differential job used `to_gene`, the returned identifiers are translated through the platform annotation tables before being returned. Because annotation coverage differs across platforms, agents should inspect the returned values rather than assume a single universal gene label.

## down_genes
Return significant hits in the negative direction of a differential result

This task is the counterpart of `up_genes`. It depends on `differential`, applies the same significance threshold, and returns the identifiers that are significant in the negative direction of the signed differential result.

As with `up_genes`, gene-level translation is applied when the underlying differential computation was run with `to_gene`. This makes `down_genes` a convenient source of input for signature matching or pathway analysis.

## barcode
Compute barcode-style activity calls from an expression matrix

This task runs barcode transformation on a dataset matrix and returns a TSV representing discretized activity calls. The matrix can be kept in probe space or translated to genes first through `to_gene`.

The `standard_deviations` parameter controls how far above the lower mode a feature must be to be called active. This task is useful when a downstream method expects a simplified active or inactive representation instead of continuous expression values.

## signatures
Build ranked signatures for the comparisons available in a dataset

This task creates a collection of ranked differential-expression signatures for a dataset. It determines the comparisons to run from the dataset subsets. When no custom comparison file is present, the workflow automatically generates comparisons and gives preference to values that look like controls, such as labels containing `control`, `wild`, or `none`. For each comparison it runs `differential` and stores the resulting identifiers ordered by t-value.

By default `to_gene` is enabled, so signatures are usually gene-based rather than probe-based. The output is a TSV keyed by signature name, typically written as `<subset>: <main> => <contrast>`, with each value containing the ordered list of identifiers for that comparison. If a custom dataset comparison file is installed for a dataset, it overrides the automatic subset-based comparisons and can define comparisons explicitly with lines of the form `subset: main vs contrast`, including explicit sample lists in brackets.

## rank_query
Rank the signatures of one dataset against user-provided gene sets

This task depends on `signatures` and evaluates how well a query gene set matches each ranked signature of the dataset. Provide `up_genes` and optionally `down_genes`; with only `up_genes` the task performs a one-sided query, and with both lists it performs an up-versus-down query. The scoring uses permutation-based p-values over the ordered signature.

The result is a TSV keyed by signature name with at least `P-value` and `Hits` fields. The task also writes hit plots for each signature into the job files directory, which can be useful when an agent needs a visual summary of where the query genes fall in the ranking.

## rank_query_batch
Run rank_query across several datasets and merge the results

This task is the batch version of `rank_query`. It accepts a `datasets` array, launches one `rank_query` job per dataset, and merges all partial outputs into a single TSV.

Each merged key is prefixed with the dataset accession so that signatures from different studies remain distinguishable. This task is the most convenient entry point when an agent wants to screen a query against multiple GEO datasets in one call.
