# SPEC_DIFFERENTIAL_PEAK_CALLING.md - second version

Status: Draft
Version: 1.0
Date: 2026-01-04
Branch target: differential_peak_calling
Applies to: nf-core/cutandrun (DSL2) fork

Table of contents
	1.	Purpose and scope
	2.	Definitions and conventions
	3.	User-facing behavior
	4.	Experimental design rules and validation
	5.	Parameters
	6.	Manifests (integrated + posthoc reproducibility)
	7.	Workflow architecture (DSL2)
	8.	Method specs
	•	8.1 DiffBind
	•	8.2 ChIPBinner
	•	8.3 SPAN / OmniPeak differential
	9.	Shared annotation layer
	10.	MultiQC integration
	11.	Output layout
	12.	Resources and performance
	13.	Testing strategy
	14.	Acceptance criteria
	15.	Out of scope (v1)

⸻

1. Purpose and scope

Add an opt‑in differential enrichment / peak analysis capability to the pipeline, supporting three complementary methods:
	•	DiffBind: peak-centric differential binding, run per (group × caller).
	•	ChIPBinner: genome-wide binning + clustering + ROTS differential (broad/global changes), run per group. ChIPbinner is designed for broad histone mark analysis via uniform binning.
	•	SPAN/OmniPeak: native compare differential when supported by the jar; otherwise a deterministic fallback producing a comparable differential table. SPAN documents an experimental compare command with treatment/control replicates and parameters such as --bin, --gap, and --fdr.

Key goals:
	•	Respect existing sample metadata: group, condition, replicate, caller, sample_id/id.
	•	Keep upstream peak calling outputs immutable.
	•	Ensure the same differential codepath runs in both:
	•	Integrated mode (during a normal run, using channels)
	•	Posthoc differential-only mode (reconstructing channels from manifests)
	•	Provide deterministic outputs and MultiQC summaries for reproducibility.

⸻

2. Definitions and conventions

2.1 Metadata fields

The pipeline already uses a condition-aware metadata model. Differential analysis assumes each sample has:
	•	sample_id (unique sample key; may be meta.id internally)
	•	group (mark / target / antibody group)
	•	condition (biological condition label)
	•	replicate (biological replicate identifier)
	•	caller (peak caller variant used to generate the sample’s peaks)

2.2 Contrast direction

All methods must report:
	•	log2FC = log2(treated / control)

The user provides the mapping via:
	•	--differential_contrast "treated,control"

2.3 Spike-in scaling convention

If spike-in normalization is enabled upstream, the pipeline may emit per-sample scale factors used for coverage scaling. When used in statistical models that expect size factors (e.g., DESeq2/edgeR), the differential module must convert scale factors to a consistent modeling convention:
	•	size_factor = 1 / scale_factor

The module must export the exact values used for transparency.

DiffBind supports providing user-defined normalization factors via dba.normalize().

⸻

3. User-facing behavior

3.1 Execution modes

Integrated mode (default when enabled)
Triggered when any of:
	•	--run_diffbind
	•	--run_chipbinner
	•	--run_span_diff

and --differential_from_run is not set.

Behavior:
	•	Differential subworkflow consumes in-memory channels emitted by upstream steps (final BAMs, peaks, optional bigWigs, scale factors).
	•	It also writes deterministic manifests under the current --outdir to enable future posthoc runs.

Posthoc differential-only mode
Triggered when:
	•	--differential_from_run <prior_outdir> is provided
	•	and Nextflow entrypoint -entry DIFFERENTIAL_ONLY is used.

Behavior:
	•	Skip all upstream computation.
	•	Load manifests from <prior_outdir>/03_peak_calling/08_differential/00_manifests/.
	•	Reconstruct the same logical channels as integrated mode.
	•	Write new results to the current --outdir (never mutating <prior_outdir>).

3.2 Typical usage patterns
	•	Integrated:
	•	Enable one or more methods and provide a contrast.
	•	Posthoc:
	•	Use the prior run’s manifests and optionally change method parameters (e.g., FDR/lfc cutoffs).

3.3 Selective execution (performance)

The module supports restricting work to subsets:
	•	--differential_groups "H3K27me3,H3K4me3" (comma-separated)
	•	--differential_callers "SEACR,MACS2" (comma-separated caller IDs matching meta.caller)

If unset, run all eligible groups/callers.

⸻

4. Experimental design rules and validation

4.1 Global requirements

If any differential method is enabled:
	•	--differential_contrast is required and must contain exactly two comma-separated labels.

If contrast labels don’t exist in any sample’s condition:
	•	Fail immediately with a clear message.

4.2 Per-group eligibility

For each group:
	•	Eligible samples are those where condition ∈ {treated, control}.
	•	Other conditions may exist in the dataset; they are ignored for the differential run, and recorded as ignored_conditions in the design manifest.

A group is eligible if:
	•	Both contrast conditions exist for the group, and
	•	Each condition has ≥ --differential_min_replicates samples.

If a group lacks one contrast condition:
	•	Mark as SKIPPED (not an error) and record the reason.

If a group has both conditions but insufficient replicates:
	•	Default behavior: FAIL (to avoid silently producing statistically invalid results).
	•	If --differential_allow_partial is true: mark as SKIPPED and record the reason.

4.3 Per-caller eligibility (DiffBind)

DiffBind is evaluated per (group, caller):
	•	The same group-level replicate rules apply.
	•	Additionally, the peaks file for each sample must exist for that caller.
	•	If peaks are missing for some samples for a given caller:
	•	Mark (group,caller) as skipped (or fail in strict mode).

⸻

5. Parameters

All parameters must be added to:
	•	nextflow_schema.json
	•	nextflow.config defaults
	•	docs/usage.md + docs/output.md (fork docs)
	•	Module meta.yml (if following nf-core module conventions)

5.1 Global toggles and mode

Parameter	Type	Default	Notes
--run_diffbind	bool	false	Enable DiffBind method.
--run_chipbinner	bool	false	Enable ChIPBinner method.
--run_span_diff	bool	false	Enable SPAN/OmniPeak differential method.
--differential_contrast	string	(required if any run_ true)*	"treated,control", defines log2FC direction.
--differential_min_replicates	int	2	Minimum replicates per condition per group.
--differential_use_spikein	enum	auto	`auto
--differential_allow_partial	bool	false	If true, skip invalid group/caller comparisons instead of failing the run.
--differential_from_run	path	null	Enables posthoc mode when using -entry DIFFERENTIAL_ONLY.
--differential_publish_manifest_only	bool	false	Generate manifests + design tables only (no heavy method runs).
--differential_groups	string	null	Comma-separated allowlist of groups to run.
--differential_callers	string	null	Comma-separated allowlist of callers (DiffBind scope).

5.2 Normalization inputs

Parameter	Type	Default	Notes
--chipbinner_ms_coeffs	path	null	Optional MS-based coefficients for ChIPBinner scaling.
--dump_scale_factors	bool	(existing)	If upstream already supports dumping spike-in factors, manifests must reference published paths.

5.3 DiffBind parameters

Parameter	Type	Default	Notes
--diffbind_fdr	float	0.05	FDR cutoff.
--diffbind_lfc	float	1.0	Absolute log2FC threshold for “up/down”.
--diffbind_min_overlap	int	2	Min overlap for consensus counting.
--diffbind_backend	enum	DESeq2	`DESeq2
--diffbind_recenter_peaks	bool	false	If true, use summit-based recentering.
--diffbind_summits	int	0	0 = full peaks; otherwise summit window width.
--diffbind_norm_method	enum	native	`native
--diffbind_extra_params	path	null	JSON/YAML passthrough for advanced runner args.
--export_diffbind_sheets	bool	true	Write sample sheets to output.

5.4 ChIPBinner parameters

Parameter	Type	Default	Notes
--chipbinner_bin_size	int	10000	Allowed: 1000,5000,10000,25000.
--chipbinner_windows_dir	path	null	Prebuilt BED windows; otherwise generate.
--chipbinner_use_input	bool	false	Optional background usage (never treated as a condition).
--chipbinner_pseudocount	int	1	Added before log transforms.
--chipbinner_hdbscan_grid_minpts	string	"100,200,500,1000"	Grid values (comma-separated).
--chipbinner_hdbscan_grid_minsamps	string	"100,200,500,1000"	Grid values (comma-separated).
--chipbinner_fdr	float	0.05	FDR cutoff for ROTS results.
--chipbinner_lfc	float	1.0	Absolute log2FC cutoff for labeling.
--chipbinner_bootstrap	int	1000	ROTS bootstrap iterations.
--chipbinner_k_value	int	100000	ROTS K parameter.
--chipbinner_functional_db	path	null	Optional enrichment resources (LOLA).

5.5 SPAN / OmniPeak differential parameters

Parameter	Type	Default	Notes
--span_diff_mode	enum	auto	`auto
--span_diff_fdr	float	0.05	FDR cutoff.
--span_diff_gap	int	5	Gap parameter (native).
--span_diff_bin	int	200	Bin size (native).
--span_diff_java_heap	string	8G	Java heap for jar.
--span_fallback_backend	enum	DESeq2	`DESeq2
--omnipeaks_jar	path	(required if run_span_diff)	Jar used for SPAN/OmniPeak operations.


⸻

6. Manifests

6.1 Location

All manifests are written in integrated mode to:

03_peak_calling/08_differential/00_manifests/

and are consumed in posthoc mode from the prior run’s outdir.

6.2 Determinism requirements
	•	Always sorted deterministically (lexicographic by group, caller, condition, replicate, sample_id).
	•	Use published output paths (paths that exist under --outdir), never work directory paths.
	•	Always written if any differential method is enabled, even if all comparisons are skipped.

6.3 Manifest files and schemas

differential_manifest.samples.tsv
One row per target sample.

Columns:
	•	sample_id
	•	group
	•	condition
	•	replicate
	•	final_bam
	•	final_bai
	•	normalisation_mode
	•	spikein_scale_factor (numeric or NA)
	•	ms_coeff (numeric or NA)
	•	bigwig_path (path or NA)

differential_manifest.peaks.tsv
One row per (sample_id × caller) peak file.

Columns:
	•	sample_id
	•	group
	•	condition
	•	caller
	•	peaks_path
	•	peaks_format (e.g., bed|narrowPeak|broadPeak|peak)

differential_manifest.design.tsv
One row per analysis unit, capturing eligibility and reasons:

Columns:
	•	group
	•	caller (for DiffBind rows; NA for group-only methods)
	•	treated_condition
	•	control_condition
	•	n_treated
	•	n_control
	•	eligible_diffbind (true|false)
	•	eligible_chipbinner (true|false)
	•	eligible_span (true|false)
	•	status (RUN|SKIP|FAIL)
	•	reason (free text; stable short codes preferred)
	•	ignored_conditions (comma-separated list or NA)

span_diff_target_pooling.tsv (optional)
Only emitted if the jar does not support replicate lists and the pipeline had to pool replicates.

Columns:
	•	group
	•	condition
	•	pooled_bam
	•	source_bams (comma-separated)

⸻

7. Workflow architecture (DSL2)

7.1 New subworkflow

Create:

subworkflows/local/differential_peak_calling.nf

This subworkflow exposes a single entry function used by both modes:

DIFFERENTIAL_PEAK_CALLING(...)

7.2 New pipeline entrypoint

Add a new entrypoint workflow:

workflow DIFFERENTIAL_ONLY { ... }

Responsibilities:
	•	Load manifests from params.differential_from_run.
	•	Rebuild channels matching integrated mode.
	•	Call DIFFERENTIAL_PEAK_CALLING(...).
	•	Publish outputs to the current --outdir.

7.3 Internal phases

Phase A: Build or load design
	•	BUILD_DIFFERENTIAL_DESIGN (integrated): derive eligibility table + manifest rows from channels.
	•	LOAD_DIFFERENTIAL_DESIGN (posthoc): read manifests into channels and re-derive eligibility checks.

Phase B: Run methods
	•	DIFFBIND_ANALYSIS (per caller × group)
	•	CHIPBINNER_ANALYSIS (per group)
	•	SPAN_DIFFERENTIAL (per group; native or fallback)

Phase C: Shared post-processing
	•	ANNOTATE_REGIONS (shared module used by all methods)
	•	BUILD_DIFFERENTIAL_SUMMARIES (aggregates method summaries + skip/fail reasons)

⸻

8. Method specifications

8.1 DiffBind (per caller × group)

8.1.1 Inputs
For each eligible (group, caller):
	•	Deduplicated/final BAM + BAI for each target sample.
	•	Peak file for each sample from the same caller variant.
	•	Spike-in scale factors (optional).
	•	Contrast definition (treated, control).

Note: nf-core/cutandrun supports running multiple peak callers, where the first caller is primary and additional callers are run and written to results.  Differential must run for all caller variants present in the manifest.

8.1.2 Sample sheet generation
Path (if enabled):
03_peak_calling/08_differential/01_diffbind/00_samplesheets/<caller>/<group>.csv

Columns:
	•	SampleID (sample_id)
	•	Factor (group)
	•	Condition (condition)
	•	Replicate (replicate)
	•	bamReads (final_bam)
	•	Peaks (peaks_path)
	•	PeakCaller (caller)
	•	Tissue (optional constant e.g. CUTRUN)

8.1.3 R runner behavior
Implement:
	•	modules/local/diffbind_run.nf
	•	bin/diffbind_run.R

Steps:
	1.	Load sample sheet into dba().
	2.	Count reads: dba.count(...)
	•	If diffbind_recenter_peaks=false, keep full peaks (summits=0).
	•	Else use summits=params.diffbind_summits.
	3.	Normalization:
	•	If spike-in enabled and factors present:
	•	compute norm = 1 / spikein_scale_factor
	•	apply via dba.normalize(..., normalize=norm, ...) (or equivalent supported interface). DiffBind supports specifying normalization factors directly.
	•	Else:
	•	use DiffBind-native normalization (optionally guided by diffbind_norm_method).
	4.	Define contrast using --differential_contrast ordering.
	5.	Analyze using selected backend (DESeq2 or edgeR).
	6.	Export:
	•	full results table
	•	significant subsets by FDR and |LFC|
	•	standard plots (PCA, correlation heatmap, MA, volcano, etc.)

8.1.4 Outputs
03_peak_calling/08_differential/01_diffbind/<caller>/<group>/

Required files:
	•	diffbind.results.tsv
	•	diffbind.results.annotated.tsv
	•	diffbind.significant.bed
	•	diffbind.significant_up.bed
	•	diffbind.significant_down.bed
	•	diffbind.summary.tsv (one-row summary for MultiQC)
	•	plots/ (PCA, correlation heatmap, MA, volcano; optional others)

Recommended reproducibility files:
	•	diffbind.dba.rds (serialized DBA object)
	•	diffbind.normalization_factors.tsv

8.1.5 Error handling
	•	Design failures are handled before launching the runner.
	•	Runtime failures:
	•	Default: fail the run.
	•	If --differential_allow_partial: record (group,caller) failure in summaries and continue other units.

⸻

8.2 ChIPBinner (per group)

8.2.1 Inputs
Per eligible group:
	•	final BAM/BAI for all target samples in treated/control conditions
	•	chrom sizes
	•	optional blacklist BED (to exclude problematic bins)
	•	optional MS coeffs file (per sample)
	•	optional pooled input/IgG if chipbinner_use_input=true (implementation must match controls deterministically)

ChIPbinner divides the genome into uniform bins to analyze broad marks and differential enrichment patterns.

8.2.2 Steps
	1.	Prepare windows
	•	If chipbinner_windows_dir provided: load <genome>.<bin>.windows.bed.
	•	Else: generate windows from chrom sizes at chipbinner_bin_size.
	•	Subtract blacklist regions if available.
	•	Cache windows per (genome, bin_size, blacklist_hash) to avoid regeneration.
	2.	Quantify
	•	Produce per-sample binned counts matrix (rows=bins, cols=samples).
	•	Persist raw count matrix.
	3.	Normalize
	•	Apply scaling in this order:
	1.	spike-in scaling (if enabled and present)
	2.	MS coefficient scaling (if provided)
	3.	pseudocount addition
	4.	optional library-size normalization (if required by ChIPbinner routines)
	4.	QC
	•	PCA and correlation plots from normalized matrix.
	5.	HDBSCAN grid search
	•	Iterate over Cartesian product of minPts × minSamps.
	•	Persist cluster assignments for each grid point.
	•	Compute stability metrics and pick a preferred parameter set deterministically (document selection rule).
	6.	Cluster extraction
	•	Extract standardized 2- and 3-cluster solutions (when possible).
	7.	Differential testing
	•	Run ROTS on bins (replicates preserved), reporting:
	•	log2FC (treated/control)
	•	p-value
	•	FDR
	•	cluster assignment
	8.	Optional enrichment
	•	If chipbinner_functional_db set: run enrichment (e.g., LOLA) per cluster.
	•	Else: emit a “not run” stub file.

8.2.3 Outputs
03_peak_calling/08_differential/02_chipbinner/<group>/

Required files:
	•	chipbinner.samplesheet.csv
	•	chipbinner.windows.bed
	•	chipbinner.bin_counts.tsv
	•	chipbinner.normalized_matrix.tsv
	•	chipbinner.hdbscan_grid_summary.tsv
	•	chipbinner.clusters.tsv
	•	chipbinner.differential.tsv
	•	chipbinner.differential.annotated.tsv
	•	chipbinner.summary.tsv
	•	plots/ (PCA, correlation, density scatter, genic/intergenic if produced)

Optional:
	•	enrichment/ outputs when enabled

8.2.4 Error handling
	•	If clustering fails for some grid points, continue; choose best among successful points.
	•	If none succeed:
	•	Default: fail.
	•	With --differential_allow_partial: record failure for the group and continue.

⸻

8.3 SPAN / OmniPeak differential (per group)

8.3.1 Capability detection
Implement a lightweight “help probe” step that determines:
	•	Whether the jar supports a native differential subcommand (e.g., compare).
	•	Which signature is supported:
	•	SPAN-like compare -t <treated_reps> -c <control_reps> ... (documented)
	•	Any alternative signature exposed by the jar (detected via --help parsing)

Decision logic:
	•	If span_diff_mode=fallback: always run fallback.
	•	If span_diff_mode=native: require native support, else fail.
	•	If span_diff_mode=auto: native if supported, else fallback.

8.3.2 Native path
For each eligible group:
	•	Build comma-separated replicate lists for treated and control.
	•	Run the jar compare with:
	•	--chrom.sizes / --cs
	•	--bin (span_diff_bin)
	•	--gap (span_diff_gap)
	•	--fdr (span_diff_fdr)
	•	output peaks file (--peaks / -p)

SPAN documents:
	•	treatment replicates via comma-separated list
	•	control replicates via comma-separated list
	•	--bin, --gap, --fdr options for compare

If replicates are not supported by the jar signature:
	•	Pool replicates per condition (samtools merge), emit span_diff_target_pooling.tsv.

8.3.3 Fallback path
Goal: produce a deterministic differential table even without native support.

Fallback region set:
	•	Prefer union of SPAN/OmniPeak peaks from the SPAN caller variant (if present).
	•	If SPAN peaks are unavailable but run_span_diff is enabled:
	•	Fail with a clear message (v1).
	•	(Future enhancement could allow using primary caller peaks.)

Steps:
	1.	Build union region BED (sorted, merged).
	2.	Count reads per region per sample.
	3.	Run DE (DESeq2 or edgeR via span_fallback_backend).
	4.	Apply spike-in derived size factors when enabled.
	5.	Output consistent tables and BED splits (up/down).

8.3.4 Outputs
03_peak_calling/08_differential/03_span/<group>/

Required files:
	•	span.differential.tsv
	•	span.differential.peaks.bed (or native peaks format, plus a converted BED)
	•	span.differential.annotated.tsv
	•	span.up.bed
	•	span.down.bed
	•	span.summary.tsv
	•	span.mode.txt (native|fallback + signature note)

Optional:
	•	span_diff_target_pooling.tsv if pooling was required

⸻

9. Shared annotation layer

9.1 Inputs
	•	Regions as BED (DiffBind significant peaks, ChIPbinner differential bins, SPAN differential regions)
	•	Gene annotation:
	•	Prefer existing pipeline gene BED (if available)
	•	Else derive from GTF (gene bodies and/or TSS) and cache

9.2 Method
	•	Use bedtools closest to annotate nearest gene/feature and distance.
	•	Append standardized columns:
	•	nearest_feature_id
	•	nearest_gene_name (if available)
	•	distance_to_feature

9.3 Output contract

All *.annotated.tsv outputs must include:
	•	original region coordinates
	•	method statistics (log2FC, p-value, FDR)
	•	annotation columns above

⸻

10. MultiQC integration

10.1 Summary TSVs

Each method produces a one-row-per-unit summary TSV:
	•	DiffBind: diffbind.summary.tsv
	•	caller,group,treated,control,n_tested,n_fdr_pass,n_up,n_down,status,reason
	•	ChIPBinner: chipbinner.summary.tsv
	•	group,treated,control,n_bins_tested,n_fdr_pass,n_up,n_down,n_clusters,chosen_minPts,chosen_minSamps,status,reason
	•	SPAN: span.summary.tsv
	•	group,treated,control,n_tested,n_fdr_pass,n_up,n_down,mode,status,reason

10.2 Aggregated MultiQC table

Create:
03_peak_calling/08_differential/multiqc/differential_summary_mqc.tsv

This is a concatenation of the per-method summaries plus design manifest status, ensuring MultiQC has stable content even when analyses were skipped.

10.3 MultiQC config

Update:
	•	assets/multiqc_config.yml

Add a new section “Differential analysis” containing:
	•	a design/eligibility table (from differential_manifest.design.tsv)
	•	per-method summary tables
	•	optional links (paths) to key result files (if MultiQC supports)

10.4 Optional differential-only report (opt-in)

If --differential_multiqc_report=true (new optional param), produce:

03_peak_calling/08_differential/multiqc/differential_multiqc_report.html

This report is focused solely on differential outputs and does not modify the pipeline’s primary MultiQC report.

⸻

11. Output layout

All differential outputs live under:

03_peak_calling/08_differential/

03_peak_calling/08_differential/
  00_manifests/
    differential_manifest.samples.tsv
    differential_manifest.peaks.tsv
    differential_manifest.design.tsv
    span_diff_target_pooling.tsv (optional)
  01_diffbind/
    00_samplesheets/<caller>/<group>.csv
    <caller>/<group>/
      diffbind.results.tsv
      diffbind.results.annotated.tsv
      diffbind.significant*.bed
      diffbind.summary.tsv
      plots/
  02_chipbinner/
    <group>/
      chipbinner.*.tsv
      chipbinner.*.csv
      chipbinner.summary.tsv
      plots/
      enrichment/ (optional)
  03_span/
    <group>/
      span.differential.tsv
      span.differential.peaks.bed
      span.differential.annotated.tsv
      span.up.bed
      span.down.bed
      span.summary.tsv
      span.mode.txt
  multiqc/
    differential_summary_mqc.tsv
    differential_multiqc_report.html (optional)


⸻

12. Resources and performance

12.1 Default resource labels (recommended)

Provide defaults in conf/base.config using labels:
	•	DiffBind analyze: 4 CPU / 16 GB
	•	ChIPBinner clustering: 8 CPU / 32 GB
	•	ChIPBinner differential: 4 CPU / 16 GB
	•	SPAN native compare: 2 CPU / 8 GB
	•	SPAN fallback DE: 4 CPU / 16 GB

12.2 Performance considerations
	•	Design gating must run before launching heavy processes.
	•	Windows generation is cached per genome/bin size.
	•	Prefer streaming and minimal intermediate files where possible, but publish key matrices for reproducibility.

⸻

13. Testing strategy

13.1 CI-fast tests
	•	Validation:
	•	missing --differential_contrast when run flags are on
	•	contrast labels absent
	•	insufficient replicates (strict vs allow_partial)
	•	Manifest emission:
	•	determinism (sorted output)
	•	required columns present
	•	Sample sheet emission:
	•	DiffBind sheets created for eligible (group,caller)
	•	ChIPbinner samplesheet created per eligible group

Use --differential_publish_manifest_only and/or a stub profile to avoid heavy compute.

13.2 Optional integration tests (non-default)
	•	Single small dataset producing a real DiffBind result TSV.
	•	SPAN signature detection test using a mocked help output.
	•	ChIPBinner smoke test with reduced bootstrap/grid sizes.

⸻

14. Acceptance criteria
	•	Integrated runs produce manifests and method outputs for eligible units.
	•	Posthoc -entry DIFFERENTIAL_ONLY --differential_from_run <outdir> reproduces the same eligibility decisions and produces outputs without rerunning upstream steps.
	•	Log2FC sign is consistent (treated/control) across all methods.
	•	Spike-in factors propagate to DiffBind normalization factors when enabled.
	•	MultiQC tables are stable and show “skipped” reasons where applicable.

⸻

15. Out of scope (v1)
	•	Multiple contrasts per run (timecourse / all pairwise).
	•	Batch covariates in statistical models.
	•	Deep cross-method concordance (beyond lightweight summary tables).
	•	Automatic peak caller benchmarking.

