# SPEC_DIFFERENTIAL_PEAK_CALLING.md
⸻

Differential Peak Calling Specification

Status: Draft (implementation-ready)
Version: 1.0
Date: 2026-01-04
Target: dhusmann/cutandrun (branch: differential_peak_calling)

⸻

Table of Contents
	1.	Overview￼
	2.	Goals and Non-goals￼
	3.	Design Principles￼
	4.	Pipeline Architecture￼
	1.	Integrated vs posthoc execution￼
	2.	Design manifest￼
	3.	Subworkflow layout￼
	4.	Data model and contrast semantics￼
	5.	Validation and skip/strict behavior￼
	5.	Shared Components￼
	1.	Peak set and region annotation￼
	2.	Summary + MultiQC inputs￼
	3.	Cross-method and cross-caller overlaps￼
	6.	Method 1: DiffBind￼
	7.	Method 2: ChIPBinner￼
	8.	Method 3: SPAN / OmniPeaks differential￼
	9.	Output Structure￼
	10.	User-facing Parameters￼
	11.	Implementation Plan￼
	12.	Testing Strategy and Acceptance Criteria￼
	13.	Resource Labels and Containers￼

⸻

Overview

This spec adds differential peak calling / differential enrichment as a first-class (opt-in) feature layer for the pipeline, using three complementary approaches:

Method	Best for	Core idea	Primary outputs
DiffBind	Peak-centric differential binding for punctate/focal marks and TFs	Counts reads in a peak consensus/union and tests with DESeq2/edgeR	Differential region TSV + up/down BED + diagnostic plots
ChIPBinner	Broad marks and global redistribution / genome-wide shifts	Genome-wide fixed binning + clustering + bin-level differential testing	Bin-level differential TSV + cluster BEDs + diagnostic plots (+ optional enrichment)
SPAN / OmniPeaks differential	SPAN users who want SPAN-consistent comparison	Use OmniPeaks differential mode if supported, otherwise explicit fallback	Differential regions BED/TSV (+ optional pooling manifest)

Key requirement: differential analysis must run for every peak caller variant specified in --peakcaller (primary + secondary), producing parallel DiffBind outputs per caller. ChIPBinner and SPAN-diff run per group (mark), not per caller (unless SPAN fallback explicitly uses SPAN-derived peaks).

⸻

Goals and Non-goals

Goals
	1.	Add opt-in differential analysis to the pipeline:
	•	DiffBind (per caller × group)
	•	ChIPBinner (per group)
	•	SPAN/OmniPeaks differential (per group; native if possible, explicit fallback otherwise)
	2.	Support two execution modes:
	•	Integrated mode (normal run): differential executes after peak calling outputs exist.
	•	Posthoc mode: run differential-only from an existing --outdir without recomputing upstream steps.
	3.	Require an explicit contrast definition to ensure consistent log2FC direction:
	•	log2FC = log2(treated / control) across all methods.
	4.	If spike-in normalization is present and enabled, incorporate external scale factors into differential analysis.
	5.	Emit stable, structured outputs including:
	•	raw result tables (TSV)
	•	annotated result tables (TSV)
	•	BED-like significant sets (up/down)
	•	diagnostic plots (PCA/correlation + method-specific)
	•	MultiQC-consumable summary tables

Non-goals (v1)
	•	Automatic modeling of batch covariates (batch-aware design matrices).
	•	“Compare normalization methods” / “compare peak callers” decision frameworks (beyond overlap metrics and summary tables).
	•	A full custom MultiQC plugin/module (v1 uses MultiQC custom-content tables and plot discovery).

⸻

Design Principles
	1.	Layered architecture: differential analysis is a separate layer that:
	•	does not mutate upstream peak calling outputs,
	•	can run integrated or posthoc,
	•	consumes published artifacts only in posthoc mode (no work/ assumptions).
	2.	Per-caller parity: any caller listed in --peakcaller is treated as a first-class source for DiffBind.
	3.	Deterministic + resumable: manifests and samplesheets must be stable-sorted to support -resume and reproducibility.
	4.	Explicit contrast direction: do not guess treated vs control from labels.
	5.	Fail/skip is user-controlled: default behavior is to skip invalid comparisons with clear warnings, but a strict mode is provided to fail fast.
	6.	Performance-aware design: avoid repeated expensive operations across methods; keep per-group computation parallelizable.

⸻

Pipeline Architecture

Integrated vs posthoc execution

Integrated mode
	•	Differential runs when any --run_diffbind, --run_chipbinner, or --run_span_diff is enabled.
	•	The workflow consumes in-memory channels produced by upstream alignment + peak calling.

Posthoc mode
	•	Set --differential_from_run <old_outdir> to run differential-only:
	•	A dedicated entry workflow DIFFERENTIAL_ONLY reads manifest(s) from <old_outdir>,
	•	reconstructs the necessary channels (BAM/BAI, per-caller peaks, optional scale factors),
	•	runs the same DIFFERENTIAL_PEAK_CALLING subworkflow.

Design manifest

To make posthoc differential reliable, the main pipeline writes deterministic manifests that contain published paths under --outdir.

Emit location (new):
03_peak_calling/06_differential/00_manifests/

Files (new):
	1.	differential_manifest.samples.tsv
One row per target sample (not IgG), minimum columns:

	•	sample_id
	•	group
	•	condition
	•	replicate
	•	final_bam
	•	final_bai
	•	normalisation_mode
	•	spikein_scale_factor (or NA)
	•	bigwig (optional; NA if absent)

	2.	differential_manifest.peaks.tsv
One row per (sample_id, caller):

	•	sample_id
	•	caller
	•	peaks_path
	•	peaks_format (bed/narrowPeak/broadPeak/other)
	•	caller_role (primary|secondary; optional)

	3.	differential_manifest.run_meta.json
Captures:

	•	pipeline version/hash
	•	reference genome id (if any)
	•	the caller list used
	•	the default contrast if provided (optional)

Manifest rules
	•	Sorted by stable keys (group, condition, replicate, sample_id, caller) before writing.
	•	Contains only published paths under --outdir.

Subworkflow layout

Add a new subworkflow:

subworkflows/local/differential_peak_calling.nf

It contains three gated method blocks and shared utilities:

DIFFERENTIAL_PEAK_CALLING
├── WRITE_DIFFERENTIAL_MANIFESTS           (integrated mode only)
├── (optional) VALIDATE_DIFFERENTIAL_DESIGN
├── DIFFBIND (per group × caller)
│   ├── MAKE_DIFFBIND_SAMPLESHEETS
│   ├── RUN_DIFFBIND_R
│   ├── ANNOTATE_REGIONS (shared)
│   └── DIFFBIND_SUMMARY
├── CHIPBINNER (per group)
│   ├── MAKE_BINS (or LOAD_BINS)
│   ├── BIN_COUNTS
│   ├── CHIPBINNER_CLUSTER_GRID (HDBSCAN grid search)
│   ├── CHIPBINNER_SELECT_MODELS
│   ├── CHIPBINNER_DIFFERENTIAL (ROTS)
│   ├── ANNOTATE_REGIONS (shared)
│   ├── (optional) CHIPBINNER_ENRICHMENT (LOLA)
│   └── CHIPBINNER_SUMMARY
├── SPAN_DIFF (per group)
│   ├── OMNIPEAKS_CAPABILITY_PROBE
│   ├── SPAN_COMPARE (native) OR SPAN_FALLBACK_DIFF
│   ├── ANNOTATE_REGIONS (shared)
│   └── SPAN_SUMMARY
└── CROSS_COMPARISON (optional)
    ├── OVERLAP_CALLERS (DiffBind across callers)
    └── OVERLAP_METHODS (DiffBind vs SPAN vs ChIPBinner)

Add a new DSL2 entry workflow:
	•	workflow DIFFERENTIAL_ONLY (new entrypoint)
	•	Reads manifests from --differential_from_run
	•	Runs DIFFERENTIAL_PEAK_CALLING only

Data model and contrast semantics

Existing metadata (already in fork per existing work):
	•	meta.group (mark/target)
	•	meta.condition
	•	meta.replicate (int)
	•	meta.sample_id / meta.id (unique)
	•	meta.caller (peak caller id)

Contrast definition (required for all methods):
	•	--differential_contrast "TREATED,CONTROL"
	•	log2FC sign convention: log2(TREATED / CONTROL)

Filtering rule (v1):
	•	For each group, differential compares only the two conditions in --differential_contrast.
	•	If additional conditions exist for that group:
	•	default: ignore them with a warning
	•	strict mode: error

Validation and skip/strict behavior

Validation is performed per group (and per caller where applicable):

Checks:
	•	both contrast conditions exist for the group
	•	each condition has ≥ --differential_min_replicates
	•	required inputs exist on disk (BAM/BAI + peak files)

Behavior:
	•	default: invalid comparisons are skipped with a clear warning and included in a skipped.tsv report
	•	--differential_strict true: any invalid comparison causes the workflow to fail

⸻

Shared Components

Peak set and region annotation

All methods must produce a BED-like region set that can be annotated consistently.

Shared annotation module:
	•	modules/local/annotate_regions.nf
	•	Uses bedtools closest against one of:
	•	--gene_bed if present, else
	•	a derived TSS BED from --gtf if present, else
	•	skip annotation with warning (still emit unannotated outputs)

Annotation fields appended to result TSV (minimum):
	•	nearest_feature_id (or gene id/name if available)
	•	distance_to_feature

Summary + MultiQC inputs

Each method writes a small summary TSV (one row per group/caller comparison) with standard columns:

Standard summary columns:
	•	method (diffbind|chipbinner|span)
	•	group
	•	caller (diffbind only; else NA)
	•	treated_condition
	•	control_condition
	•	n_tested
	•	n_fdr_pass
	•	n_up
	•	n_down

Additional method-specific columns may be appended (e.g., chosen HDBSCAN parameters).

MultiQC integration (v1):
	•	MultiQC custom content tables discover *.summary.tsv under differential output tree.
	•	No custom plugin required for v1.

Cross-method and cross-caller overlaps

To provide immediate utility and “best-of-both” value, emit overlap metrics:
	1.	Across peak callers (DiffBind significant sets)

	•	For each group:
	•	compute pairwise overlap counts and Jaccard between caller outputs

	2.	Across methods

	•	For each group:
	•	overlap DiffBind significant BED (per caller) with:
	•	SPAN significant BED (if present)
	•	ChIPBinner cluster-derived differential BEDs (if present)
	•	Emit both summary TSV and optional plots (UpSet if feasible; v1 can be summary-only)

Implementation:
	•	bedtools-based overlap (fast, deterministic)
	•	03_peak_calling/06_differential/04_cross_comparison/

⸻

Method 1: DiffBind

Purpose

Run classic differential binding using called peaks as inputs. Runs separately for each group × caller.

Inputs

Per sample:
	•	final BAM + BAI
	•	called peaks for the selected caller variant

Per comparison:
	•	explicit contrast (treated/control)
	•	optional spike-in scale factor per sample (if enabled)

Samplesheet generation

Write one DiffBind samplesheet CSV per (group, caller):

Path:
03_peak_calling/06_differential/01_diffbind/00_samplesheets/<caller>/<group>.diffbind.csv

Columns (minimum DiffBind-compatible):
	•	SampleID = sample id
	•	Factor = group
	•	Condition = condition
	•	Replicate = replicate
	•	bamReads = final bam path
	•	Peaks = peak path

Recommended extra columns:
	•	Tissue (constant, e.g. CUTRUN)
	•	Caller / PeakCaller (for provenance)
	•	SpikeinScaleFactor (optional; used by runner if enabled)

Execution

Create a new module:
	•	modules/local/diffbind_run.nf
	•	bin/diffbind_run.R

Runner responsibilities:
	1.	Validate design for the specific group×caller subset.
	2.	Construct DiffBind object from samplesheet.
	3.	Count reads in peaks:
	•	default: no summit recentering (use full peak widths)
	•	optional: if --diffbind_summits > 0, use summit-centered widths
	4.	Normalization:
	•	if external scaling enabled and scale factors exist:
	•	inject user-provided scaling into the downstream DE test (preferred)
	•	else:
	•	use DiffBind’s standard normalization/backend
	5.	Contrast:
	•	enforce treated vs control and consistent sign convention
	6.	Report:
	•	full TSV
	•	significant BED sets:
	•	up (log2FC > 0 & FDR pass)
	•	down (log2FC < 0 & FDR pass)
	7.	Plots (minimum):
	•	PCA
	•	correlation heatmap
	•	MA
	•	volcano
	•	heatmap (binding affinity)
	•	(optional) Venn/upset across conditions (if supported)

Spike-in scaling integration (v1 implementation contract)

If --differential_use_spikein resolves to true:
	•	Runner must attempt to use external sample scaling so differential log2FC reflects spike-in normalization.
	•	Implementation must be explicit and testable:
	•	either by setting normalization factors in DiffBind/DESeq2/edgeR, or
	•	by extracting the count matrix and running DESeq2/edgeR with user-provided size factors.

Outputs

Per (caller, group):

03_peak_calling/06_differential/01_diffbind/<caller>/<group>/

Minimum artifacts:
	•	diffbind.results.tsv
	•	diffbind.results.annotated.tsv (if annotation available)
	•	diffbind.significant.bed
	•	diffbind.significant_up.bed
	•	diffbind.significant_down.bed
	•	plots/PCA.pdf
	•	plots/correlation_heatmap.pdf
	•	plots/MA.pdf
	•	plots/volcano.pdf
	•	plots/heatmap.pdf
	•	diffbind.summary.tsv

⸻

Method 2: ChIPBinner

Purpose

Broad-mark differential enrichment with robustness to global changes:
	•	fixed genomic bins (default 10kb)
	•	normalization (external scaling when available)
	•	clustering (HDBSCAN grid)
	•	bin-level differential testing (ROTS)
	•	optional enrichment (LOLA)

Inputs

Per sample:
	•	final target BAM + BAI
	•	group / condition / replicate metadata
	•	optional external scale factors:
	•	spike-in scale factors (preferred when present)
	•	MS coefficients YAML (optional)

Reference:
	•	chrom sizes (*.sizes)
	•	optional blacklist BED

Workflow (per group)
	1.	Bin definition
	•	default: generate bins from chrom sizes
	•	optional: load precomputed windows from --chipbinner_windows_dir
	•	optional: subtract blacklist regions
	2.	Quantification
	•	compute bin-level counts matrix across samples
	3.	Normalization
	•	apply external scaling factors if enabled
	•	add pseudocount (default 1) before log transforms for clustering/plots
	4.	QC
	•	PCA plot
	•	correlation heatmap
	5.	Clustering (grid search)
	•	run HDBSCAN over bins using a parameter grid:
	•	min_cluster_size values
	•	min_samples/minPts values
	•	emit:
	•	chipbinner.hdbscan_grid_summary.tsv ranking solutions by:
	•	number of clusters (excluding noise)
	•	fraction of bins assigned (non-noise)
	•	mean cluster stability/persistence (if available)
	•	replicate consistency proxy (optional metric)
	6.	Model selection
	•	choose:
	•	best overall solution
	•	best solution yielding ~2 clusters (if any)
	•	best solution yielding ~3 clusters (if any)
	7.	Differential testing
	•	run ROTS on bins (replicates separate)
	•	produce bin-level log2FC/FDR with consistent sign
	8.	Cluster labeling
	•	label clusters by behavior based on per-cluster mean signals:
	•	stable (similar signal)
	•	control_enriched
	•	treated_enriched
	•	noise (if applicable)
	9.	Annotation
	•	annotate bins and/or cluster BEDs to nearest features (shared module)
	10.	(Optional) Enrichment

	•	if a LOLA DB is supplied, run enrichment per cluster

Outputs (per group)

03_peak_calling/06_differential/02_chipbinner/<group>/

Minimum artifacts:
	•	bins/<group>.<bin_size>.windows.bed
	•	counts/chipbinner.bin_counts.tsv
	•	counts/chipbinner.normalized_matrix.tsv
	•	clustering/chipbinner.hdbscan_grid_summary.tsv
	•	clustering/chipbinner.clusters.best.tsv
	•	differential/chipbinner.differential.tsv
	•	differential/chipbinner.differential.annotated.tsv (if annotation available)
	•	bed/chipbinner.significant_up.bed
	•	bed/chipbinner.significant_down.bed
	•	plots/PCA.pdf
	•	plots/correlation_heatmap.pdf
	•	plots/density_scatter.pdf (treated vs control means, colored by cluster)
	•	chipbinner.summary.tsv
Optional enrichment (if enabled):
	•	enrichment/<cluster>_lola.tsv
	•	enrichment/<cluster>_lola_plot.pdf

⸻

Method 3: SPAN / OmniPeaks differential

Purpose

Provide a SPAN-consistent differential workflow:
	•	Prefer OmniPeaks native differential mode if supported.
	•	Otherwise, run an explicit fallback and label it clearly.

Capability detection

Add a probe step:
	•	OMNIPEAKS_CAPABILITY_PROBE
	•	Runs java -jar omnipeaks.jar --help (or equivalent) and parses supported subcommands.
	•	Emits omnipeaks_capabilities.json for audit.

Native mode (if supported)

If a native differential subcommand exists (e.g., compare):
	•	Run it per group using the two contrast conditions.
	•	Replicate handling:
	•	if tool supports multiple BAMs per condition: pass them
	•	else: pool BAMs per condition and write pooling manifest

Fallback mode (explicit)

If no native differential exists:
	•	Do not silently pretend it does.
	•	Implement “SPAN peaks + DESeq2/edgeR” fallback:
	1.	Require SPAN peak outputs to exist (SPAN caller present in --peakcaller or peaks provided in manifests).
	2.	Build union peak set across samples in the contrast.
	3.	Count reads per union region from final BAMs.
	4.	Run DESeq2/edgeR with external scaling if enabled.
	5.	Emit outputs under span_fallback/ with clear naming.

Outputs (per group)

03_peak_calling/06_differential/03_span/<group>/

Native:
	•	span.differential.tsv
	•	span.differential.bed
	•	span.differential.annotated.tsv (if annotation available)
	•	span.significant_up.bed
	•	span.significant_down.bed
	•	span.summary.tsv
	•	00_manifests/span_pooling.tsv (only if pooling used)

Fallback:
	•	same outputs but with span_fallback.* prefix and a span_fallback.readme.txt describing the method.

⸻

Output Structure

All outputs live under:

03_peak_calling/06_differential/

03_peak_calling/06_differential/
├── 00_manifests/
│   ├── differential_manifest.samples.tsv
│   ├── differential_manifest.peaks.tsv
│   └── differential_manifest.run_meta.json
├── 01_diffbind/
│   ├── 00_samplesheets/<caller>/<group>.diffbind.csv
│   ├── <caller>/<group>/
│   │   ├── diffbind.results.tsv
│   │   ├── diffbind.results.annotated.tsv
│   │   ├── diffbind.significant.bed
│   │   ├── diffbind.significant_up.bed
│   │   ├── diffbind.significant_down.bed
│   │   ├── diffbind.summary.tsv
│   │   └── plots/*.pdf
├── 02_chipbinner/
│   └── <group>/
│       ├── bins/*.bed
│       ├── counts/*.tsv
│       ├── clustering/*.tsv
│       ├── differential/*.tsv
│       ├── bed/*.bed
│       ├── plots/*.pdf
│       ├── enrichment/ (optional)
│       └── chipbinner.summary.tsv
├── 03_span/
│   └── <group>/
│       ├── span.differential.tsv
│       ├── span.differential.bed
│       ├── span.differential.annotated.tsv
│       ├── span.significant_up.bed
│       ├── span.significant_down.bed
│       ├── span.summary.tsv
│       └── 00_manifests/ (optional)
├── 04_cross_comparison/ (optional)
│   ├── overlap_callers.tsv
│   ├── overlap_methods.tsv
│   └── plots/ (optional)
└── multiqc/
    ├── differential.summary.all_methods.tsv
    ├── differential.skipped.tsv
    └── (optional) differential_multiqc_report.html


⸻

User-facing Parameters

Top-level toggles
	•	--run_diffbind (bool; default: false)
	•	--run_chipbinner (bool; default: false)
	•	--run_span_diff (bool; default: false)

Execution mode
	•	--differential_from_run <path> (string; default: null)
	•	if set: run differential-only using manifests from that outdir

Contrast and validation behavior
	•	--differential_contrast "TREATED,CONTROL" (string; required if any method enabled)
	•	--differential_min_replicates (int; default: 2)
	•	--differential_strict (bool; default: false)
	•	--differential_use_spikein (string; default: "auto")
	•	allowed: auto|true|false
	•	auto resolves to true iff pipeline normalization mode is spike-in and scale factors exist

Shared outputs / QC
	•	--differential_annotate (bool; default: true)
	•	--differential_cross_compare (bool; default: true)
	•	--differential_run_multiqc (bool; default: true)

DiffBind parameters
	•	--diffbind_fdr (float; default: 0.05)
	•	--diffbind_lfc (float; default: 1.0)
	•	--diffbind_backend (string; default: DESeq2)
	•	allowed: DESeq2|edgeR
	•	--diffbind_min_overlap (int; default: 2)
	•	--diffbind_summits (int; default: 0)
	•	0 = do not recenter; >0 = summit-based fixed width
	•	--diffbind_extra_params (path; default: null)
	•	JSON/YAML override file for advanced DiffBind runner options

ChIPBinner parameters
	•	--chipbinner_bin_size (int; default: 10000)
	•	--chipbinner_windows_dir (path; default: null)
	•	if set, load bins from this dir; else generate
	•	--chipbinner_blacklist (path; default: null)
	•	--chipbinner_use_input (bool; default: false)
	•	--chipbinner_pseudocount (int; default: 1)
	•	--chipbinner_fdr (float; default: 0.05)
	•	--chipbinner_lfc (float; default: 1.0)
	•	--chipbinner_bootstrap (int; default: 1000)
	•	--chipbinner_k_value (int; default: 100000)
	•	--chipbinner_hdbscan_grid_min_cluster_size (string; default: "100,200,500,1000")
	•	--chipbinner_hdbscan_grid_min_samples (string; default: "100,200,500,1000")
	•	--chipbinner_ms_coeffs (path; default: null)
	•	--chipbinner_lola_db (path; default: null)
	•	--chipbinner_run_lola (bool; default: false)

SPAN / OmniPeaks parameters
	•	--omnipeaks_jar (path; required if --run_span_diff)
	•	--span_diff_mode (string; default: auto)
	•	auto|native|fallback
	•	--span_diff_fdr (float; default: 0.05)
	•	--span_diff_gap (int; default: 5)
	•	--span_diff_bin (int; default: 200)
	•	--span_diff_java_heap (string; default: 8G)

⸻

Implementation Plan

New files

Subworkflows:
	•	subworkflows/local/differential_peak_calling.nf
	•	workflows/differential_only.nf (or add new entrypoint in existing workflow file)

Modules:
	•	modules/local/write_differential_manifests.nf
	•	modules/local/diffbind_run.nf
	•	modules/local/chipbinner_bins.nf
	•	modules/local/chipbinner_counts.nf
	•	modules/local/chipbinner_hdbscan_grid.nf
	•	modules/local/chipbinner_rots.nf
	•	modules/local/span_capability_probe.nf
	•	modules/local/span_compare.nf
	•	modules/local/span_fallback_diff.nf
	•	modules/local/annotate_regions.nf
	•	modules/local/differential_overlap.nf

Scripts:
	•	bin/diffbind_run.R
	•	bin/chipbinner_rots.R
	•	bin/chipbinner_select_models.R (or embedded in rots runner)
	•	bin/span_fallback_deseq2.R (or shared DE runner)
	•	bin/hdbscan_grid.py (if clustering done in python)
	•	bin/annotate_regions.sh (optional wrapper around bedtools)

Modified files
	•	workflows/cutandrun.nf (call DIFFERENTIAL_PEAK_CALLING in integrated mode)
	•	nextflow.config (params defaults)
	•	nextflow_schema.json (new params + help text)
	•	conf/base.config (resource labels for new processes)
	•	assets/multiqc_config.yml (custom content tables ordering + file cleanup)
	•	docs/usage.md (new params & examples)
	•	docs/output.md (new output tree section)
	•	tests/ (new minimal tests and/or profile behavior)

⸻

Testing Strategy and Acceptance Criteria

Minimal CI-friendly tests (recommended)
	1.	Parameter validation tests:

	•	missing --differential_contrast when methods enabled → fails
	•	only one condition present → skipped (default) or fails (strict)
	•	<2 replicates/condition → skipped (default) or fails (strict)

	2.	Manifest generation test:

	•	after a small peak-calling test run, confirm:
	•	differential_manifest.samples.tsv exists with expected columns
	•	differential_manifest.peaks.tsv exists and includes all callers

	3.	DiffBind samplesheet generation test:

	•	enabling --run_diffbind should create per-caller per-group samplesheets

Optional heavier integration tests
	•	Run a small dataset with 2×2 design and confirm:
	•	DiffBind outputs exist for at least one caller and group
	•	ChIPBinner outputs exist for one group
	•	SPAN diff produces either native or fallback outputs (depending on tool availability)

Acceptance criteria
	•	DiffBind runs per group × caller and respects log2FC direction.
	•	Posthoc differential-only mode runs without upstream recompute (using manifests).
	•	Spike-in scaling is applied when enabled and scale factors exist.
	•	Each method emits: result TSV, annotated TSV (if possible), significant BEDs, plots, summary TSV.
	•	MultiQC custom content tables can load the summaries.

⸻

Resource Labels and Containers

Resource labels (conf/base.config)

Add/extend labels:
	•	process_single (default)
	•	process_medium (DiffBind, annotation)
	•	process_high (ChIPBinner binning/clustering, DE tests)
	•	process_long (SPAN compare, heavy clustering)

Containers / environments

DiffBind:
	•	R (>=4.2), Bioconductor DiffBind, DESeq2, edgeR, ggplot2, rtracklayer

ChIPBinner:
	•	R (>=4.2), ROTS, (optional) LOLA, plus clustering dependencies
	•	If python used: python3 + hdbscan + numpy/pandas

SPAN:
	•	Java runtime and omnipeaks.jar

All modules must emit versions.yml for nf-core provenance.
