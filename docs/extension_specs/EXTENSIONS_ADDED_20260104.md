# What we added to the fork (branch adding_new_peak_callers)
 
  - Samplesheet now supports condition (new header: group,condition,replicate,fastq_1,fastq_2,control). Legacy
  sheets without condition still work (condition defaults to NA). This enables group × condition-aware
  behavior throughout.
  - New grouping metadata: every sample gets group_condition = "${group}_${condition}" and filenames/IDs avoid
  collisions across conditions.
  - Group/condition-aware normalisation controls
      - --normalisation_scope {all,group,group_condition}: for --normalisation_mode Spikein, spike-in scale
  factors can be computed globally (legacy) or separately per group / per group+condition.
      - --dump_scale_factors true: writes spike-in scale-factor TSVs under
  03_peak_calling/00_normalisation_factors/.
      - --igg_scale_scope {legacy,group_condition,sample}: for IgG/control scaling when using read-count
  normalisation modes (e.g. BPM/CPM/RPKM).
  - Condition-aware pooled IgG controls for callers that require it:
      - Pooled BAMs per (control_group, condition) under 03_peak_calling/01_pooled_controls/.
      - Fallbacks (if no exact condition-match control exists) are recorded in 03_peak_calling/07_qc_tables/
  control_pooling_fallbacks.tsv.
  - New peak callers + variants (in --peakcaller, comma-separated, case-insensitive):
      - MACS2_NARROW, MACS2_BROAD (no-control variants)
      - GOPEAKS_NARROW, GOPEAKS_BROAD (emits *_gopeaks.json for MultiQC)
      - EPIC2_200BP, EPIC2_150BP, EPIC2_25BP (requires pooled controls; needs --epic2_genome if not inferable)
      - SPAN_DEFAULT, SPAN_STRINGENT (requires pooled controls; requires --omnipeaks_jar)
      - Plus --peakcaller_preset {standard,extended} (extended = all new callers/variants)
  - Consensus peaks can be condition-aware
      - --consensus_grouping {group,group_condition} (defaults to group_condition when condition exists in the
  samplesheet).
  - Primary-caller semantics are enforced
      - The first entry in --peakcaller (or the first in the preset list) is the “primary” used for downstream
  consensus/QC; others are output-only.
      - EPIC2_* and SPAN_* error early if --use_control false.

# DONE (from docs/extension_specs/TODO.md)
 
  - #1 Add condition to samplesheet: Implemented (condition column supported; legacy sheets still work with
  default condition=NA; metadata now includes group_condition).
  - #2 Group-specific handling (Normalization + Peak-calling): Implemented (new scoping knobs like
  --normalisation_scope and --igg_scale_scope, condition-aware control pooling for control-required callers,
  and condition-aware consensus via --consensus_grouping).
  - #3 Support additional peakcallers: Implemented (new caller variants incl. GoPeaks, epic2, SPAN/OMNIPEAKS,
  MACS2 narrow/broad; plus --peakcaller_preset extended and primary-caller semantics).

# Peak Caller Variants integrated into our peak-calling pipeline;
  - seacr (SEACR; supports --use_control true/false)
  - macs2 (legacy MACS2; supports --use_control true/false, with broad/narrow controlled by existing MACS2
    params)
  - macs2_narrow (MACS2 narrow, no control)
  - macs2_broad (MACS2 broad, no control)
  - gopeaks_narrow (GoPeaks narrow; emits peaks + *_gopeaks.json for MultiQC)
  - gopeaks_broad (GoPeaks broad; emits peaks + *_gopeaks.json for MultiQC)
  - epic2_200bp, epic2_150bp, epic2_25bp (epic2 with different binning; requires pooled controls →
    --use_control true)
  - span_default, span_stringent (SPAN/OMNIPEAKS; requires pooled controls → --use_control true)

# How to run your modified pipeline (modeled on your 20251120_CNR scripts)
 
  1. Make a condition-aware samplesheet (this is what actually “turns on” the condition/group additions).
      - Your current samplesheet.csv encodes condition inside group (e.g. H3K4me3_sgControl). Convert to
  something like:
 
  group,condition,replicate,fastq_1,fastq_2,control
  IgG,sgControl,1,fastq/.../C1_1_..._1.fq.gz,fastq/.../C1_1_..._2.fq.gz,
  IgG,sgControl,2,fastq/.../C2_1_..._1.fq.gz,fastq/.../C2_1_..._2.fq.gz,
  IgG,sgNSD3,1,fastq/.../S1_1_..._1.fq.gz,fastq/.../S1_1_..._2.fq.gz,
  IgG,sgNSD3,2,fastq/.../S2_1_..._1.fq.gz,fastq/.../S2_1_..._2.fq.gz,
  H3K4me3,sgControl,1,fastq/.../C1_2_..._1.fq.gz,fastq/.../C1_2_..._2.fq.gz,IgG
  H3K4me3,sgControl,2,fastq/.../C2_2_..._1.fq.gz,fastq/.../C2_2_..._2.fq.gz,IgG
  H3K4me3,sgNSD3,1,fastq/.../S1_2_..._1.fq.gz,fastq/.../S1_2_..._2.fq.gz,IgG
  ...
 
  2. Create a new params JSON (copy your existing nf-params.narrow.json / nf-params.spikein.json and add the
  new knobs). Examples:
 
  - BPM run using condition-aware consensus + condition-aware IgG scaling + extended callers
      - Set at least:
          - "input": "samplesheet.with_condition.csv"
          - "consensus_grouping": "group_condition" (or omit; defaults this way when condition exists)
          - "igg_scale_scope": "group_condition" (new)
          - Either "peakcaller_preset": "extended" or explicitly choose/order callers via "peakcaller": "..."
      - If you include SPAN/epic2, add:
          - "omnipeaks_jar": "/path/to/omnipeaks.jar"
          - "epic2_genome": "mm39" (important for your GRCm39 resources unless you set genome in a way that
  can be inferred)
  - Spike-in run using group-scoped spike-in scaling + scale-factor TSVs
      - Add:
          - "normalisation_mode": "Spikein"
          - "normalisation_scope": "group" (or "group_condition")
          - "dump_scale_factors": true
 
  3. Run your fork instead of nf-core/cutandrun
 
  - Replace:
      - nextflow run nf-core/cutandrun -r 3.2.1 ...
  - With either:
      - nextflow run dhusmann/cutandrun -r adding_new_peak_callers -profile singularity -work-dir
  <scratch_workdir> -params-file <your_new_json> -resume
      - or run the local clone directly:
          - nextflow run /home/users/dhusmann/.nextflow/assets/nf-core/cutandrun_dev/cutandrun -profile
  singularity -work-dir <scratch_workdir> -params-file <your_new_json> -resume
 
  That will produce the new outputs under 03_peak_calling/ (pooled controls, per-caller peak folders,
  per-(group,condition) consensus peaks, optional normalisation-factor TSVs, and control fallback QC).

# REMAINS TODO
  - 4 Differential peak calling subworkflow (DiffBind samplesheets + standalone ChIPBinner + standalone
  SPAN): Not implemented yet (only design/docs scaffolding exists).
  - 5 --compare_norm_methods mode: Not implemented.
  - 6 --compare_peak_callers mode: Not implemented.
  - 7 Differential peak/enrichment characterization (ChromHMM, cCRE, LOLA/GSEA sets from group×cond, etc.):
  Not implemented.

