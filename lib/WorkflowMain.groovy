//
// This file holds several functions specific to the main.nf workflow in the nf-core/cutandrun pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "  https://doi.org/10.5281/zenodo.5653535\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log, args) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)
        // Check that the profile doesn't contain spaces and doesn't end with a trailing comma
        checkProfile(workflow.profile, args, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided unless running differential-only
        def differential_only = isDifferentialOnly(workflow)
        if (!params.input && !(differential_only && params.differential_from_run)) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

    //
    // Detect whether the workflow is running the differential-only entrypoint
    //
    public static boolean isDifferentialOnly(workflow) {
        def entry = getWorkflowEntryName(workflow)
        return entry?.toString()?.toUpperCase() == 'DIFFERENTIAL_ONLY'
    }

    //
    // Get the entrypoint name in a Nextflow-version-safe way
    //
    private static String getWorkflowEntryName(workflow) {
        if (!workflow) {
            return null
        }
        if (workflow.hasProperty('entry')) {
            return workflow.entry
        }
        if (workflow.hasProperty('entryName')) {
            return workflow.entryName
        }
        if (workflow.manifest && workflow.manifest.hasProperty('entryName')) {
            return workflow.manifest.entryName
        }
        if (workflow.commandLine) {
            def matcher = (workflow.commandLine =~ /(?:^|\\s)-entry(?:=|\\s+)(\\S+)/)
            if (matcher) {
                return matcher[0][1]
            }
        }
        return null
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttributeSpikeIn(params, attribute) {
        def val = ''
        if (params.genomes && params.spikein_genome && params.genomes.containsKey(params.spikein_genome)) {
            if (params.genomes[ params.spikein_genome ].containsKey(attribute)) {
                val = params.genomes[ params.spikein_genome ][ attribute ]
            }
        }
        return val
    }

    //
    // Exit pipeline if --profile contains spaces
    //
    private static void checkProfile(profile, args, log) {
        if (profile.endsWith(',')) {
            Nextflow.error "Profile cannot end with a trailing comma. Please remove the comma from the end of the profile string.\nHint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
        if (args[0]) {
            log.warn "nf-core pipelines do not accept positional arguments. The positional argument `${args[0]}` has been detected.\n      Hint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
    }
}
