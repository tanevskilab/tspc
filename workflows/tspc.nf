/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_tspc_pipeline'

include { BACKSUB } from '../modules/nf-core/backsub/main'
include { MAX_PROJECTION } from '../modules/local/max_projection/main'
include { PRINT_TEST } from '../modules/local/print_test/main'
include { CELLPOSESAM } from '../modules/local/cellposesam/main'
include { MCQUANT } from '../modules/nf-core/mcquant/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TSPC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'tspc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //ch_samplesheet.view()

    image_tuple = ch_samplesheet
            .map {meta, image ->
                [[id: meta.id], image[0]]}

    //    image_tuple.view()

    markersheet_tuple = ch_samplesheet
        .map {meta, image ->
            [[id: meta.id], image[1]]}

    //    markersheet_tuple.view()

    //channels_projection = ch_samplesheet
    //    .map {meta, image -> 
    //        [[id: meta.id], image[2]]}

    //channels_projection.view()
    //
    // MODULE: BACKSUB
    //
    // Sometimes the microscope already does background subtraction on the data, so it can be skiped:

    if (params.do_backsub) {
        BACKSUB(image_tuple, markersheet_tuple)
    } else {
        println "Skipping background subtraction as per user request (params.do_backsub = false)"
    }

    //BACKSUB.out.backsub_tif.view()

    PRINT_TEST(params.channels_projection_list)

    //
    // MODULE: MAX_PROJECTION
    //
    if (params.do_backsub) { 
        MAX_PROJECTION(BACKSUB.out.backsub_tif, params.channels_projection_list)
    } else {
        MAX_PROJECTION(image_tuple, params.channels_projection_list)
    }

    //
    // Cellpose-SAM segmentation
    //
    CELLPOSESAM(MAX_PROJECTION.out.max_proj_img)

    //
    // Mcquant
    //
    if (params.do_backsub) {
        MCQUANT(BACKSUB.out.backsub_tif, CELLPOSESAM.out.mask, markersheet_tuple)
    } else {
        MCQUANT(image_tuple, CELLPOSESAM.out.mask, markersheet_tuple)
    }
    
    // mc_quant_ch = MCQUANT(backsub_img[0], cellpose_ch, markers_nsclc_ch)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
