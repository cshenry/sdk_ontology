/*
A KBase module: sdk_ontology
This module convert given KBase annotations of a genome to GO terms.
*/

module sdk_ontology {
	/*
        A binary boolean
    */
    typedef int bool;
    /*
        A string representing a Genome id.
    */
    typedef string genome_id;
    /*
        A string representing a workspace name.
    */
    typedef string workspace_name;
	/* 
        The workspace ID for a Report object
        @id ws KBaseReport.Report
    */
	typedef string ws_report_id;
	
    /*
        workspace - the name of the workspace for input/output
        input_genome - reference to the input genome object
        ontology_translation - optional reference to user specified ontology translation map
        output_genome - the name of the mapped genome annotation object

        @optional ontology_translation
    */
    typedef structure {
        string workspace;
        string input_genome;
        string ontology_translation;
        string translation_behavior;
 	    string custom_translation;
        string clear_existing;
        string output_genome;
    } ElectronicAnnotationParams;

    typedef structure {
        string report_name;
        string report_ref;
        string output_genome_ref;
        int n_total_features;
        int n_features_mapped;
    } ElectronicAnnotationResults;

    funcdef annotationtogo(ElectronicAnnotationParams params) returns (ElectronicAnnotationResults output)
        authentication required;
    
    typedef structure {
        string report_name;
        ws_report_id report_ref;
    } StandardFunctionOutput;
    
    typedef structure {
        workspace_name workspace;
        workspace_name genome_workspace;
		genome_id genome_id;
		genome_id genome_output_id;
    } AnnotateGenomeWithInterproPipelineParams;
    
    funcdef annotate_genome_with_interpro_pipeline(AnnotateGenomeWithInterproPipelineParams params) returns (StandardFunctionOutput output)
        authentication required;
};
