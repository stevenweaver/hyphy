RequireVersion("2.3.4");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("libv3/models/rate_variation.bf");

LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein/REV.bf");
LoadFunctionLibrary("libv3/models/protein.bf");
LoadFunctionLibrary("ProteinGTRFit_helper.ibf");


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);
utility.ToggleEnvVariable ("PRODUCE_OPTIMIZATION_LOG", 1); 

utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 1); // Uncomment for testing to make it all run faster.
//utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 0.001);

protein_gtr.analysis_banner = {
    terms.io.info: "Fit a general time reversible model to a collection
    of training protein sequence alignments. Generate substitution
    and scoring matrices following the procedures described in Nickle et al 2007",
    terms.io.version: "0.01",
    terms.io.reference: "Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL (2007) HIV-Specific Probabilistic Models of Protein Evolution. PLoS ONE 2(6): e503. doi:10.1371/journal.pone.0000503",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu"
};
io.DisplayAnalysisBanner(protein_gtr.analysis_banner);


protein_gtr.filename_to_index = terms.data.filename_to_index;
protein_gtr.logl = terms.fit.log_likelihood;
protein_gtr.phase = terms.fit.phase;
protein_gtr.json.information = "information";
protein_gtr.baseline_phase = "Baseline Phase";
protein_gtr.final_phase    = "REV-Final";
protein_gtr.rev_phase_prefix = "REV-Phase-";
protein_gtr.bl_phase_prefix = "REV-local-Phase-";
protein_gtr.options.convergence_type = "convergence type";
protein_gtr.options.tolerance        = "tolerance";
protein_gtr.options.baseline_model   = "baseline model";
protein_gtr.options.rate_variation   = "use rate variation";

protein_gtr.analysis_results = {terms.json.analysis: protein_gtr.analysis_banner,
                                terms.json.input: {},
                                terms.json.timers: {}};
                                
protein_gtr.timers = {};




/********************************************** MENU PROMPTS ********************************************************/
/********************************************************************************************************************/

// Load file containing paths to alignments for fitting and assess whether to start from scratch or resume a cached analysis

SetDialogPrompt ("Specify a multiple sequence alignment file");
protein_gtr.alignment_info  = alignments.ReadNucleotideDataSet ("protein_gtr.dataset", None);

protein_gtr.name_mapping = protein_gtr.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];
if (None == protein_gtr.name_mapping) {
    protein_gtr.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("protein_gtr.dataset"), "_value_", "`&protein_gtr.name_mapping`[_value_] = _value_");
}


protein_gtr.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (protein_gtr.alignment_info[utility.getGlobalValue("terms.data.partitions")], protein_gtr.name_mapping);
protein_gtr.partition_count = Abs (protein_gtr.partitions_and_trees);
io.CheckAssertion ("protein_gtr.partition_count==1", "This analysis can only handle a single partition");
protein_gtr.filter_specification = alignments.DefineFiltersForPartitions (protein_gtr.partitions_and_trees, "protein_gtr.dataset" , "protein_gtr.filter.", protein_gtr.alignment_info);
protein_gtr.trees = utility.Map (protein_gtr.partitions_and_trees, "_value_", '_value_[terms.data.tree]');
protein_gtr.data_filter = utility.Map (protein_gtr.filter_specification, "_value_", "_value_[terms.data.name]");



io.ReportProgressMessageMD ("relative_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_rates", "Data", "Loaded **" +
                            protein_gtr.alignment_info [terms.data.sequences] + "** sequences, **" +
                            protein_gtr.alignment_info [terms.data.sites] + "** sites, and **" + protein_gtr.partition_count + "** partitions from \`" + protein_gtr.alignment_info [terms.data.file] + "\`");

protein_gtr.file = protein_gtr.alignment_info[terms.data.file];
protein_gtr.json_file = protein_gtr.file  + ".json";
protein_gtr.final_likelihood_function = protein_gtr.file  + "_Final-Phase-LF.nex";
protein_gtr.model_file = protein_gtr.file  + ".fitted_model";

protein_gtr.baseline_model = "JTT";
protein_gtr.use_rate_variation = "Gamma"; 

protein_gtr.analysis_results[utility.getGlobalValue("terms.json.options")] = {utility.getGlobalValue("protein_gtr.options.convergence_type"): protein_gtr.convergence_type,
                                                                              utility.getGlobalValue("protein_gtr.options.tolerance"): protein_gtr.tolerance,
                                                                              utility.getGlobalValue("protein_gtr.options.baseline_model"): protein_gtr.baseline_model,
                                                                              utility.getGlobalValue("protein_gtr.options.rate_variation") : protein_gtr.use_rate_variation};
protein_gtr.analysis_results[utility.getGlobalValue("terms.json.input")] = { 
                                 utility.getGlobalValue("terms.json.file"): protein_gtr.file,
                                 utility.getGlobalValue("terms.json.sequences"): protein_gtr.file_info[utility.getGlobalValue("terms.data.sequences")],
                                 utility.getGlobalValue("terms.json.sites"): protein_gtr.file_info[utility.getGlobalValue("terms.data.sites")],
                                 utility.getGlobalValue("terms.json.trees"): (protein_gtr.tree["0"])[utility.getGlobalValue("terms.trees.newick_with_lengths")],
                               };


protein_gtr.baseline_model_name = protein_gtr.baseline_model + "+F, with 4 category Gamma rates";
protein_gtr.baseline_model_desc = "protein_gtr.Baseline.ModelDescription.withGamma";
protein_gtr.rev_model           = "models.protein.REVML.ModelDescription.withGamma";
/********************************************************************************************************************/



protein_gtr.startTimer (protein_gtr.timers, "Total time");


protein_gtr.queue = mpi.CreateQueue ({  utility.getGlobalValue("terms.mpi.Headers")   : utility.GetListOfLoadedModules ("libv3/") ,
                                        utility.getGlobalValue("terms.mpi.Functions") :
                                        {
                                            {"models.protein.REV.ModelDescription.withGamma",
                                             "models.protein.REV.ModelDescription.withGDD4",
                                             "protein_gtr.REV.ModelDescription",
                                             "protein_gtr.REV.ModelDescription.withGamma",
                                             "protein_gtr.REV.ModelDescription.withGDD4",
                                             "protein_gtr.REV.ModelDescription.freqs",
                                             "protein_gtr.Baseline.ModelDescription.withGamma",
                                             "protein_gtr.Baseline.ModelDescription.withGDD4",
                                             "protein_gtr.Baseline.ModelDescription",
                                             "protein_gtr.fitBaselineToFile"
                                            }
                                        },
                                        utility.getGlobalValue("terms.mpi.Variables") : {{
                                            "protein_gtr.shared_EFV",
                                            "protein_gtr.baseline_model_desc",
                                            "protein_gtr.rev_model",
                                            "protein_gtr.baseline_model",
                                            "protein_gtr.analysis_results",
                                            "protein_gtr.baseline_phase",
                                            "protein_gtr.file",
                                            "protein_gtr.trees",
                                            "protein_gtr.data_filter"
                                        }}
                                     });
  

                        
io.ReportProgressMessageMD ("Protein GTR Fitter", "Initial branch length fit", "Initial branch length fit");

protein_gtr.fit_phase = 0;
protein_gtr.scores = {};

/*************************** STEP ONE ***************************
Perform an initial fit of Baseline model+F(+/-4G) to the data
*****************************************************************/
console.log("\n\n Performing initial branch length optimization using " + protein_gtr.baseline_model);

protein_gtr.startTimer (protein_gtr.timers, protein_gtr.baseline_phase);
protein_gtr.timer_count +=1; 

utility.EnsureKey(protein_gtr.analysis_results, protein_gtr.baseline_phase);
protein_gtr.baseline_mle = estimators.FitSingleModel_Ext(protein_gtr.data_filter,
                                                        protein_gtr.trees,
                                                        protein_gtr.baseline_model_desc,
                                                        None,
                                                        None);
protein_gtr.baseline_mle - terms.global; // delete empty key
protein_gtr.analysis_results[protein_gtr.baseline_phase] = protein_gtr.baseline_mle;

protein_gtr.stopTimer (protein_gtr.timers, protein_gtr.baseline_phase);


console.log("\n\n Optimizing model in full, single pass.");
protein_gtr.startTimer (protein_gtr.timers, protein_gtr.final_phase);

protein_gtr.initial_values = {terms.global : {}, terms.branch_length : {}};

for (l1 = 0; l1 < 20; l1 += 1) {
    for (l2 = l1 + 1; l2 < 20; l2 += 1) {
        (protein_gtr.initial_values[terms.global]) [terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2])] =  {terms.fit.MLE : 0.1};
    }
}
protein_gtr.initial_values[terms.branch_length] = (protein_gtr.baseline_mle[terms.branch_length])[0];

utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);
utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", 0);

    
protein_gtr.trees2 = {"0": {terms.trees.newick :  (protein_gtr.baseline_mle[terms.fit.trees])[0]}};
protein_gtr.rev_mle = estimators.FitSingleModel_Ext (protein_gtr.data_filter,
                                                        protein_gtr.trees2,
                                                        protein_gtr.rev_model_final,
                                                        protein_gtr.initial_values,
                                         {terms.run_options.retain_lf_object : TRUE}
                                   );
protein_gtr.stopTimer (protein_gtr.timers, protein_gtr.final_phase);

                                   
console.log (""); // clear past the optimization progress line
utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
lf_id = protein_gtr.rev.mle[terms.likelihood_function];
Export(protein_gtr.finalphase_LF, ^lf_id);
protein_gtr.rev.mle - terms.likelihood_function;
fprintf(protein_gtr.final_likelihood_function, protein_gtr.finalphase_LF);                   
                              
                                                                                              


/* save, write custom model */
fprintf(protein_gtr.model_file, "custom_Rij = " + protein_gtr.extract_rates(protein_gtr.rev_mle) + ";");
fprintf(protein_gtr.model_file, "\n\n\n");
fprintf(protein_gtr.model_file, "custom_EFV = " + protein_gtr.extract_efv(protein_gtr.rev_mle) + ";");


/* Save the JSON */
protein_gtr.stopTimer (protein_gtr.timers, "Total time");
protein_gtr.analysis_results[terms.json.timers] = protein_gtr.timers;
io.SpoolJSON(protein_gtr.analysis_results, protein_gtr.json_file);
