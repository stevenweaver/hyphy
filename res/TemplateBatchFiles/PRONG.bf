RequireVersion("2.3.9");
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
LoadFunctionLibrary("PRONG_helper.ibf"); // Functions, model definitions used for this batchfile.


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);
utility.ToggleEnvVariable ("PRODUCE_OPTIMIZATION_LOG", 1); 
//utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 0.1);
prong.analysis_banner = {
    terms.io.info: "PronG, `PROteiN GTR fitter`: Fit a general time reversible (GTR) model to a collection of training protein sequence alignments.",
    terms.io.version: "0.01",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu",
    terms.io.requirements: "All alignments must be in HyPhy-format: Each file must contain a protein multiple sequence alignment and newick phylogeny. NEXUS input is not accepted."
};
io.DisplayAnalysisBanner(prong.analysis_banner);

prong.baseline_phase   = "Baseline Fit";
prong.final_phase      = "GTR Fit";

prong.options.frequency_type = "frequency estimation";
prong.ml_freq                = "ML";
prong.emp_freq               = "+F";

prong.single   = "Single";
prong.multiple = "Multiple";

prong.options.baseline_model   = "baseline model";

prong.output_hyphy = "HyPhy";
prong.output_paml  = "PAML";
prong.output_raxml = "RAxML";
prong.output_all   = "All";
prong.hyphy_model_ext = ".fitted_model";
prong.paml_model_ext  = ".paml";
prong.raxml_model_ext = ".raxml";

prong.analysis_results = {terms.json.analysis: prong.analysis_banner,
                                terms.json.input: {},
                                terms.json.timers: {}};
prong.timers = {};




/********************************************** MENU PROMPTS ********************************************************/
/********************************************************************************************************************/


// Prompt for number of files to analyze, and read file list accordingly //
prong.one_or_many  = io.SelectAnOption ({{prong.multiple, "Infer a protein model from multiple training datasets (this is more common)."}, 
                                               {prong.single, "Infer a protein model from a single training datasets."}}, 
                                                "How many datasets will be used to fit the protein model?");
prong.listfile = "";                                     
if (prong.one_or_many == prong.single)
{
    prong.alignment_file = io.PromptUserForString ("Provide the filename of the alignment to analyze");
    prong.listfile            = prong.alignment_file + ".list";
    fprintf(prong.listfile, CLEAR_FILE, prong.alignment_file);
    prong.file_list           = {{prong.alignment_file}};
    prong.output_model_prefix = prong.alignment_file;
    prong.json_file           = prong.alignment_file  + ".json";

}

if (prong.one_or_many == prong.multiple)
{
    SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
    fscanf (PROMPT_FOR_FILE, "Lines", prong.file_list);
    prong.listfile            = utility.getGlobalValue("LAST_FILE_PATH");
    prong.output_model_prefix = prong.listfile;
    prong.json_file           = prong.listfile  + ".json";
}

prong.file_list         = io.validate_a_list_of_files (prong.file_list);
prong.file_list_count   = Abs (prong.file_list);
prong.index_to_filename = utility.SwapKeysAndValues(prong.file_list);

// Prompt for baseline AA model //
prong.baseline_model  = io.SelectAnOption (models.protein.empirical_models,
                                                "Select an empirical protein model to use for optimizing the provided branch lengths:");

// Prompt for F inference //
prong.frequency_type  = io.SelectAnOption ({{prong.emp_freq, "+F (Empirical)"}, {prong.ml_freq, "Maximum likelihood"}},
                                                "Select an frequency specification:");
                     
// Prompt for output format //
prong.output_format  = io.SelectAnOption ({
                                                  {prong.output_hyphy, "HyPhy-formatted model (extension `.fitted_model`)"},
                                                  {prong.output_paml, "PAML-formatted model (extension `.paml`)"},
                                                  {prong.output_raxml, "RAXML-formatted model (extension `.raxml`)"},
                                                  {prong.output_all, "Output all file formats"}},
                                                 "Select an output format for the fitted model:");
prong.use_rate_variation = "Gamma"; 

prong.save_options();

prong.baseline_model_name = prong.baseline_model + "+F, with 4 category Gamma rates";
prong.baseline_model_desc = "prong.Baseline.ModelDescription.withGamma";
prong.initial_rates       = Eval("models.protein." + prong.baseline_model + ".Rij");

if (prong.frequency_type == prong.emp_freq){
    prong.rev_model = "models.protein.REV.ModelDescription.withGamma";
}
if (prong.frequency_type == prong.ml_freq){
    prong.rev_model = "models.protein.REVML.ModelDescription.withGamma";
}



/********************************************************************************************************************/
/********************************************* ANALYSIS BEGINS HERE *************************************************/
/********************************************************************************************************************/


prong.startTimer (prong.timers, "Total time");


prong.queue = mpi.CreateQueue ({  utility.getGlobalValue("terms.mpi.Headers")   : utility.GetListOfLoadedModules ("libv3/") ,
                                        utility.getGlobalValue("terms.mpi.Functions") :
                                        {
                                            {"models.protein.REV.ModelDescription.withGamma",
                                             "models.protein.REV.ModelDescription.withGDD4",
                                             "prong.REV.ModelDescription",
                                             "prong.REV.ModelDescription.withGamma",
                                             "prong.REV.ModelDescription.freqs",
                                             "prong.Baseline.ModelDescription.withGamma",
                                             "prong.Baseline.ModelDescription",
                                             "prong.fitBaselineToFile"
                                            }
                                        },
                                        utility.getGlobalValue("terms.mpi.Variables") : {{
                                            "prong.baseline_model_desc",
                                            "prong.rev_model",
                                            "prong.baseline_model",
                                            "prong.index_to_filename",
                                            "prong.analysis_results",
                                            "prong.baseline_phase",
                                            "prong.shared_EFV"
                                        }}
                                     });
  


/******************************************* STEP ONE *******************************************************
        Perform an initial fit of the Baseline model+F+4G to the each dataset independently
*************************************************************************************************************/
console.log("\n\n[PHASE 1] Performing initial branch length optimization using " + prong.baseline_model);

prong.startTimer (prong.timers, prong.baseline_phase);
prong.timer_count +=1; 

for (file_index = 0; file_index < prong.file_list_count; file_index += 1) {

     io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                                     "Dispatching file '" + prong.file_list[file_index]);
    mpi.QueueJob (prong.queue, "prong.fitBaselineToFile", {"0" : prong.file_list[file_index]},
                                                            "prong.handle_baseline_callback");
}
mpi.QueueComplete (prong.queue);
prong.stopTimer (prong.timers, prong.baseline_phase);

prong.baseline_fit_logL = math.Sum (utility.Map (utility.Filter (prong.analysis_results, "_value_", "_value_/prong.baseline_phase"), "_value_", "(_value_[prong.baseline_phase])[terms.fit.log_likelihood]"));
io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                            "Overall Log(L) = " + prong.baseline_fit_logL);
/*************************************************************************************************************/
/*************************************************************************************************************/




/******************************************* STEP TWO *******************************************************
        Fit a full GTR model to all dataset(s) jointly, using the baseline model as initial rates
*************************************************************************************************************/
console.log("\n\n[PHASE 2] Optimizing protein model");

prong.startTimer (prong.timers, prong.final_phase);

prong.baseline_fit = utility.Map (utility.Filter (prong.analysis_results, "_value_", "_value_/'" + prong.baseline_phase + "'"), "_value_", "_value_['" + prong.baseline_phase + "']");
prong.gtr_fit = prong.fitGTR(prong.baseline_fit);
                                                                                                                              
prong.stopTimer (prong.timers, prong.final_phase);
/*************************************************************************************************************/
/*************************************************************************************************************/



/*********************** Save custom model to file(s) as specified **************************/
prong.final_rij = prong.extract_rates();
prong.final_efv = prong.extract_efv();
prong.write_model_to_file();


/************************************* Save analysis JSON ***********************************/
prong.stopTimer (prong.timers, "Total time");
prong.analysis_results[terms.json.timers] = prong.timers;
io.SpoolJSON(prong.analysis_results, prong.json_file);

console.log("\n\nAnalysis complete!");

