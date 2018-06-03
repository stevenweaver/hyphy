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
LoadFunctionLibrary("PogoFit_helper.ibf"); // Functions, model definitions used for this batchfile.


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);
//utility.ToggleEnvVariable ("PRODUCE_OPTIMIZATION_LOG", 1); 
//utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 0.1);
pogofit.analysis_banner = {
    terms.io.info: "PogoFit, *P*rotein *G*TR *Fit*ter: Fit a general time reversible (GTR) model to a collection of training protein sequence alignments.",
    terms.io.version: "0.01",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu",
    terms.io.requirements: "All alignments must be in HyPhy-format: Each file must contain a protein multiple sequence alignment and newick phylogeny. NEXUS input is not accepted."
};
io.DisplayAnalysisBanner(pogofit.analysis_banner);

pogofit.baseline_phase   = "Baseline Fit";
pogofit.final_phase      = "GTR Fit";

pogofit.options.frequency_type = "frequency estimation";
pogofit.ml_freq                = "ML";
pogofit.emp_freq               = "Emp";

pogofit.single   = "Single";
pogofit.multiple = "Multiple";

pogofit.options.baseline_model   = "baseline model";

pogofit.output_hyphy = "HyPhy";
pogofit.output_paml  = "PAML";
pogofit.output_raxml = "RAxML";
pogofit.output_all   = "All";
pogofit.hyphy_model_ext = ".fitted_model";
pogofit.paml_model_ext  = ".paml";
pogofit.raxml_model_ext = ".raxml";

pogofit.analysis_results = {terms.json.analysis: pogofit.analysis_banner,
                                terms.json.input: {},
                                terms.json.timers: {}};
pogofit.timers = {};




/********************************************** MENU PROMPTS ********************************************************/
/********************************************************************************************************************/


// Prompt for number of files to analyze, and read file list accordingly //
pogofit.one_or_many  = io.SelectAnOption ({{pogofit.multiple, "Infer a protein model from multiple training datasets (this is more common)."}, 
                                               {pogofit.single, "Infer a protein model from a single training datasets."}}, 
                                                "How many datasets will be used to fit the protein model?");
pogofit.listfile = "";                                     
if (pogofit.one_or_many == pogofit.single)
{
    pogofit.alignment_file = io.PromptUserForString ("Provide the filename of the alignment to analyze");
    pogofit.listfile            = pogofit.alignment_file + ".list";
    fprintf(pogofit.listfile, CLEAR_FILE, pogofit.alignment_file);
    pogofit.file_list           = {{pogofit.alignment_file}};
    pogofit.output_model_prefix = pogofit.alignment_file;
    pogofit.json_file           = pogofit.alignment_file  + ".json";

}

if (pogofit.one_or_many == pogofit.multiple)
{
    SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
    fscanf (PROMPT_FOR_FILE, "Lines", pogofit.file_list);
    pogofit.listfile            = utility.getGlobalValue("LAST_FILE_PATH");
    pogofit.output_model_prefix = pogofit.listfile;
    pogofit.json_file           = pogofit.listfile  + ".json";
}

pogofit.file_list         = io.validate_a_list_of_files (pogofit.file_list);
pogofit.file_list_count   = Abs (pogofit.file_list);
pogofit.index_to_filename = utility.SwapKeysAndValues(pogofit.file_list);

// Prompt for baseline AA model //
pogofit.baseline_model  = io.SelectAnOption (models.protein.empirical_models,
                                                "Select an empirical protein model to use for optimizing the provided branch lengths:");

// Prompt for F inference //
pogofit.frequency_type  = io.SelectAnOption ({{pogofit.emp_freq, "+F (Empirical)"}, {pogofit.ml_freq, "Maximum likelihood"}},
                                                "Select an frequency specification:");
                     
// Prompt for output format //
pogofit.output_format  = io.SelectAnOption ({
                                                  {pogofit.output_hyphy, "HyPhy-formatted model (extension `.fitted_model`)"},
                                                  {pogofit.output_paml, "PAML-formatted model (extension `.paml`)"},
                                                  {pogofit.output_raxml, "RAXML-formatted model (extension `.raxml`)"},
                                                  {pogofit.output_all, "Output all file formats"}},
                                                 "Select an output format for the fitted model:");
pogofit.use_rate_variation = "Gamma"; 

pogofit.save_options();

pogofit.baseline_model_name = pogofit.baseline_model + "+F, with 4 category Gamma rates";
pogofit.baseline_model_desc = "pogofit.Baseline.ModelDescription.withGamma";
pogofit.initial_rates       = Eval("models.protein." + pogofit.baseline_model + ".Rij");

if (pogofit.frequency_type == pogofit.emp_freq){
    pogofit.rev_model = "models.protein.REV.ModelDescription.withGamma";
}
if (pogofit.frequency_type == pogofit.ml_freq){
    pogofit.rev_model = "models.protein.REVML.ModelDescription.withGamma";
}



/********************************************************************************************************************/
/********************************************* ANALYSIS BEGINS HERE *************************************************/
/********************************************************************************************************************/


pogofit.startTimer (pogofit.timers, "Total time");


pogofit.queue = mpi.CreateQueue ({  utility.getGlobalValue("terms.mpi.Headers")   : utility.GetListOfLoadedModules ("libv3/") ,
                                        utility.getGlobalValue("terms.mpi.Functions") :
                                        {
                                            {"models.protein.REV.ModelDescription.withGamma",
                                             "models.protein.REV.ModelDescription.withGDD4",
                                             "pogofit.REV.ModelDescription",
                                             "pogofit.REV.ModelDescription.withGamma",
                                             "pogofit.REV.ModelDescription.freqs",
                                             "pogofit.Baseline.ModelDescription.withGamma",
                                             "pogofit.Baseline.ModelDescription",
                                             "pogofit.fitBaselineToFile"
                                            }
                                        },
                                        utility.getGlobalValue("terms.mpi.Variables") : {{
                                            "pogofit.baseline_model_desc",
                                            "pogofit.rev_model",
                                            "pogofit.baseline_model",
                                            "pogofit.index_to_filename",
                                            "pogofit.analysis_results",
                                            "pogofit.baseline_phase",
                                            "pogofit.shared_EFV"
                                        }}
                                     });
  


/******************************************* STEP ONE *******************************************************
        Perform an initial fit of the Baseline model+F+4G to the each dataset independently
*************************************************************************************************************/
console.log("\n\n[PHASE 1] Performing initial branch length optimization using " + pogofit.baseline_model);

pogofit.startTimer (pogofit.timers, pogofit.baseline_phase);
pogofit.timer_count +=1; 

for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {

     io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                                     "Dispatching file '" + pogofit.file_list[file_index]);
    mpi.QueueJob (pogofit.queue, "pogofit.fitBaselineToFile", {"0" : pogofit.file_list[file_index]},
                                                            "pogofit.handle_baseline_callback");
}
mpi.QueueComplete (pogofit.queue);
pogofit.stopTimer (pogofit.timers, pogofit.baseline_phase);

pogofit.baseline_fit_logL = math.Sum (utility.Map (utility.Filter (pogofit.analysis_results, "_value_", "_value_/pogofit.baseline_phase"), "_value_", "(_value_[pogofit.baseline_phase])[terms.fit.log_likelihood]"));
io.ReportProgressMessageMD ("Protein GTR Fitter", " * Initial branch length fit",
                            "Overall Log(L) = " + pogofit.baseline_fit_logL);
/*************************************************************************************************************/
/*************************************************************************************************************/




/******************************************* STEP TWO *******************************************************
        Fit a full GTR model to all dataset(s) jointly, using the baseline model as initial rates
*************************************************************************************************************/
console.log("\n\n[PHASE 2] Optimizing protein model");

pogofit.startTimer (pogofit.timers, pogofit.final_phase);

pogofit.baseline_fit = utility.Map (utility.Filter (pogofit.analysis_results, "_value_", "_value_/'" + pogofit.baseline_phase + "'"), "_value_", "_value_['" + pogofit.baseline_phase + "']");
pogofit.gtr_fit = pogofit.fitGTR(pogofit.baseline_fit);
                                                                                                                              
pogofit.stopTimer (pogofit.timers, pogofit.final_phase);
/*************************************************************************************************************/
/*************************************************************************************************************/



/*********************** Save custom model to file(s) as specified **************************/
pogofit.final_rij = pogofit.extract_rates();
pogofit.final_efv = pogofit.extract_efv();
pogofit.write_model_to_file();


/************************************* Save analysis JSON ***********************************/
pogofit.stopTimer (pogofit.timers, "Total time");
pogofit.analysis_results[terms.json.timers] = pogofit.timers;
io.SpoolJSON(pogofit.analysis_results, pogofit.json_file);

console.log("\n\nAnalysis complete!");

