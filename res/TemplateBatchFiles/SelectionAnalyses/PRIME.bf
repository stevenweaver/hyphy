RequireVersion("2.3");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");

/* Display analysis information */

prime.analysis_description = {
    terms.io.info: "PRIME",
    terms.io.version: "1.00",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};

io.DisplayAnalysisBanner(prime.analysis_description);

/* Environment Setup */
utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/* Globals */
prime.site_alpha_0 = "Property 0 Site relative synonymous rate";
prime.site_alpha_1 = "Property 1 Site relative synonymous rate";
prime.site_alpha_2 = "Property 2 Site relative synonymous rate";
prime.site_alpha_3 = "Property 3 Site relative synonymous rate";
prime.site_alpha_4 = "Property 4 Site relative synonymous rate";

// default cutoff for printing to screen
prime.p_value = 0.1;
prime.scaler_prefix = "PRIME.scaler";

// The dictionary of results to be written to JSON at the end of the run
prime.json = {
    terms.json.analysis: prime.analysis_description,
    terms.json.input: {},
    terms.json.fits: {},
    terms.json.timers: {},
};

selection.io.startTimer (prime.json [terms.json.timers], "Total time", 0);

prime.table_headers = {
                         {"alpha_0", "Synonymous substitution rate at a site"}
                         {"alpha_0_pval", "The rate estimate under the neutral model"}
                         {"alpha_1", "Synonymous substitution rate at a site"}
                         {"alpha_1_pval", "The rate estimate under the neutral model"}
                         {"alpha_2", "Synonymous substitution rate at a site"}
                         {"alpha_2_pval", "The rate estimate under the neutral model"}
                         {"alpha_3", "Synonymous substitution rate at a site"}
                         {"alpha_3_pval", "The rate estimate under the neutral model"}
                         {"alpha_4", "Synonymous substitution rate at a site"}
                         {"alpha_4_pval", "The rate estimate under the neutral model"}
                     };


/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/
prime.table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "LRT", "Selection detected?"}};
prime.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};

// Load all the file information
namespace prime {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("prime");
}


selection.io.startTimer (prime.json [terms.json.timers], "Model fitting",1);

namespace prime {
    doGTR ("prime");
}

estimators.fixSubsetOfEstimates(prime.gtr_results, prime.gtr_results[terms.global]);

io.ReportProgressMessageMD ("prime", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

namespace prime {
    doPartitionedMG ("prime", FALSE);
}

prime.final_partitioned_mg_results = estimators.FitMGREV(prime.filter_names, prime.trees, prime.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: prime.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, prime.partitioned_mg_results);

io.ReportProgressMessageMD("prime", "codon-refit", "* Log(L) = " + Format(prime.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));

prime.global_dnds = selection.io.extract_global_MLE_re (prime.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (prime.global_dnds, "_value_", 'io.ReportProgressMessageMD ("PRIME", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

estimators.fixSubsetOfEstimates(prime.final_partitioned_mg_results, prime.final_partitioned_mg_results[terms.global]);

//Store MG94 to JSON
selection.io.json_store_lf_GTR_MG94 (prime.json,
                            terms.json.global_mg94xrev,
                            prime.final_partitioned_mg_results[terms.fit.log_likelihood],
                            prime.final_partitioned_mg_results[terms.parameters],
                            prime.sample_size,
                            utility.ArrayToDict (utility.Map (prime.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (prime.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            prime.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (prime.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(prime.json, terms.json.global_mg94xrev, terms.branch_length, prime.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((prime.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');




selection.io.stopTimer (prime.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

// TODO: This is where we need to implement the QCAP model

prime.site.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "prime_mg", {
            "0": parameters.Quote(terms.local),
            "1": prime.codon_data_info[terms.code]
        },
        prime.filter_names,
        None);


prime.site_model_mapping = {"prime_mg" : prime.site.mg_rev};

/* set up the local constraint model */

// TODO : will need to update
prime.alpha = model.generic.GetLocalParameter (prime.site.mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
prime.beta = model.generic.GetLocalParameter (prime.site.mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
io.CheckAssertion ("None!=prime.alpha && None!=prime.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

selection.io.startTimer (prime.json [terms.json.timers], "prime analysis", 2);


// Full Model
function prime.set_global_parameters() {

    model.generic.AddGlobal(prim.site.mg_rev, "prime.alpha_0", prime.parameter_site_alpha);
    parameters.DeclareGlobal ("prime.alpha_0", {});

    model.generic.AddGlobal(prime.site.mg_rev, "prime.alpha_1", prime.parameter_site_alpha);
    parameters.DeclareGlobal ("prime.alpha_1", {});

    model.generic.AddGlobal (prime.site.mg_rev, "prime.alpha_scaler", prime.site_alpha);
    model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_test", prime.site_beta);
    model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_nuisance", prime.site_beta_nuisance);

    parameters.DeclareGlobal (prime.scalers, {});


}

lfunction prime.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping) {


    GetString   (lfInfo, ^lf_fel,-1);

    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    GetString (lfInfo, ^lf_bsrel,-1);
    __make_filter ((lfInfo["Datafilters"])[0]);

    bsrel_tree_id = (lfInfo["Trees"])[0];

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    ^"prime.site_alpha" = 1;
    ^"prime.site_beta_plus"  = 1;
    ^"prime.site_beta_nuisance"  = 1;
    
    Optimize (results, ^lf_fel);

    fel = estimators.ExtractMLEs (lf_fel, model_mapping);
    fel[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];


     ^"meme.site_mixture_weight" = 0.75;
     if (^"meme.site_alpha" > 0) {
         ^"meme.site_omega_minus" = 1;
     } else {
         ^"meme.site_omega_minus" = ^"meme.site_beta_plus" / Max (^"meme.site_alpha", 1e-6);
         /* avoid 0/0 by making the denominator non-zero*/
     }

    Optimize (results, ^lf_bsrel);

    alternative = estimators.ExtractMLEs (lf_bsrel, model_mapping);
    alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];


    ancestral_info = ancestral.build (lf_bsrel,0,FALSE);

    //TODO
    branch_substitution_information = selection.substitution_mapper (ancestral_info ["MATRIX"],
                                                      ancestral_info ["TREE_AVL"],
                                                      ancestral_info ["AMBIGS"],
                                                      ^"meme.pairwise_counts",
                                                      ancestral_info ["MAPPING"],
                                                      (^"meme.codon_data_info")[utility.getGlobalValue("terms.code")]);


    DeleteObject (ancestral_info);

    branch_ebf       = {};
    branch_posterior = {};

    if (^"meme.site_beta_plus" > ^"meme.site_alpha" && ^"meme.site_mixture_weight" > 0) {

        LFCompute (^lf_bsrel,LF_START_COMPUTE);
        LFCompute (^lf_bsrel,baseline);

        utility.ForEach (^bsrel_tree_id, "_node_name_",
        '
            if ((meme.selected_branches [^"`&partition_index`"])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                _node_name_res_ = meme.compute_branch_EBF (^"`&lf_bsrel`", ^"`&bsrel_tree_id`", _node_name_, ^"`&baseline`");
                (^"`&branch_ebf`")[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.empirical_bayes_factor")];
                (^"`&branch_posterior`")[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.posterior")];
            } else {
                (^"`&branch_ebf`")[_node_name_] = None;
                (^"`&branch_posterior`")[_node_name_] = None;
            }
        '
        );

        LFCompute (^lf_bsrel,LF_DONE_COMPUTE);

        ^"meme.site_beta_plus" := ^"meme.site_alpha";
        Optimize (results, ^lf_bsrel);

        null = estimators.ExtractMLEs (lf_bsrel, model_mapping);
        null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];



    } else {
        null = alternative;
        utility.ForEach (^bsrel_tree_id, "_node_name_",
        '
            if ((meme.selected_branches [^"`&partition_index`"])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                (^"`&branch_ebf`")[_node_name_] = 1.0;
                (^"`&branch_posterior`")[_node_name_] = 0.0;
            } else {
                (^"`&branch_ebf`")[_node_name_] = None;
                (^"`&branch_posterior`")[_node_name_] = None;
            }
        '
        );
    }

    return {"fel" : fel,
            utility.getGlobalValue("terms.alternative") : alternative,
            utility.getGlobalValue("terms.posterior") : branch_posterior,
            utility.getGlobalValue("terms.empirical_bayes_factor") : branch_ebf,
            utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information, //TODO: keep this attr?
            utility.getGlobalValue("terms.null"): null};
}

