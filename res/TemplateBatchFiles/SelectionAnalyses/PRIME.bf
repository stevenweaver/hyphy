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

utility.ToggleEnvVariable ("OPTIMIZATION_PRECISION", 1);
utility.ToggleEnvVariable ("OPTIMIZATION_TIME_HARD_LIMIT", 1);


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

prime.table_headers = {
                        {"alpha_0", "Synonymous substitution rate at a site"}
                        {"alpha_1", "Synonymous substitution rate at a site"}
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
prime.scalers = {{"prime.alpha_scaler", "prime.beta_scaler_test", "prime.beta_scaler_nuisance"}};
model.generic.AddGlobal (prime.site.mg_rev, "prime.alpha_scaler", prime.site_alpha);
model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_test", prime.site_beta);
model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_nuisance", prime.site_beta_nuisance);
parameters.DeclareGlobal (prime.scalers, {});

prime.alpha = model.generic.GetLocalParameter (prime.site.mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
prime.beta = model.generic.GetLocalParameter (prime.site.mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));

io.CheckAssertion ("None!=prime.alpha && None!=prime.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

selection.io.startTimer (prime.json [terms.json.timers], "PRIME analysis", 2);

prime.report.count = {{0}};

prime.report.positive_site = {{"" + (1+((prime.filter_specification[prime.report.partition])[terms.data.coverage])[prime.report.site]),
                                    prime.report.partition + 1,
                                    Format(prime.report.row[0],7,3),
                                    Format(prime.report.row[3],7,3),
                                    Format(prime.report.row[4],7,3),
                                    Format(prime.report.row[5],7,3),
                                    "Yes, p = " + Format(prime.report.row[6],7,4),
                                    Format(prime.report.row[7],0,0)
}};

prime.site_results = {};

for (prime.partition_index = 0; prime.partition_index < prime.partition_count; prime.partition_index += 1) {

    prime.report.header_done = FALSE;
    prime.table_output_options[utility.getGlobalValue("terms.table_options.header")] = TRUE;

    model.ApplyModelToTree( "prime.site_tree_fel", prime.trees[prime.partition_index], {terms.default : prime.site.mg_rev}, None);

    prime.site_patterns = alignments.Extract_site_patterns ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")]);

    utility.ForEach (prime.case_respecting_node_names, "_node_",
            '_node_class_ = (prime.selected_branches[prime.partition_index])[_node_];
             if (_node_class_ == terms.tree_attributes.test) {
                _beta_scaler = prime.scalers[1];
             } else {
                _beta_scaler = prime.scalers[2];
             }
             prime.apply_proportional_site_constraint ("prime.site_tree", _node_, prime.alpha, prime.beta, prime.scalers[0], _beta_scaler, (( prime.final_partitioned_mg_results[terms.branch_length])[prime.partition_index])[_node_]);
        ');

    // create the likelihood function for this site
    ExecuteCommands (alignments.serialize_site_filter
                                       ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")],
                                       ((prime.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
                     ));

    __make_filter ("prime.site_filter");

    LikelihoodFunction prime.site_likelihood = (prime.site_filter, prime.site_tree_fel);

    estimators.ApplyExistingEstimates ("prime.site_likelihood", prime.site_model_mapping, prime.final_partitioned_mg_results,
                                        terms.globals_only);

    prime.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{"prime.site_likelihood"}},
                                   terms.mpi.Models : {{"prime.site.mg_rev"}},
                                   terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
                                   terms.mpi.Variables : {{"prime.selected_branches","prime.codon_data_info"}},
                                   terms.mpi.Functions : {{"prime.compute_branch_EBF"}}
                                 });


    // We need to queue a new job for each alpha
    constraints = {
        "0" : {"name": "alpha_0", "constraint" : "alpha_0 := 0", "cleanup" : "alpha_0 = 1"},
        "1" : {"name": "alpha_1", "constraint" : "alpha_1 := 0", "cleanup" : "alpha_1 = 1"},
        "2" : {"name": "alpha_2", "constraint" : "alpha_2 := 0", "cleanup" : "alpha_2 = 1"},
        "3" : {"name": "alpha_3", "constraint" : "alpha_3 := 0", "cleanup" : "alpha_3 = 1"},
        "4" : {"name": "alpha_4", "constraint" : "alpha_4 := 0", "cleanup" : "alpha_4 = 1"}
    };

    utility.ForEachPair (prime.site_patterns, "_pattern_", "_pattern_info_",
        '

            for (i = 0; i < 5; i += 1) {
                if (_pattern_info_[terms.data.is_constant]) {
                    prime.store_results (-1,None,{"0" : "prime.site_likelihood",
                                                 "1" : None,
                                                 "2" : prime.partition_index,
                                                 "3" : _pattern_info_,
                                                 "4" : prime.site_model_mapping,
                                                 "5" : constraints[i]
                                         });
                } else {
                    mpi.QueueJob (prime.queue, "prime.handle_a_site", {"0" : "prime.site_likelihood",
                                                                     "1" : alignments.serialize_site_filter
                                                                       ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")],
                                                                       (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                     "2" : prime.partition_index,
                                                                     "3" : _pattern_info_,
                                                                     "4" : prime.site_model_mapping,
                                                                     "5" : constraints[i]
                                                                        },
                                                                        "prime.store_results");
                }
            }
        '
    );

    mpi.QueueComplete (prime.queue);
    prime.partition_matrix = {Abs (prime.site_results[prime.partition_index]), Rows (prime.table_headers)};

    fprintf(stdout, prime.site_results);

    utility.ForEachPair (prime.site_results[prime.partition_index], "_key_", "_value_",
    '
        for (prime.index = 0; prime.index < Rows (prime.table_headers); prime.index += 1) {
            prime.partition_matrix [0+_key_][prime.index] = _value_["0"][prime.index];
        }
    '
    );

    prime.site_results[prime.partition_index] = prime.partition_matrix;

}

prime.json [terms.json.MLE] = {terms.json.headers   : prime.table_headers,
                               terms.json.content : prime.site_results };

// conduct prime post processing

io.ReportProgressMessageMD ("prime", "results", "** Found _" + prime.report.count[0] + "_ sites under episodic diversifying positive selection at p <= " + prime.pvalue + "**");

selection.io.stopTimer (prime.json [terms.json.timers], "Total time");
selection.io.stopTimer (prime.json [terms.json.timers], "PRIME analysis");

io.SpoolJSON (prime.json, prime.codon_data_info[terms.json.json]);

function prime.apply_proportional_site_constraint(tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    prime.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * prime.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * prime.branch_length__;
    ");
}

// Full Model
function prime.set_global_parameters() {

    //model.generic.AddGlobal(prim.site.mg_rev, "prime.alpha_0", prime.parameter_site_alpha);
    //parameters.DeclareGlobal ("prime.alpha_0", {});

    //model.generic.AddGlobal(prime.site.mg_rev, "prime.alpha_1", prime.parameter_site_alpha);
    //parameters.DeclareGlobal ("prime.alpha_1", {});

    //model.generic.AddGlobal (prime.site.mg_rev, "prime.alpha_scaler", prime.site_alpha);
    //model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_test", prime.site_beta);
    //model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_nuisance", prime.site_beta_nuisance);

    //parameters.DeclareGlobal (prime.scalers, {});

    prime.scalers = {{"prime.alpha_scaler", "prime.beta_scaler_test", "prime.beta_scaler_nuisance"}};
    model.generic.AddGlobal (prime.site.mg_rev, "prime.alpha_scaler", prime.site_alpha);
    model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_test", prime.site_beta);
    model.generic.AddGlobal (prime.site.mg_rev, "prime.beta_scaler_nuisance", prime.site_beta_nuisance);
    parameters.DeclareGlobal (prime.scalers, {});

}

//----------------------------------------------------------------------------------------
lfunction prime.compute_branch_EBF (lf_id, tree_name, branch_name, baseline) {

    parameter_name = "`tree_name`.`branch_name`." + ^"prime.branch_mixture";
    ^parameter_name = 1;

    LFCompute (^lf_id,LOGL0);

    utility.ExecuteInGlobalNamespace (parameter_name + ":= prime.site_mixture_weight");

    if (^"prime.site_mixture_weight" != 1 && ^"prime.site_mixture_weight" != 0) {
        _priorOdds = (1-^"prime.site_mixture_weight")/^"prime.site_mixture_weight";
    } else {
        _priorOdds = 0;
    }

    normalizer  = -Max (LOGL0,baseline);


    p1 = Exp(LOGL0+normalizer) * ^"prime.site_mixture_weight";
    p2 = (Exp(baseline+normalizer) - p1);

    _posteriorProb = {{p1,p2}};

    _posteriorProb = _posteriorProb * (1/(+_posteriorProb));
    if ( _priorOdds != 0) {
        eBF = _posteriorProb[1] / (1 - _posteriorProb[1]) / _priorOdds;
    } else {
        eBF = 1;
    }
    return {utility.getGlobalValue("terms.empirical_bayes_factor") : eBF__, utility.getGlobalValue("terms.posterior") : _posteriorProb__[1]};
}


// Performs MLE optimization per site. This function is executed in an embarrassingly parallel way.
lfunction prime.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping, test_condition) {

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    // declares constraint
    ExecuteCommands(test_condition["constraint"]);

    Optimize (results, ^lf);

    alpha_fit = estimators.ExtractMLEs (lf, model_mapping);
    alpha_fit[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];

    fit = {
        "name" : test_condition["name"],
        "fit" : alpha_fit
    };

    ExecuteCommands(test_condition["cleanup"]);

    return fit;
    
}

//
function prime.report.echo (prime.report.site, prime.report.partition, prime.report.row) {

    prime.print_row = None;

    if (prime.report.row [4] < prime.pvalue) {
        if (prime.report.row[0] < prime.report.row[1]) {
            prime.print_row = prime.report.positive_site;
            prime.report.counts[0] += 1;
        } else {
            prime.print_row = prime.report.negative_site;
            prime.report.counts [1] += 1;
        }
    }

    if (None != prime.print_row) {

        if (!prime.report.header_done) {
            io.ReportProgressMessageMD("prime", "" + prime.report.partition, "For partition " + (prime.report.partition+1) + " these sites are significant at p <=" + prime.pvalue + "\n");
            fprintf (stdout,
                io.FormatTableRow (prime.table_screen_output,prime.table_output_options));
            prime.report.header_done = TRUE;
            prime.table_output_options[terms.table_options.header] = FALSE;
        }

        fprintf (stdout,
            io.FormatTableRow (prime.print_row,prime.table_output_options));

    }

}

lfunction prime.store_results (node, result, arguments) {

    partition_index = arguments [2];
    pattern_info = arguments [3];

    fit_name = result["name"];
    fit_results = result["fit"];

    result_row = {{ 0, 0 }};

    // not a constant site
    if (None != result) { 
        result_row [0] = estimators.GetGlobalMLE (fit_results, ^"prime.site_alpha");
    }

    utility.EnsureKey (^"prime.site_results", partition_index);

    utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_prime_result_",
        '
            utility.EnsureKey (prime.site_results[`&partition_index`], _prime_result_);
            ((prime.site_results[`&partition_index`])[_prime_result_])[`&fit_name`] = `&result_row`;
        '
    );

    // Only report to screen when all tests have been completed
    //prime.report.echo (_prime_result_, `&partition_index`, `&result_row`);

}

// Called after all optimizations across sites takes place.
lfunction prime.compute_likelihoods(site_results) {


}


