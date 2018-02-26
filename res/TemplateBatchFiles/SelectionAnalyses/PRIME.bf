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


prime.setDefaultVariables();
io.DisplayAnalysisBanner(prime.analysis_description);
selection.io.startTimer (prime.json [terms.json.timers], "Total time", 0);

// Load all the file information
namespace prime {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("prime");
}

selection.io.startTimer (prime.json [terms.json.timers], "Model fitting",1);

// Initial global fit with GTR
namespace prime {
    doGTR ("prime");
}


// MG94 
estimators.fixSubsetOfEstimates(prime.gtr_results, prime.gtr_results[terms.global]);

namespace prime {
    doPartitionedMG ("prime", FALSE);
}

io.ReportProgressMessageMD ("prime", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

prime.final_partitioned_mg_results = estimators.FitMGREV(prime.filter_names, prime.trees, prime.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: prime.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, prime.partitioned_mg_results);

io.ReportProgressMessageMD("prime", "codon-refit", "* Log(L) = " + Format(prime.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));

prime.global_dnds = selection.io.regexExtractGlobalMLE(prime.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (prime.global_dnds, "_value_", 'io.ReportProgressMessageMD ("PRIME", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

estimators.fixSubsetOfEstimates(prime.final_partitioned_mg_results, prime.final_partitioned_mg_results[terms.global]);

//Store MG94 to JSON
selection.io.formatMG94ResultsToJSON(prime.json,
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

// Initial model fitting completed
// PRIME model set-up

// TODO: This is where we need to implement the PRIME model
prime.site.model = model.generic.DefineModel("models.codon.QCAP.ModelDescription",
        "prime_model", {
            "0": parameters.Quote(terms.local),
            "1": prime.codon_data_info[terms.code]
        },
        prime.filter_names,
        None);

// set global parameters 
prime.setGlobalParameters();

prime.site_model_mapping = {"prime_model" : prime.site.model};

// Report variable initialization
selection.io.startTimer (prime.json [terms.json.timers], "PRIME analysis", 2);

// Iterate through each partition and spawn optimization routines
for (prime.partition_index = 0; prime.partition_index < prime.partition_count; prime.partition_index += 1) {

    prime.report.header_done = FALSE;
    prime.table_output_options[utility.getGlobalValue("terms.table_options.header")] = TRUE;

    model.ApplyModelToTree("prime.site_tree_fel", prime.trees[prime.partition_index], {terms.default : prime.site.mg_rev}, None);
    prime.site_patterns = alignments.Extract_site_patterns ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")]);


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
                                   terms.mpi.Functions : {}
                                 });


    //// SW 20180221 : Need to include entire model
    //constraints = {
    //    "0" : {"name": "alpha_0", "constraint" : "alpha_0 := 0", "cleanup" : "alpha_0 = 1"},
    //    "1" : {"name": "alpha_1", "constraint" : "alpha_1 := 0", "cleanup" : "alpha_1 = 1"},
    //    "2" : {"name": "alpha_2", "constraint" : "alpha_2 := 0", "cleanup" : "alpha_2 = 1"},
    //    "3" : {"name": "alpha_3", "constraint" : "alpha_3 := 0", "cleanup" : "alpha_3 = 1"},
    //    "4" : {"name": "alpha_4", "constraint" : "alpha_4 := 0", "cleanup" : "alpha_4 = 1"}
    //};

    utility.ForEachPair (prime.site_patterns, "_pattern_", "_pattern_info_",
        '

            for (i = 0; i < 5; i += 1) {
                if (_pattern_info_[terms.data.is_constant]) {
                } else {
                    mpi.QueueJob (prime.queue, "prime.handleSite", {"0" : "prime.site_likelihood",
                                                                     "1" : alignments.serialize_site_filter
                                                                       ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")],
                                                                       (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                     "2" : prime.partition_index,
                                                                     "3" : _pattern_info_,
                                                                     "4" : prime.site_model_mapping,
                                                                     "5" : constraints[i]
                                                                        },
                                                                        "prime.storeResults");
                }
            }
        '
    );

    mpi.QueueComplete (prime.queue);

    prime.partition_matrix = {Abs (prime.site_results[prime.partition_index]), Rows (prime.table_headers)};

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


function prime.setDefaultVariables() {

    /* Display analysis information */

    prime.analysis_description = {
        terms.io.info: "PRIME",
        terms.io.version: "1.00",
        terms.io.reference: "TBD",
        terms.io.authors: "Sergei L Kosakovsky Pond",
        terms.io.contact: "spond@temple.edu",
        terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
    };


    /* Environment Setup */
    utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

    /* Globals */
    prime.site_alpha = "Site relative synonymous rate";
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


    prime.site_results = {};
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

}

function prime.setGlobalParameters() {

    // TODO: More descriptive names
    prime.site_alpha_0 = "alpha_0";
    prime.site_alpha_1 = "alpha_1";
    prime.site_alpha_2 = "alpha_2";
    prime.site_alpha_3 = "alpha_3";
    prime.site_alpha_4 = "alpha_4";

    model.generic.AddGlobal (prime.site.model, "prime.alpha_0", prime.site_alpha_0);
    model.generic.AddGlobal (prime.site.model, "prime.alpha_1", prime.site_alpha_1);
    model.generic.AddGlobal (prime.site.model, "prime.alpha_2", prime.site_alpha_2);
    model.generic.AddGlobal (prime.site.model, "prime.alpha_3", prime.site_alpha_3);
    model.generic.AddGlobal (prime.site.model, "prime.alpha_4", prime.site_alpha_4);
    
    model.generic.AddGlobal (prime.scaler, "prime.scaler", prime.site_alpha);
    parameters.DeclareGlobal (prime.scalers, {"prime.scaler"});

    // Applies prime.scaler constraint to all branches
    prime.applyProportianalSiteConstraints();

}

function prime.applyProportianalSiteConstraints () {

    // apply constraints to the site tree
    // alpha = alpha_scaler * branch_length
    utility.ForEach (prime.case_respecting_node_names, "_node_",
        '
            // Get branch length from mg_rev fit
            branch_length = (( prime.final_partitioned_mg_results[terms.branch_length])[prime.partition_index])[_node_];
            prime.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];
            node_name = "prime.site_tree." + node_name;

            prime_scaler_label = "prime.scaler";

            ExecuteCommands ("
                `node_name`.`prime_scaler_label` := (`prime.scaler`) * prime.branch_length__;
            ");

        ');

}

// Performs MLE optimization per site. This function is executed in an embarrassingly parallel way.
lfunction prime.handleSite (lf, filter_data, partition_index, pattern_info, model_mapping, test_condition) {

    //GetString (lfInfo, ^lf,-1);
    //ExecuteCommands (filter_data);

    //__make_filter ((lfInfo["Datafilters"])[0]);

    //utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    //// declares constraint
    //ExecuteCommands(test_condition["constraint"]);

    //Optimize (results, ^lf);

    //alpha_fit = estimators.ExtractMLEs (lf, model_mapping);
    //alpha_fit[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];

    //fit = {
    //    "name" : test_condition["name"],
    //    "fit" : alpha_fit
    //};

    //ExecuteCommands(test_condition["cleanup"]);

    //return fit;
    
}

lfunction prime.storeResults(node, result, arguments) {

    partition_index = arguments [2];
    pattern_info = arguments [3];

    fit_name = result["name"];
    fit_results = result["fit"];

    result_row = {{ 0, 0 }};

    // not a constant site
    if (None != result) { 
        result_row[0] = estimators.GetGlobalMLE(fit_results, ^"prime.site_alpha");
    }

    utility.EnsureKey(^"prime.site_results", partition_index);

    utility.ForEach(pattern_info[utility.getGlobalValue("terms.data.sites")], "_prime_result_",
        '
            (prime.site_results[`&partition_index`])[_prime_result_] = `&result_row`;
        '
    );

    // Only report to screen when all tests have been completed
    prime.reportResults(prime.site_results);

}

function prime.reportResults(prime.site_results) {

    //utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_fel_result_",
    //    '
    //        (fel.site_results[`&partition_index`])[_fel_result_] = `&result_row`;
    //        fel.report.echo (_fel_result_, `&partition_index`, `&result_row`);
    //    '
    //);

    //prime.print_row = None;

    //if (prime.report.row [4] < prime.pvalue) {
    //    if (prime.report.row[0] < prime.report.row[1]) {
    //        prime.print_row = prime.report.positive_site;
    //        prime.report.counts[0] += 1;
    //    } else {
    //        prime.print_row = prime.report.negative_site;
    //        prime.report.counts [1] += 1;
    //    }
    //}

    //if (None != prime.print_row) {
    //    if (!prime.report.header_done) {
    //        io.ReportProgressMessageMD("prime", "" + prime.report.partition, "For partition " + (prime.report.partition+1) + " these sites are significant at p <=" + prime.pvalue + "\n");
    //        fprintf (stdout,
    //            io.FormatTableRow (prime.table_screen_output,prime.table_output_options));
    //        prime.report.header_done = TRUE;
    //        prime.table_output_options[terms.table_options.header] = FALSE;
    //    }
    //    fprintf (stdout,
    //        io.FormatTableRow (prime.print_row,prime.table_output_options));
    //}

}

