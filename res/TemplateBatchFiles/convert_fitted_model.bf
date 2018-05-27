/*
    Convert a .fitted_model file, output by the protfitter, into  RAxML and/or PAML-formatted models.
    Users can therefore use this output to infer phylogenies with external programs.
*/
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/models/protein.bf");

convert_fitted_model.external_order       = {{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}};
convert_fitted_model.hyphy_order_dict = utility.MatrixToDict(models.protein.alphabet); 


namespace convert_fitted_model {
    custom_model_file = io.PromptUserForString("Provide name of model file returned by protfitter (extension `.fitted_model`)");
    LoadFunctionLibrary(custom_model_file); // custom_Rij; custom_EFV
}

convert_fitted_model.convert = io.SelectAnOption ({{"PAML", "PAML-format for use with IQTree, PhyML, PAML"},
                                     {"RAxML", "RAxML-format for use with RAxML, ExaML"},
                                     {"All", "Output both PAML- and RAxML-formatted models"}},
                                    "Select an output format:");

if (convert_fitted_model.convert == "PAML")
{
    convert_fitted_model.convert_to_paml();
}
if (convert_fitted_model.convert == "RAxML")
{
    convert_fitted_model.convert_to_raxml();
}
if (convert_fitted_model.convert == "All")
{
    convert_fitted_model.convert_to_raxml();
    convert_fitted_model.convert_to_paml();
}







function convert_fitted_model.convert_to_raxml(){

    convert_fitted_model.raxml_output = "";
    convert_fitted_model.raxml_file = convert_fitted_model.custom_model_file + ".RAXML";
    for (i = 0; i < 20; i +=1)
    {
        for (j = 0; j < 20; j += 1)
        {
        
            if (i == j)
            {
                convert_fitted_model.raxml_output += "0.0\n";
            }
            else
            {
                convert_fitted_model.aa1 = convert_fitted_model.external_order[i];
                convert_fitted_model.aa2 = convert_fitted_model.external_order[j];
        
                if (convert_fitted_model.aa1 < convert_fitted_model.aa2)
                {
                    convert_fitted_model.rate = (convert_fitted_model.custom_Rij[convert_fitted_model.aa1])[convert_fitted_model.aa2];
                }
                else
                {
                    convert_fitted_model.rate = (convert_fitted_model.custom_Rij[convert_fitted_model.aa2])[convert_fitted_model.aa1];
                }
                convert_fitted_model.raxml_output += convert_fitted_model.rate;
                convert_fitted_model.raxml_output += "\n"; //WTF WHY WONT YOU WORK
            }
        }    
    }
    for (i = 0; i < 20; i += 1)
    {
        convert_fitted_model.raxml_output += convert_fitted_model.custom_EFV[ convert_fitted_model.hyphy_order_dict[convert_fitted_model.external_order[i]] ];
        // strip hack
        if (i <= 18){
            convert_fitted_model.raxml_output += "\n";
        }
    }

    fprintf(convert_fitted_model.raxml_file, CLEAR_FILE, convert_fitted_model.raxml_output);
}




function convert_fitted_model.convert_to_paml(){

    convert_fitted_model.paml_output = "";
    convert_fitted_model.paml_file = convert_fitted_model.custom_model_file + ".PAML";

    for (i = 1; i < 20; i +=1)
    {
        convert_fitted_model.row = "";
        for (j = 0; j < i; j += 1)
        {
        
            convert_fitted_model.aa1 = convert_fitted_model.external_order[i];
            convert_fitted_model.aa2 = convert_fitted_model.external_order[j];
        
            if (convert_fitted_model.aa1 < convert_fitted_model.aa2){
                convert_fitted_model.rate = (convert_fitted_model.custom_Rij[convert_fitted_model.aa1])[convert_fitted_model.aa2];
            }
            else{
                convert_fitted_model.rate = (convert_fitted_model.custom_Rij[convert_fitted_model.aa2])[convert_fitted_model.aa1];
            }
            
            convert_fitted_model.row += convert_fitted_model.rate;
            // strip hack
            if (j != (i-1)){
                convert_fitted_model.row += " ";
            }
        }
    
        convert_fitted_model.paml_output += convert_fitted_model.row + "\n";
    }

    convert_fitted_model.paml_output += "\n";
    for (i = 0; i < 20; i += 1)
    {
        convert_fitted_model.paml_output += convert_fitted_model.custom_EFV[ convert_fitted_model.hyphy_order_dict[convert_fitted_model.external_order[i]] ];
        // strip hack
        if (i <= 18){
            convert_fitted_model.paml_output += " ";
        }
    }

    fprintf(convert_fitted_model.paml_file, CLEAR_FILE, convert_fitted_model.paml_output);
}
      