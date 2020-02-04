__precompile__()

"""
This module orchestrates a GWAS analysis.
"""
module MendelGWAS
#
# Required OpenMendel packages and modules.
#
using MendelBase
using MendelPlots
using SnpArrays
#
# Required external modules.
#
using CSV
using DataFrames
using Distributions
using GLM
using LinearAlgebra
using Missings
using Printf
using StatsBase
using StatsModels

export GWAS

"""
This is the wrapper function for the GWAS analysis option.
"""
function GWAS(control_file = ""; args...)
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("      GWAS analysis option")
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["distribution"] = "" # Binomial(), Gamma(), Normal(), Poisson(), etc.
  keyword["link"] = ""         # LogitLink(), IdentityLink(), LogLink(), etc.
  keyword["lrt_threshold"] = 5e-8
  keyword["maf_threshold"] = 0.01
  keyword["min_success_rate_per_sample"] = 0.98
  keyword["min_success_rate_per_snp"] = 0.98
  keyword["manhattan_plot_file"] = ""
  keyword["output_table"] = "" # Table of all SNPs and their p-values
  keyword["qq_plot_file"] = ""
  keyword["regression"] = ""   # linear, logistic, or poisson
  keyword["regression_formula"] = ""
  keyword["pcs"] = 0      # Number of Principal Components to include in model.
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "gwas")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "GWAS"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, person_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Check if SNP data were read.
  #
  if snpdata.snps == 0
    println(" \n\nERROR: This analysis requires SNP data and none were read!\n")
  else
  #
  # If principal components are requested to be included in the model,
  # then call a function that will use the PCA routine in SnpArrays,
  # and will add these PCs to the person_frame. Once in the pedigree frame
  # these PCs will be included in the regression procedure.
  # Each PC will be named PCn where n is its number, e.g., "PC1".
  #
##!!  if keyword["pcs"] > 0
##!!    add_pcs!(person_frame, snpdata, keyword)
##!!  end
  #
  # Execute the specified analysis.
  #
    println(" \nAnalyzing the data.\n")
    execution_error = gwas_option(person, snpdata, person_frame, keyword)
    if execution_error
      println(" \n \nERROR: Mendel terminated prematurely!\n")
    else
      println(" \n \nMendel's analysis is finished.\n")
    end
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function GWAS

"""
This function performs GWAS on a set of traits.
"""
function gwas_option(person::Person, snpdata::SnpDataStruct,
  person_frame::DataFrame, keyword::Dict{AbstractString, Any})

  TAB_CHAR :: Char = Char(9)

  people = person.people
  snps = snpdata.snps
  io = keyword["output_unit"]
  lrt_threshold = keyword["lrt_threshold"]
  maf_threshold = keyword["maf_threshold"]
  min_success_rate_per_sample = keyword["min_success_rate_per_sample"]
  min_success_rate_per_snp = keyword["min_success_rate_per_snp"]
  #
  # Recognize the three basic GWAS regression models:
  # linear, logistic, and poisson.
  # For these three we will use fast internal regression code,
  # unless an interaction term is detected.
  #
  regression_type = lowercase(keyword["regression"])
  distribution_family = keyword["distribution"]
  link = keyword["link"]
  if regression_type == "linear" ||
     (distribution_family == Normal() && link == IdentityLink())
    keyword["distribution"] = Normal()
    keyword["link"] = IdentityLink()
    keyword["regression"] = "linear"
    fast_method = true
  elseif regression_type == "logistic" ||
     (distribution_family == Binomial() && link == LogitLink())
    keyword["distribution"] = Binomial()
    keyword["link"] = LogitLink()
    keyword["regression"] = "logistic"
    fast_method = true
  elseif regression_type == "poisson" ||
     (distribution_family == Poisson() && link == LogLink())
    keyword["distribution"] = Poisson()
    keyword["link"] = LogLink()
    keyword["regression"] = "poisson"
    fast_method = true
  elseif regression_type == "" && distribution_family == "" && link == ""
    throw(ArgumentError(
      "No regression type has been defined for this analysis.\n" *
      "Set the keyword regression to either:\n" *
      "  linear (for usual quantitative traits),\n" *
      "  logistic (for usual qualitative traits),\n" *
      "  poisson, or blank (for unusual analyses).\n" *
      "If keyword regression is not assigned, then the keywords\n" *
      "'distribution' and 'link' must be assigned names of functions.\n \n"))
  else
    fast_method = false
  end
  regression_type = lowercase(keyword["regression"])
  distribution_family = keyword["distribution"]
  link = keyword["link"]
  if regression_type != "" && regression_type != "linear" &&
     regression_type != "logistic" && regression_type != "poisson"
    throw(ArgumentError(
     "The keyword regression (currently: $regression_type) must be assigned\n" *
     "the value 'linear', 'logistic', 'poisson', or blank.\n \n"))
  end
  #
  # Retrieve and error check the regression formula.
  #
  regression_formula = keyword["regression_formula"]
  if regression_formula == ""
    throw(ArgumentError(
      "The keyword regression_formula appears blank.\n" *
      "A regression formula must be provided.\n \n"))
  end
  if !occursin("~", regression_formula)
    throw(ArgumentError(
      "The value of the keyword regression_formula ('$regression_formula')\n" *
      "does not contain the required '~' that should separate\n" *
      "the trait and the predictors.\n \n"))
  end
  side = split(regression_formula, "~")
  side[1] = strip(side[1])
  side[2] = strip(side[2])
  if side[1] == ""
    throw(ArgumentError(
      "The left hand side of the formula specified in\n" *
      "the keyword regression_formula appears blank.\n" *
      "The left hand side should contain only the name of the trait field\n" *
      "in the Pedigree data. If the trait field is unnamed, use 'Trait'.\n \n"))
  end
  entry_separators = [' ', TAB_CHAR, ',', ';', '+', '*', '&']
  if something(findfirst((in)(entry_separators), side[1]), 0) > 0
    lhs_string = side[1]
    throw(ArgumentError(
      "The left hand side ('$lhs_string') of the formula specified in\n" *
      "the keyword regression_formula appears to have multiple entries.\n" *
      "The left hand side should contain only the name of the trait field\n" *
      "in the Pedigree data. If the trait field is unnamed, use 'Trait'.\n \n"))
  end
  if side[2] == ""; side[2] = "1"; end
  lhs = Meta.parse(side[1])
  rhs = Meta.parse(side[2])
  if !(lhs in names(person_frame))
    lhs_string = string(lhs)
    throw(ArgumentError(
      "The field named on the left hand side of the formula specified in\n" *
      "the keyword regression_formula (currently: '$lhs_string')\n" *
      "is not in the Pedigree data file.\n" *
      "The left hand side should contain only the name of the trait field\n" *
      "in the Pedigree data. If the trait field is unnamed, use 'Trait'.\n \n"))
  end
  #
  # If the regression formula includes an interaction term,
  # do not use the fast internal regression code.
  #
  interaction_symbols = ['*', '&']
  if something(findfirst((in)(interaction_symbols), string(rhs)), 0) > 0
    fast_method = false
  end
  #
  # Create a list of the names of the columns used in the formula.
  # NB: Since this procedure only looks for string matches
  # it will be flummoxed by, for example, a column named "l",
  # if the "log" function is used in the formula.
  #
  namecolumns = names(person_frame)
  totalcolumns = size(namecolumns, 1)
  indicat = falses(totalcolumns)
  for i = 1:totalcolumns
    namestring = string(namecolumns[i])
    if namestring == lhs ||
       first(something(findfirst(namestring, string(rhs)), 0:-1)) > 0
      indicat[i] = true
    end
  end
  names_in_formula = fill(:placeholder, sum(indicat))
  j = 0
  for i = 1:totalcolumns
    if indicat[i]
      j = j + 1
      names_in_formula[j] = namecolumns[i]
    end
  end
  if j != sum(indicat)
    throw(ArgumentError("Inconsistency in formula column count./n /n"))
  end
  #
  # Change sex designations to 1.0 (for females) and -1.0 (for males).
  # Since the field :Sex may have type String,
  # create a new field of type Float64 that replaces :Sex.
  # Note that unrecognized sexes are labeled female.
  #
  if first(something(findfirst("Sex", string(rhs)), 0:-1)) > 0 &&
     in(:Sex, names(person_frame))
    person_frame[!, :NumericSex] = ones(people)
    for i = 1:people
      s = person_frame[i, :Sex]
      if !isa(Meta.parse(string(s), raise=false), Number); s = lowercase(s); end
      if !(s in keyword["female"]); person_frame[i, :NumericSex] = -1.0; end
    end
    names_list = names(person_frame)
    deleteat!(names_list, findall((in)([:Sex]), names_list))
    select!(person_frame, names_list)
    rename!(person_frame, :NumericSex => :Sex)
  end
  #
  # For Logistic regression make sure the cases are 1.0,
  # non-cases are 0.0, and missing data is "missing".
  # Again, since the trait field may be of type String,
  # create a new field of type Float64 that replaces it.
  #
  case_label = keyword["affected_designator"]
  if regression_type == "logistic" && case_label != ""
    person_frame[!, :NumericTrait] = zeros(people)
    for i = 1:people
      s = strip(string(person_frame[i, lhs]))
      if s == "" || s == "NaN" || s == "missing" || s == "NA"
        person_frame[i, :NumericTrait] = missing
      elseif s == case_label
        person_frame[i, :NumericTrait] = 1.0
      end
    end
    names_list = names(person_frame)
    deleteat!(names_list, findall((in)([lhs]), names_list))
    select!(person_frame, names_list)
    rename!(person_frame, :NumericTrait => lhs)
  end
  #
  # To filter the SNP data, first find the SNPs and samples
  # that surpass the requested success rates and MAF.
  # Note that the sample_mask is on the people in entry order!
  #
  sample_mask, snp_mask = SnpArrays.filter(snpdata.snpmatrix,
    min_success_rate_per_row = min_success_rate_per_snp,
    min_success_rate_per_col = min_success_rate_per_sample,
    min_maf = maf_threshold)
  #
  # Create a model frame, reordering the individuals (i.e., the rows)
  # by their entry order to align with the corresponding SNP data.
  #
  model_frame = deepcopy(person_frame)
  sort!(model_frame, :EntryOrder)
  #
  # Restrict the model_frame to the columns used in the regression formula.
  # (Perhaps wait on this until we know it's working well.)
  #
##  select!(model_frame, names_in_formula)
  #
  # Remove the individuals with low genotyping success rates
  # using the sample_mask from SnpArrays.
  # Note that now rows = people - too.few.genotypes.
  #
  model_frame = model_frame[sample_mask, :]
  #
  # Next, drop all individuals with missing values among the variables
  # used in the regression formula. Keep track of the corresponding mask
  # to use with the SNP data. Thus, model_frame has no missing data.
  # Note that the completeness_mask should be applied after the sample_mask.
  # Note that now rows = people - too.few.genotypes - incomplete.predictors.
  #
  completeness_mask = completecases(model_frame, names_in_formula)
  dropmissing!(model_frame, names_in_formula, disallowmissing=true)
  cases = size(model_frame, 1) # also = sum(completeness_mask)
##
##   model = ModelFrame(fm, person_frame[sample_mask, :])
##   #
##   # To ensure that the trait and SNPs occur in the same order, sort the
##   # the model dataframe by the entry-order column of the pedigree frame.
##   #
##   model.df[!, :EntryOrder] = person_frame[sample_mask, :EntryOrder]
##   sort!(model.df, :EntryOrder)
##   names_list = names(model.df)
##   deleteat!(names_list, findall((in)([:EntryOrder]), names_list))
##   select!(model.df, names_list)
##   #
##   # Find which predictors have complete data.
##   # Note that the completeness_mask should be applied after the sample_mask.
##   #
##   completeness_mask = completecases(model.df)
  #
  # First consider the base model with no SNPs included.
  # If using one of the standard regression models without interactions,
  # then use the fast regression code in OpenMendel's general utilities file.
  # Otherwise, for general distributions, use the code in the GLM package.
  #
  if fast_method
    #
    # Create a model matrix using the regression formula and the model_frame.
    # The model_frame is just a subset and rearrangement of the person_frame
    # (a subset due to minimum genotype success rates, and a rearrangement
    # due to the need to match the original SNP Array data order).
    # The model matrix has columns for the intercept, any interaction terms,
    # and any variable functions (such as "log(BMI)").
    # The functions and macros used here are in the StatsModels package.
    # Here y is the vector of responses and X is the matrix of predictors.
    #
    fm = @eval(@formula($lhs ~ $rhs))
##    fm = @eval(@formula($side[1] ~ $side[2]))
    schem = schema(fm, model_frame)
    fs = apply_schema(fm, schem)
    predictor_names = coefnames(fs)[2]
    response, X = modelcols(fs, model_frame)
    predictors = size(X, 2)
    y = zeros(cases)
    y[1:end] = model_frame[!, lhs]
    if predictors != size(predictor_names, 1)
      throw(ArgumentError( "Inconsistency in number of predictors.\n \n"))
    end
##    #
##    # Create the model matrix, initially still with any missing predictor values.
##    # Copy only the complete rows into design matrix X and response vector y,
##    # thus they have size = people - too.few.genotypes - incomplete.predictors.
##    #
##    modelmatrx = ModelMatrix(model)
##    cases = sum(completeness_mask)
##    predictors = size(modelmatrx.m, 2)
##    X = zeros(cases, predictors)
##    base_estimate = zeros(predictors)
##    for j = 1:predictors
##      X[:, j] = modelmatrx.m[completeness_mask, j]
##    end
##    y = zeros(cases)
##    y[1:end] = model.df[completeness_mask, lhs]
    #
    # Estimate parameters under the base model.
    # Then, for linear regression, record the vector of residuals.
    #
    base_estimate = zeros(predictors)
    (base_estimate, base_loglikelihood) = fast_regress(X, y, regression_type)
    if regression_type == "linear"
      base_residual = y - (X * base_estimate)
    end
    #
    # Output the results of the base model.
    #
    println(io, "\nResults for Base Model:\n  ", fm, "\n")
    println(io, "Regression model: ", regression_type)
    println(io, "Link function: ", "canonical")
    println(io, "Base components' effect estimates: ")
    for j = 1:predictors
        println(io, "   ", string(predictor_names[j]), " : ",
          round(base_estimate[j], sigdigits = 6))
    end
    println(io, "Base model loglikelihood: ",
      round(base_loglikelihood, sigdigits = 8))
    println(io, " ")
  else
    #
    # Let the GLM package estimate parameters. Output results.
    # Note that model_frame only contains complete cases.
    #
    fm = @eval(@formula($lhs ~ $rhs))
    base_model = glm(fm, model_frame, distribution_family, link)
    println(io, "\nResults for Base Model:\n  ", base_model.mf.f, "\n")
    print(io, coeftable(base_model))
    println(io, "\n \n")
  end
  #
  # Now consider the alternative model with a SNP included.
  # Add a column to the design matrix to hold the SNP dosages.
  # Change the regression formula to include the SNP.
  #
  if fast_method; X = [X zeros(cases)]; end
  rhs = Meta.parse(side[2] * " + SNP")
  fm = @eval(@formula($lhs ~ $rhs))
  #
  # Analyze the SNP predictors one by one.
  #
  dosage = zeros(sum(sample_mask))
  pvalue = ones(snps)
  if fast_method
    alt_estimate = zeros(predictors+1)
    if regression_type == "linear"
      alt_residual = zeros(cases)
    end
  end
  skipped_snps = 0
  signif_snps = 0
  signif_threshold =  0.05 / snps
  chr_max_bp = zeros(Integer, 50) # Assumes number of chromosomes <= 50.
  current_chr_number = 1
  current_chr = snpdata.chromosome[1]
  for snp = 1:snps
    #
    # Find maximum basepair listed for each chromosome.
    # Note that this section assumes SNPs on the same chromosome
    # are listed in a contigous block.
    #
    if snpdata.chromosome[snp] != current_chr
      current_chr = snpdata.chromosome[snp]
      current_chr_number = current_chr_number + 1
    end
    if !ismissing(snpdata.basepairs[snp])
      chr_max_bp[current_chr_number] =
        max(chr_max_bp[current_chr_number], snpdata.basepairs[snp])
    end
    #
    # Skip analysis of SNPs with MAF or genotype success rate
    # below the specified thresholds.
    #
    if !snp_mask[snp]
      skipped_snps = skipped_snps + 1
      continue
    end
    #
    # Copy the filtered SNP genotypes into a dosage vector,
    # allowing missing genotypes to be simplistically imputed based on MAF.
    #
    copyto!(dosage, view(snpdata.snpmatrix, sample_mask, snp); impute = true)
    #
    # For the three basic regression types, analyze the alternative model
    # using internal score test code. If the score test p-value
    # is below the specified threshold, carry out a likelihood ratio test.
    #
    if fast_method
##      copyto!(X[:, end], dosage[completeness_mask])
##      copyto!(X[:, end], view(dosage, completeness_mask))
##      X[:, end] = deepcopy(dosage[completeness_mask])
##      X[:, end] .= view(dosage, completeness_mask)
      X[:, end] = dosage[completeness_mask]
##if snp == 1; println(" X[:, end] = ", typeof(X), " :: ", X[:, end]); end
      alt_estimate = [base_estimate; 0.0]
      score_test = fast_score_test(X, y, alt_estimate, regression_type)
      pvalue[snp] = ccdf(Chisq(1), score_test)
      if pvalue[snp] < lrt_threshold
        (alt_estimate, alt_loglikelihood) = fast_regress(X, y, regression_type)
        lrt = 2.0 * (alt_loglikelihood - base_loglikelihood)
        pvalue[snp] = ccdf(Chisq(1), lrt)
      end
      #
      # Record the vector of residuals under this alternative model:
      # alt_residual = y - (X * alt_estimate)
      #
      if regression_type == "linear"
        mul!(alt_residual, X, alt_estimate) # alt_residual = X * alt_estimate
        alt_residual .= y .- alt_residual
      end
    #
    # For other distributions, analyze the alternative model using GLM package.
    # Add the (abridges) SNP data to the design matrix and then fit the model.
    #
    else
      model_frame[!, :SNP] = dosage[completeness_mask]
      snp_model = fit(GeneralizedLinearModel, fm, model_frame,
        distribution_family, link)
      #
      # Find which row in the output holds the SNP predictor results
      # and then store the p-value.
      #
      sm_predictor_names = coeftable(snp_model).rownms
      sm_predictors = size(sm_predictor_names, 1)
      snp_predictor = 0
      for i = 1:sm_predictors
        if sm_predictor_names[i] == "SNP"
          snp_predictor = i
          continue
        end
      end
      if snp_predictor == 0
        throw(ArgumentError( "Predictor named SNP not found in model.\n \n"))
      end
      pvalue[snp] = coeftable(snp_model).cols[4][snp_predictor]
    end
    #
    # Output regression results for potentially significant SNPs.
    #
    if pvalue[snp] < signif_threshold
      signif_snps = signif_snps + 1
      if signif_snps == 1
        println(io, "\n\nSummary Results for Alternative Model:")
        if fast_method
          println(io, "  ", fm)
        else
          println(io, "  ", snp_model.mf.f)
        end
        println(io, "Using significance threshold ", signif_threshold)
      end
      println(io, "\nResults for SNP ", snpdata.snpid[snp])
      println(io, " on chromosome ", snpdata.chromosome[snp],
        " at basepair ", snpdata.basepairs[snp])
      println(io, "SNP p-value: ", round(pvalue[snp], sigdigits = 4))
      println(io, "Minor allele frequency: ",
        round(snpdata.maf[snp], sigdigits = 4))
      if uppercase(snpdata.chromosome[snp]) == "X"
        hw = xlinked_hardy_weinberg_test(dosage, person.male)
      else
        hw = hardy_weinberg_test(dosage)
      end
      println(io, "Hardy-Weinberg p-value: ", round(hw, sigdigits = 4))
      if fast_method && pvalue[snp] < lrt_threshold
        println(io, "SNP effect estimate: ",
          round(alt_estimate[end], sigdigits = 4))
        println(io, "SNP model loglikelihood: ",
          round(alt_loglikelihood, sigdigits = 8))
      elseif fast_method
        println(io, "SNP effect estimate: ",
          round(alt_estimate[end], sigdigits = 4))
      else
        println(io, "\n", coeftable(snp_model), "\n")
      end
      #
      # For linear models, output the proportion of the base model's variance
      # that is explained by including this SNP in the model.
      #
      if fast_method && regression_type == "linear"
        variance_explained = 1 - (norm(alt_residual)^2 / norm(base_residual)^2)
        println(io, "Proportion of base-model variance explained: ",
          round(variance_explained, sigdigits = 4))
      end
    end
##if snp == 1; return; end
  end
  #
  # Report if there were no significant SNPs.
  #
  if signif_snps == 0
    println(io, "\n\nSummary Results for Alternative Model:")
    if fast_method
      println(io, "  ", fm)
    else
      println(io, "  ", snp_model.mf.f)
    end
    println(io, "Using significance threshold ", signif_threshold)
    println(io, "\nNo SNPs passed the significance threshold.\n")
  end
  #
  # Output false discovery rates.
  #
  fdr = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
  (number_passing, threshold) =
    simes_fdr(pvalue[snp_mask], fdr, snps - skipped_snps)
  println(io, " \n \n ")
  println(io, "        P-value   Number of Passing")
  println(io, "FDR    Threshold     Predictors \n")
  for i = 1:length(fdr)
   @printf(io, "%4.2f   %8.5f   %9i\n", fdr[i], threshold[i], number_passing[i])
  end
  println(io, " ")
  #
  # Create a new dataframe that will hold the data to output.
  #
  output_frame = DataFrame(SNP = snpdata.snpid,
    Chromosome = snpdata.chromosome,
    BasePair = snpdata.basepairs,
    Pvalue = pvalue,
    NegLog10Pvalue = -log10.(pvalue))
  #
  # If requested, output a full listing of the p-values in csv format.
  #
  table_file = string(keyword["output_table"])
  if table_file != ""
    CSV.write(table_file, output_frame;
      writeheader = true, delim = keyword["output_field_separator"],
      missingstring = keyword["output_missing_value"])
    println(" \nTable of p-values for all SNPs is in file: ", table_file)
  end
  #
  # If requested, output a Manhattan plot or QQ plot.
  #
  manhat_plot_file = string(keyword["manhattan_plot_file"])
  qq_plot_file = string(keyword["qq_plot_file"])
  if manhat_plot_file != "" || qq_plot_file != ""
    println(" \nPlots being drawn using MendelPlots package.")
  end
  if manhat_plot_file != ""
    manhattan(output_frame, outfile = manhat_plot_file,
##      linecolor = colorant"black",
      chrvar = "Chromosome", posvar = "BasePair", pvalvar = "Pvalue")
    println(" \nManhattan plot of GWAS results is in file: ", manhat_plot_file)
  end
  if qq_plot_file != ""
    qq(output_frame, outfile = qq_plot_file, pvalvar = "Pvalue")
    println(" \nQQ plot of GWAS results is in file: ", qq_plot_file)
  end
  return execution_error = false
end # function gwas_option

##!!"""
##!!Add principal components into a pedigree frame and regression formula.
##!!The PCs are calculated via pca() from the SnpArrays module.
##!!"""
##!!function add_pcs!(person_frame::DataFrame, snpdata::SnpDataStruct,
##!!  keyword::Dict{AbstractString, Any})
##!!
##!!  pcs = keyword["pcs"]
##!!  regress_form = keyword["regression_formula"]
##!!  #
##!!  # Perform PCA on the SNP data.
##!!  #
##!!  pcscore, pcloading, pcvariance = pca(snpdata.snpmatrix, pcs)
##!!  #
##!!  # Include the new covariates in the pedigree dataframe.
##!!  # NB: this process is *slow* for a large number of PCs!
##!!  #
##!!  for i = 1:pcs
##!!    person_frame[!, Symbol("PC$i")] = zscore(pcscore[:,i])
##!!    regress_form = regress_form * " + PC$i"
##!!  end
##!!  #
##!!  # Save the expanded regression formula to the keyword dictionary.
##!!  #
##!!  keyword["regression_formula"] = regress_form
##!!  return nothing
##!!end # function add_pcs!
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelGWAS
