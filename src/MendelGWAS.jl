"""
This module orchestrates a GWAS analysis.
"""
module MendelGWAS
#
# Required OpenMendel packages and modules.
#
using MendelBase
# using DataStructures                  # Now in MendelBase.
# using GeneticUtilities                # Now in MendelBase.
# using GeneralUtilities                # Now in MendelBase.
using SnpArrays
#
# Required external modules.
#
using DataFrames                        # From package DataFrames.
using Distributions                     # From package Distributions.
using GLM                               # From package GLM.
using StatsBase                         # From package StatsBase.

import Plots

export GWAS

"""
This is the wrapper function for the GWAS analysis option.
"""
function GWAS(control_file = ""; args...)

  const GWAS_VERSION :: VersionNumber = v"0.1.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("      GWAS analysis option")
  println("        version ", GWAS_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["distribution"] = "" # Binomial(), Gamma(), Normal(), Poisson(), etc.
  keyword["link"] = "" # LogitLink(), IdentityLink(), LogLink(), etc.
  keyword["lrt_threshold"] = 5e-8
  keyword["maf_threshold"] = 0.01
  keyword["manhattan_plot_file"] = ""
  keyword["regression"] = ""   # Linear, Logistic, or Poisson
  keyword["regression_formula"] = ""
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
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = gwas_option(person, snpdata, pedigree_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
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
function gwas_option(person::Person, snpdata::SnpData,
  pedigree_frame::DataFrame, keyword::Dict{ASCIIString, Any})

  people = person.people
  snps = snpdata.snps
  io = keyword["output_unit"]
  lrt_threshold = keyword["lrt_threshold"]
  maf_threshold = keyword["maf_threshold"]
  #
  # Recognize the three basic GWAS regression models:
  # Linear, Logistic, and Poisson.
  # For these three we will use fast internal regression code.
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
  # Retrieve the regression formula and create a model frame.
  #
  regression_formula = keyword["regression_formula"]
  if regression_formula == ""
    throw(ArgumentError(
      "The keyword regression_formula appears blank.\n" *
      "A regression formula must be provided.\n \n"))
  end
  if !contains(regression_formula, "~")
    throw(ArgumentError(
      "The value of the keyword regression_formula ($regression_formula)\n" *
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
      "This should be the name of the trait field in the Pedigree file.\n \n"))
  end
  if side[2] == ""; side[2] = "1"; end
  lhs = parse(side[1])
  rhs = parse(side[2])
  if !(lhs in names(pedigree_frame))
    lhs_string = string(lhs)
    throw(ArgumentError(
      "The field named on the left hand side of the formula specified in\n" *
      "the keyword regression_formula (currently: $lhs_string)\n" *
      "is not in the Pedigree data file.\n \n"))
  end
  fm = Formula(lhs, rhs)
  model = ModelFrame(fm, pedigree_frame)
  #
  # Change sex designations to -1.0 (females) and +1.0 (males).
  # Since the field :sex may have type string,
  # create a new field of type Float64 that will replace :sex.
  #
  if searchindex(string(rhs), "Sex") > 0 && in(:Sex, names(model.df))
    model.df[:NewSex] = ones(person.people)
    for i = 1:person.people
      s = model.df[i, :Sex]
      if !isa(parse(string(s), raise=false), Number); s = lowercase(s); end
      if !(s in keyword["male"]); model.df[i, :NewSex] = -1.0; end
    end
    names_list = names(model.df)
    deleteat!(names_list, findin(names_list, [:Sex]))
    model.df = model.df[:, names_list]
    rename!(model.df, :NewSex, :Sex)
  end
  #
  # For Logistic regression make sure the cases are 1.0,
  # non-cases are 0.0, and missing data is NaN.
  # Again, since the trait field may be of type string,
  # create a new field of type Float64 that will replace it.
  #
  case_label = keyword["affected_designator"]
  if regression_type == "logistic" && case_label != ""
    model.df[:NewTrait] = zeros(person.people)
    for i = 1:person.people
      s = string(model.df[i, lhs])
      if s == ""
        model.df[i, :NewTrait] = NaN
      elseif s == case_label
        model.df[i, :NewTrait] = 1.0
      end
    end
    names_list = names(model.df)
    deleteat!(names_list, findin(names_list, [lhs]))
    model.df = model.df[:, names_list]
    rename!(model.df, :NewTrait, lhs)
  end
  #
  # To ensure that the trait and SNPs occur in the same order, sort the
  # the model dataframe by the entry-order column of the pedigree frame.
  #
  model.df[:EntryOrder] = pedigree_frame[:EntryOrder]
  sort!(model.df, cols = [:EntryOrder])
  names_list = names(model.df)
  deleteat!(names_list, findin(names_list, [:EntryOrder]))
  model.df = model.df[:, names_list]
  #
  # First consider the base model with no SNPs included.
  # If using one of the standard regression models,
  # then use the fast regression code in OpenMendel's general utilities file.
  # Otherwise, for general distributions, use the code in the GLM package.
  #
  if fast_method
    #
    # Copy the complete rows into a design matrix X and a response vector y.
    #
    mm = ModelMatrix(model) # Extract the model matrix with missing values.
    complete = complete_cases(model.df)
    cases = sum(complete)
    predictors = size(mm.m, 2)
    X = zeros(cases, predictors)
    for j = 1:predictors
      X[:, j] = mm.m[complete, j]
    end
    y = zeros(cases)
    y[1:end] = model.df[complete, lhs]
    #
    # Estimate parameters under the base model. Output results.
    #
    (base_estimate, base_loglikelihood) = regress(X, y, regression_type)
    println(io, " ")
    println(io, "Summary for Base Model with ", fm)
    println(io, "Regression Model: ", regression_type)
    println(io, "Link Function: ", "canonical")
    names_list = names(model.df)
    model_names = size(names_list,1)
    outcome_index = findin(names_list, [lhs])[1]
    println(io, "Base Components Effect Estimates: ")
    println(io, "   (Intercept) : ", signif(base_estimate[outcome_index], 6))
    for j = 1:model_names
      if j != outcome_index
        println(io, "   ", names_list[j], " : ", signif(base_estimate[j], 6))
      end
    end
    println(io, "Base Loglikelihood: ", signif(base_loglikelihood, 8))
    println(io, " ")
  else
    #
    # Let the GLM package estimate parameters. Output results.
    #
    base_model = glm(fm, model.df, distribution_family, link)
    println(io, "Summary for Base Model:\n ")
    print(io, base_model)
    println(io, " ")
  end
  #
  # Now consider the alternative model with a SNP included.
  # Add a column to the design matrix to hold the SNP dosages.
  # Change the regression formula to include the SNP.
  #
  X = [X zeros(cases)]
  if side[2] == ""
    rhs = parse("SNP")
  else
    rhs = parse(side[2] * " + " * "SNP")
  end
  fm = Formula(lhs, rhs)
  #
  # Analyze the SNP predictors one by one.
  #
  dosage = zeros(people)
  pvalue = ones(snps)
  for snp = 1:snps
    #
    # Ignore SNPs with MAF below the specified threshold.
    #
    if snpdata.maf[snp] <= maf_threshold; continue; end
    #
    # Copy the current SNP genotypes into a dosage vector.
    #
    copy!(dosage, slice(snpdata.snpmatrix, :, snp); impute = false)
    #
    # For the three basic regression types, analyze the alternative model
    # using internal score test code. If the score test p-value
    # is below the specified threshold, carryout a likelihood ratio test.
    #
    if fast_method
      X[:, end] = dosage[complete]
      estimate = [base_estimate; 0.0]
      score_test = glm_score_test(X, y, estimate, regression_type)
      pvalue[snp] = ccdf(Chisq(1), score_test)
      if pvalue[snp] < lrt_threshold
        (estimate, loglikelihood) = regress(X, y, regression_type)
        lrt = 2.0 * (loglikelihood - base_loglikelihood)
        pvalue[snp] = ccdf(Chisq(1), lrt)
      end
    #
    # For other distributions analyze the alternative model
    # using the GLM package.
    #
    else
      model.df[:SNP] = dosage
      snp_model = fit(GeneralizedLinearModel, fm, model.df,
        distribution_family, link)
      pvalue[snp] = coeftable(snp_model).cols[end][end].v
    end
    #
    # Output regression results for potentially significant SNPs.
    #
    if pvalue[snp] < 0.05 / snps
      println(io, " \n")
      println(io, "Summary for SNP ", snpdata.snpid[snp])
      println(io, " on chromosome ", snpdata.chromosome[snp],
        " at basepair ", snpdata.basepairs[snp])
      println(io, "SNP p-value: ", signif(pvalue[snp], 6))
      println(io, "Minor Allele Frequency: ", round(snpdata.maf[snp], 4))
      if uppercase(snpdata.chromosome[snp]) == "X"
        hw = xlinked_hardy_weinberg_test(dosage, person.male)
      else
        hw = hardy_weinberg_test(dosage)
      end
      println(io, "Hardy-Weinberg p-value: ", round(hw, 4))
      if fast_method
        println(io, "SNP Effect Estimate: ", signif(estimate[end], 4))
        println(io, "SNP Effect Loglikelihood: ", signif(loglikelihood, 8))
      else
        println(io, "")
        println(io, snp_model)
      end
    end
  end
  #
  # Output false discovery rates.
  #
  fdr = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
  (number_passing, threshold) = simes_fdr(pvalue, fdr, snps)
  println(io, " \n \n ")
  println(io, "        P-value   Number of Passing")
  println(io, "FDR    Threshold     Predictors \n")
  for i = 1:length(fdr)
    @printf(io,"%4.2f   %8.5f   %9i\n", fdr[i], threshold[i], number_passing[i])
  end
  println(io, " ")
  #
  # If requested, output a Manhattan Plot in .png format.
  # First, create a new dataframe that will hold the SNP data to plot.
  # Next, sort the data by chromosome and basepair location.
  #
  plot_file = keyword["manhattan_plot_file"]
  if plot_file != ""
    if !contains(plot_file, ".png"); string(plot_file, ".png"); end
    # generate dataframe for plotting
    df = DataFrame(
        NegativeLogPvalue = -log10(pvalue),
        # Basepairs = snpdata.basepairs,
        Chromosome = snpdata.chromosome,
        X = 1:length(pvalue)
    )
    # Get tick marks
    xticks = by(df, :Chromosome, df -> mean(df[:X]))
    # initialize plot
    plt = Plots.scatter(
        df[:NegativeLogPvalue], xlabel = "Chromosome", ylabel = "\$log_{10}(p-value)\$",
        group = df[:Chromosome], markersize = 3, markerstrokewidth = 0,
        legend = false, palette = :viridis,
        xticks = (sort(xticks[:x1].data)[1:2:end], 1:2:size(xticks, 1))
    )
    # add horizontal line for Bonferonni Correction
    Plots.abline!(plt, 0, -log10(.05 / size(df, 1)), color = :black, line = :dash)

    Plots.savefig(plt, plot_file)
    display(plt)
  end
  close(io)
  return execution_error = false
end # function gwas_option

end # module MendelGWAS
