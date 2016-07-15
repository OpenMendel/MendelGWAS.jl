"""
This module orchestrates a GWAS analysis.
"""
module MendelGWAS
#
# Other OpenMendel modules.
#
# using DataStructures
# using GeneticUtilities
using MendelBase
using SnpArrays
#
# External modules.
#
using DataFrames    # From package DataFrames.
using Distributions # From package Distributions.
using GLM           # From package GLM.
using PyPlot        # From package PyPlot.
using StatsBase     # From package StatsBase.

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
  # Keywords unique to this analysis may be defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = value
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
  #
  # Retrieve the regression formula and create a model frame.
  #
  regression_formula = keyword["regression_formula"]
  side = split(regression_formula, "~")
  lhs = parse(side[1])
  rhs = parse(side[2])
  fm = Formula(lhs, rhs)
  model = ModelFrame(fm, pedigree_frame)
  #
  # Change sex designations to -1 (females) and  +1 (males).
  #
  if searchindex(string(rhs), "Sex") > 0 && in(:Sex, names(model.df))
    for i = 1:person.people
      s = model.df[i, :Sex]
      if !isa(parse(string(s), raise=false), Number); s = lowercase(s); end
      if s in keyword["male"]
        model.df[i, :Sex] = 1.0
      else
        model.df[i, :Sex] = -1.0
      end
    end
  end
  #
  # Regress on the non-SNP predictors.
  #
  base_model = glm(fm, model.df, Normal(), IdentityLink())
  println(io, "Summary for base Model:")
  print(io, base_model)
  #
  # To ensure that the trait and SNPs occur in the same order, sort the
  # the model dataframe by the entry-order column of the pedigree frame.
  #
  model.df[:EntryOrder] = pedigree_frame[:EntryOrder]
  sort!(model.df, cols = [:EntryOrder]) 
  #
  # Add a column to the Pedigree frame to hold the SNP dosages.
  # Change the regression formula to include the current SNP.
  #
  model.df[:SNP] = @data(zeros(people))
  rhs = parse(side[2] * " + " * "SNP")
  fm = Formula(lhs, rhs)
  #
  # Analyze the SNP predictors one by one. These will be inserted
  # into the extra column of the Pedigree frame.
  #
  dosage = zeros(people)
  pvalue = ones(snps)
  for snp = 1:snps
    if snpdata.maf[snp] <= 0.01; continue; end
    #
    # Copy the current SNP genotypes into a dosage vector.
    #
    copy!(dosage, slice(snpdata.snpmatrix, :, snp); impute = false)
    #
    # Estimate parameters for the SNP model. Typical distributions for
    # use with glm and their canonical link functions are Binomial (LogitLink),
    # Gamma (InverseLink), Normal (IdentityLink), and Poisson (LogLink).
    # The currently available (as of 2015) Link types are CauchitLink, 
    # CloglogLink, IdentityLink, InverseLink, LogitLink, LogLink, ProbitLink, 
    # and SqrtLink.
    #
    model.df[:SNP] = dosage
    if keyword["distribution"] == Normal() && keyword["link"] == IdentityLink()
      snp_model = fit(LinearModel, fm, model.df)
    else
      snp_model = fit(GeneralizedLinearModel, fm, model.df,
        keyword["distribution"], keyword["link"])
    end
    pvalue[snp] = coeftable(snp_model).cols[end][end].v
    #
    # Output regression results for potentially significant SNPs.
    #
    if pvalue[snp] < 0.05 / snps
      println(io, "")
      println(io, "Summary for SNP ", snpdata.snpid[snp], ":")
      println(io, "Minor Allele Frequency: ", snpdata.maf)
      if uppercase(snpdata.chromosome[1]) == "X"
        hw = xlinked_hardy_weinberg_test(dosage, person.male)
      else
        hw = hardy_weinberg_test(dosage)
      end
      println(io, "Hardy-Weinberg Pvalue: ", round(hw, 5))
      println(io, "")
      println(io, snp_model)
    end
  end
  #
  # Output false discovery rates.
  #
  fdr = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
  (number_passing, threshold) = simes_fdr(pvalue, fdr, snps)
  println(io, "")
  println(io, )
  println(io, "             P-value    Number of Passing")
  println(io, "    FDR       Threshold     Predictors")
  for i = 1:length(fdr)
    println(io, "    ", round(fdr[i], 6), "      ", round(threshold[i], 6),
      "      ", number_passing[i])
  end
  #
  # Sort the SNPs by chromosome and basepair location.
  #
  manhattan_frame = DataFrame()
  manhattan_frame[:NegativeLogPvalue] = -log10(pvalue)
  manhattan_frame[:Basepairs] = snpdata.basepairs
  manhattan_frame[:Chromosome] = snpdata.chromosome
  sort!(manhattan_frame::DataFrame, cols = [:Chromosome, :Basepairs])
  #
  # Print a Manhattan plot of the pvalues.
  #
  x = collect(1:snps)
  y = manhattan_frame[:NegativeLogPvalue]
  plot(x, y, ".")
  title("Manhattan Plot of Negative Log P-values")
  xlabel("Chromosome Location")
  ylabel("-log p-value")
  plot_file = keyword["manhattan_plot_file"]
  plot_format = split(plot_file, ".")[2]
  savefig(plot_file, format = plot_format)
  close()
  return execution_error = false
end # function gwas_option

end # module MendelGWAS

