"""
This module orchestrates a GWAS analysis.
"""
module MendelGWAS
#
# Required OpenMendel packages and modules.
#
using MendelBase
using SnpArrays
#
# Required external modules.
#
using Compat
import Compat: view

using DataFrames                        # From package DataFrames.
using Distributions                     # From package Distributions.
using GLM                               # From package GLM.
using StatsBase                         # From package StatsBase.
#
# Use Plots as the plotting frontend and select a backend.
#
using Plots                             # From package Plots.
gr()
## pyplot()

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
  keyword["regression"] = ""   # Linear, Logistic, or Poisson
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
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # If principal components are requested to be included in the model,
  # then call a function that will use the PCA routine in SnpArrays,
  # and will add these PCs to the pedigree_frame. Once in the pedigree frame
  # these PCs will be included in the regression procedure.
  # Each PC will be named PCn where n is its number, e.g., "PC1".
  #
  if keyword["pcs"] > 0
    add_pcs!(pedigree_frame, snpdata, keyword)
  end
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
  pedigree_frame::DataFrame, keyword::Dict{AbstractString, Any})

  const TAB_CHAR :: Char = Char(9)

  people = person.people
  snps = snpdata.snps
  io = keyword["output_unit"]
  lrt_threshold = keyword["lrt_threshold"]
  maf_threshold = keyword["maf_threshold"]
  min_success_rate_per_sample = keyword["min_success_rate_per_sample"]
  min_success_rate_per_snp = keyword["min_success_rate_per_snp"]
  #
  # Recognize the three basic GWAS regression models:
  # Linear, Logistic, and Poisson.
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
      "  Poisson, or blank (for unusual analyses).\n" *
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
     "the value 'linear', 'logistic', 'Poisson', or blank.\n \n"))
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
  if !contains(regression_formula, "~")
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
      "in the Pedigree file. If the trait field is unnamed, use 'Trait'.\n \n"))
  end
  if search(side[1], [' ', TAB_CHAR, ',', ';', '+', '*', '&']) > 0
    lhs_string = side[1]
    throw(ArgumentError(
      "The left hand side ('$lhs_string') of the formula specified in\n" *
      "the keyword regression_formula appears to have multiple entries.\n" *
      "The left hand side should contain only the name of the trait field\n" *
      "in the Pedigree file. If the trait field is unnamed, use 'Trait'.\n \n"))
  end
  if side[2] == ""; side[2] = "1"; end
  lhs = parse(side[1])
  rhs = parse(side[2])
  if !(lhs in names(pedigree_frame))
    lhs_string = string(lhs)
    throw(ArgumentError(
      "The field named on the left hand side of the formula specified in\n" *
      "the keyword regression_formula (currently: '$lhs_string')\n" *
      "is not in the Pedigree data file.\n" *
      "The left hand side should contain only the name of the trait field\n" *
      "in the Pedigree file. If the trait field is unnamed, use 'Trait'.\n \n"))
  end
  #
  # If the regression formula includes an interaction term,
  # do not use the fast internal regression code.
  #
  if search(string(rhs), ['*', '&']) > 0; fast_method = false; end
  #
  # Change sex designations to 1.0 (females) and -1.0 (males).
  # Since the field :Sex may have type String,
  # create a new field of type Float64 that replaces :Sex.
  #
  if searchindex(string(rhs), "Sex") > 0 && in(:Sex, names(pedigree_frame))
    pedigree_frame[:NumericSex] = ones(people)
    for i = 1:people
      s = pedigree_frame[i, :Sex]
      if !isa(parse(string(s), raise=false), Number); s = lowercase(s); end
      if !(s in keyword["female"]); pedigree_frame[i, :NumericSex] = -1.0; end
    end
    names_list = names(pedigree_frame)
    deleteat!(names_list, findin(names_list, [:Sex]))
    pedigree_frame = pedigree_frame[:, names_list]
    rename!(pedigree_frame, :NumericSex => :Sex)
  end
  #
  # For Logistic regression make sure the cases are 1.0,
  # non-cases are 0.0, and missing data is "missing".
  # Again, since the trait field may be of type String,
  # create a new field of type Float64 that replaces it.
  #
  case_label = keyword["affected_designator"]
  if regression_type == "logistic" && case_label != ""
    pedigree_frame[:NumericTrait] = zeros(people)
    for i = 1:people
      s = strip(string(pedigree_frame[i, lhs]))
      if s == "" || s == "NaN" || s == "missing" || s == "NA"
        pedigree_frame[i, :NumericTrait] = missing
      elseif s == case_label
        pedigree_frame[i, :NumericTrait] = 1.0
      end
    end
    names_list = names(pedigree_frame)
    deleteat!(names_list, findin(names_list, [lhs]))
    pedigree_frame = pedigree_frame[:, names_list]
    rename!(pedigree_frame, :NumericTrait => lhs)
  end
  #
  # To filter the SNP data, first find the SNPs and samples
  # that surpass the requested success rates.
  #
  snp_mask, sample_mask = SnpArrays.filter(snpdata.snpmatrix,
    min_success_rate_per_snp, min_success_rate_per_sample)
  #
  # Create a model frame.
  # Note that the model data collections invoke the sample_mask,
  # thus have size = people - too.few.genotypes.
  #
  fm = Formula(lhs, rhs)
  model = ModelFrame(fm, pedigree_frame[sample_mask, :])
  #
  # To ensure that the trait and SNPs occur in the same order, sort the
  # the model dataframe by the entry-order column of the pedigree frame.
  #
  model.df[:EntryOrder] = pedigree_frame[sample_mask, :EntryOrder]
  sort!(model.df, cols = [:EntryOrder])
  names_list = names(model.df)
  deleteat!(names_list, findin(names_list, [:EntryOrder]))
  model.df = model.df[:, names_list]
  #
  # First consider the base model with no SNPs included.
  # If using one of the standard regression models without interactions,
  # then use the fast regression code in OpenMendel's general utilities file.
  # Otherwise, for general distributions, use the code in the GLM package.
  #
  if fast_method
    #
    # Create the model matrix, initially still with any missing predictor values.
    # Copy only the complete rows into design matrix X and response vector y,
    # thus they have size = people - too.few.genotypes - incomplete.predictors.
    # Note that the completeness_mask should be applied after the sample_mask.
    #
    modelmatrx = ModelMatrix(model)
    completeness_mask = completecases(model.df)
    cases = sum(completeness_mask)
    predictors = size(modelmatrx.m, 2)
    X = zeros(cases, predictors)
    base_estimate = zeros(predictors)
    for j = 1:predictors
      X[:, j] = modelmatrx.m[completeness_mask, j]
    end
    y = zeros(cases)
    y[1:end] = model.df[completeness_mask, lhs]
    #
    # Estimate parameters under the base model.
    # Then, for linear regression, record the vector of residuals.
    #
    (base_estimate, base_loglikelihood) = regress(X, y, regression_type)
    if regression_type == "linear"
      base_residual = y - (X * base_estimate)
    end
    #
    # Output the results of the base model.
    #
    println(io, " ")
    println(io, "Summary for Base Model with ", fm)
    println(io, "Regression model: ", regression_type)
    println(io, "Link function: ", "canonical")
    names_list = names(model.df)
    model_names = size(names_list,1)
    outcome_index = findin(names_list, [lhs])[1]
    println(io, "Base components' effect estimates: ")
    println(io, "   (Intercept) : ", signif(base_estimate[outcome_index], 6))
    for j = 1:model_names
      if j != outcome_index
        println(io, "   ", names_list[j], " : ", signif(base_estimate[j], 6))
      end
    end
    println(io, "Base model loglikelihood: ", signif(base_loglikelihood, 8))
    println(io, " ")
  else
    #
    # Let the GLM package estimate parameters. Output results.
    # ?? Assumming the glm estimates only use complete cases??
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
  if fast_method; X = [X zeros(cases)]; end
  rhs = parse(side[2] * " + SNP")
  fm = Formula(lhs, rhs)
  #
  # Analyze the SNP predictors one by one.
  #
##  dosage = zeros(people)
  dosage = zeros(sum(sample_mask))
  pvalue = ones(snps)
  alt_estimate = zeros(predictors+1)
  if fast_method && regression_type == "linear"
    alt_residual = zeros(cases)
  end
  skipped_snps = 0
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
    if snpdata.maf[snp] <= maf_threshold || !snp_mask[snp]
      skipped_snps = skipped_snps + 1
      continue
    end
    #
    # Copy the filtered SNP genotypes into a dosage vector,
    # allowing missing genotypes to be simplistically imputed based on MAF.
    #
##    copy!(dosage, view(snpdata.snpmatrix, :, snp); impute = true)
    copy!(dosage, view(snpdata.snpmatrix, sample_mask, snp); impute = true)
    #
    # For the three basic regression types, analyze the alternative model
    # using internal score test code. If the score test p-value
    # is below the specified threshold, carry out a likelihood ratio test.
    #
    if fast_method
##      copy!(X[:, end], dosage[1:cases])
##      @views X[:, end] = copy(dosage[1:cases])
##      copy!(X[:, end], dosage[completeness_mask])
##      copy!(X[:, end], view(dosage, completeness_mask))
##      BLAS.blascopy!(cases, dosage[1:cases], 1, X[:, end], 1)
##      X[:, end] .= view(dosage, completeness_mask)
      X[:, end] = dosage[completeness_mask]
##if snp == 1; println(" X[:, end] = ", typeof(X), " :: ", X[:, end]); end
      alt_estimate = [base_estimate; 0.0]
      score_test = glm_score_test(X, y, alt_estimate, regression_type)
      pvalue[snp] = ccdf(Chisq(1), score_test)
      if pvalue[snp] < lrt_threshold
        (alt_estimate, alt_loglikelihood) = regress(X, y, regression_type)
        lrt = 2.0 * (alt_loglikelihood - base_loglikelihood)
        pvalue[snp] = ccdf(Chisq(1), lrt)
      end
      #
      # Record the vector of residuals under this alternative model:
      # alt_residual = y - (X * alt_estimate)
      #
      if regression_type == "linear"
        A_mul_B!(alt_residual, X, alt_estimate)
        alt_residual .= y .- alt_residual
      end
    #
    # For other distributions analyze the alternative model
    # using the GLM package.
    #
    else
      model.df[:SNP] = dosage[completeness_mask]
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
      println(io, "SNP p-value: ", signif(pvalue[snp], 4))
      println(io, "Minor allele frequency: ", round(snpdata.maf[snp], 4))
      if uppercase(snpdata.chromosome[snp]) == "X"
        hw = xlinked_hardy_weinberg_test(dosage, person.male)
      else
        hw = hardy_weinberg_test(dosage)
      end
      println(io, "Hardy-Weinberg p-value: ", round(hw, 4))
      if fast_method && pvalue[snp] < lrt_threshold
        println(io, "SNP effect estimate: ", signif(alt_estimate[end], 4))
        println(io, "SNP model loglikelihood: ", signif(alt_loglikelihood, 8))
      elseif fast_method
        println(io, "SNP effect estimate: ", signif(alt_estimate[end], 4))
      else
        println(io, "")
        println(io, snp_model)
      end
      #
      # For linear models, output the proportion of the base model's variance
      # that is explained by including this SNP in the model.
      #
      if regression_type == "linear"
        variance_explained = 1 - (norm(alt_residual)^2 / norm(base_residual)^2)
        println(io, "Proportion of base-model variance explained: ",
          round(variance_explained, 4))
      end
    end
##if snp == 1; return; end
  end
  #
  # Output false discovery rates.
  #
  fdr = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
  (number_passing, threshold) = simes_fdr(pvalue, fdr, snps - skipped_snps)
  println(io, " \n \n ")
  println(io, "        P-value   Number of Passing")
  println(io, "FDR    Threshold     Predictors \n")
  for i = 1:length(fdr)
    @printf(io, "%4.2f   %8.5f   %9i\n", fdr[i], threshold[i], number_passing[i])
  end
  println(io, " ")
  #
  # If requested, output a Manhattan Plot in .png format.
  #
  plot_file = keyword["manhattan_plot_file"]
  if plot_file != ""
    println(" \nCreating a Manhattan plot from the GWAS results.\n")
    if !contains(plot_file, ".png"); string(plot_file, ".png"); end
    #
    # Determine if any non-zero basepairs are present in the data set.
    # If so, use basepair position as the x-axis for the Manhattan plot.
    # Otherwise, use SNP number.
    #
    using_basepairs =  (maximum(chr_max_bp) > 0)
    #
    # Find adjusted basepair position for each SNP,
    # which is the position from the start of the data.
    #
    if using_basepairs
      current_chr_number = 1
      running_bp_level = 0
      adj_bp = copy(snpdata.basepairs)
      for snp = 2:snps
        if snpdata.chromosome[snp] != snpdata.chromosome[snp-1]
          running_bp_level = running_bp_level + chr_max_bp[current_chr_number]
          current_chr_number = current_chr_number + 1
        end
        if current_chr_number > 1
          adj_bp[snp] = adj_bp[snp] + running_bp_level
        end
      end
    end
    #
    # Create a new dataframe that will hold the data to plot:
    # (X) adjusted basepair positions or SNP numbers and (Y) -log10(p-values).
    #
    # Create the scatter plot of the -log10(p-values) grouped by chromosome.
    # Set the size, shape, and color of the plotted elements.
    #
    # Specify the x-axis tick marks to be at the center of the chromosomes.
    #
    if using_basepairs
      plot_frame = DataFrame(AdjBasepairs = adj_bp,
        Chromosome = snpdata.chromosome,
        NegativeLogPvalue = -log10.(pvalue))

      plt = scatter(plot_frame[:AdjBasepairs], plot_frame[:NegativeLogPvalue],
        group = plot_frame[:Chromosome],
        markersize = 3, markerstrokewidth = 0, color_palette = :rainbow)

      xticks = by(plot_frame, :Chromosome,
        plot_frame -> mean(plot_frame[:AdjBasepairs]))
    else
      plot_frame = DataFrame(SNPnumber = 1:snps,
        Chromosome = snpdata.chromosome,
        NegativeLogPvalue = -log10.(pvalue))

      plt = scatter(plot_frame[:SNPnumber], plot_frame[:NegativeLogPvalue],
        group = plot_frame[:Chromosome],
        markersize = 3, markerstrokewidth = 0, color_palette = :rainbow)

      xticks = by(plot_frame, :Chromosome,
        plot_frame -> mean(plot_frame[:SNPnumber]))
    end
    #
    # Use x-axis tick marks only in the odd numbered chromosomes.
    # Also, label the x-axis.
    #
    xaxis!(plt, xticks = (sort(xticks[:x1])[1:2:end], 1:2:size(xticks, 1)))
    xaxis!(plt, xlabel = "Chromosome")
    #
    # Add the y-axis information.
    #
    yaxis!(plt, ylabel = "-log10(p-value)")
    #
    # Use a grey grid and remove the legend.
    #
    plot!(plt, gridcolor = :lightgrey, legend = false)
    #
    # Add an overall title.
    #
    plot!(plt, title = "Manhattan Plot", window_title = "Manhattan Plot")
    #
    # Add a dashed horizontal line that indicates the Bonferonni threshold.
    #
    Plots.abline!(plt, 0, -log10(.05 / length(pvalue)), color = :black,
        line = :dash)
##    hline!(plt, -log10(0.05 / length(pvalue)), color = :black, line = :dash)
    #
    # Display the plot and then save the plot to a file.
    #
    ## display(plt)
    savefig(plt, plot_file)
  end
  return execution_error = false
end # function gwas_option

"""
Add principal components into a pedigree frame and regression formula.
The PCs are calculated via pca() from the SnpArrays module.
"""
function add_pcs!(pedigree_frame::DataFrame, snpdata::SnpData,
  keyword::Dict{AbstractString, Any})

  pcs = keyword["pcs"]
  regress_form = keyword["regression_formula"]
  #
  # Perform PCA on the SNP data.
  #
  pcscore, pcloading, pcvariance = pca(snpdata.snpmatrix, pcs)
  #
  # Include the new covariates in the pedigree dataframe.
  # NB: this process is *slow* for a large number of PCs!
  #
  for i = 1:pcs
    pedigree_frame[Symbol("PC$i")] = zscore(pcscore[:,i])
    regress_form = regress_form * " + PC$i"
  end
  #
  # Save the expanded regression formula to the keyword dictionary.
  #
  keyword["regression_formula"] = regress_form
  return nothing
end # function add_pcs!

end # module MendelGWAS
