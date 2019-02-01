### Overview
Mendel GWAS is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. This analysis option performs a standard Genome-Wide Association Study (GWAS) to assess common genetic variants in unrelated individuals.

### Appropriate Problems and Data Sets
Mendel GWAS analysis input data is unrelated individuals genotyped at a large number of autosomal or X-linked SNPs. Mendel GWAS uses the compressed SNP data files. 

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [MendelSearch](https://openmendel.github.io/MendelSearch.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelGWAS:

   pkg> add https://github.com/OpenMendel/MendelGWAS.jl.git

This package supports Julia v1.0+

### Input Files
The MendelGWAS analysis package uses the following input files. Example input files can be found in the [data](https://github.com/OpenMendel/MendelGWAS.jl/tree/master/data) subfolder of the MendelGWAS project.

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Pedigree File](https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, parental information, and sex.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, and allele names.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set and must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file, you must also have a SNP definition file.

<a id="control-file"></a>
### Control file
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run GWAS:

	#
	# Input and Output files.
	#
	plink_input_basename = gwas 1 data
	output_file = gwas 1 Output.txt
	manhattan_plot_file = gwas 1 Manhattan Plot Output.png
	#
	# Analysis parameters for GWAS option.
	#
	regression = linear
	regression_formula = Trait ~ Sex

In the example above, there are five keywords. The keyword *plink_input_basename* tells MendelGWAS that the input data files will comprise three [PLINK format](http://zzz.bwh.harvard.edu/plink) data files: *gwas 1 data.fam*, *gwas 1 data.bim*, and *gwas 1 data.bed*. The next two keywords specify the output files: *gwas 1 Output.txt* the results file - and *gwas 1 Manhattan Plot Output.png* - a plot of the results. The last two keywords specify analysis parameters. The text after the "=" are the keyword values.

### Keywords<a id="keywords-table"></a>
This is a list of OpenMendel keywords relevant to GWAS. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword *values* may be case sensitive.)

Keyword          |   Default Value    | Allowed Values   |  Short Description       
---------------- |  ----------------  | ---------------- |  -----------------
affected_designator | affected | | Under logistic regression, the trait label assigned to cases
distribution |  | Binomial(), Gamma(), Normal(), Poisson(), etc. | Name of distribution function
link |  | LogitLink(), IdentityLink(), LogLink(), etc. |  Name of link function
lrt_threshold | 5e-8 | | Threshold for the score test p-value below which a likelihood ratio test is performed
maf_threshold | 0.01 | | Threshold for the minor allele frequency below which SNPs are not analyzed
manhattan_plot_file | | | Name of file to hold Manhattan plot
regression | |   Linear, Logistic, or Poisson  | Type of standard regression to perform; if left blank, the keywords distribution and link must be assigned values
regression_formula | | | Defines the regression model to analyze under the null hypothesis. See below for a description of the syntax to use

The value for the keyword regression_formula takes the following form: the trait variable is separated from the predictors by "~". For example, `regression_formula = Case_Control ~ Sex + BMI` means the trait labeled Case_Control will be modeled as a grand mean (intercept term) along with the covariates of Sex and BMI. Case_Control and BMI must be field names used in the pedigree file. (If you use a PLINK style FAM file, which does not use a header row, then the trait value should be referred to as "Trait".) Interactions between predictors are represented in the regression_formula using an "&", for example, `Case_Control ~ Sex + BMI + Sex&BMI` would add the interaction term to the previous main effects. (A product term, for example, Sex*BMI, is a short cut for main effects plus their interaction.) If the right hand side (after the "~") is blank, then the model uses only the grand mean. Do not include a term for the SNPs as that is added automatically for each alternative model.

Under logistic regression, the cases are those individuals whose value at the trait field is the same as the label assigned to the keyword affected_designator. The trait field is the field listed on the left hand side (before the "~") in the regression_formula. The controls are those individuals with non-missing values at the trait that are not cases. Of course individuals with missing values at the trait are neither cases nor controls. 

### Data Files
GWAS requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file) and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data is provided in a [SNP data file](https://openmendel.github.io/MendelBase.jl/#snp-data-file), with a [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file) describing the SNPs. OpenMendel will also accept [PLINK format](http://zzz.bwh.harvard.edu/plink) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the GWAS [data](https://github.com/OpenMendel/MendelGWAS.jl/tree/master/data) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelGWAS

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in your Control file, for example, Control_file.txt, use the command:

     julia> GWAS("Control_file.txt")

*Note: The package is called* MendelGWAS *but the analysis function is called simply* GWAS.

<!--- ### Interpreting the results --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
