using MendelGWAS
using Base.Test

# write your own tests here
@test 1 == 1

MendelGWAS.GWAS("../docs/gwas 1 Control.txt")
MendelGWAS.GWAS("../docs/gwas 2 Control.txt")
MendelGWAS.GWAS("../docs/gwas 3 Control.txt")
MendelGWAS.GWAS("../docs/gwas 4 Control.txt")
