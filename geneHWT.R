# License: GPL-2


library("LDlinkR")


#------------------------------ A function to perform HWT
# Input:
#	AA, Aa and aa: Vector of number of individuals with genotype AA, Aa and aa
#
# Output: 
#	z: Vector of z from HWT
 
HWz <- function(AA,Aa,aa){
    N <- AA + Aa + aa
    p1 <- ((2*AA)+Aa)/(2*N)
    p2 <-  1-p1
    d <- (Aa - 2*N*p1*p2)

    if(min(p1^2*N,p2^2*N) <= 5){
	z <- (d-0.5*sign(d))/(2*p1*p2*sqrt(N))
	} else {
    	z <- (d)/(2*p1*p2*sqrt(N))
    }
    return(z)}

#------------------------------ A function to perform HWT
# Arguments:
# - dsin: A dataframe with the following variables.
#   - gene : Gene name
#   - rsID: The rsID of the genetic variant.
#   - AA: The number of individuals with genotype AA.
#   - Aa: The number of individuals with genotype Aa.
#   - aa: The number of individuals with genotype aa.
# - token: A personal access token for LDlinkR package. See LDlinkR package documentation for details.
# - pop: A 1000 Genomes Project population (e.g., YRI or CEU). See LDlinkR package documentation for details.
# - genome_build: Choose between one of three options: grch37 for genome build GRCh37 (hg19), grch38 for GRCh38 (hg38), or grch38_high_coverage for GRCh38 High Coverage (hg38) 1000 Genome Project data sets. Default is GRCh37 (hg19). See LDlinkR package documentation for details.
# - n_r2: Number of attempts to retrieve r2 using LDmatrix. If the number of attempts exceeds n_r2, the function will abort.
#
# Outputs: A list containing the following variables:
# - gene: Gene name.
# - n_vari: Number of variants used in the analysis (number of rsIDs).
# - HWz_com: Z-value for geneHWT.
# - HWp_com: P-value for geneHWT.

geneHWT <- function(dsin, token, pop, genome_build="grch37",n_r2=10){
	y=NA
	dsin$HWz <- HWz(dsin$AA, dsin$Aa, dsin$aa)
	if(length(dsin[,1])==1) {#----------------------The case of one rsID
		ldmat <- matrix(1, ncol=2)
	} else if(length(dsin[,1])>1000) {#----------------------The case of a number of rsIDs greater than 1000.
		return(list(gene=dsin$gene[1], n_vari=length(dsin$rsID),HWz_com=NA,HWp_com=NA))
	} else {                 
		for(i in 1:(n_r2=+1)){#----------------------Retry up to n_r2 times when LDmatrix malfunctions.
			y<-try(ldmat <- LDmatrix(dsin$rsID, pop=pop, r2d="r2", token=token, genome_build=genome_build))
			if (class(y) == "try-error") {
				if (i==(n_r2+1)) return(list(gene=dsin$gene[1], n_vari=length(dsin$rsID),HWz_com=NA,HWp_com=NA))
				Sys.sleep(3)
				next
			}
			break # success
		}
	}
	if(length(dsin[,1])>1) {
		# remove NA cols and rows (when LDmatrix returns NA)
		if(sum(is.na(ldmat))>0) ldmat <- ldmat[-which(is.na(ldmat[,2])),-which(is.na(ldmat[1,]))] 
		if(length(ldmat)==0) return(list(gene=dsin$gene[1], n_vari=length(dsin$rsID),HWz_com=NA,HWp_com=NA)) 
		dsin <- dsin[dsin$rsID %in% ldmat[,1],] # limit rsIDs to ones having r2 values .
	}
	HWz_com <- sum(dsin$HWz)/sqrt(sum(ldmat[,-1]))
	HWp_com <- 1-pchisq(HWz_com^2,1)
	return(list(gene=dsin$gene[1], n_vari=length(dsin$rsID),HWz_com=HWz_com,HWp_com=HWp_com))
}

#------------------------------ A function to perform HWT
# Arguments:
# - dsin: A dataframe with the following variables.
#   - gene : Gene name
#   - rsID: The rsID of the genetic variant.
#   - AA: The number of individuals with genotype AA.
#   - Aa: The number of individuals with genotype Aa.
#   - aa: The number of individuals with genotype aa.
#
# - ldmat: Matrix of linkage disequilibrium coefficients (r2) for rsIDs. It must be ordered according to the input rsIDs.
#
# Outputs: A list containing the following variables:
# - gene: Gene name.
# - n_vari: Number of variants used in the analysis (number of rsIDs).
# - HWz_com: Z-value for geneHWT.
# - HWp_com: P-value for geneHWT.

#dsin: required vars, rsID,  AA, Aa,  aa, ldmat:matrix of r2
geneHWT2 <- function(dsin,ldmat){
	dsin$HWz <- HWz(dsin$AA, dsin$Aa, dsin$aa)
	HWz_com <- sum(dsin$HWz)/sqrt(sum(ldmat))
	HWp_com <- 1-pchisq(HWz_com^2,1)
	return(list(gene=dsin$gene[1], n_vari=length(dsin$rsID),HWz_com=HWz_com,HWp_com=HWp_com))
}



#------------------------------ Test data (test_dt)

gene <- c(rep("AK8", 20), rep("NPNT", 6), rep("PARP1", 12), rep("QRFPR", 9), rep("ZNF680", 2))

rsID <- c("rs1000012", "rs10117989", "rs10901214", "rs11243900", "rs12000803", "rs1410969", 
	 "rs1609597", "rs17149822", "rs215156", "rs215165", "rs215192", "rs2771994", 
	"rs2772007", "rs415831", "rs4962211", "rs4962215", "rs4962218", "rs929293", 
	"rs932887", "rs954390", "rs10000111", "rs11945032", "rs13125117", "rs17036413",
	"rs17036425", "rs6811135", "rs1000033", "rs1805415", "rs2271347", "rs2377312", 
	"rs3219058", "rs3219070", "rs3219094", "rs3219126", "rs3219128", "rs3754370",
	"rs747657", "rs747658", "rs10857080", "rs11732033", "rs1513708", "rs17372524",
	"rs2013332", "rs6813974", "rs691182", "rs691712", "rs7696277", "rs10243954", 
	"rs17540510")

AA <- c(50, 93, 90, 77, 3, 163, 171, 171, 6, 41, 7, 10, 4, 40, 77, 7, 85, 92, 6, 174,
        142, 99, 62, 129, 15, 97, 101, 37, 165, 51, 101, 156, 166, 166, 166, 101, 120,
	101, 36, 139, 163, 147, 78, 16, 11, 12, 158, 41, 76)

Aa <- c(95, 85, 79, 77, 41, 22, 22, 22, 39, 96, 46, 62, 48, 97, 77, 45, 94, 85, 43, 20,
        46, 83, 95, 60, 15, 90, 78, 94, 28, 95, 78, 38, 28, 28, 29, 78, 68, 78, 85, 44,
	22, 34, 85, 58, 25, 46, 34, 69, 77)

aa <- c(43, 17, 26, 41, 151, 2, 2, 2, 150, 55, 141, 123, 143, 58, 41, 143, 16, 18, 146, 1,
        5, 13, 38, 6, 111, 7, 16, 64, 2, 49, 16, 1, 1, 1, 0, 16, 7, 16, 74, 12, 5, 14, 32, 
	121, 159, 137, 0, 85, 39)

test_dt <- data.frame(gene = gene, rsID = rsID, AA = AA, Aa = Aa, aa = aa)

