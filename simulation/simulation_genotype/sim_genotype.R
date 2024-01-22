
#Modified version of the CARMA simulations script
#change working directory for reading the files
ref.table<-read.csv('Simulation loci.csv')

num.snps=1000
n=100

# see https://zenodo.org/records/7772462 for the original version
#Population   from 1000G, can be changed to EAS, AFR, etc..
pop.names<-sample ( c('EUR','EAS','AFR','SAS','AMR'),size=1)
pop.names
#####1000G vcf file can be found at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#####Here we use the first locus listed in the 'Simulation loci.csv' file to demonstrate the simulation.
###Due to copyright, we can not re-distribute the complete 1000G vcf file. 
##Therefore, we used PLINK to extract the variants from the first locus listed in the 'Simulation loci.csv' file.
##The vcf file here only contains the variants in the first locus listed in the 'Simulation loci.csv' file.

a.vcf<-readVCF(paste0(ref.table$chr[1],'_',
                      as.character(ref.table$region_start [1]),'_',ref.table$region_end[1],'.vcf' ),
               maxNumberOfVariants = num.snps,min_maf=0.05,max_maf = 0.99)
###1000G sample file; the sample file is the complete list of 2504 individuals in 1000G;
sample.file<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)
### match the population in 1000G
pop.index<-which(sample.file$super_pop==pop.names)
readGeneticMapFromFile(paste0("genetic_map_GRCh37_",ref.table$chr[1],".txt.gz"))

startSimulation(a.vcf, totalNumberOfIndividuals =n)

ids = generateUnrelatedIndividuals(n)
genotype = retrieveGenotypes(ids)
image (cor(genotype))
