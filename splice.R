library("IsoformSwitchAnalyzeR")

args = commandArgs(trailingOnly=TRUE)


salmonQuant <- importIsoformExpression( parentDir ="/rds/project/rds-O11U8YqSuCk/SM/SLX-21015_2/")
summary(salmonQuant)
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)
myDesign 

aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "gencode.v40.chr_patch_hapl_scaff.annotation.gtf",
    isoformNtFasta       = "gencode.v40.transcripts.fa",
    showProgress = TRUE,
    ignoreAfterBar = TRUE
)

head(aSwitchList$isoformFeatures,2)
write.csv(aSwitchList$isoformFeatures, "aSwitchList_isoforms.csv") 

head(aSwitchList$exons,2)
write.csv(aSwitchList$exons, "aSwitchList_exons.csv")

head(aSwitchList$ntSequence,2)
write.csv(aSwitchList$ntSequence, "aSwitchList_ntSeq.csv")


##Prefilter 

aSwitchListFilteredStrict <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 5,
    isoformExpressionCutoff = 3,
    removeSingleIsoformGenes = TRUE
)
summary(aSwitchListFilteredStrict)
head(aSwitchListFilteredStrict) 
##Part one 

subsetSwitchList <- subsetSwitchAnalyzeRlist(
    switchAnalyzeRlist = aSwitchListFilteredStrict,
    subset = abs(aSwitchList$isoformFeatures$dIF) > 0.1
)

subsetSwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist   = subsetSwitchList,
    pathToOutput = '/rds/project/rds-O11U8YqSuCk/SM/SLX-21015_2',
    outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
    prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)

head(subsetSwitchList)
write.csv(subsetSwitchList$isoformFeatures, "subset_isoformsfeatures.csv")
write.csv(subsetSwitchList$exons, "subset_exons.csv")
write.csv(subsetSwitchList$ntSequence, "subset_ntsequence.csv") 

extractSwitchSummary( subsetSwitchList )

DEXseqSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = subsetSwitchList,
    reduceToSwitchingGenes=TRUE
)

# Summarize switching features
extractSwitchSummary(DEXseqSwitchListAnalyzed)

write.csv(DEXseqSwitchListAnalyzed$isoformFeatures, "dexseq_isoformsfeatures.csv")
write.csv(DEXseqSwitchListAnalyzed$exons, "dexseq_exons.csv")
write.csv(DEXseqSwitchListAnalyzed$ntSequence, "dexseq_ntseq.csv")
#
#subsetSwitchList <- isoformSwitchAnalysisPart2(
#  switchAnalyzeRlist        = subsetSwitchList, 
#  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
#  removeNoncodinORFs        = TRUE,
#  pathToCPC2resultFile      = system.file("extdata/cpc2_result.txt"         , package = "IsoformSwitchAnalyzeR"),
#  pathToPFAMresultFile      = system.file("extdata/pfam_results.txt"        , package = "IsoformSwitchAnalyzeR"),
#  pathToIUPred2AresultFile  = system.file("extdata/iupred2a_result.txt.gz"  , package = "IsoformSwitchAnalyzeR"),
#  pathToSignalPresultFile   = system.file("extdata/signalP_results.txt"     , package = "IsoformSwitchAnalyzeR"),
#  outputPlots               = FALSE  # keeps the function from outputting the plots from this example
#)

