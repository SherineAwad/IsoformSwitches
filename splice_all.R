library("IsoformSwitchAnalyzeR")

args = commandArgs(trailingOnly=TRUE)


salmonQuant <- importIsoformExpression( parentDir ="/rds/project/rds-O11U8YqSuCk/SM/SLX-21015_2")
summary(salmonQuant)
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = c("E601K","E601K","E601K","RNAb","RNAb","RNAb","RNAb","RNAb","WT","WT","WT")
)


myDesign
#condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
#condition = c("E601K","E601K","E601K","RNAb","RNAb","RNAb","RNAb","RNAb","WT","WT","WT"),     

mySwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "gencode.v40.chr_patch_hapl_scaff.annotation.gtf",
    isoformNtFasta       = "gencode.v40.transcripts.fa",
    showProgress = TRUE,
    ignoreAfterBar = TRUE
)

head(mySwitchList$isoformFeatures,2)
write.csv(mySwitchList$isoformFeatures, "mySwitchList_isoforms.csv")

head(mySwitchList$exons,2)
write.csv(mySwitchList$exons, "mySwitchList_exons.csv")

head(mySwitchList$ntSequence,2)
write.csv(mySwitchList$ntSequence, "mySwitchList_ntSeq.csv")

mySwitchList <- subsetSwitchAnalyzeRlist(
    mySwitchList,
    subset = abs(mySwitchList$isoformFeatures$dIF) > 0.4
)

##Prefilter 

mySwitchList <- preFilter(
    mySwitchList,
    geneExpressionCutoff = 5,
    isoformExpressionCutoff = 3,
    removeSingleIsoformGenes = TRUE
)
summary(mySwitchList)
head(mySwitchList)



mySwitchList <- isoformSwitchTestDEXSeq(
    mySwitchList,
    reduceToSwitchingGenes=TRUE
)

#Summarize switching features
extractSwitchSummary(mySwitchList)



mySwitchList <- extractSequence(
    mySwitchList,
    pathToOutput = '/rds/project/rds-O11U8YqSuCk/SM/SLX-21015_2',
    writeToFile=TRUE

)

extractSwitchSummary(mySwitchList)


write.csv(mySwitchList$isoformFeatures, "part1_isoformsfeatures.csv")
write.csv(mySwitchList$exons, "part1_exons.csv")
write.csv(mySwitchList$ntSequence, "part1_ntseq.csv")


mySwitchList <- analyzeIntronRetention(
    mySwitchList,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)

summary(mySwitchList)


mySwitchList <- analyzeAlternativeSplicing(
   mySwitchList,
    onlySwitchingGenes=TRUE,
    alpha=0.05,
    dIFcutoff = 0.1,
    showProgress=TRUE,
    quiet=TRUE
)


mySwitchList <- analyzeCPC2( mySwitchList, 'cpc2output.txt',
    removeNoncodinORFs = TRUE,
    codingCutoff = 0.725,
    quiet=FALSE
 )


mySwitchList <- analyzeIntronRetention(
    mySwitchList,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)


#consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')
consequencesOfInterest <- c('intron_retention','coding_potential', 'ORF_seq_similarity')

mySwitchList <- analyzeSwitchConsequences(
    mySwitchList,
    consequencesToAnalyze = consequencesOfInterest,
    dIFcutoff = 0.1, # very high cutoff for fast runtimes - you should use the default (0.1)
    showProgress=FALSE
)

pdf(file = 'ConsequenceSummary.pdf', onefile = TRUE, height=6, width = 9)
extractConsequenceSummary(
    mySwitchList,
    consequencesToAnalyze='all',
    plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()

pdf(file = 'ConsequenceEnrichment.pdf', onefile = TRUE, height=6, width = 9)
extractConsequenceEnrichment(
    mySwitchList,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    countGenes = TRUE,
    analysisOppositeConsequence=TRUE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 12),
    minEventsForPlotting = 10,
)
dev.off()

pdf(file = 'ConsequenceEnrichmentComparison.pdf', onefile = TRUE, height=6, width = 9)
extractConsequenceEnrichmentComparison(
    mySwitchList,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    countGenes = TRUE,
    analysisOppositeConsequence=TRUE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 14),
    minEventsForPlotting = 10,
)
dev.off()

pdf(file = 'SplicingEnrichment.pdf', onefile = TRUE, height=6, width = 9)
extractSplicingEnrichmentComparison(
    mySwitchList,
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    countGenes = TRUE,
    plot = TRUE,
    localTheme = theme_bw(base_size = 14),
    minEventsForPlotting = 10,
    returnResult=FALSE
)
dev.off()

#https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#predicting-alternative-splicing


pdf(file = 'SwitchOverlap.pdf', onefile = TRUE, height=6, width = 9)
extractSwitchOverlap(
    mySwitchList,
    filterForConsequences=TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    scaleVennIfPossible=TRUE,
    plotIsoforms = TRUE,
    plotSwitches = TRUE,
    plotGenes = TRUE
)
dev.off()

