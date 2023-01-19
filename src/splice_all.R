library("IsoformSwitchAnalyzeR")

args = commandArgs(trailingOnly=TRUE)
control_condition = args[1]
treat1_condition = args[2] 
treat2_condition = args[3] 
Nc = args[4] 
Nt1 = args[5] 
Nt2 =args[6] 
workdir = getwd()
salmonQuant <- importIsoformExpression( parentDir =workdir) 
summary(salmonQuant)
myDesign <- data.frame(
   sampleID = colnames(salmonQuant$abundance)[-1],
   condition = c( rep(treat1_condition, Nt1), rep(treat2_condition, Nt2), rep(control_condition, Nc) ) 
)
 
myDesign 

conditions <- data.frame(condition_1 =c(control_condition, control_condition), condition_2= c(treat1_condition, treat2_condition) )
summary(conditions) 

mySwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "gencode.v40.chr_patch_hapl_scaff.annotation.gtf",
    isoformNtFasta       = "gencode.v40.transcripts.fa",
    showProgress = TRUE,
    ignoreAfterBar = TRUE,
    comparisonsToMake=conditions 
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
    pathToOutput = workdir,
    removeLongAAseq=TRUE, 
    alsoSplitFastaFile=TRUE,
    removeShortAAseq=TRUE, 
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

mySwitchList <-  analyzePFAM( mySwitchList, 'pfam.txt')

mySwitchList <- analyzeSignalP( mySwitchList, 'protein_type.txt')

consequencesOfInterest <- c('intron_retention','coding_potential', 'ORF_seq_similarity') 
#consequencesOfInterest <- c('isoform_seq_similarity','isoform_length','intron_retention','coding_potential','NMD_status','domain_length','domains_identified','ORF_length','ORF_seq_similarity', 'signal_peptide_identified')
mySwitchList <- analyzeSwitchConsequences(
    mySwitchList,
    consequencesToAnalyze = consequencesOfInterest,
    dIFcutoff = 0.1, # very high cutoff for fast runtimes - you should use the default (0.1)
    showProgress=FALSE
)

pdf(file = 'ConsequenceSummary.pdf', onefile = TRUE, height=6, width = 14, pointsize =7)
consq_summary <- extractConsequenceSummary(
    mySwitchList,
    consequencesToAnalyze= consequencesOfInterest,
    plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE,      # enables analysis of fraction of significant features
    returnResult = TRUE
)
dev.off()
write.csv(consq_summary, "consequences_summary.csv")

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
    minEventsForPlotting = 5,
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


pdf(file = 'SwitchOverlap.pdf', onefile = TRUE, height=6, width = 9)
switchoverlap <- extractSwitchOverlap(
    mySwitchList,
    filterForConsequences=TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    scaleVennIfPossible=TRUE,
    plotIsoforms = TRUE,
    plotSwitches = TRUE,
    plotGenes = TRUE
)
head(switchoverlap)
write.csv(switchoverlap, "switchoverlap.csv")
dev.off()

pdf(file = 'splicingsummary.pdf', onefile = FALSE, height=6, width = 9)
extractSplicingSummary( mySwitchList,  asFractionTotal = FALSE,
    plotGenes=TRUE )
dev.off()

pdf(file = 'GeneEnrichment.pdf', onefile = FALSE, height=6, width = 9)
extractSplicingEnrichment(
    mySwitchList,
    minEventsForPlotting=5,
    plot = TRUE
)
dev.off()



pdf(file = 'GenomeWide.pdf', onefile = FALSE, height=6, width = 9)
myIsoforms <- extractSplicingGenomeWide(
    mySwitchList,
    featureToExtract = 'isoformUsage',
    splicingToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
)
dev.off()

mytoplist <- extractTopSwitches(
    mySwitchList,
    filterForConsequences = TRUE,
    n =50,
    sortByQvals = TRUE
)

head(mytoplist)
write.csv(mytoplist, "toplist_all.csv")

switchPlotTopSwitches( mySwitchList, n=20, filterForConsequences = TRUE )

