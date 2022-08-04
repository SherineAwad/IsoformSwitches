library("IsoformSwitchAnalyzeR")
packageDescription("IsoformSwitchAnalyzeR")$Version


args = commandArgs(trailingOnly=TRUE)


salmonQuant <- importIsoformExpression( parentDir ="/rds/project/rds-O11U8YqSuCk/SM/SLX-21015_2/PL")
summary(salmonQuant)
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = c("WT", "WT","WT","RNAb","RNAb","RNAb","RNAb","RNAb") ) 

myDesign 

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

#mySwitchList <- analyzeNovelIsoformORF( mySwitchList, analysisAllIsoformsWithoutORF, genomeObject = NULL)


mySwitchList <- extractSequence(
    mySwitchList, 
    pathToOutput = '/rds/project/rds-O11U8YqSuCk/SM/SLX-21015_2/PL',
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


mySwitchList <- analyzeCPC2( mySwitchList, 'cpc2output.txt',
    removeNoncodinORFs = TRUE,
    codingCutoff = 0.725,
    quiet=FALSE
 )

#mySwitchList <- analyzePFAM(
#    mySwitchList,
#    'pfam_results',
#    showProgress=TRUE,
#    quiet=FALSE
#)

#mySwitchList <- analyzeSignalP(
#    mySwitchList ,
#    'output_protein_type.txt')
#)

#mySwitchList <- analyzeIUPred2A(
#    mySwitchList,
#    'iupred2a_results.txt',
#    showProgress = FALSE
#)

#consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')
mySwitchList <- analyzeAlternativeSplicing(
   mySwitchList,
    onlySwitchingGenes=TRUE,
    alpha=0.05,
    dIFcutoff = 0.1,
    showProgress=TRUE,
    quiet=TRUE
)
summary(mySwitchList$AlternativeSplicingAnalysis)

### overview of number of intron retentions (IR)
write.csv(mySwitchList$AlternativeSplicingAnalysis, "AlternativeSplicing_PL.csv") 


mySwitchList <- analyzeIntronRetention(
    mySwitchList,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)

## Analyse consequences
mySwitchList <- analyzeSwitchConsequences(
    mySwitchList,
    consequencesToAnalyze=c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity'
    ),
    showProgress=TRUE,
    quiet=FALSE
)

write.csv(mySwitchList$switchConsequence, "consequence.csv")

mytoplist <- extractTopSwitches(
    mySwitchList,
    filterForConsequences = TRUE,
    n =50,
    sortByQvals = TRUE
)
head(mytoplist)
write.csv(mytoplist, "toplist_PL.csv")

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
extractSplicingGenomeWide(
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


### Visual analysis
# Indiviudal switches
switchPlotTopSwitches( mySwitchList, n=20, filterForConsequences = TRUE )

pdf(file = 'AMER3.pdf', onefile = TRUE, height=6, width = 9)
switchPlot(mySwitchList, gene='AMER3')
dev.off()


### Summary
extractSwitchSummary(mySwitchList, filterForConsequences = TRUE)

pdf(file = 'VolcanoPlots.pdf', onefile = TRUE, height=6, width = 9)

ggplot(data=mySwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    #facet_wrap( ~ condition_2) +
    facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()

dev.off()


pdf(file = 'OverviewPlot.pdf', onefile = TRUE, height=6, width = 9)
ggplot(data=mySwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    #facet_wrap(~ condition_2) +
    facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
dev.off()


subset(
    extractTopSwitches(
        mySwitchList,
        filterForConsequences = TRUE,
        n=10,
        inEachComparison = TRUE
    )[,c('gene_name','condition_1','condition_2','gene_switch_q_value','Rank')],
    gene_name == 'AMER3'
)

