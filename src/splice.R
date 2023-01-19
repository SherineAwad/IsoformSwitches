library("IsoformSwitchAnalyzeR")
packageDescription("IsoformSwitchAnalyzeR")$Version


args = commandArgs(trailingOnly=TRUE)
control_condition = args[1]
treat_condition = args[2]
Nc = args[3] 
Nt = args[4] 

workdir = getwd()
salmonQuant <- importIsoformExpression( parentDir =workdir)
summary(salmonQuant)
myDesign <- data.frame(
   sampleID = colnames(salmonQuant$abundance)[-1],
   condition =  c(rep(treat_condition, Nt), rep(control_condition,Nc) ) )

myDesign 

conditions <- data.frame(condition_1 = control_condition, condition_2=treat_condition)


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

out = paste(treat_condition,"mySwitchList_isoforms.csv", sep="_")
head(mySwitchList$isoformFeatures,2)
write.csv(mySwitchList$isoformFeatures, out) 

head(mySwitchList$exons,2)
out = paste(treat_condition, "mySwitchList_exons.csv", sep ="_") 
write.csv(mySwitchList$exons, out)

head(mySwitchList$ntSequence,2)
out = paste(treat_condition, "mySwitchList_ntSeq.csv", sep = "_")
write.csv(mySwitchList$ntSequence, out)



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

out = paste(treat_condition, "part1_isoformsfeatures.csv", sep ="_") 
write.csv(mySwitchList$isoformFeatures, out)
out = paste(treat_condition, "part1_exons.csv", sep ="_")
write.csv(mySwitchList$exons, out)
out = paste(treat_condition, "part1_ntseq.csv", sep= "_")
write.csv(mySwitchList$ntSequence, out)


mySwitchList <- analyzeIntronRetention(
    mySwitchList,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)

summary(mySwitchList) 


out =paste(treat_condition, "cpc2output.txt", sep ="_")
out
mySwitchList <- analyzeCPC2( mySwitchList, out, 
    removeNoncodinORFs = TRUE,
    codingCutoff = 0.725,
    quiet=FALSE
 )



out = paste(treat_condition, 'pfam.txt', sep = "_")
mySwitchList <-  analyzePFAM( mySwitchList, out)

out = paste(treat_condition, "protein_type.txt", sep ="_") 
mySwitchList <- analyzeSignalP( mySwitchList, out)

consequencesOfInterest <- c('intron_retention','coding_potential', 'ORF_seq_similarity')
#consequencesOfInterest <- c('isoform_seq_similarity','isoform_length','intron_retention','coding_potential','NMD_status','domain_length','domains_identified','ORF_length','ORF_seq_similarity', 'signal_peptide_identified')
mySwitchList <- analyzeSwitchConsequences(
    mySwitchList,
    consequencesToAnalyze = consequencesOfInterest,
    dIFcutoff = 0.1, # very high cutoff for fast runtimes - you should use the default (0.1)
    showProgress=FALSE
)

out = paste(treat_condition, "ConsequenceSummary.pdf", sep ="_") 
pdf(file = out, onefile = TRUE, height=6, width = 14, pointsize =7)
consq_summary <- extractConsequenceSummary(
    mySwitchList,
    consequencesToAnalyze= consequencesOfInterest,
    plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE,      # enables analysis of fraction of significant features
    returnResult = TRUE
)
dev.off()

out = paste(treat_condition, "consequence_summary.csv", sep ="_") 
write.csv(consq_summary, out)


mytoplist <- extractTopSwitches(
    mySwitchList,
    filterForConsequences = TRUE,
    n =50,
    sortByQvals = TRUE
)

head(mytoplist)
out = paste(treat_condition, "toplist.csv", sep ="_")
write.csv(mytoplist, out)

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

out = paste(treat_condition, "AlternativeSplicing.csv", sep ="_") 
write.csv(mySwitchList$AlternativeSplicingAnalysis, out) 


mySwitchList <- analyzeIntronRetention(
    mySwitchList,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)
out = paste(treat_condition, "consequence.csv", sep ="_") 
write.csv(mySwitchList$switchConsequence, out)


out = paste(treat_condition, "splicingsummary.pdf", sep="_")
pdf(file = out, onefile = FALSE, height=6, width = 9)
extractSplicingSummary( mySwitchList,  asFractionTotal = FALSE,
    plotGenes=TRUE )
dev.off()

out = paste(treat_condition, "GeneEnrichment.pdf", sep="_") 
pdf(file = out, onefile = FALSE, height=6, width = 9)
extractSplicingEnrichment(
    mySwitchList,
    minEventsForPlotting=5,
    plot = TRUE
)
dev.off()

out = paste(treat_condition, "GenomeWide.pdf", sep ="_") 
pdf(file = out, onefile = FALSE, height=6, width = 9)
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

### Summary
extractSwitchSummary(mySwitchList, filterForConsequences = TRUE)


out = paste(treat_condition, "volcanoplots.pdf", sep ="_") 
pdf(file = out, onefile = TRUE, height=6, width = 9)

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

out = paste(treat_condition, "overview.pdf", sep="_") 
pdf(file = out, onefile = TRUE, height=6, width = 9)
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



out = paste(treat_condition, "ConsequenceEnrichment.pdf", sep ="_")
pdf(file = out, onefile = TRUE, height=6, width = 9)
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




