# IDEAS:
	# - conducting a second temporal signal analysis with the BETS approach of Duchene et al. (2020)?
	# - completing Figure 2 with some boxes/graphical elements summarising known ASFV dispersal history?

library(ape)
library(diagram)
library(fields)
library(igraph)
library(lubridate)
library(maptools)
library(MetBrewer)
library(RColorBrewer)
library(raster)
library(rgeos)
library(seraphim)
library(vegan)

# 1. Preparing the different alignment files for the global data set
# 2. Generating an overall sampling map for all the samples
# 3. Generating an haplotype network for the global data set (Network)
# 4. Testing for a signal of recombination within the global alignment
# 5. Inferring a ML phylogeny for the alignment without MW306192
# 6. Investigating the temporal signal associated with the ML tree
# 7. Extracting the spatio-temporal information embedded in posterior trees
# 8. Estimating dispersal statistics for the ASFV lineages (RRW analysis)
# 9. Mapping the inferred dispersal history of ASFV lineages (RRW analysis)
# 10. Mapping the inferred dispersal history of ASFV lineages (DTA analysis)

writingFiles = FALSE; savingPlots = FALSE

different_countries = c("Belgium","Czechia","Georgia","Germany","Hungary","Italy","Lithuania","Moldova","Poland","Russia","Ukraine")
tab = read.csv("Alignment_140823.csv", head=T, sep=";"); different_colours = colorRampPalette(brewer.pal(11,"Spectral"))(11)

e_Palearctic = extent(-13, 150, 35, 73); e_Europe = extent(-13, 55, 35, 73)
countries1 = crop(gBuffer(shapefile("World_countries_shps/World_countries_shapefile.shp") , byid=T, width=0), e_Palearctic)
borders1 = crop(shapefile("International_borders/Only_international_borders.shp"), e_Palearctic)
coasts1 = crop(shapefile("Coast_lines_borders/Only_coast_lines_borders.shp"), e_Palearctic)
countries2 = crop(countries1, e_Europe); borders2 = crop(borders1, e_Europe); coasts2 = crop(coasts1, e_Europe)

# 1. Preparing the different alignment files for the global data set

txt1 = scan(paste0("Alignment_140823.fasta"), what="", sep="\n", quiet=T)
tab[,1] = gsub(", complete genome","",tab[,1]); IDs = c(); txt2 = c()
for (i in 1:length(txt1))
	{
		if (grepl(">",txt1[i]))
			{
				txt2 = c(txt2, gsub(", complete genome","",txt1[i]))
			}	else	{
				if (grepl(">",txt2[length(txt2)]))
					{
						txt2 = c(txt2, txt1[i])
					}	else		{
						txt2[length(txt2)] = paste0(txt2[length(txt2)], txt1[i])
					}
			}
	}
IDs = txt2[grepl(">",txt2)]; IDs = gsub(">","",IDs)
collectionDates = rep(NA, dim(tab)[1])
for (i in 1:length(collectionDates))
	{
		if (length(unlist(strsplit("\\/",tab[i,"sampling_date"]))) == 3)
			{
				collectionDates[i] = as.character(dmy(gsub("\\/","-",tab[i,"sampling_date"])))
			}	else	{
				collectionDates[i] = as.character(tab[i,"sampling_date"])
			}
	}
for (i in 1:length(IDs))
	{
		ID = gsub("_","",gsub("Ukraine","Ukra",IDs[i])); ID = substr(ID,nchar(ID)-5,nchar(ID))
		if (nchar(ID) == 5) ID = paste0(substr(ID,1,4),"0",substr(ID,5,5))
		IDs[i] = ID # Network only allows 6 characters for the sequence IDs
	}
rowNames = rep(NA, dim(tab)[1])
for (i in 1:dim(tab)[1])
	{
		ID = gsub("_","",gsub("Ukraine","Ukra",tab[i,"seqID"])); ID = substr(ID,nchar(ID)-5,nchar(ID))
		if (nchar(ID) == 5) ID = paste0(substr(ID,1,4),"0",substr(ID,5,5))
		rowNames[i] = ID # Network only allows 6 characters for the sequence IDs
	}
row.names(tab) = rowNames
if (length(unique(IDs)) != length(IDs)) cat("Warning: some shorter sequence IDs are identical")
seqs1 = txt2[!grepl(">",txt2)]; seqs2 = matrix(nrow=length(seqs1), ncol=nchar(seqs1)[1])
for (i in 1:dim(seqs2)[1])
	{
		seqs2[i,] = unlist(strsplit(seqs1[i],""))
	}	# n.b.: if this loop creates continuous error messages, kill the AppleShed process with "kill -9 [PID]"
txt3 = paste0(length(IDs)," ",dim(seqs2)[2])
for (i in 1:length(IDs))
	{
		txt3 = c(txt3, paste0(IDs[i],"   ",paste(seqs2[i,],collapse="")))
	}
if (writingFiles) write(txt3, "ASFV_alignment_1.phy")
seqs3 = seqs2
for (i in 1:dim(seqs3)[1])
	{
		if (seqs3[i,1] == "-") seqs3[i,1] = "£"
		if (seqs3[i,dim(seqs3)[2]] == "-") seqs3[i,dim(seqs3)[2]] = "£"
		indices = which(seqs3[i,] == "-")
		for (j in 1:length(indices))
			{
				if (seqs3[i,indices[j]-1] == "£") seqs3[i,indices[j]] = "£"
				if (seqs3[i,indices[j]+1] == "£") seqs3[i,indices[j]] = "£"
			}
	}
columnsToKeep = c()
for (i in 1:dim(seqs3)[2])
	{
		if (sum(seqs3[,i]!="£") == dim(seqs3)[1])
			{
				columnsToKeep = c(columnsToKeep, i)
			}
	}
seqs4 = seqs3[,columnsToKeep]
txt3 = paste0(length(IDs)," ",dim(seqs4)[2])
for (i in 1:length(IDs))
	{
		txt3 = c(txt3, paste0(IDs[i],"    ",paste(seqs4[i,],collapse="")))
	}
if (writingFiles) write(txt3, "ASFV_alignment_2.phy") # alignment after 5'/3' trimming
columnsToKeep = c()
for (i in 1:dim(seqs2)[2])
	{
		if (sum(seqs2[,i]%in%c("a","c","t","g")) == dim(seqs2)[1])
			{
				columnsToKeep = c(columnsToKeep, i)
			}
	}
seqs3 = seqs2[,columnsToKeep]
txt3 = paste0(length(IDs)," ",dim(seqs3)[2])
for (i in 1:length(IDs))
	{
		txt3 = c(txt3, paste0(IDs[i],"   ",paste(seqs3[i,],collapse="")))
	}
if (writingFiles) write(txt3, "ASFV_alignment_3.phy") # only with the "a/c/t/g" sites
columnsToKeep = c()
for (i in 1:dim(seqs3)[2])
	{
		if (length(unique(seqs3[,i])) != 1)
			{
				columnsToKeep = c(columnsToKeep, i)
			}
	}
seqs4 = seqs3[,columnsToKeep]
txt4 = paste0(length(IDs)," ",dim(seqs4)[2])
for (i in 1:length(IDs))
	{
		txt4 = c(txt4, paste0(IDs[i],"   ",paste(seqs4[i,],collapse="")))
	}
if (writingFiles) write(txt4, "ASFV_alignment_4.phy") # only with the polymorphic sites
metadata = matrix(nrow=length(IDs), ncol=5); metadata[,1] = IDs
for (i in 1:length(IDs))
	{
		index1 = which(row.names(tab)==IDs[i])
		metadata[i,2] = tab[index1,"admin0"]
		index2 = which(different_countries==metadata[i,2])
		metadata[i,3] = different_colours[index2]
	}
decimalDates = rep(NA, length(collectionDates))
for (i in 1:length(decimalDates))
	{
		if (length(unlist(strsplit(collectionDates[i],"\\/"))) == 3)
			{
				decimalDates[i] = decimal_date(dmy(collectionDates[i]))
			}
		if (length(unlist(strsplit(collectionDates[i],"\\/"))) == 2)
			{
				decimalDates[i] = decimal_date(dmy(paste0(collectionDates[i],"-15")))
			}
		if (length(unlist(strsplit(collectionDates[i],"\\/"))) == 1)
			{
				decimalDates[i] = as.numeric(collectionDates[i])+0.5
			}
	}
names(decimalDates) = IDs; colourScale = colorRampPalette(brewer.pal(11,"BrBG"))(101)
colIndices = (((decimalDates-min(decimalDates,na.rm=T))/(max(decimalDates,na.rm=T)-min(decimalDates,na.rm=T)))*100)+1
if (savingPlots)
	{
		rast = raster(as.matrix(c(min(decimalDates,na.rm=T),max(decimalDates,na.rm=T)))); par(lwd=0.2, col="gray30")
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.500,0.512,0.40,0.750),
			 legend.args=list(text="", cex=0.8, line=0.3, col="gray30"), horizontal=F,
		     axis.args=list(cex.axis=0.8, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.5,0)))
	}
cols = colourScale[colIndices]; metadata[,4] = decimalDates; metadata[,5] = cols
populations = cbind(different_countries, different_colours)
if (writingFiles) write.table(metadata, "ASFV_alignment_4.csv", col.names=F, row.names=F, quote=F, sep=";")
if (writingFiles) write.table(metadata[,1:2], "ASFV_Network_pops.csv", col.names=F, row.names=F, quote=F, sep=";")
if (writingFiles) write.table(populations, "ASFV_Network_cols.csv", col.names=F, row.names=F, quote=F, sep=";")

# 2. Generating an overall sampling map for all the samples

all_coordinates = rep(NA, dim(tab)[1])
for (i in 1:dim(tab)[1]) all_coordinates[i] = paste0(tab[i,"latitude"],", ",tab[i,"longitude"])
unique_coordinates = unique(all_coordinates); sampling = as.data.frame(matrix(nrow=length(unique_coordinates), ncol=5))
colnames(sampling) = c("longitude","latitude","country","colour","samples")
for (i in 1:dim(sampling)[1])
	{
		sampling[i,1] = as.numeric(unlist(strsplit(unique_coordinates[i],", "))[2])
		sampling[i,2] = as.numeric(unlist(strsplit(unique_coordinates[i],", "))[1])
		indices = which(all_coordinates==unique_coordinates[i])
		sampling[i,3] = tab[indices[1],"admin0"]; sampling[i,5] = length(indices)
		index = which(different_countries==tab[indices[1],"admin0"])
		sampling[i,4] = different_colours[index]
	}
sampling = sampling[which(!is.na(sampling[,1])),]
pdf("ASFV_sampling_m1.pdf", width=8, height=3.5)
par(oma=c(0,0,0,0), mar=c(1,1,1,0.5), lwd=0.2, col="gray30")
plot(countries1, col="gray90", border=NA, ann=F, axes=F)
plot(borders1, col="white", lwd=0.3, add=T)
plot(coasts1, col="gray70", lwd=0.5, add=T)
for (i in 1:dim(sampling)[1])
	{
		cex = (((sqrt(sampling[i,5]/pi))*2)/((sqrt(1/pi))*2))*0.75
		points(sampling[i,1:2], col=sampling[i,4], cex=cex, pch=16)
		points(sampling[i,1:2], cex=cex, pch=1, lwd=0.5, col="gray30")
	}
rect(-13, 35, 150, 73, lwd=0.2, border="gray30")
dev.off()
pdf("ASFV_sampling_m2.pdf", width=5.3, height=5)
par(oma=c(0,0,0,0), mar=c(1,1,1,0.4), lwd=0.2, col="gray30")
plot(countries2, col="gray90", border=NA, ann=F, axes=F)
plot(borders2, col="white", lwd=0.3, add=T)
plot(coasts2, col="gray70", lwd=0.5, add=T)
for (i in 1:dim(sampling)[1])
	{
		cex = (((sqrt(sampling[i,5]/pi))*2)/((sqrt(1/pi))*2))*0.75
		points(sampling[i,1:2], col=sampling[i,4], cex=cex, pch=16)
		points(sampling[i,1:2], cex=cex, pch=1, lwd=0.5, col="gray30")
	}
rect(-13, 35, 55, 73, lwd=0.2, border="gray30")
dev.off()

# 3. Generating an haplotype network for the global data set (Network)

source("networkGraph_SPADS.r")
network_out = "ASFV_alignment_4.out"
populations = "ASFV_Network_pops.csv"
colour_code = "ASFV_Network_cols.csv"
scaleFactor1 = 4; scaleFactor2 = 1
networkGraph_SPADS(network_out, populations, colour_code, scaleFactor1, scaleFactor2)

# 4. Testing for a signal of recombination within the global alignment

	# --> phi test on the overall aligment: the phi test did find significant evidence for recombination (p = 0.03994)
	# --> analyses conducted with RDP4 highlight a particular sequence associated with a large recombinant fragment: MW306192
	# --> phi test after having excluded that sequence: did not find significant evidence for recombination (p = 0.07103)

	# TO DO: re-running RDP4 on the alignment without MW306192 ??

# 5. Inferring a ML phylogeny for the alignment without MW306192

system("IQTREE_1.6.12_MacOS/iqtree -s ASFV_1_wo_306192.fas -m MFP -mem 10Go -mset GTR -nt 4 -b 100")

tre = readAnnotatedNexus("ASFV_1_wo_306192.tre")
tre$tip.label = paste0("  ",gsub("'","",tre$tip.label))
for (i in 1:length(tre$tip.label))
	{
		country = tab[which(tab[,"seqID"]==gsub("  ","",tre$tip.label[i])),"admin0"]
		tre$tip.label[i] = paste0(tre$tip.label[i]," (",country,")")
	}
dev.new(width=9, height=7); par(mar=c(0.7,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o"); indices = c()
plot(tre, show.tip.label=T, show.node.label=F, edge.width=0.75, cex=0.5, col="gray30", edge.color="gray30")
tre = readAnnotatedNexus("ASFV_1_wo_306192.tre"); tre$tip.label = paste0("  ",gsub("'","",tre$tip.label))
for (i in 1:dim(tre$edge)[1])
	{
		if (!tre$edge[i,2]%in%tre$edge[,1])
			{
				country = tab[which(tab[,"seqID"]==gsub("  ","",tre$tip.label[tre$edge[i,2]])),"admin0"]
				colour = different_colours[which(different_countries==country)]
				nodelabels(node=tre$edge[i,2], pch=16, cex=0.70, col=colour)
				nodelabels(node=tre$edge[i,2], pch=1, cex=0.70, col="gray30", lwd=0.4)
				# tiplabels(tab[index,"seqID"], tip=tre$edge[i,2], col="gray30", cex=0.50, frame="none", adj=-0.2)
			}	else	{
				if (tre$annotations[[i]]$label >= 70)
					{
						nodelabels(node=tre$edge[i,2], pch=16, cex=0.50, col="gray30")
					}
			}
	}
add.scale.bar(x=0.000185, y=1.1, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.6)

# 6. Investigating the temporal signal associated with the ML tree

tre = read.nexus("ASFV_1_wo_306192.tre")
tre$tip.label = gsub("'","",tre$tip.label)
mat = matrix(nrow=length(tre$tip.label), ncol=5)
colnames(mat) = c("trait","date","latitude","longitude","location")
mat[,"trait"] = tre$tip.label
for (i in 1:dim(mat)[1])
	{
		index = which(tab[,"seqID"]==tre$tip.label[i])
		mat[i,"latitude"] = tab[index,"latitude"]
		mat[i,"longitude"] = tab[index,"longitude"]
		mat[i,"location"] = tab[index,"admin0"]
		if (tab[index,"admin0"] == "Russia")
			{
				mat[i,"location"] = tab[index,"admin1"]
			}
		txt = tab[index,"sampling_date"]
		if (length(unlist(strsplit(txt,"\\/"))) == 1)
			{
				date = txt
			}
		if (length(unlist(strsplit(txt,"\\/"))) == 2)
			{
				year = unlist(strsplit(txt,"\\/"))[2]
				month = unlist(strsplit(txt,"\\/"))[1]
				date = paste(year,month,sep="\\/")
			}
		if (length(unlist(strsplit(txt,"\\/"))) == 3)
			{
				year = unlist(strsplit(txt,"\\/"))[3]
				month = unlist(strsplit(txt,"\\/"))[2]
				day = unlist(strsplit(txt,"\\/"))[1]
				date = paste(year,month,day,sep="-")
			}		
		mat[i,"date"] = date
	}
write.table(mat, "ASFV_1_wo_306192.txt", row.names=F, quote=F, sep="\t")
write.tree(drop.tip(tre, c("MG939584","Lith_1")), "ASFV_1_wo_3_tips.tre")

	# For the phylogeographic reconstructions, 4 sequences have been discarded:
		# - MW306192 (recombinant sequence identified with the RDP4 analyses)
		# - MG939584 and Lith_1 because identified as outliers by TempEst
		# - Ukraine_24 because not associated with a known collection date

# 7. Extracting the spatio-temporal information embedded in posterior trees

nberOfTreesToSample = 1000; burnIn = 101 # nberOfTreesToSample = 100; burnIn = 101
log = scan(paste0("ASFV_1_wo_4tips_RRW/ASFV_1_wo_4_tips.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = 6+burnIn; index2 = length(log); interval = round((index2-index1)/nberOfTreesToSample)
indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
write(log[c(5,indices)], paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.log"))
trees = scan(paste0("ASFV_1_wo_4tips_RRW/ASFV_1_wo_4_tips.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]
index2 = index1 + burnIn + 1
indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
selected_trees = c(trees[c(1:index1,indices)],"End;")
write(selected_trees, paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.trees"))

	# To do: getting and annotating the MCC tree with TreeAnnotator, using the 1,000 selected trees as an input
	# Tracer: ucld.mean = 6.68E-6, 95% HPD = [3.67E-6, 1.016E-5]; root age = 2005.59, 95% HPD = [2001.10-2007.94]

source("Tree_data_extraction1.r") # for the MCC tree
source("Tree_data_extraction2.r") # for the posterior trees
mostRecentSamplingDatum = 2022 # ?? (unprecise...)
localTreesDirectory = paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected_ext")
mcc_tre = readAnnotatedNexus(paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.tree"))
mcc_tab = Tree_data_extraction1(mcc_tre, mostRecentSamplingDatum)
	# not used as leading to intermediate internal node coordinates that did not really make sense
mcc_tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_664.csv"), head=T)
	# posterior tree whose topology was selected by TreeAnnotator to annotate the MCC tree
write.csv(mcc_tab, paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.csv"), row.names=F, quote=F)
allTrees = readAnnotatedNexus(paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.trees"))
for (i in 1:length(allTrees))
	{
		csv = Tree_data_extraction2(allTrees[[i]], mostRecentSamplingDatum)
		write.csv(csv, paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

# 8. Estimating dispersal statistics for the ASFV lineages (RRW analysis)

localTreesDirectory = paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected_ext"); extractionFiles = list.files(localTreeDirectory)
nberOfExtractionFiles = length(extractionFiles[which(grepl("TreeExtraction",extractionFiles))])
timeSlices = 100; onlyTipBranches = F; showingPlots = F
outputName = paste0("ASFV_1_wo_4tips_RRW/All_dispersal_statistics/",localTreesDirectory); nberOfCores = 10; slidingWindow = 1
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)

# 9. Mapping the inferred dispersal history of ASFV lineages (RRW analysis)

localTreesDirectory = paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected_ext")
nberOfExtractionFiles = length((list.files(localTreesDirectory)))
mostRecentSamplingDatum = 2022 # ?? (unprecise...)
mcc = read.csv(paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.csv"), head=T)
mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
mcc_tre = readAnnotatedNexus(paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.tree")); mcc_tre$tip.label = gsub("'","",mcc_tre$tip.label)
rootHeight = max(nodeHeights(mcc_tre)); root_time = mostRecentSamplingDatum-rootHeight
minYear = mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
colour_scale = met.brewer(name="Hiroshige", n=111, type="continuous")[1:101]

pdf(paste0("ASFV_1_wo_4_tips_NEW1.pdf"), width=3.5, height=4.05); # dev.new(width=3.5, height=4.05)
par(mar=c(1,0,0,0.5), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o", col="gray30"); plottingRootNode = TRUE
plot(mcc_tre, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, 
	 x.lim=c(minYear-(maxYear-max(nodeHeights(mcc_tre))), max(nodeHeights(mcc_tre))), col="gray30", edge.color="gray30")
mcc_tre_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); rootBarPlotted = FALSE
for (j in 1:dim(mcc_tre$edge)[1])
	{
		endYear = root_time+nodeHeights(mcc_tre)[j,2]
		endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		endYear_colour = colour_scale[endYear_index]
		if ((mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])&(length(mcc_tre$annotations[[j]]$`height_95%_HPD`) > 1))
			{
				x1 = (mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(mcc_tre_obj$yy[mcc_tre$edge[j,2]],2), lwd=3.5, lend=0, col=paste0(endYear_colour,"40"))
			}
		if ((rootBarPlotted == FALSE)&&(!mcc_tre$edge[j,1]%in%mcc_tre$edge[,2]))
			{
				endYear = root_time+nodeHeights(mcc_tre)[j,1]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = colour_scale[endYear_index]
				x1 = (mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(mcc_tre_obj$yy[mcc_tre$edge[j,1]],2), lwd=3.5, lend=0, col=paste0(endYear_colour,"40"))
				rootBarPlotted = TRUE
			}				
	}
for (j in 1:dim(mcc_tre$edge)[1])
	{
		endYear = root_time+nodeHeights(mcc_tre)[j,2]
		endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		endYear_colour = colour_scale[endYear_index]
		if ((mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])&&(mcc_tre$annotations[[j]]$posterior >= 0.95))
			{
				nodelabels(node=mcc_tre$edge[j,2], pch=16, cex=0.50, col=endYear_colour)
				nodelabels(node=mcc_tre$edge[j,2], pch=1, cex=0.50, col="gray30", lwd=0.2)
			}
		if (!mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])
			{
				nodelabels(node=mcc_tre$edge[j,2], pch=16, cex=0.50, col=endYear_colour)
				nodelabels(node=mcc_tre$edge[j,2], pch=1, cex=0.50, col="gray30", lwd=0.2)
			}
		if ((plottingRootNode == TRUE)&&(!mcc_tre$edge[j,1]%in%mcc_tre$edge[,2]))
			{
				endYear = root_time+nodeHeights(mcc_tre)[j,1]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = colour_scale[endYear_index]; plottingRootNode = TRUE		
				nodelabels(node=mcc_tre$edge[j,1], pch=16, cex=0.50, col=endYear_colour)
				nodelabels(node=mcc_tre$edge[j,1], pch=1, cex=0.50, col="gray30", lwd=0.2)
			}
	}
for (j in 1:dim(mcc_tre$edge)[1])
	{
		if (!mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])
			{
				tipLabel1 = mcc_tre$tip.label[mcc_tre$edge[j,2]]
				tipLabel2 = paste0("       ",tab[which(tab[,"seqID"]==tipLabel1),"outbreak"])
				tiplabels(tipLabel2, tip=mcc_tre$edge[j,2], col="gray30", frame="none", bg=NULL, cex=0.4)
			}
	}
selectedDates = seq(1990,2020,5); selectedLabels = selectedDates
selectedDates = c(minYear, selectedDates, maxYear); selectedLabels = c("", selectedLabels, "")
axis(lwd=0.3, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.60, mgp=c(0,-0.10,-0.3), lwd.tick=0.3, 
	 col.lab="gray30", col="gray30", tck=-0.010, side=1)
dev.off()

e_Palearctic = extent(-13, 150, 35, 73); e_Europe = extent(-13, 55, 35, 73)
countries1 = crop(gBuffer(shapefile("World_countries_shps/World_countries_shapefile.shp") , byid=T, width=0), e_Palearctic)
borders1 = crop(shapefile("International_borders/Only_international_borders.shp"), e_Palearctic)
coasts1 = crop(shapefile("Coast_lines_borders/Only_coast_lines_borders.shp"), e_Palearctic)
countries2 = crop(countries1, e_Europe); borders2 = crop(borders1, e_Europe); coasts2 = crop(coasts1, e_Europe)
	
ancestors_BE_clade = c() # monophyletic BE clade only in 822 out of 1000 posterior trees
for (i in 1:nberOfExtractionFiles)
	{
		csv = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
		node1a = csv[which(csv[,"tipLabel"]=="LR536725"),"node1"] # 1° Belgian sequence
		node1b = csv[which(csv[,"tipLabel"]=="MK543947"),"node1"] # 2° Belgian sequence
		if (node1a == node1b)
			{
				node2 = csv[which(csv[,"node2"]==node1a),"node1"]
				coords = csv[which(csv[,"node2"]==node2),c("endLon","endLat")]
				ancestors_BE_clade = rbind(ancestors_BE_clade, coords)
			}
	}
prob = 0.95; H = Hpi(cbind(ancestors_BE_clade[,1],ancestors_BE_clade[,2]))
kde = kde(cbind(ancestors_BE_clade[,1],ancestors_BE_clade[,2]), H=H, compute.cont=T, gridsize=c(1000,1000))
contourLevel = contourLevels(kde, prob=(1-prob)); polygons = list()
contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
for (j in 1:length(contourLines)) polygons[[j]] = Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
ps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps)); # plot(contourPolygons, add=T)
ancestors_BE_clade_95HPD = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))

prob = 0.95; precision = 1; startDatum = minYear
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		polygons_colours[i] = paste0(colour_scale[polygon_index],"20")
	}

cutOffs = c(maxYear); croppingPolygons = FALSE
for (h in 1:length(cutOffs))
	{
		pdf(paste0("ASFV_1_wo_4_tips_NEW2.pdf"), width=8, height=3.5)
		par(oma=c(0,0,0,0), mar=c(1,1,1,0.5), mgp=c(0,0.1,0), lwd=0.2, bty="o")
		plot(countries1, col="gray90", border=NA, ann=F, axes=F)
		plot(borders1, col="white", lwd=0.3, add=T)
		plot(coasts1, col="gray70", lwd=0.5, add=T)
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
		if (h == length(cutOffs))
			{
				plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.70,0.65,0.66),
	 				 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
					 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-1.0, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
			 }
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) <= cutOffs[h])
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(countries1)
						if (croppingPolygons == TRUE) pol = crop(pol, countries1)
						plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
				  	}
			}
		for (i in dim(mcc)[1]:1)
			{
				if ((mcc[i,"startYear"] <= cutOffs[h])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.6)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.6)
					}
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.6)
					}
			}
		rect(-13, 35, 150, 73, lwd=0.2, border="gray30")
		dev.off()
	}
cutOffs = c(maxYear); croppingPolygons = FALSE; plottingAllNodes = FALSE
for (h in 1:length(cutOffs))
	{
		pdf(paste0("ASFV_1_wo_4_tips_NEW3.pdf"), width=5.3, height=5)
		par(oma=c(0,0,0,0), mar=c(1,1,1,0.4), mgp=c(0,0.1,0), lwd=0.2, bty="o")
		plot(countries2, col="gray90", border=NA, ann=F, axes=F)
		plot(borders2, col="white", lwd=0.3, add=T)
		plot(coasts2, col="gray70", lwd=0.5, add=T)
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
		if (plottingAllNodes == TRUE)
			{
				for (i in 1:nberOfExtractionFiles)
					{
						csv = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
						internalNodes = which(csv[,"node2"]%in%csv[,"node1"]); tipNodes = which(!csv[,"node2"]%in%csv[,"node1"])
						points(csv[internalNodes,c("endLon","endLat")], cex=0.5, pch=1, lwd=0.1, col="gray30")
						if (i == nberOfExtractionFiles)
							{
								points(csv[tipNodes,c("endLon","endLat")], cex=0.3, pch=16, col="red")
							}
					}
				internalNodes = which(mcc[,"node2"]%in%mcc[,"node1"])
				points(mcc[internalNodes,c("endLon","endLat")], cex=0.3, pch=16, col="black")
			}
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) <= cutOffs[h])
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(countries2)
						if (croppingPolygons == TRUE) pol = crop(pol, countries2)
						plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
				  	}
			}
		for (i in dim(mcc)[1]:1)
			{
				if ((mcc[i,"startYear"] <= cutOffs[h])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.5)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.5)
					}
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.5)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.5)
					}
			}
		rect(-13, 35, 55, 73, lwd=0.2, border="gray30")
		dev.off()
	}
pdf(paste0("ASFV_1_wo_4_tips_NEW4.pdf"), width=5.3, height=5)
par(oma=c(0,0,0,0), mar=c(1,1,1,0.4), mgp=c(0,0.1,0), lwd=0.2, bty="o")
plot(countries2, col="gray90", border=NA, ann=F, axes=F)
plot(borders2, col="white", lwd=0.3, add=T)
plot(coasts2, col="gray70", lwd=0.5, add=T)
for (i in 1:length(ancestors_BE_clade_95HPD@polygons))
	{
		ancestors_BE_clade_95HPD@polygons[[i]] = checkPolygonsHoles(ancestors_BE_clade_95HPD@polygons[[i]])
	}
pol = ancestors_BE_clade_95HPD; crs(pol) = crs(countries2)
if (croppingPolygons == TRUE) pol = crop(pol, countries2)
plot(pol, axes=F, col=rgb(222,67,39,150,maxColorValue=255), add=T, border=NA)
rect(-13, 35, 55, 73, lwd=0.2, border="gray30")
dev.off()

# 10. Mapping the inferred dispersal history of ASFV lineages (DTA analysis)

mcc_tre = readAnnotatedNexus("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1.tree")
mostRecentSamplingDatum = 2022 # ?? (unprecise...)
mcc = read.csv(paste0("ASFV_1_wo_4tips_RRW/ASFV_all_selected.csv"), head=T)
mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
rootHeight = max(nodeHeights(mcc_tre)); root_time = mostRecentSamplingDatum-rootHeight
minYear = mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]]
maxYear = mostRecentSamplingDatum; countries = c(); nodes_cols = c()
for (i in 1:length(mcc_tre$annotations)) countries = c(countries, mcc_tre$annotations[[i]]$location)
countries = unique(countries); countries = countries[order(countries)]
countries = countries[c(2,3,4,5,6,7,10,11,12,1,8,9,13,14)]
cols = different_colours[c(1:9,10,10,10,10,11)]
root_countries_prob = matrix(0,nrow=1, ncol=length(countries))
countries_prob = matrix(0,nrow=length(mcc_tre$annotations), ncol=length(countries))
colnames(root_countries_prob) = countries; colnames(countries_prob) = countries
index = which(different_countries == mcc_tre$root.annotation$location); root_node_col = cols[index]
for (i in 1:length(mcc_tre$root.annotation$location.set))
	{
		root_countries_prob[1,mcc_tre$root.annotation$location.set[[i]][1]] = mcc_tre$root.annotation$location.set.prob[[i]][1]
	}
for (i in 1:length(mcc_tre$annotations))
	{
		index = which(countries == mcc_tre$annotations[[i]]$location); nodes_cols = c(nodes_cols, cols[index])
		for (j in 1:length(mcc_tre$annotations[[i]]$location.set))
			{
				countries_prob[i,mcc_tre$annotations[[i]]$location.set[[j]][1]] = mcc_tre$annotations[[i]]$location.set.prob[[j]][1]
			}
	}
pieCharts = TRUE; tipNodes = c()
for (i in 1:dim(mcc_tre$edge)[1])
	{
		if (!mcc_tre$edge[i,2]%in%mcc_tre$edge[,1]) tipNodes = c(tipNodes, TRUE)
		if (mcc_tre$edge[i,2]%in%mcc_tre$edge[,1]) tipNodes = c(tipNodes, FALSE)
	}
root_countries_prob[root_countries_prob[]<0.05] = 0; countries_prob[countries_prob[]<0.05] = 0

pdf(paste0("ASFV_1_wo_4_tips_NEW5.pdf"), width=5.0, height=4.05); # dev.new(width=3.5, height=4.05)
par(mar=c(1,0,0,0.5), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o", col="gray30"); plottingRootNode = TRUE
plot(mcc_tre, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, 
	 x.lim=c(minYear-(maxYear-max(nodeHeights(mcc_tre))), max(nodeHeights(mcc_tre))), col="gray30", edge.color="gray30")
mcc_tre_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); rootBarPlotted = FALSE
gray90_transparent = rgb(229, 229, 229, 150, names=NULL, maxColorValue=255)
for (j in 1:dim(mcc_tre$edge)[1])
	{
		# endYear = root_time+nodeHeights(mcc_tre)[j,2]
		# endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		# endYear_colour = colour_scale[endYear_index]
		if ((mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])&(length(mcc_tre$annotations[[j]]$`height_95%_HPD`) > 1))
			{
				x1 = (mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(mcc_tre_obj$yy[mcc_tre$edge[j,2]],2), lwd=3.5, lend=0, col=gray90_transparent) # paste0(endYear_colour,"40")
			}
		if ((rootBarPlotted == FALSE)&&(!mcc_tre$edge[j,1]%in%mcc_tre$edge[,2]))
			{
				# endYear = root_time+nodeHeights(mcc_tre)[j,1]
				# endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				# endYear_colour = colour_scale[endYear_index]
				x1 = (mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(mcc_tre_obj$yy[mcc_tre$edge[j,1]],2), lwd=3.5, lend=0, col=gray90_transparent) # paste0(endYear_colour,"40")
				rootBarPlotted = TRUE
			}
	}
for (i in 1:dim(mcc_tre$edge)[1])
	{
		if (pieCharts == FALSE)
			{
				nodelabels(node=mcc_tre$edge[i,2], pch=16, cex=0.8, col=nodes_cols[i])
				nodelabels(node=mcc_tre$edge[i,2], pch=1, cex=0.8, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
			}	else		{
				if ((sum(countries_prob[i,]==1) > 0)|(sum(countries_prob[i,]==0) == 13))
					{
						nodelabels(node=mcc_tre$edge[i,2], cex=0.20, pie=t(countries_prob[i,]), piecol=cols)
					}	else		{
						nodelabels(node=mcc_tre$edge[i,2], cex=0.60, pie=t(countries_prob[i,]), piecol=cols)
					}	
			}
	}
selectedDates = seq(1990,2020,5); selectedLabels = selectedDates
selectedDates = c(minYear, selectedDates, maxYear); selectedLabels = c("", selectedLabels, "")
axis(lwd=0.3, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.60, mgp=c(0,-0.10,-0.3), lwd.tick=0.3, 
	 col.lab="gray30", col="gray30", tck=-0.010, side=1)
dev.off()

source("Tree_data_extraction3.r")
dir.create(file.path("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1_ext"), showWarnings=F)
trees = readAnnotatedNexus("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1.trees")[102:1101]
for (j in 1:length(trees))
	{
		tab = Tree_data_extraction3(trees[[j]], mostRecentSamplingDatum)	
		write.csv(tab, paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1_ext/TreeExtractions_",j,".csv"), row.names=F, quote=F)
	}
nberOfExtractionFiles = 1000; matrices = list()
for (i in 1:nberOfExtractionFiles)
	{
		mat = matrix(0, nrow=length(countries), ncol=length(countries))
		row.names(mat) = countries; colnames(mat) = countries
		tab = read.csv(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1_ext/TreeExtractions_",i,".csv"), head=T)
		for (j in 1:dim(tab)[1])
			{
				index1 = which(countries==tab[j,"startLoc"])
				index2 = which(countries==tab[j,"endLoc"])
				mat[index1,index2] = mat[index1,index2]+1
			}
		matrices[[i]] = mat
	}
saveRDS(matrices, "ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1.rds")

log1 = scan(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
write(log1[which(!grepl("# ",log1))], paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_4.log"))
log1 = read.table(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_4.log"), header=T, sep="\t"); log1 = log1[102:1001,]
log2 = scan(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_TSW_1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
write(log2[which(!grepl("# ",log2))], paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_TSW_4.log"))
log2 = read.table(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_TSW_4.log"), header=T, sep="\t"); log2 = log2[102:1001,]
BFs1 = matrix(nrow=length(countries), ncol=length(countries))
BFs2 = matrix(nrow=length(countries), ncol=length(countries))
row.names(BFs1) = countries; colnames(BFs1) = countries
row.names(BFs2) = countries; colnames(BFs2) = countries
for (i in 1:length(countries))
	{
		for (j in 1:length(countries))
			{
				if (i != j)
					{
						colName = paste0("location.indicators.",gsub(" ",".",countries[i]),".",gsub(" ",".",countries[j]))
						index1 = which(colnames(log1)==colName); index2 = which(colnames(log2)==colName)
						p = sum(log1[,index1]==1)/dim(log1)[1]
						K = length(countries)
						q = (log(2)+K-1)/(K*(K-1))
						BFs1[i,j] = (p/(1-p))/(q/(1-q))
						p1 = sum(log1[,index1]==1)/dim(log1)[1]
						p2 = sum(log2[,index2]==1)/dim(log2)[1]
						BFs2[i,j] = (p1/(1-p1))/(p2/(1-p2))
					}
			}
	}
write.table(round(BFs1,1), paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_4.csv"), sep=",", quote=F)
write.table(round(BFs2,1), paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_TSW_4.csv"), sep=",", quote=F)

countries1@data$NAME[which(countries1@data$NAME=="Czech Republic")] = "Czechia"
countries1@data$NAME[which(countries1@data$NAME=="Republic of Moldova")] = "Moldova"
selected_countries = subset(countries1, countries1@data$NAME%in%countries)
centroids = coordinates(selected_countries); row.names(centroids) = gsub(" ","",selected_countries@data$NAME)
centroids = centroids[countries[which(countries%in%row.names(centroids))],]
buffer = rbind(cbind(128.13,52.43), cbind(21.41,54.75), cbind(43.42,43.45), cbind(133.21,45.09))
row.names(buffer) = c("Armur","Kaliningrad","Karbardino-Balkaria","Primorsky")
centroids = rbind(centroids, buffer); centroids = centroids[c(1:9,11:14,10),]; matrix_mean_list = list()
matrices = readRDS(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_1.rds"))
matrix_mean = matrix(0, nrow=length(countries), ncol=length(countries)); nberOfExtractionFiles = 1000
for (i in 1:nberOfExtractionFiles) matrix_mean = matrix_mean+matrices[[i]]
matrix_mean = matrix_mean/nberOfExtractionFiles
minVals1 = min(diag(matrix_mean)); maxVals1 = max(diag(matrix_mean))
mat = matrix_mean; diag(mat) = NA; minVals2 = min(mat, na.rm=T); maxVals2 = max(mat, na.rm=T)

pdf("ASFV_1_wo_4_tips_NEW6.pdf", width=5.3, height=5); par(oma=c(0,0,0,0), mar=c(1,1,1,0.4), lwd=0.2, col="gray30")
plot(countries2, col="gray90", border=NA, ann=F, axes=F); adjustedBFs = TRUE
plot(borders2, col="white", lwd=0.3, add=T); plot(coasts2, col="gray70", lwd=0.5, add=T)
mat = matrix_mean; multiplier1 = 200; multiplier2 = 2; multiplier3 = 0.1
if (adjustedBFs == TRUE) BFs = read.csv(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_TSW_4.csv"), head=T)
if (adjustedBFs == FALSE) BFs = read.csv(paste0("ASFV_1_wo_4tips_DTA/ASFV_1_wo_4_tips_DTA_4.csv"), head=T)
points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
for (i in 1:dim(centroids)[1])
	{
		for (j in 1:dim(centroids)[1])
			{
				if ((i!=j)&(mat[i,j]<=0.75))
					{
						LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.02
						# curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=2,
									# lcol="gray70", arr.col="gray70", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")								
					}
			}
	}
for (i in 1:dim(centroids)[1])
	{
		for (j in 1:dim(centroids)[1])
			{
				if ((i!=j)&(mat[i,j]>0.75))
					{
						LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.02
						if ((!is.na(BFs[i,j]))&&(BFs[i,j]<3))
							{
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="gray70", arr.col="gray70", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")								
							}
					}
			}
	}
for (i in 1:dim(centroids)[1])
	{
		for (j in 1:dim(centroids)[1])
			{
				if ((i!=j)&(mat[i,j]>0.75))
					{
						LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.02
						if ((!is.na(BFs[i,j]))&&(BFs[i,j]>3))
							{
								cat("\t",countries[i]," - ",countries[j],": ",round(mat[i,j],2),"\n",sep="")
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
							}
					}
			}
	}
points(cbind(rep(-2,3),rep(68,3)), cex=sqrt((multiplier1*((c(1,10,30)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
rect(-13, 35, 55, 73, lwd=0.2, border="gray30")
dev.off()

