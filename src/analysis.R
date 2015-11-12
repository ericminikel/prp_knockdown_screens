options(stringsAsFactors=FALSE)

require(rcdk)
require(fingerprint)
require(sqldf)
source('../miscellanea/r_helper.r')

setwd('~/d/sci/src/prp_knockdown_screens')

####
# Preface: reading in data
####

axis_col = '#FF2015'

# pubchem data have the annoying property that the header row with the
# column names is followed by several rows of extra metadata that is a
# different data type (usually strings) than the actual data (usually numeric)
# so I need to read the header and data in separately, with this function
read_pubchem = function(path, skip) {
  header = colnames(read.table(path,sep=',',header=TRUE))
  data = read.table(path,sep=',',skip=skip,header=FALSE)
  colnames(data) = fix_colnames(header)
  return (data)
}

# function to quickly just look at a structure in R interactive mode
plotmol = function(molecule,width=500,height=500) {
  par(mar=c(0,0,0,0)) # set margins to zero since this isn't a real plot
  temp1 = view.image.2d(molecule,width,height) # get Java representation into an image matrix. set number of pixels you want horiz and vertical
  plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='') # create an empty plot
  rasterImage(temp1,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
}

# PRNP 5'UTR primary screen data
prnp_5utr_primary = read_pubchem('data/AID_488862_datatable_all.csv',skip=4)

# SNCA 5'UTR primary screen data
snca_5utr_primary = read_pubchem('data/AID_1813_datatable_all.csv',skip=4)

# APP 5'UTR primary screen data
app_5utr_primary = read_pubchem('data/AID_1285_datatable_all.csv',skip=4)

# APP 5'UTR counterscreen - actives should be removed from hits
app_counterscreen = read_pubchem('data/AID_540343_datatable_all.csv',skip=6)

# PRNP dose-response #1 - inactives should be removed from hits
prnp_dose_1 = read_pubchem('data/AID_504932_datatable_all.csv',skip=6)

# PRNP dose-response #2 - inactives should be removed from hits
prnp_dose_2 = read_pubchem('data/AID_504539_datatable_all.csv',skip=5)

# PRNP cherry-pick re-test - inactives should be removed from hits
prnp_cherrypick = read_pubchem('data/AID_504592_datatable_all.csv',skip=5)

# Luciferase inhibitor counterscreen - actives should be removed from hits
luciferase_inhibitors = read_pubchem('data/AID_588451_datatable_all.csv',skip=6)

# GFP counterscreen - inactives should be removed from hits
gfp_counterscreen = read_pubchem('data/AID_588507_datatable_all.csv',skip=6)


# PrP-FEHTA primary screen data
fehta_primary = read_pubchem('data/AID_720596_datatable_all.csv',skip=5)

# generate list of all pubchem CIDs to go get SMILES for them and later compute fingerprints
cids = unique(c(prnp_5utr_primary$pubchem_cid, snca_5utr_primary$pubchem_cid, fehta_primary$pubchem_cid))
cids = cids[!is.na(cids)]
write.table(data.frame(cid=cids), 'data/cids_to_look_up.tsv', row.names=FALSE, col.names=FALSE, quote=FALSE)
# upload this file to https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi and download a CID\tSMILES .txt.gz file
cid_smiles = read.table('data/4128562969890747401.txt',header=FALSE,sep='\t',quote='',comment.char='')
colnames(cid_smiles) = c('pubchem_cid','smiles')

####
# Chapter 1: Hit calling
####

# Plot PRNP vs. SNCA 5'UTR hit data to call hits specific to PRNP
plot_data = sqldf("
select   p.pubchem_cid,
         (p.replicate_a_activity_score + p.replicate_b_activity_score)/2 prnp_mean,
         (s.replicate_a_activity_score + s.replicate_b_activity_score)/2 snca_mean
from     prnp_5utr_primary p, snca_5utr_primary s
where    p.pubchem_cid = s.pubchem_cid
;")

par(mar=c(6,6,6,6))
plot(-plot_data$prnp_mean, -plot_data$snca_mean, xlim=c(-100,100), ylim=c(-150,150), axes=FALSE, xlab='', ylab='', pch=20, cex=.2)
abline(h=0, col=axis_col, lwd=2)
abline(v=0, col=axis_col, lwd=2)
mtext(side=2, at=0, las=2, line=0, text=expression(italic('PRNP')~'down'), col=axis_col)
mtext(side=4, at=0, las=2, line=0, text=expression(italic('PRNP')~'up'), col=axis_col)
mtext(side=1, at=0, line=0, text=expression(italic('SNCA')~'down'), col=axis_col)
mtext(side=3, at=0, line=0, text=expression(italic('SNCA')~'up'), col=axis_col)
rect(xleft=-100,xright=-30,ybottom=-25,ytop=25,border=axis_col,lwd=4)

n_hits = sum(-plot_data$prnp_mean < -50 & abs(plot_data$snca_mean < 25), na.rm=TRUE)
n_hits


# Table counting hits and how they are progressively eliminated
hit_table = data.frame(step=character(0), count=integer(0))
prnp = prnp_5utr_primary[!duplicated(cid_smiles$smiles[match(prnp$pubchem_cid,cid_smiles$pubchem_cid)]),]
prnp$smiles_avail = prnp$pubchem_cid %in% cid_smiles$cid[!is.na(cid_smiles$smiles)] # all TRUE
hit_table = rbind(hit_table,cbind('Unique compounds tested',nrow(prnp)))
prnp$hit = prnp$pubchem_cid %in% plot_data$pubchem_cid[-plot_data$prnp_mean < -50 & abs(plot_data$snca_mean < 25)]
hit_table = rbind(hit_table,cbind('Hits - PRNP active and SNCA inactive',sum(prnp$hit)))
prnp$hit = prnp$hit & prnp$reproducibility_cosine_transform > .95 
hit_table = rbind(hit_table,cbind('Inter-replicate reproducibility > 0.95',sum(prnp$hit)))
prnp$hit = prnp$hit & !(prnp$pubchem_cid %in% app_counterscreen$pubchem_cid[app_counterscreen$pubchem_activity_outcome=='Active'])
hit_table = rbind(hit_table,cbind('Inactive in APP counterscreen',sum(prnp$hit)))
prnp$hit = prnp$hit & !(prnp$pubchem_cid %in% gfp_counterscreen$pubchem_cid[gfp_counterscreen$pubchem_activity_outcome=='Inactive'])
hit_table = rbind(hit_table,cbind('Active in GFP confirmatory screen',sum(prnp$hit)))
prnp$hit = prnp$hit & !(prnp$pubchem_cid %in% luciferase_inhibitors$pubchem_cid[luciferase_inhibitors$pubchem_activity_outcome=='Active'])
hit_table = rbind(hit_table,cbind('Not luciferase inhibitors',sum(prnp$hit)))
prnp$hit = prnp$hit & !(prnp$pubchem_cid %in% prnp_dose_1$pubchem_cid[prnp_dose_1$pubchem_activity_outcome=='Inactive'])
prnp$hit = prnp$hit & !(prnp$pubchem_cid %in% prnp_dose_2$pubchem_cid[prnp_dose_2$pubchem_activity_outcome=='Inactive'])
hit_table = rbind(hit_table,cbind('Not found inactive in PRNP dose-response',sum(prnp$hit)))
prnp$hit = prnp$hit & !(prnp$pubchem_cid %in% prnp_cherrypick$pubchem_cid[prnp_cherrypick$pubchem_activity_outcome=='Inactive'])
hit_table = rbind(hit_table,cbind('Not found inactive in cherry pick re-testing',sum(prnp$hit)))

hit_table

# final hit rate
sum(prnp$hit) / nrow(prnp)

# generate ECFP6 fingerprints for each compound
fps = list(nrow(cid_smiles))
for (i in 1:nrow(cid_smiles)) {
  mol = parse.smiles(cid_smiles$smiles[i])[[1]]
  fp = get.fingerprint(mol, type='circular', fp.mode = 'bit', depth=6, size=1024, verbose=FALSE)
  fps[i] = fp
}

# to save compute time:
# cluster only the hits, and hits vs. non-hits, but don't cluster the non-hits with each other
cid_smiles$screened = cid_smiles$pubchem_cid %in% prnp$pubchem_cid
cid_smiles$hit = cid_smiles$pubchem_cid %in% prnp$pubchem_cid[prnp$hit]
hit_fps = fps[cid_smiles$hit]
screened_fps = fps[cid_smiles$screened]
screened_fps_test = screened_fps[1:1000]

hit_vs_all_sim = fp.sim.matrix(fplist=hit_fps[1:10], fplist2=screened_fps_test, method='tanimoto')
hit_vs_all_dist = 1 - hit_vs_all_sim
haclust = hclust(as.dist(hit_vs_all_dist))
plot(haclust)

hit_dist = 1 - fp.sim.matrix(hit_fps, method='tanimoto')
clustering = hclust(as.dist(hit_dist))
pdf('~/to_delete/cluster.pdf',height=40,width=60)
plot(clustering)
dev.off()
names(clustering)
nrow(clustering$merge)
head(clustering$height)
head(clustering$order)
head(clustering$labels)
clustering$height
length(clustering$height)
length(unique(clustering$height))
cid_smiles$smiles[cid_smiles$hit][378]
plotmol(parse.smiles(cid_smiles$smiles[cid_smiles$hit][75])[[1]])
cutree(clustering, 4)


cid_smiles$pubchem_cid[cid_smiles$hit]
clusters = cutree(clustering, h=.75)
head(clusters)
length(unique(clusters))
