options(stringsAsFactors=FALSE)

require(rcdk)
require(fingerprint)
require(sqldf)
source('../miscellanea/r_helper.r')

setwd('~/d/sci/src/prp_knockdown_screens')
outdir = '~/d/j/cureffi/media/2015/11'

####
# Chapter 1: reading in data
####

axis_col = '#000000'
hl_col = '#660000'
pt_col = '#FF9912'

# pubchem data have the annoying property that the header row with the
# column names is followed by several rows of extra metadata that is a
# different data type (usually strings) than the actual data (usually numeric)
# so I need to read the header and data in separately, with this function
read_pubchem = function(path, skip) {
  header = colnames(read.table(path,sep=',',header=TRUE))
  data = read.table(path,sep=',',skip=skip,header=FALSE)
  colnames(data) = fix_colnames(header)
  data = data[!is.na(data$pubchem_cid),] # remove rows without a PubChem CID
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

# PrP-FEHTA confirmaton screen data
fehta_conf = read_pubchem('data/AID_743200_datatable_all.csv',skip=5)

# IEC-6 cytotoxicity
iec6_cytotox = read_pubchem('data/AID_1825_datatable_all.csv',skip=5)

# Ba/F3 cytotoxicity
baf3_cytotox = read_pubchem('data/AID_1486_datatable_all.csv',skip=5)

# Silber 2014 (PMID: 24530226) validated hits, correspondence to PubChem CIDs from 3413826967951477504.txt.gz
silber2014 = read_pubchem('data/AID_1072100_datatable_all.csv',skip=5)

silber2014 = read.table('data/silber-2014-validated-hits.tsv',sep='\t',header=TRUE,allowEscapes=FALSE,comment.char='',quote='')
silber2014$pubchem_cid = as.integer(silber2014$pubchem_cid)
silber2014 = silber2014[!is.na(silber2014$pubchem_cid),] # for convenience, remove the 1 w/o a CID

# astemizole
astemizole_smiles = 'COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F'
astemizole_cid = compounds$pubchem_cid[compounds$smiles==astemizole]

# generate list of all pubchem CIDs to go get SMILES for them and later compute fingerprints
cids = unique(c(prnp_5utr_primary$pubchem_cid, snca_5utr_primary$pubchem_cid, fehta_primary$pubchem_cid))
cids = cids[!is.na(cids)]
write.table(data.frame(cid=cids), 'data/cids_to_look_up.tsv', row.names=FALSE, col.names=FALSE, quote=FALSE)
# upload this file to https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi and download a CID\tSMILES .txt.gz file
cid_smiles = read.table('data/4128562969890747401.txt',header=FALSE,sep='\t',quote='',comment.char='')
colnames(cid_smiles) = c('pubchem_cid','smiles')

compounds = unique(rbind(cid_smiles, silber2014[,c('pubchem_cid','smiles')]))

####
# Chapter 2: Pre-process some data
####

# generate ECFP6 fingerprints for each compound. NOTE: this takes >10 minutes
fps = list(nrow(compounds))
for (i in 1:nrow(compounds)) {
  mol = parse.smiles(compounds$smiles[i])[[1]]
  fp = get.fingerprint(mol, type='circular', fp.mode = 'bit', depth=6, size=1024, verbose=FALSE)
  fps[i] = fp
}

# merge the cytotox data
# what fraction of compounds from the cytotox screens are in the CID set?
mean(iec6_cytotox$pubchem_cid %in% compounds$pubchem_cid)
mean(baf3_cytotox$pubchem_cid %in% compounds$pubchem_cid)

min(iec6_cytotox$inhibition[iec6_cytotox$pubchem_activity_outcome=='Active'])
min(baf3_cytotox$inhibition[baf3_cytotox$pubchem_activity_outcome=='Active'])

compounds$iec6_tox = iec6_cytotox$inhibition[match(compounds$pubchem_cid,iec6_cytotox$pubchem_cid)]
compounds$baf3_tox = baf3_cytotox$inhibition[match(compounds$pubchem_cid,baf3_cytotox$pubchem_cid)]
compounds$either_tox = compounds$iec6_tox > 20 | compounds$baf3_tox > 20
mean(compounds$either_tox, na.rm=TRUE)

cor.test(compounds$iec6_tox, compounds$baf3_tox, method='pearson')

compounds$max_tox = pmax(pmin(pmax(scale(compounds$iec6_tox), scale(compounds$baf3_tox)),5),-5)





####
# Chapter 4: 5'UTR hit calling & analysis
####

prnp_5utr_primary$mean = (prnp_5utr_primary$replicate_a_activity_score + prnp_5utr_primary$replicate_b_activity_score)/2
compounds$prnp_z = pmax(pmin(scale(prnp_5utr_primary$mean[match(compounds$pubchem_cid,prnp_5utr_primary$pubchem_cid)]),5),-5)
snca_5utr_primary$mean = (snca_5utr_primary$replicate_a_activity_score + snca_5utr_primary$replicate_b_activity_score)/2
compounds$snca_z = pmax(pmin(scale(snca_5utr_primary$mean[match(compounds$pubchem_cid,snca_5utr_primary$pubchem_cid)]),5),-5)
compounds$app_z = pmax(pmin(scale(app_5utr_primary$ps_inh_2um[match(compounds$pubchem_cid, app_5utr_primary$pubchem_cid)]),5),-5)



png(paste(outdir,'prnp-vs-snca-app.png',sep='/'),width=1200,height=600)
par(mar=c(7,7,2,2), mfrow=c(1,2))
# SNCA
plot(compounds$prnp_z, compounds$snca_z, xlim=c(-5,5.1), ylim=c(-5,5), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='', pch=20, cex=.2, col=pt_col)
axis(side=1, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, cex.axis=1.5)
axis(side=2, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, las=2, cex.axis=1.5)
mtext(side=1, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=1, line=5, text='PRNP activity (Z score)', cex=2)
mtext(side=2, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=2, line=5, text='SNCA activity (Z score)', cex=2)
text(x=3.75,y=-1,pos=1,labels='hits',font=2,cex=1.5,col=hl_col)
rect(xleft=2.5, xright=5, ybottom=-1, ytop=1, border=hl_col, lwd=4)
# APP
plot(compounds$prnp_z, compounds$app_z, xlim=c(-5,5.1), ylim=c(-5,5), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='', pch=20, cex=.2, col=pt_col)
axis(side=1, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, cex.axis=1.5)
axis(side=2, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, las=2, cex.axis=1.5)
mtext(side=1, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=1, line=5, text='PRNP activity (Z score)', cex=2)
mtext(side=2, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=2, line=5, text='APP activity (Z score)', cex=2)
text(x=3.75,y=-1,pos=1,labels='hits',font=2,cex=1.5,col=hl_col)
rect(xleft=2.5, xright=5, ybottom=-1, ytop=1, border=hl_col, lwd=4)
dev.off()

mean(!is.na(compounds$snca_z[!is.na(compounds$prnp_z)]))
mean(!is.na(compounds$app_z[!is.na(compounds$prnp_z)]))
mean(!is.na(compounds$app_z[!is.na(compounds$prnp_z)]) | !is.na(compounds$snca_z[!is.na(compounds$prnp_z)]))

# Table counting hits and how they are progressively eliminated
hit_table = data.frame(step=character(0), count=integer(0))
hit_table = rbind(hit_table,cbind("Unique compounds tested against PRNP 5'UTR",nrow(prnp)))
compounds$prnp_hit = !is.na(compounds$prnp_z) & (!is.na(compounds$snca_z) | !is.na(compounds$app_z))
hit_table = rbind(hit_table,cbind("Also tested against SNCA and/or APP 5'UTR",sum(compounds$prnp_hit)))
compounds$prnp_hit = !is.na(compounds$prnp_z) & compounds$prnp_z > 2.5 & (!is.na(compounds$snca_z) | !is.na(compounds$app_z)) & (is.na(compounds$snca_z) | abs(compounds$snca_z) < 2) & (is.na(compounds$app_z) | abs(compounds$app_z) < 2)
hit_table = rbind(hit_table,cbind('Hits - PRNP active but SNCA and/or APP inactive',sum(compounds$prnp_hit)))
compounds$prnp_hit = compounds$prnp_hit & compounds$pubchem_cid %in% prnp_5utr_primary$pubchem_cid[prnp_5utr_primary$reproducibility_cosine_transform > .95]
hit_table = rbind(hit_table,cbind('Inter-replicate reproducibility > 0.95',sum(compounds$prnp_hit)))
compounds$prnp_hit = compounds$prnp_hit & ((!compounds$either_tox) | is.na(compounds$either_tox))
hit_table = rbind(hit_table,cbind('Non-toxic (if tested) in Ba/F3 & IEC6 cells',sum(compounds$prnp_hit)))
compounds$prnp_hit = compounds$prnp_hit & !(compounds$pubchem_cid %in% gfp_counterscreen$pubchem_cid[gfp_counterscreen$pubchem_activity_outcome=='Inactive'])
compounds$prnp_hit = compounds$prnp_hit & !(compounds$pubchem_cid %in% luciferase_inhibitors$pubchem_cid[luciferase_inhibitors$pubchem_activity_outcome=='Active'])
hit_table = rbind(hit_table,cbind('Not found to be luciferase inhibitors',sum(compounds$prnp_hit)))
compounds$prnp_hit = compounds$prnp_hit & !(compounds$pubchem_cid %in% prnp_dose_1$pubchem_cid[prnp_dose_1$pubchem_activity_outcome=='Inactive'])
compounds$prnp_hit = compounds$prnp_hit & !(compounds$pubchem_cid %in% prnp_dose_2$pubchem_cid[prnp_dose_2$pubchem_activity_outcome=='Inactive'])
hit_table = rbind(hit_table,cbind('Active (if tested) in PRNP dose-response',sum(compounds$prnp_hit)))
compounds$prnp_hit = compounds$prnp_hit & !(compounds$pubchem_cid %in% prnp_cherrypick$pubchem_cid[prnp_cherrypick$pubchem_activity_outcome=='Inactive'])
hit_table = rbind(hit_table,cbind('Active (if tested) in 2.5 uM cherry pick',sum(compounds$prnp_hit)))
compounds$prnp_hit = compounds$prnp_hit & compounds$fehta_z > 2 & !is.na(compounds$fehta_z)
hit_table = rbind(hit_table,cbind('Also active in Scripps FRET screen',sum(compounds$prnp_hit)))
hit_table


mean(!is.na(compounds$fehta_z[!is.na(compounds$prnp_z)]))
mean(!is.na(compounds$fehta_z[compounds$prnp_hit]))

png(paste(outdir,'utr-vs-cell-surface.png',sep='/'),width=600,height=600)
par(mar=c(7,7,2,2))
plot(compounds$prnp_z, compounds$fehta_z, xlim=c(-5,5.2), ylim=c(-5,5.2), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='', pch=20, cex=.2, col=pt_col)
axis(side=1, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, cex.axis=1.5)
axis(side=2, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, las=2, cex.axis=1.5)
mtext(side=1, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=1, line=5, text="5'UTR assay (Z score)", cex=2)
mtext(side=2, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=2, line=5, text='Cell surface PrP assay (Z score)', cex=2)
rect(xleft=2.5, xright=5.1, ybottom=2, ytop=5.1, border=hl_col, lwd=4)
text(x=3.75,y=2,pos=1,labels='shockingly few',font=2,cex=1.5,col=hl_col)
dev.off()

compounds$utr_nominal = compounds$pubchem_cid %in% prnp_5utr_primary$pubchem_cid[prnp_5utr_primary$pubchem_activity_outcome=='Active']
screened_subset = compounds[compounds$pubchem_cid %in% prnp_5utr_primary$pubchem_cid & compounds$pubchem_cid %in% fehta_primary$pubchem_cid,]
contingency_table = table(screened_subset[,c('utr_nominal','fehta_nominal')])
fisher.test(as.matrix(contingency_table))
contingency_table


compounds$both_nominal = compounds$utr_nominal & compounds$fehta_nominal
sum(compounds$both_nominal)
sum(compounds$both_nominal & !is.na(compounds$snca_z) & !is.na(compounds$app_z))
sum(compounds$both_nominal & !is.na(compounds$snca_z) & !is.na(compounds$app_z) & (compounds$snca_z > 2.5 | compounds$app_z > 2.5))

cat(paste(compounds$smiles[compounds$prnp_hit],collapse='.')) # to paste into ChemDraw

hit_fps = fps[compounds$prnp_hit]
hit_dist = 1 - fp.sim.matrix(fplist=hit_fps, method='tanimoto')
hit_clustering = hclust(as.dist(hit_dist))
png(paste(outdir,'prnp_5utr_clusters.png',sep='/'),width=600,height=600)
par(mar=c(4,4,4,4))
plot(hit_clustering, main="PRNP 5'UTR hit clustering", xlab='', ylab='', yaxt='n', col=pt_col)
axis(side=2, at=(0:5)/5, labels=formatC((0:5)/5,digits=1,format='f'), las=2)
mtext(side=2, line=2.5, text='height')
dev.off()

hits_export = compounds[compounds$prnp_hit,c('smiles','pubchem_cid')]
hits_export$cluster = cutree(hit_clustering, h=.6)
hits_export = hits_export[hits_export$cluster %in% which(table(hits_export$cluster) > 1),]
write.table(hits_export[,c('cluster','pubchem_cid','smiles')],'results/prnp_5utr_hits.tsv',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
sqldf("
select   cluster, group_concat(smiles, '.')
from     hits_export
group by 1
order by 1
;")
# here, I just copied the smiles and pasted structures into ChemDraw

plotmol(parse.smiles(compounds$smiles[compounds$prnp_hit][8])[[1]])

screened_fps = fps[!is.na(compounds$prnp_z)]



####
# Chapter 4: Cell surface PrP hit calling & analysis
####


# specific compounds

astemizole_cid = 2247
amcinonide_cid = 443958
#tunicamycin_a_cid = 11104835
#tunicamycin_b_cid = 56927836
#fk506_cid = 445643

fehta_primary[fehta_primary$pubchem_cid==astemizole_cid,]
fehta_primary[fehta_primary$pubchem_cid==amcinonide_cid,]
scale(fehta_primary$inhibition_at_13_8_um)[fehta_primary$pubchem_cid==astemizole_cid]
rank(fehta_primary$inhibition_at_13_8_um)[fehta_primary$pubchem_cid==astemizole_cid]/nrow(fehta_primary)

fehta_silber_overlap = fehta_primary[fehta_primary$pubchem_cid %in% compounds$pubchem_cid[compounds$silber_hit],]
fehta_silber_overlap
silber2014$smiles[silber2014$pubchem_cid %in% fehta_primary$pubchem_cid[fehta_primary$pubchem_activity_score==29]]

# for these EC50s, the CIDs are: 2964934, 2997481, 669092, 684756, 2558390, 1069778, 1473035
fehta_silber_overlap$silber_ec50 = c(0.24, 0.68, 4.13, 0.53, 0.44, 0.69, 1.32)
plot(log10(fehta_silber_overlap$silber_ec50), fehta_silber_overlap$inhibition_at_13_8_um, pch=20)

compounds[compounds$silber_hit & !is.na(compounds$fehta_z),]
sum(compounds$silber_hit & !is.na(compounds$fehta_z))
sum(compounds$silber_hit)

# the hits they called
compounds$fehta_nominal = compounds$pubchem_cid %in% fehta_primary$pubchem_cid[fehta_primary$pubchem_activity_outcome=='Active']
# z score of the inhibition metric for new hit calling
compounds$fehta_z = pmax(pmin(scale(fehta_primary$inhibition_at_13_8_um[match(compounds$pubchem_cid, fehta_primary$pubchem_cid)]),5),-5)
png(paste(outdir,'cell-surface-prp-vs-cytotox.png',sep='/'),width=600,height=600)
par(mar=c(7,7,2,2))
plot(compounds$fehta_z, compounds$max_tox, xlim=c(-5,5.2), ylim=c(-5,5.2), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='', pch=20, cex=.2, col=pt_col)
axis(side=1, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, cex.axis=1.5)
axis(side=2, at=(-2:2)*2.5, col=axis_col, col.axis=axis_col, las=2, cex.axis=1.5)
mtext(side=1, line=3.5, text='inactive \u2190                                                \u2192 active', cex=1.5)
mtext(side=1, line=5, text='cell surface PrP reduction (Z score)', cex=2)
mtext(side=2, line=3.5, text='non-toxic \u2190                                               \u2192 toxic', cex=1.5)
mtext(side=2, line=5, text='cytotoxicity (max Z score)', cex=2)
rect(xleft=3, xright=5.1, ybottom=-2, ytop=2, border=hl_col, lwd=4)
text(x=4,y=-2,pos=1,labels='hits',font=2,cex=1.5,col=hl_col)
dev.off()


hit_table = data.frame(step=character(0), count=integer(0))
hit_table = rbind(hit_table,cbind("Unique compounds tested against cell surface PrP",sum(!is.na(compounds$fehta_z))))
compounds$fehta_hit = !is.na(compounds$fehta_z) & compounds$fehta_z > 2.5 & !is.na(compounds$max_tox) & compounds$max_tox > -2 & compounds$max_tox < 1
hit_table = rbind(hit_table,cbind("Active and non-toxic",sum(compounds$fehta_hit)))
compounds$fehta_hit = compounds$fehta_hit & !(compounds$pubchem_cid %in% fehta_conf$pubchem_cid[fehta_conf$pubchem_activity_outcome=='Inactive'])
hit_table = rbind(hit_table,cbind("Active (if tested) in confirmation screen",sum(compounds$fehta_hit)))
hit_table

compounds$silber_hit = compounds$pubchem_cid %in% silber2014$pubchem_cid

hit_fps = fps[compounds$fehta_hit | compounds$silber_hit]
hit_dist_matrix = 1 - fp.sim.matrix(fplist=hit_fps, method='tanimoto')
hit_clustering = hclust(as.dist(hit_dist_matrix))
plot(hit_clustering, col=pt_col, labels=FALSE, axes=FALSE, main='Cell surface & total PrP hit clustering', xlab='')
axis(side=2, at=(0:5)/5, labels=formatC((0:5)/5,format='f',digits=1),las=2)
abline(h=c(.6,.8), lty=3, lwd=3, col=axis_col)

hits_export = compounds[compounds$fehta_hit | compounds$silber_hit,c('smiles','pubchem_cid','fehta_hit','silber_hit')]

# if you cluster them at height 0.6, there are no clusters with compounds
# shared between screens
hits_export$cluster = cutree(hit_clustering, h=.6)
sqldf('
select   cluster,
      sum(case when fehta_hit then 1 else 0 end) fehta_hits,
      sum(case when silber_hit then 1 else 0 end) silber_hits
from     hits_export h
group by 1
having   fehta_hits > 0 and silber_hits > 0
order by 1
;')
cat(paste(hits_export$smiles[hits_export$cluster==204],collapse='.'))

# if you cluster at height 0.8, there are a few clusters with compounds from both screens
hits_export$cluster = cutree(hit_clustering, h=.8)
sqldf('
select   cluster,
      sum(case when fehta_hit then 1 else 0 end) fehta_hits,
      sum(case when silber_hit then 1 else 0 end) silber_hits
from     hits_export h
group by 1
having   fehta_hits > 0 and silber_hits > 0
order by 1
;')
# cluster 109
paste(hits_export$smiles[hits_export$cluster==109],collapse='.')
# cluster 201
paste(hits_export$smiles[hits_export$cluster==201],collapse='.')


hits_export$smiles[hits_export$cluster==109 & hits_export$silber_hit]

hits_export = hits_export[hits_export$cluster %in% which(table(hits_export$cluster) > 1),]

# write them out to use the Shoichet tool
write.table(hits_export[,c('smiles','pubchem_cid')],'results/hits.smi',sep=' ',row.names=FALSE,col.names=TRUE)


# final hit rate
mean(compounds$prnp_hit[!is.na(compounds$prnp_z)])




# screened_fps_test = screened_fps[1:1000]
# hit_vs_all_sim = fp.sim.matrix(fplist=hit_fps[1:10], fplist2=screened_fps_test, method='tanimoto')
# hit_vs_all_dist = 1 - hit_vs_all_sim
# haclust = hclust(as.dist(hit_vs_all_dist))
# plot(haclust)

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



