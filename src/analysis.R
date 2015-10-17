require(rcdk)
require(fingerprint)
require(sqldf)

axis_col = '#FF2015'

header = colnames(read.table('data/AID_488862_datatable_all.csv',sep=',',header=TRUE))
prnp_5utr_primary = read.table('data/AID_488862_datatable_all.csv',sep=',',skip=4,header=FALSE)
colnames(prnp_5utr_primary) = tolower(header)
head(prnp_5utr_primary)

header = colnames(read.table('data/AID_1813_datatable_all.csv',sep=',',header=TRUE))
snca_5utr_primary = read.table('data/AID_1813_datatable_all.csv',sep=',',skip=4,header=FALSE)
colnames(snca_5utr_primary) = tolower(header)
head(snca_5utr_primary)

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

n_hits = sum(-plot_data$prnp_mean < -30 & abs(plot_data$snca_mean < 25), na.rm=TRUE)
n_hits
