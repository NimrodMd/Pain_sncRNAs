library(edgeR) ; library(ggplot2) ; library(matrixStats) ; library(reshape) ; library(ggfortify)  ; library(corrplot) ; library(tidyr)

# setwd('C:/Users/nmisr/Downloads/')
hnni<-read.csv('hanani.csv') ; hnni$group<-paste(hnni$gngl,hnni$trt,sep=' ')
trfs<-read.csv('tRNA_Exclusive_Combined_data.csv') ; rownames(trfs)<-trfs$X ; trfs<-trfs[,-1]
colnames(trfs)<-unlist(lapply(colnames(trfs),function(x) strsplit(as.character(x),'X')[[1]][2]))
cpm_trfs<-as.data.frame(cpm(trfs))
meta<-read.csv('trf_meta.csv')[,2:5] ; colnames(meta)<-c('trf','seq','type','details')
meta$len<-unlist(lapply(meta$trf,function(x) strsplit(as.character(x),'-')[[1]][2]))
meta$trna<-unlist(lapply(meta$details,function(x) substr(strsplit(as.character(x),'_')[[1]][2],1,3)))
meta$codon<-unlist(lapply(meta$details,function(x) substr(strsplit(as.character(x),'_')[[1]][2],4,6)))
meta$gene<-unlist(lapply(meta$details,function(x) strsplit(as.character(x),'_')[[1]][1]))
meta$origin<-factor(meta$gene %in% c('MT','trnalookalike8'),labels=c('Nuclear','Mitochondrial'))

mirs<-read.csv('unnormcounts_compiled.csv') ; rownames(mirs)<-mirs$miRNA ; mirs<-mirs[,-c(1,2)]
colnames(mirs)<-unlist(lapply(colnames(mirs),function(x) strsplit(as.character(x),'X')[[1]][2]))
cpm_mirs<-as.data.frame(cpm(mirs))

sgs<-c('mmu-miR-144-5p','mmu-miR-183-5p','mmu-miR-143-3p','mmu-miR-128-3p','mmu-miR-144-3p','mmu-miR-451a','mmu-miR-96-5p','mmu-miR-29b-3p',
       'mmu-miR-1249-3p','mmu-miR-221-5p','mmu-miR-146a-5p','mmu-miR-155-5p','mmu-miR-122-5p','tRF-17-W96KM8N','tRF-24-7SHRMFWRE2',
       'tRF-30-86V8WPMN1E8Y','tRF-28-86V8WPMN1E0J','tRF-35-9EYK46S9Y81H93','tRF-30-86J8WPMN1E8Y','tRF-29-86V8WPMN1E3J')
dn<-c('mmu-miR-221-5p','tRF-30-86V8WPMN1E8Y','tRF-30-86J8WPMN1E8Y')
up<-c('mmu-miR-183-5p','mmu-miR-143-3p','mmu-miR-128-3p','mmu-miR-96-5p','mmu-miR-29b-3p',
       'mmu-miR-1249-3p','mmu-miR-146a-5p','tRF-24-7SHRMFWRE2')

# Check data all sncRNAs
tmp1<-rbind(trfs,mirs) ; tmp1<-tmp1[,colnames(tmp1) %in% hnni$id] ; tmp1<-tmp1[rowSums(tmp1)>0,] 
tmp1<-tmp1[,colnames(tmp1) %in% hnni$id] ; tmp1<-tmp1[rowSums(tmp1)>0,]
pcaData1<-as.data.frame(t(tmp1)) ; pcaData1$id<-rownames(pcaData1) ; pcaData1<-merge(pcaData1,hnni,by='id')
rownames(pcaData1)<-pcaData1$id ; pca_res1<-prcomp(pcaData1[,2:nrow(tmp1)],scale. = T)
c='trt' # tissue timepoint hits_to_head
autoplot(pca_res1, x=1,y=2 , data = pcaData1, colour = c,label=F,size=3)+facet_wrap(gngl~sex,scales='free')
# ggsave('PCA Hanani no filter.svg',width=4,height=3.5)

# Check data all sncRNAs
tmp1<-rbind(trfs[rownames(trfs) %in% sgs,],mirs[rownames(mirs) %in% sgs,])
tmp1<-tmp1[,colnames(tmp1) %in% hnni$id] ; tmp1<-tmp1[rowSums(tmp1)>0,]
pcaData1<-as.data.frame(t(tmp1)) ; pcaData1$id<-rownames(pcaData1) ; pcaData1<-merge(pcaData1,hnni,by='id')
rownames(pcaData1)<-pcaData1$id ; pca_res1<-prcomp(pcaData1[,2:nrow(tmp1)],scale. = T)
c='trt' # tissue timepoint hits_to_head
autoplot(pca_res1, x=1,y=2 , data = pcaData1, colour = c,label=F,size=3)+facet_wrap(gngl~sex,scales='free')
# ggsave('PCA Hanani sg sncRNAs.svg',width=4,height=3.5)

# tgg7d<-cpm_trfs[,colnames(cpm_trfs) %in% subset(hnni$id,hnni$gngl=='T' & hnni$hours==168)]
# tgg7d<-tgg7d[rowMedians(as.matrix(tgg7d))>10,]
# 
# pvTab<-as.data.frame(t(data.frame(row.names=c('tRF','log2FC','p'))))
# for(t in rownames(tgg7d)){
#   tmp<-data.frame(id=colnames(tgg7d),exp=as.numeric(tgg7d[t,])) ; tmp<-merge(tmp,hnni,by='id')
#   tt<-t.test(tmp$exp~tmp$trt)
#   pvTab<-rbind(pvTab,data.frame('tRF'=t,'log2FC'=log2(tt$estimate[2]/tt$estimate[1]),'p'=tt$p.value))
# } ; pvTab$fdr<-p.adjust(pvTab$p,'fdr')
# 
# ggplot(pvTab,aes(log2FC,-log10(p),col=p<0.05))+theme_classic()+geom_point()
# 
# 
# scg7d<-cpm_trfs[,colnames(cpm_trfs) %in% subset(hnni$id,hnni$gngl=='S' & hnni$hours==168)]
# scg7d<-scg7d[rowMedians(as.matrix(scg7d))>10,]
# 
# pv2Tab<-as.data.frame(t(data.frame(row.names=c('tRF','log2FC','p'))))
# for(t in rownames(scg7d)){
#   tmp<-data.frame(id=colnames(scg7d),exp=as.numeric(scg7d[t,])) ; tmp<-merge(tmp,hnni,by='id')
#   tt<-t.test(tmp$exp~tmp$trt)
#   pv2Tab<-rbind(pv2Tab,data.frame('tRF'=t,'log2FC'=log2(tt$estimate[2]/tt$estimate[1]),'p'=tt$p.value))
# } ; pv2Tab$fdr<-p.adjust(pv2Tab$p,'fdr')
# 
# ggplot(pv2Tab,aes(log2FC,-log10(p),col=p<0.05))+theme_classic()+geom_point()

cmb<-rbind(cpm_mirs,cpm_trfs)
scr<-data.frame(id=colnames(cmb),d=colSums(cmb[rownames(cmb) %in% dn,]),u=colSums(cmb[rownames(cmb) %in% up,]))
scr$Score<-scr$u/scr$d ; scr<-merge(scr,hnni,by='id')
scr$time<-factor(scr$hours,labels=c('4h','1d','1w'))

ggplot(scr,aes(time,Score,col=trt))+theme_classic()+facet_wrap(~gngl,scales='fixed')+geom_hline(yintercept=200)+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitterdodge(jitter.width=0.15,jitter.height=0))
# ggsave('Score changes hanani.svg',width=5,height=2)
TukeyHSD(aov(scr$Score~scr$time*scr$gngl*scr$sex*scr$trt))

pvtab<-as.data.frame(t(data.frame(row.names=c('gng','time','p'))))
for(g in unique(scr$gngl)){for(t in unique(scr$time)){
  tmp<-subset(scr,scr$gngl==g & scr$time==t) ; tt<-t.test(tmp$Score~tmp$trt)
  pvtab<-rbind(pvtab,data.frame('gng'=g,'time'=t,'p'=tt$p.value))
}} ; pvtab$FDR<-p.adjust(pvtab$p,'fdr')

tmp<-subset(scr,scr$gngl==g) ; TukeyHSD(aov(tmp$Score~tmp$time*tmp$trt))

##### edgeR #####
p1<-as.data.frame(t(data.frame(row.names=c('Ganglion','Timepoint','RNA','FoldChange')))) ; CutOff<-0.01
s1<-c()

for(g in c('S','T')){
  for(t in c(4,24,168)){
    n<-paste('sg_',g,t,'_Genes',sep='')
    cts<-trfs[,colnames(trfs) %in% subset(hnni$id,hnni$gngl==g & hnni$hours==t)] 
    cld<-subset(hnni,hnni$gngl==g & hnni$hours==t) ; cts<-cts[,as.character(cld$id)] 
    cts<-cts[,order(cld$id)] ; cld<-cld[order(cld$id),] ; nrow(cld)==sum(cld$id==colnames(cts))
    y <- DGEList(counts=cts,group=cld$trt) ; keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE]
    y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)
    
    cts1<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix
    dsgn <- model.matrix(~sex+trt, data = cld) ; y <- estimateDisp(y, dsgn, robust = T) ; head(dsgn,3)
    fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 3) 
    
    # toRemove<-c() ; for(g in rownames(cts1)){
    #   if(quantile(as.numeric(cts1[g,]),prob=0.85)<mean(as.numeric(cts1[g,]))){toRemove<-append(toRemove,g)}}
    sgGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1)))
    sgGens$log2FC<-log2(2.718^sgGens$logFC)
    sgGens$isSg<-'N.S.' ; if(sum(sgGens$FDR<0.05)>=1){sgGens$isSg<-factor(as.numeric(sgGens$FDR<0.051),labels = c('N.S.','Sgnificant'))}
    sgGens$trf<-rownames(sgGens) ; sgGens<-merge(sgGens,meta,by='trf')
    # assign(n,sgGens) ; write.csv(sgGens,paste0(n,'.csv'))
    s1<-append(s1,subset(sgGens$trf,sgGens$FDR<0.05 & abs(sgGens$log2FC)>1))
    p<-ggplot(sgGens,aes(log2FC,-log10(PValue),shape=FDR<0.05,alpha=FDR<0.05,col=type))+theme_classic()+
      geom_point(size=3)+scale_alpha_manual(values=c(0.5,1))
    # ggsave(plot=p,paste0(n,'.svg'),width=4.5,height=3.5)
    if(length(subset(sgGens$log2FC,sgGens$PValue<CutOff))>0){
      tmp<-data.frame('Ganglion'=g,'Timepoint'=t,'RNA'='tRFs','FoldChange'=subset(sgGens$log2FC,sgGens$PValue<CutOff))
    } else{tmp<-data.frame('Ganglion'=g,'Timepoint'=t,'RNA'='tRFs','FoldChange'=0)} ; p1<-rbind(p1,tmp)
  }
}

for(g in c('S','T')){
  for(t in c(4,24,168)){
    n<-paste('sg_',g,t,'_Mirs',sep='')
    cts<-mirs[,colnames(mirs) %in% subset(hnni$id,hnni$gngl==g & hnni$hours==t)] 
    cld<-subset(hnni,hnni$gngl==g & hnni$hours==t) ; cts<-cts[,as.character(cld$id)] 
    cts<-cts[,order(cld$id)] ; cld<-cld[order(cld$id),] ; nrow(cld)==sum(cld$id==colnames(cts))
    y <- DGEList(counts=cts,group=cld$trt) ; keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE]
    y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)
    
    cts1<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix
    dsgn <- model.matrix(~sex+trt, data = cld) ; y <- estimateDisp(y, dsgn, robust = T) ; head(dsgn,3)
    fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 3) 
    
    # toRemove<-c() ; for(g in rownames(cts1)){
    #   if(quantile(as.numeric(cts1[g,]),prob=0.85)<mean(as.numeric(cts1[g,]))){toRemove<-append(toRemove,g)}}
    sgGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1)))
    sgGens$log2FC<-log2(2.718^sgGens$logFC)
    sgGens$isSg<-'N.S.' ; if(sum(sgGens$FDR<0.05)>=1){sgGens$isSg<-factor(as.numeric(sgGens$FDR<0.051),labels = c('N.S.','Sgnificant'))}
    sgGens$mir<-rownames(sgGens)
    assign(n,sgGens) ; #write.csv(sgGens,paste0(n,'.csv'))
    s1<-append(s1,subset(sgGens$mir,sgGens$FDR<0.05 & abs(sgGens$log2FC)>1))
    p<-ggplot(sgGens,aes(log2FC,-log10(PValue),col=FDR<0.05,alpha=FDR<0.05))+theme_classic()+
      geom_point(size=3)+scale_alpha_manual(values=c(0.5,1))
    # ggsave(plot=p,paste0(n,'.svg'),width=4.5,height=3.5)
    if(length(subset(sgGens$log2FC,sgGens$PValue<CutOff))>0){
      tmp<-data.frame('Ganglion'=g,'Timepoint'=t,'RNA'='miRs','FoldChange'=subset(sgGens$log2FC,sgGens$PValue<CutOff))
    } else{tmp<-data.frame('Ganglion'=g,'Timepoint'=t,'RNA'='miRs','FoldChange'=0)} ; p1<-rbind(p1,tmp)
  }
} ; p1$Ganglion<-factor(p1$Ganglion,levels=c('T','S'),labels=c('Trigeminal','Superior cervical')) ; s1<-unique(s1)

ggplot(tmp,aes(Timepoint,FoldChange,col=RNA,group=RNA))+theme_classic()+facet_wrap(~Ganglion,nrow=2)+
  geom_point(position=position_jitterdodge(jitter.width=0.5,jitter.height=0,dodge.width=1))+
  geom_smooth(method='loess',formula='y~x')+geom_hline(yintercept=0)+xlab('Time (hours)')
# ggsave('Foldchange over time.svg',width=3.5,height=4.5)

tmp<-append(subset(sg_S168_Mirs$mir,sg_S168_Mirs$logFC>0 & sg_S168_Mirs$FDR<0.05),
            subset(sg_S24_Mirs$mir,sg_S24_Mirs$logFC>0 & sg_S24_Mirs$FDR<0.05))
tmp<-append(tmp,subset(sg_S4_Mirs$mir,sg_S4_Mirs$logFC>0 & sg_S4_Mirs$FDR<0.05))
tmp<-append(tmp,subset(sg_T168_Mirs$mir,sg_T168_Mirs$logFC>0 & sg_T168_Mirs$FDR<0.05))
tmp<-append(tmp,subset(sg_T24_Mirs$mir,sg_T24_Mirs$logFC>0 & sg_T24_Mirs$FDR<0.05))
tmp<-append(tmp,subset(sg_T4_Mirs$mir,sg_T4_Mirs$logFC>0 & sg_T4_Mirs$FDR<0.05))
tmp<-data.frame(table(tmp))

tmp1<-rbind(cpm_trfs[rownames(cpm_trfs) %in% sgs,],cpm_mirs[rownames(cpm_mirs) %in% sgs,])
tmp<-cor(t(tmp1)) ; round(tmp,2) ; corrplot(tmp,type="upper",order="hclust",tl.col="black",tl.srt=45)
ng<-c(c("tRF-28-86V8WPMN1E0J","tRF-30-86J8WPMN1E8Y","tRF-30-86V8WPMN1E8Y","mmu-miR-221-5p",
        "mmu-miR-146a-5p","mmu-miR-155-5p","mmu-miR-143-3p","mmu-miR-29b-3p"))

p_sgs<-data.frame(long=colSums(tmp1[rownames(tmp1) %in% ng,]),shrt=colSums(tmp1[!rownames(tmp1) %in% ng,]))
p_sgs$score<-p_sgs$shrt/(1+p_sgs$long) ; p_sgs$id<-rownames(p_sgs) ; p_sgs<-merge(p_sgs,hnni,by='id')

ggplot(p_sgs,aes(trt,score,col=trt))+theme_classic()+scale_x_discrete(guide=guide_axis(angle=45))+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.15,height=0))+
  theme(axis.title.x=element_blank(),legend.title=element_blank())+facet_wrap(~gngl,scales='free_y')
# ggsave('PainScore Hannani.svg',width=2.7,height=3)

tmp<-subset(p_sgs,p_sgs$gngl=='T') ; wilcox.test(tmp$score~tmp$trt)
tmp<-subset(p_sgs,p_sgs$gngl=='S') ; wilcox.test(tmp$score~tmp$trt)
p.adjust(c(0.1782,0.02418),'fdr') # 0.17820 0.04836

TukeyHSD(aov(p_sgs$score~p_sgs$trt*p_sgs$gngl))

sgSNCs<-as.data.frame(t(data.frame(row.names=c('rna','fc','gngl','Time'))))
for(t in unique(hnni$hours)){for(g in unique(hnni$gngl)){
  tmp<-cmb[s1,subset(hnni$id,hnni$hours==t & hnni$gngl==g)]
  t2mp<-data.frame(rna=rownames(tmp),
                   lps=rowSums(tmp[,colnames(tmp) %in% subset(hnni$id,hnni$trt=='LPS')]),
                   cnt=rowSums(tmp[,colnames(tmp) %in% subset(hnni$id,hnni$trt=='CNT')]))
  t2mp$fc<-t2mp$lps/(t2mp$cnt+1) ; t2mp<-t2mp[,c('rna','fc')] ; t2mp$gngl<-g ; t2mp$Time<-t ; sgSNCs<-rbind(sgSNCs,t2mp)}}
sgSNCs$rna<-gsub('mmu-','',sgSNCs$rna) ; sgSNCs$type<-unlist(lapply(sgSNCs$rna, function(x) strsplit(x,'-')[[1]][1]))
sgSNCs$gngl<-factor(sgSNCs$gngl,levels=c('T','S'),labels=c('Trigeminal','Superior cervical'))

ggplot(sgSNCs,aes(Time,log2(fc),col=type))+theme_classic()+facet_grid(gngl~type)+
  geom_line(aes(group=rna))+ geom_smooth(method='loess',formula='y~x',se=T,span=0.5,col='black')+
  geom_hline(yintercept=0,lwd=0.5,col='grey33')+geom_vline(xintercept=c(4,24,168),lwd=0.35,col='grey33',linetype='dashed')
# ggsave('Foldchange over time improved.svg',width=4.5,height=3)

ggplot(sg_T168_Genes,aes(log2FC,-log10(PValue),shape=FDR<0.05,alpha=FDR<0.05,col=type))+theme_classic()+geom_point()+
  scale_alpha_manual(values=c(0.5,1))

g<-'tRF-18-W9RKM80E'
tmp<-data.frame(id=colnames(cpm_trfs),exp=as.numeric(cpm_trfs[g,])) ; tmp<-merge(tmp,hnni,by='id')
ggplot(tmp,aes(as.factor(hours),exp,col=trt))+theme_classic()+facet_wrap(~gngl)+
  geom_point(position=position_jitterdodge(jitter.width=0.2,jitter.height=0))


# S ganglion
cts<-trfs[,colnames(trfs) %in% subset(hnni$id,hnni$gngl=='S' & hnni$hours==168)] 
cld<-subset(hnni,hnni$gngl=='S' & hnni$hours==168) ; cts<-cts[,as.character(cld$id)] 
cts<-cts[,order(cld$id)] ; cld<-cld[order(cld$id),] ; nrow(cld)==sum(cld$id==colnames(cts))
y <- DGEList(counts=cts,group=cld$trt) ; keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE]
y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)

cts1<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix
dsgn <- model.matrix(~sex+trt, data = cld) ; y <- estimateDisp(y, dsgn, robust = T) ; head(dsgn,3)
fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 3) 

toRemove<-c() ; for(g in rownames(cts1)){
  if(quantile(as.numeric(cts1[g,]),prob=0.85)<mean(as.numeric(cts1[g,]))){toRemove<-append(toRemove,g)}}
sg_sGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1)))
sg_sGens$log2FC<-log2(2.718^sg_sGens$logFC)
sg_sGens$isSg<-factor(as.numeric(sg_sGens$FDR<0.051),labels = c('N.S.','Sgnificant'))
sg_sGens$trf<-rownames(sg_sGens) ; sg_sGens<-merge(sg_sGens,meta,by='trf')
# write.csv(sg_sGens,'EdgeR_S_ganglion.csv')

ggplot(sg_sGens,aes(log2FC,-log10(PValue),shape=FDR<0.05,col=type))+theme_classic()+geom_point()

g<-'tRF-18-W9RKM80E'
tmp<-data.frame(id=colnames(cpm_trfs),exp=as.numeric(cpm_trfs[g,])) ; tmp<-merge(tmp,hnni,by='id')
ggplot(tmp,aes(as.factor(hours),exp,col=trt))+theme_classic()+facet_wrap(~gngl)+
  geom_point(position=position_jitterdodge(jitter.width=0.2,jitter.height=0))

##### Kinetics #####
ks<-data.frame(trf=rownames(cpm_trfs),
               'T_4'=as.numeric((rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='LPS' & hnni$hours==4 & hnni$gngl=='T')])+1)/
                                                    (rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='CNT' & hnni$hours==4 & hnni$hours==4 & hnni$gngl=='T')])+1)),
               'T_24'=as.numeric((rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='LPS' & hnni$hours==24 & hnni$gngl=='T')])+1)/
                                                  (rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='CNT' & hnni$hours==24 & hnni$gngl=='T')])+1)),
               'T_168'=as.numeric((rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='LPS' & hnni$hours==168 & hnni$gngl=='T')])+1)/
                                                  (rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='CNT' & hnni$hours==168 & hnni$gngl=='T')])+1)),
               'S_4'=as.numeric((rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='LPS' & hnni$hours==4 & hnni$gngl=='S')])+1)/
                                                    (rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='CNT' & hnni$hours==4 & hnni$gngl=='S')])+1)),
               'S_24'=as.numeric((rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='LPS' & hnni$hours==24 & hnni$gngl=='S')])+1)/
                                                     (rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='CNT' & hnni$hours==24 & hnni$gngl=='S')])+1)),
               'S_168'=as.numeric((rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='LPS' & hnni$hours==168 & hnni$gngl=='S')])+1)/
                                                      (rowSums(cpm_trfs[,subset(hnni$id,hnni$trt=='CNT' & hnni$hours==168 & hnni$gngl=='S')])+1)))
ks<-ks[rowSums(ks[,2:7])>0,] ; ks<-ks[abs(log2(rowMedians(as.matrix(ks[,2:7]))))>0.6,] 
ks<-melt(ks,id.vars='trf',variable_name='hours')
ks$gngl<-unlist(lapply(ks$hours,function(x) strsplit(as.character(x),'_')[[1]][1]))
ks$hours<-as.numeric(unlist(lapply(ks$hours,function(x) strsplit(as.character(x),'_')[[1]][2])))
ks<-merge(ks,meta,by='trf')

ggplot(ks)+theme_classic()+facet_grid(gngl~type)+geom_line(aes(hours,log2(value),col=type,group=trf))+
  geom_smooth(aes(hours,log2(value),col=type),method='loess',formula='y~x',col='black')+ylab('log2(FC)')+xlab('Time (hours)')

