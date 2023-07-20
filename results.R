library(dplyr)
library(data.table)

### Exome-wide burden testing
dir<-commandArgs(trailingOnly=FALSE)[1]
setwd(dir)

system('head -n1 parts/chr1_0_19.txt > all.txt; for i in `ls parts/`; do tail -n+2 parts/${i} >> all.txt; done')

x<-fread('all.txt', stringsAsFactors=F, header=T, data.table=F)
y <- rbind((x[which(x$Variant_class!='synonymous/synonymous'),] %>% group_by(Variant_class, Population, LD_thinning_r2) %>% summarise(O=sum(Observed_biallelic_genotypes), E=sum(Expected_freq_biallelic_genotypes*N_probands), `O-E`=sum(Observed_biallelic_genotypes - (Expected_freq_biallelic_genotypes*N_probands)), `O-E/n`=sum(Observed_biallelic_genotypes - (Expected_freq_biallelic_genotypes*N_probands))/max(N_probands), p_poisson=poisson.test(x=sum(Observed_biallelic_genotypes), T=1, r=sum(Expected_freq_biallelic_genotypes*N_probands),alternative='greater')$p.value, lb95ci=poisson.test(x=sum(Observed_biallelic_genotypes), T=1, r=sum(Expected_freq_biallelic_genotypes*N_probands))$conf.int[1], ub95ci=poisson.test(x=sum(Observed_biallelic_genotypes), T=1, r=sum(Expected_freq_biallelic_genotypes*N_probands))$conf.int[2], N_trios=max(N_probands))),  (x[which(x$Variant_class=='synonymous/synonymous'),] %>% group_by(Variant_class, Population, LD_thinning_r2) %>% summarise(O=sum(Observed_biallelic_genotypes), E=sum(Expected_freq_biallelic_genotypes*N_probands), `O-E`=sum(Observed_biallelic_genotypes - (Expected_freq_biallelic_genotypes*N_probands)), `O-E/n`=sum(Observed_biallelic_genotypes - (Expected_freq_biallelic_genotypes*N_probands))/max(N_probands), p_poisson=poisson.test(x=sum(Observed_biallelic_genotypes), T=1, r=sum(Expected_freq_biallelic_genotypes*N_probands))$p.value, lb95ci=poisson.test(x=sum(Observed_biallelic_genotypes), T=1, r=sum(Expected_freq_biallelic_genotypes*N_probands))$conf.int[1], ub95ci=poisson.test(x=sum(Observed_biallelic_genotypes), T=1, r=sum(Expected_freq_biallelic_genotypes*N_probands))$conf.int[2], N_trios=max(N_probands))))
x[which(x$Variant_class=='lof/missense'),c('Observed_biallelic_genotypes', 'Expected_freq_biallelic_genotypes', 'corrected_Expected_freq_biallelic_genotypes')]<-x[which(x$Variant_class=='lof/missense'),c('Observed_biallelic_genotypes', 'Expected_freq_biallelic_genotypes', 'corrected_Expected_freq_biallelic_genotypes')]+x[which(x$Variant_class=='lof/lof'),c('Observed_biallelic_genotypes', 'Expected_freq_biallelic_genotypes', 'corrected_Expected_freq_biallelic_genotypes')]

tmp<-x[which(x$Variant_class=='lof/missense'|x$Variant_class=='missense/missense'),]%>%group_by(Gene, Population, LD_thinning_r2)%>%summarize(Variant_class="all", Observed_biallelic_genotypes=sum(Observed_biallelic_genotypes), N_probands=max(N_probands), N_parents=max(N_parents), Expected_freq_biallelic_genotypes=sum(Expected_freq_biallelic_genotypes), corrected_Expected_freq_biallelic_genotypes=sum(corrected_Expected_freq_biallelic_genotypes))
tmp<-tmp[, names(x)]
x<-bind_rows(x, tmp)
fwrite(x, 'all.txt', row.names=F, quote=F, sep=' ')
write.table(y, 'results.txt', row.names=F, quote=F)

#### Per-gene testing

x$Expected<-x$Expected_freq_biallelic_genotypes*x$N_probands
x[which(x$Observed_biallelic_genotypes>0 & x$Expected==0),'cExpected']<-x[which(x$Observed_biallelic_genotypes>0 & x$Expected==0),'corrected_Expected_freq_biallelic_genotypes']*x[which(x$Observed_biallelic_genotypes>0 & x$Expected==0),'N_probands']
final<- x%>%group_by(Gene, Variant_class)%>%summarise(p=ifelse(sum(Observed_biallelic_genotypes)>0 & sum(Expected)==0, ppois(sum(Observed_biallelic_genotypes)-1, lambda=sum(cExpected), lower.tail = F), ppois(sum(Observed_biallelic_genotypes)-1, lambda=sum(Expected), lower.tail = F)), O=sum(Observed_biallelic_genotypes), E=sum(Expected))
final$fdr<-1.0
final[which(final$Variant_class!='synonymous/synonymous'),]$fdr<-p.adjust(c(final[which(final$Variant_class!='synonymous/synonymous'),]$p, rep(1.0,(17320-length(unique(final$Gene)))*4)) , method='fdr')[1:nrow(final[which(final$Variant_class!='synonymous/synonymous'),])]
gn<-read.csv('data/gene_pos.txt.gz', stringsAsFactors=F, sep='\t')
gn<-gn[,c('Gene.stable.ID', 'HGNC.symbol')]
final<-left_join(final, gn, by=c('Gene'='Gene.stable.ID'))
write.table(final, 'sumpois.txt', row.names=F, quote=F)
