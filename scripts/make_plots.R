args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
w <- read.delim(file.path(base,"stats/window_variant_density_1Mb.tsv"),check.names=FALSE)
h <- read.delim(file.path(base,"stats/heterozygosity_per_sample.tsv"),check.names=FALSE)
f <- read.delim(file.path(base,"stats/allele_frequency_spectrum.tsv"),check.names=FALSE)
plot_density <- function() {
  chr_order <- c(as.character(1:19),"X")
  op <- par(mfrow=c(4,5),mar=c(2.7,2.8,2,0.8),oma=c(1,1,2,0),mgp=c(1.7,.5,0))
  on.exit(par(op))
  for (chr in chr_order) {
    d <- w[w$chrom==chr,]
    if (nrow(d)==0) next
    plot((d$window_start_1based-1)/1e6,d$snps_per_Mb,type="h",lwd=1,
      xlab="Position (Mb)",ylab="SNP/Mb",main=paste0("Chr ",chr),
      col=if(chr=="X") "#8B5CF6" else "#2563EB",cex.main=.85,cex.axis=.7,cex.lab=.75)
    abline(h=0,col="grey80")
  }
  mtext("1-Mb SNP density in autosomes and X non-PAR intervals",outer=TRUE,cex=1.05)
}
plot_het <- function() {
  g <- factor(h$group,levels=c("ICR1","ICR2","KM"))
  cols <- c("#2563EB","#DC2626","#059669")
  boxplot(h$heterozygous_proportion_of_nonmissing_diploid_GT~g,
    col=cols,outline=FALSE,xlab="Group",ylab="Heterozygous GT / nonmissing diploid autosomal GT",
    main="Per-sample autosomal heterozygosity")
  stripchart(h$heterozygous_proportion_of_nonmissing_diploid_GT~g,
    vertical=TRUE,method="jitter",add=TRUE,pch=16,cex=.65,col="#22222288")
  grid(nx=NA,ny=NULL,col="grey85")
}
plot_af <- function() {
  d <- f[f$genomic_class %in% c("autosome","X_total"),]
  bins <- unique(d$AF_bin)
  groups <- c("autosome","X_total")
  mat <- sapply(groups,function(g) {
    x <- d[d$genomic_class==g,]
    setNames(x$proportion_within_class, x$AF_bin)[bins]
  })
  mat[is.na(mat)] <- 0
  barplot(t(mat),beside=TRUE,names.arg=bins,las=2,cex.names=.65,
    col=c("#2563EB","#8B5CF6"),xlab="Alternate allele frequency bin (AC/AN)",
    ylab="Proportion within genomic class",
    main="Main-set allele-frequency spectrum by genomic class")
  legend("topright",legend=c("Autosome 1–19","X total"),fill=c("#2563EB","#8B5CF6"),bty="n",cex=.8)
}
for (ext in c("svg","png")) {
  if (ext=="svg") svg(file.path(base,"plots/window_variant_density_1Mb.svg"),width=14,height=9)
  else png(file.path(base,"plots/window_variant_density_1Mb.png"),width=2800,height=1800,res=200)
  plot_density(); dev.off()
  if (ext=="svg") svg(file.path(base,"plots/sample_heterozygosity_by_group.svg"),width=7,height=5)
  else png(file.path(base,"plots/sample_heterozygosity_by_group.png"),width=1400,height=1000,res=200)
  plot_het(); dev.off()
  if (ext=="svg") svg(file.path(base,"plots/allele_frequency_spectrum.svg"),width=10,height=6)
  else png(file.path(base,"plots/allele_frequency_spectrum.png"),width=2000,height=1200,res=200)
  plot_af(); dev.off()
}
cat("PLOTS_OK\n")
