args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop('Usage: plot_known_novel_spectrum.R <QC/annotation output directory>')
base <- normalizePath(args[[1]], mustWork=TRUE)
data <- read.delim(file.path(base, 'stats/known_candidate_novel_allele_frequency_spectrum.tsv'),
                   stringsAsFactors=FALSE, check.names=FALSE)
classes <- c('dbSNP_rs_match', 'candidate_novel_vs_dbSNP')
scopes <- c('autosome', 'X_nonPAR')
cols <- c('#276FBF', '#D1495B')
labels <- unique(data$AF_bin)

draw <- function() {
  old <- par(no.readonly=TRUE)
  on.exit(par(old))
  par(mfrow=c(1,2), mar=c(8,4.5,3.3,1), oma=c(2.5,0,2,0))
  for (scope in scopes) {
    d <- data[data$genomic_class == scope & data$variant_class %in% classes,]
    values <- d$proportion_within_class
    values <- values[is.finite(values)]
    ymax <- if (length(values)) max(values) else 1
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    plot(seq_along(labels), rep(NA_real_, length(labels)), type='n', xaxt='n',
         xlab='', ylab='Proportion within variant class', ylim=c(0,ymax*1.08),
         main=if (scope == 'autosome') 'Autosomes 1-19' else 'X non-PAR')
    axis(1, at=seq_along(labels), labels=labels, las=2, cex.axis=0.62)
    for (i in seq_along(classes)) {
      x <- d[d$variant_class == classes[[i]],]
      x <- x[match(labels, x$AF_bin),]
      y <- x$proportion_within_class
      if (length(y) && any(is.finite(y)))
        lines(seq_along(labels), y, type='b', pch=19,
              cex=0.55, lwd=1.6, col=cols[[i]])
    }
    legend('topright', legend=c('dbSNP rs match','Candidate novel vs dbSNP'),
           col=cols, lty=1, pch=19, bty='n', cex=0.78)
  }
  mtext('Main set allele frequency spectrum by VEP dbSNP allele match', outer=TRUE, cex=1.1)
  mtext('AF = AC/AN; each line is normalized within its variant class. Sites without CSQ are reported separately.',
        side=1, outer=TRUE, line=0.2, cex=0.72)
}

svg(file.path(base, 'plots/known_candidate_novel_AF_spectrum.svg'), width=12, height=7)
draw()
dev.off()
pdf(file.path(base, 'plots/known_candidate_novel_AF_spectrum.pdf'), width=12, height=7, useDingbats=FALSE)
draw()
dev.off()
