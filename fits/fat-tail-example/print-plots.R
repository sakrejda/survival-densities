
pl_quantile <- readRDS('qq-plots.rds')
pl_histogram <- readRDS('histogram-plots.rds')

pdf(file='plots.pdf', width=9.5, height=7)
print(pl_quantile$comparison)
print(pl_histogram$comparison)
print(pl_histogram$full_comparison)
dev.off()




