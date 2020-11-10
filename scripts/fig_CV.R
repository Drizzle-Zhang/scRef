library(ggplot2)

datasets <- c('Campbell', 'panc8_indrop', 'pbmcsca_10Xv2', 'Tasic')
GSE_ids <- c('GSE93374', 'GSE84133', 'GSE132044', 'GSE71585')
path.output <- '/home/zy/scRef/cross_validation/train4_test1/'

df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    sub.res <- read.delim(paste0(path.output, dataset, '/results_', dataset, '.txt'), stringsAsFactors = F)
    sub.res <- sub.res[sub.res$term == 'Accuracy',]
    sub.res$dataset <- rep(GSE_ids[i], dim(sub.res)[1])
    df.plot <- rbind(df.plot, sub.res)
}

accuracy.mean <- aggregate(df.plot$value, by = list(df.plot$method), FUN = mean)
df.plot$method <- factor(df.plot$method, 
                            levels = accuracy.mean$Group.1[order(accuracy.mean$x, decreasing = T)])

plot.bar <- ggplot(df.plot,
                   aes(x = method, y = value, fill = dataset)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    theme_bw() +
    # facet_wrap(~ Evaluator, scales = 'free_x', ncol = 2) +
    labs(title = "", y = 'Accuracy', x = 'DataSet', fill = 'Methods') + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = 'grey', size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12))
ggsave(filename = 'barplot_CV.png', 
       path = '/home/zy/scRef/figure/', plot = plot.bar,
       units = 'cm', height = 12, width = 18)
