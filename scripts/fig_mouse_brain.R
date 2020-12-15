library(ggplot2)
library(gridExtra)

ref.dataset <- c('Tasic', 'MCA')
ref.GSE_id <- c('GSE71585', 'GSE108097')

sc.dataset <- c('Tasic2018', 'Habib', 'HochgernerA', 'Mizrak')
sc.GSE_id <- c('GSE115746', 'GSE93374', 'GSE95315', 'GSE109447')

path.res <- '/home/zy/scRef/Benchmark/mouse_brain/'
path.fig <- '/home/zy/scRef/figure/mouse_brain'

# Accuracy
for (i in 1:length(ref.dataset)) {
    ref.data <- ref.dataset[i]
    ref.GSE <- ref.GSE_id[i]
    df.plot <- data.frame(stringsAsFactors = F)
    for (j in 1:length(sc.dataset)) {
        sc.data <- sc.dataset[j]
        sc.GSE <- sc.GSE_id[j]
        file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
        sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
        sub.res <- sub.res[sub.res$term == 'Accuracy',]
        acc.scRef <- sub.res[sub.res$method == 'scMAGIC', 'value']
        sub.plot <- data.frame(other = sub.res[sub.res$method != 'scMAGIC', 'value'], 
                               scRef = rep(acc.scRef, nrow(sub.res)-1), 
                               method = sub.res[sub.res$method != 'scMAGIC', 'method'],
                               dataset = rep(sc.GSE, nrow(sub.res)-1))
        df.plot <- rbind(df.plot, sub.plot)
    }
    plot.scatter <- 
        ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
        geom_point(size = 3) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + 
        ylim(0, 1) + 
        labs(x = 'Accuracy(Other method)', y = 'Accuracy(scMAGIC)', 
             shape = 'Query dataset', color = 'Method') + 
        theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
    ggsave(filename = paste0('scatter_', ref.data, '.png'), 
           path = path.fig, plot = plot.scatter,
           units = 'cm', height = 12, width = 16)
}

# Macro F1
for (i in 1:length(ref.dataset)) {
    ref.data <- ref.dataset[i]
    ref.GSE <- ref.GSE_id[i]
    df.plot <- data.frame(stringsAsFactors = F)
    for (j in 1:length(sc.dataset)) {
        sc.data <- sc.dataset[j]
        sc.GSE <- sc.GSE_id[j]
        file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
        sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
        sub.res <- sub.res[sub.res$term == 'Macro F1',]
        acc.scRef <- sub.res[sub.res$method == 'scMAGIC', 'value']
        sub.plot <- data.frame(other = sub.res[sub.res$method != 'scMAGIC', 'value'], 
                               scRef = rep(acc.scRef, nrow(sub.res)-1), 
                               method = sub.res[sub.res$method != 'scMAGIC', 'method'],
                               dataset = rep(sc.GSE, nrow(sub.res)-1))
        df.plot <- rbind(df.plot, sub.plot)
    }
    plot.scatter <- 
        ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
        geom_point(size = 3) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + 
        ylim(0, 1) + xlim(0, 1) + 
        labs(x = 'Macro F1(Other method)', y = 'Macro F1(scMAGIC)', 
             shape = 'Query dataset', color = 'Method') + 
        theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
    ggsave(filename = paste0('scatter_MacroF1_', ref.data, '.png'), 
           path = path.fig, plot = plot.scatter,
           units = 'cm', height = 12, width = 16)
}

# Accuracy
for (i in 1:length(ref.dataset)) {
    ref.data <- ref.dataset[i]
    ref.GSE <- ref.GSE_id[i]
    df.plot <- data.frame(stringsAsFactors = F)
    for (j in 1:length(sc.dataset)) {
        sc.data <- sc.dataset[j]
        sc.GSE <- sc.GSE_id[j]
        file.res <- paste0(path.res, 'results_', ref.data, '_', sc.data, '.txt')
        sub.res <- read.delim(file = file.res, row.names = 1, stringsAsFactors = F)
        sub.res <- sub.res[sub.res$term == 'Accuracy',]
        acc.scRef <- sub.res[sub.res$method == 'scMAGIC', 'value']
        sub.plot <- data.frame(other = sub.res[sub.res$method != 'scMAGIC', 'value'], 
                               scRef = rep(acc.scRef, nrow(sub.res)-1), 
                               method = sub.res[sub.res$method != 'scMAGIC', 'method'],
                               dataset = rep(sc.GSE, nrow(sub.res)-1))
        df.plot <- rbind(df.plot, sub.plot)
    }
    plot.scatter <- 
        ggplot(data = df.plot, aes(x = other, y = scRef, color = method, shape = dataset)) + 
        geom_point(size = 3) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + 
        ylim(0, 1) + 
        labs(x = 'Accuracy(Other method)', y = 'Accuracy(scMAGIC)', 
             shape = 'Query dataset', color = 'Method') + 
        theme(panel.background = element_rect(fill = 'transparent', color = 'gray'))
    ggsave(filename = paste0('scatter_', ref.data, '.png'), 
           path = path.fig, plot = plot.scatter,
           units = 'cm', height = 12, width = 16)
}


# heatmap
plot.heatmap <- function(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref) {
    rda.scRef <- paste0(path.res, ref.dataset, '_', dataset, '_', method, '.Rdata')
    pred.scRef <- readRDS(rda.scRef)
    pred.scRef[!(pred.scRef %in% names.ref)] <- 'Unassigned'
    mytable <- table(true.tags, pred.scRef)
    mydata <- data.frame(stringsAsFactors = F)
    table.true <- table(true.tags)
    for (label1 in rownames(mytable)) {
        row.sum <- table.true[label1]
        for (label2 in colnames(mytable)) {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                                count = mytable[label1, label2], 
                                               prop = mytable[label1, label2]/row.sum))
        }
    }
    mydata$origin <- factor(mydata$origin, levels = names.sc)
    mydata$annotation <- factor(mydata$annotation, levels = c(names.ref, 'Unassigned'))
    
    plot.heatmap <- 
        ggplot(data = mydata, aes(x = origin, y = annotation)) + 
        geom_tile(aes(fill = prop)) + 
        scale_fill_continuous(low = "#FFFAFA", high = "#CD5C5C") + 
        labs(fill = 'Proportion', title = method) +
        # labs(fill = 'Proportion', title = 'scMAGIC') +
        theme_bw() +
        theme(
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        ) + 
        geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5)
    ggsave(filename = paste0('heatmap_', ref.dataset, '_', dataset, '_', method, '.png'), 
           path = path.fig, plot = plot.heatmap,
           units = 'cm', height = 12, width = 22)
}

ref.dataset <- 'Tasic'
dataset <- 'Habib'
OUT <- readRDS(paste0(path.res, dataset, '.Rdata'))
label_sc <- OUT$label.filter
true.tags <- label_sc$label.unlabeled.use.cols...
names.ref <- c('Neuron', 'Oligodendrocyte', 'Oligodendrocyte Precursor Cell', 
               'Microglia', 'Endothelial Cell', 'Astrocyte')
names.sc <- c('Neurons', 'Oligodend', 'OPC', 'microglia', 'EndothelialCells', 'Astrocyte', 
              'Ependymocytes', 'Fibroblast', 'MuralCells', 'ParsTuber', 'Tanycyte')

method <- 'scRef'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'sciBet'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'scPred'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'singleCellNet'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

method <- 'scClassify'
plot.heatmap(ref.dataset, dataset, true.tags, path.res, path.fig, method, names.sc, names.ref)

