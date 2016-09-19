library(reshape2)
library(ggplot2)
library(plyr)
library(scales)
library(grid)

#{{{ HELPER FUNCTIONS

type_of_change = function(c){gsub("(.)(.)(.)>(.)", '\\2>\\4',c)}
change_5_prime = function(c){gsub("(.)(.)(.)>(.)", '\\1',c)}
change_3_prime = function(c){gsub("(.)(.)(.)>(.)", '\\3',c)}
context        = function(c){gsub("(.)(.)(.)>(.)", '\\1\\2\\3',c)}
 
get_position = function(mutations){
  sapply(mutations, function(mutation){
    mutation = strsplit(mutation, "@", fixed = T)[[1]][1]
    position = as.numeric(strsplit(mutation, ":", fixed = T)[[1]][2]);
  })
}


get_flat_position.old = function(mutations, chr_offset){
  sapply(mutations, function(mutation){
    mutation = strsplit(mutation, "@", fixed = T)[[1]][1]
    chromosome = strsplit(mutation, ":", fixed = T)[[1]][1];
    position = as.numeric(strsplit(mutation, ":", fixed = T)[[1]][2]);
    offset = 12000 #= chr_offset[chromosome]
    sum = position + offset
    sum
  })
}

get_flat_position = function(mutation, chr_offset){
    mutation = strsplit(mutation, "@", fixed = T)[[1]][1]
    chromosome = strsplit(mutation, ":", fixed = T)[[1]][1];
    position = as.numeric(strsplit(mutation, ":", fixed = T)[[1]][2]);

    offset = chr_offset[chromosome]
    sum = position + offset

    sum
}

get_distance = function(mutations){
  sapply(mutations, function(mutation){
      as.numeric(strsplit(mutation, "@", fixed = T)[[1]][2])
  })
}

#{{{ PLOTS

signature_heatmap <- function(data, channel_counts, scale = 'none'){
    signatures = melt(as.matrix(data));
    names(signatures) = c("Change", "Sample", "Counts")

    channel_proportions = channel_counts / sum(as.numeric(channel_counts))

    signatures = ddply(signatures, .(), transform, type=type_of_change(Change))
    signatures = ddply(signatures, .(), transform, context=context(Change))
    signatures = ddply(signatures, .(), transform, adjacent.5=change_5_prime(Change))
    signatures = ddply(signatures, .(), transform, adjacent.3=change_3_prime(Change))

    signatures$channel_proportion = channel_proportions[signatures$context]
    signatures = ddply(signatures, .(), transform, channel.rescale=Counts / channel_proportion )
    signatures = ddply(signatures, .(Sample), transform, sample.rescale=rescale(Counts))
    signatures = ddply(signatures, .(Sample), transform, sample.rescale.channel=rescale(channel.rescale))

    if (scale == 'none'){ signatures$value = signatures$Counts }
    if (scale == 'sample'){ signatures$value = signatures$sample.rescale }
    if (scale == 'channel'){ signatures$value = signatures$channel.rescale }
    if (scale == 'both'){ signatures$value = signatures$sample.rescale.channel }

    p <- ggplot(signatures, aes(adjacent.5, adjacent.3));

    p <- p + theme_grey() + 
      geom_tile(aes(fill = log10(value)),colour = "white") + 
      theme(panel.margin=unit(0.02 , "lines")) + 
      theme(legend.position = "none", axis.ticks = element_blank()) + 
      theme(strip.text.y = element_text(angle = 0)) +
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0));

    p <- p + facet_grid(Sample~type, margins=F)

    p <- p + theme(axis.text.y=element_text(size=8))

    p
}

rainfall_plot <- function(data, chr_size){
    distances = melt(data)

    distances$value = as.character(distances$value)

    names(distances) = c("Mutation", "Chromosome")

    distances = ddply(distances, .(), transform, position=get_position(Mutation))
    distances = ddply(distances, .(), transform, distance=get_distance(Mutation))

    chromosomes = c("1", "2", "3",  "4",  "5",  "6",  "7",  "8",  "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "M",  "X", "Y")

    chr_offset = cumsum(c(0, chr_size[chromosomes, 1])) 
    names(chr_offset) <- chromosomes

    flat_position = apply(distances, 1, function(row){ get_flat_position(row["Mutation"], chr_offset)})
    
    distances = cbind(distances, flat_position = flat_position)

    p <- qplot(jitter(flat_position), log(sapply(distance, function(d){min(d,1000)})), data=distances, color=Chromosome, size=I(1.5))
    p <- p + geom_vline(xintercept = chr_offset[chromosomes], colour="gray", linetype = "longdash")
    p
}

factor_profile_plot <- function(factor_profiles){
    factor_profiles$Change = as.factor(rownames(factor_profiles))

    m = melt(factor_profiles)
    m = ddply(m, .(), transform, Base.change=type_of_change(Change))
    m = ddply(m, .(), transform, Context=context(Change))
    m = ddply(m, .(), transform, adjacent.3=change_3_prime(Change))

    pasted = apply(m, 1, function(x){paste(x["Base.change"],x["Change"],sep=":")})
    order = unique(sapply(sort(pasted), function(x){strsplit(x,":")[[1]][[2]]}))
    m$Change = factor(as.character(m$Change), levels=order)

    #p <- qqplot(Change ~ value, data= m, fill=Base.change, color=adjacent.3 , geom='bar', stat='identity') + facet_grid(variable~.) #+ theme(axis.text.x = element_text(angle = 90))
    #p <- p + scale_colour_brewer('adjacent3', palette="Set1")
    p <- ggplot(m, aes(x=Change, y=value)) + geom_bar(stat='identity') + facet_grid(variable~.)
    p <- p + theme(axis.text.x = element_text(angle = 90, size=10, family='mono'))      

    p
}

sample_profile_plot <- function(sample_profiles){
    mat = as.matrix(sample_profiles)
    class(mat) <- "numeric"
    mat = mat %*% diag(1/colSums(mat))
    #sample_profiles.aux = as.data.frame(mat)
    #sample_profiles.aux = sample_profiles
    #sample_profiles.aux$Factor = rownames(sample_profiles)

    m = melt(mat)
    names(m) = c("Factor", "Sample", "value")
    m$value = as.numeric(m$value)

    #p <- ggplot(m, aes(Sample,fill=Factor)) + geom_bar() + theme(axis.text.x = element_text(angle = 90))      
    #p <- ggplot(m, aes(Sample,fill=Factor)) + geom_bar()
    p <- ggplot(m, aes(x=Sample, y=value, fill=Factor)) + geom_bar(stat='identity')
    p <- p + theme(axis.text.x = element_text(angle = 90, size=5, family='mono'))      
    p
}
