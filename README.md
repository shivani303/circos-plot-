# circos-plot

library(RCircos)

data("UCSC.HG38.Human.CytoBandIdeogram")

cyto.info = UCSC.HG38.Human.CytoBandIdeogram

RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)
                            
#track inside -to view something inside the circle similar track outside

RCircos.Set.Plot.Area() 

RCircos.Chromosome.Ideogram.Plot()

![image](https://user-images.githubusercontent.com/66779651/197962126-09bc0342-e5bc-47b6-bc4e-54b35403fa1c.png)

# Customising the outer band

cyto.info = UCSC.HG38.Human.CytoBandIdeogram

cyto.info$Name = NA

cyto.info$Stain = NA

RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)
                            
RCircos.Set.Plot.Area()

RCircos.Chromosome.Ideogram.Plot()

ideo = RCircos.Get.Plot.Ideogram()

ideo$BandColor = 'purple'

num = which(ideo$Chromosome == 'chrX')

ideo[num, 'BandColor'] = 'chartreuse'

num = which(ideo$Chromosome == 'chrY')

ideo[num, 'BandColor'] = 'orange'

RCircos.Reset.Plot.Ideogram(ideo)

RCircos.Set.Plot.Area()

RCircos.Chromosome.Ideogram.Plot()
![image](https://user-images.githubusercontent.com/66779651/197962582-eee060da-a15c-4e49-8b43-d76bdc58fbf5.png)

library(biomaRt)

mat= readRDS("mat.RDS")
View(mat)

m = useMart('ensembl', dataset='hsapiens_gene_ensembl')

coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(mat)),
               mart = m)
View(coords)

coords$chromosome_name = paste0('chr', coords$chromosome_name)

chr_order = unique(cyto.info$Chromosome)

coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)


num = which(is.na(coords$chromosome_name))

coords = coords[-num, ]

# upregulated genes
up = which((mat$pval < 0.01) &
             (mat$logfc2 > 1))
upmat = mat[up, ]

num = which(coords$hgnc_symbol %in% rownames(upmat))

coords1 = coords[num, ]

View(coords1)

table(coords1[,1]) #total upregulated genes in each chromosome
![image](https://user-images.githubusercontent.com/66779651/197965748-adce5708-81f3-4851-a6e1-9f88a3c6a584.png)


RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)

RCircos.Get.Gene.Name.Plot.Parameters()
![image](https://user-images.githubusercontent.com/66779651/197965869-ec570b8d-005f-4ff5-9917-8bde31190613.png)

genes = intersect(rownames(mat), coords$hgnc_symbol)

mat1 = mat[genes, ]

df = cbind.data.frame(rownames(mat1), mat1[, c(1,2,4)])

colnames(df)[1] = 'hgnc_symbol'

data = merge(coords, df, by = 'hgnc_symbol', all.x = T)

data = data[, c('chromosome_name', 'start_position',
                'end_position', 'hgnc_symbol',
                'control', 'test', 'logfc2')]
summary(mat$logfc2)
RCircos.Heatmap.Plot(data, data.col = 7, track.num = 6, side = "in",
                     min.value = -3, max.value = 2, 
                     is.sorted = F)

RC.param = RCircos.Get.Plot.Parameters()

RC.param['heatmap.color'] = "GreenWhiteRed"

RCircos.Reset.Plot.Parameters(RC.param)

RCircos.Heatmap.Plot(data, data.col = 7, track.num = 10, side = "in",
                     min.value = -3, max.value = 2,
                     is.sorted = F)

![image](https://user-images.githubusercontent.com/66779651/199909944-8c44631b-0517-4528-a1a7-deb31c722344.png)

#INTERPRETATION

Circos plots are a great way to show genomic data. Circos can illustrate genomic rearrangements, where a relationship between two elements (genomic positions) represents a structural fusion.
Above ideogram plot shows the differentially expressed upregulated genes expression patterns and the distribution of the chromosomal location where they are located, with the outer circle representing the chromosome and the location of the gene in the chromosome, and the heatmap in the inner circle representing the expression of all upregulated differentially expressed genes (DEGs)
analysis of which genes from sample dataset are present on which chromosome and their regulation i.e.. up-regulated and down-regulated is done.
Outermost circle shows the chormosomes (autosomes + sex chromosomes) , chromosome number and genes from sample data present on which chromosome.
Inner circle shows the heatmap of regulation of genes, which gene is down-regulated and up-regulated.
In above plot we can see that most of the genes are present on chromosome number 16. 


