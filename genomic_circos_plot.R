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

#Customising the outer band
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

#write.csv(coords, file = 'coords.csv', row.names = F)


#coords = read.csv('coords.csv')
coords$chromosome_name = paste0('chr', coords$chromosome_name)
chr_order = unique(cyto.info$Chromosome)
coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)

num = which(is.na(coords$chromosome_name))
coords = coords[-num, ]

#upregulated genes
up = which((mat$pval < 0.01) &
             (mat$logfc2 > 1))
upmat = mat[up, ]

num = which(coords$hgnc_symbol %in% rownames(upmat))
coords1 = coords[num, ]
View(coords1)
table(coords1[,1]) #total upregulated genes in each chromosome
RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)

RCircos.Get.Gene.Name.Plot.Parameters()

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


