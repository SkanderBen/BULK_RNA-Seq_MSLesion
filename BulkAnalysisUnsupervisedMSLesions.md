Bulk analysis of MS lesions(Unsupervised)
================

This is an analysis performed on the bulk RNA-Seq dataset

``` r
library(edgeR)
library(fgsea)
library(qusage)
library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichR)
library(VennDiagram)
library(ggpubr)
library(ggrepel)
library(PCAtools)
library(RColorBrewer)
library(ggiraph)
library(UpSetR)
library(rhandsontable)
library(nloptr)
library(ggVennDiagram)
library(BiocManager)
library(stringr)
library(knitr)
library(ggfortify)
```

### Counts list

``` r
dat = read.csv('/Users/skander/Library/CloudStorage/OneDrive-Personnel/Documents/UDEM/Master/BulkData/data/bulkExtra/Totalbulk/merged_counts.csv',sep=';', header = TRUE, row.names=1)
knitr::kable(dat[1:12, 1:12], caption ="A Knitr kable")
```

|                 | gene_name | S13_b1 | S14_b1 |    S23 |    S26 |    S29 |    S33 |     S37 |    S38 |  S3_b1 |     S43 |     S45 |
|:----------------|:----------|-------:|-------:|-------:|-------:|-------:|-------:|--------:|-------:|-------:|--------:|--------:|
| ENSG00000000003 | TSPAN6    | 131.00 |  46.00 |  83.00 |  59.00 |  79.00 |  98.00 |   56.00 |  57.00 | 220.00 |   74.00 |  725.00 |
| ENSG00000000005 | TNMD      |   1.00 |   1.00 |   0.00 |   2.00 |   2.00 |   1.00 |    0.00 |   0.00 |   5.00 |    1.00 |    3.00 |
| ENSG00000000419 | DPM1      |  81.00 | 103.00 | 121.00 | 128.00 | 141.00 |  91.00 |   61.00 |  68.00 |  66.00 |  122.00 |  112.00 |
| ENSG00000000457 | SCYL3     |  91.24 | 105.74 | 120.97 | 107.41 | 115.13 |  93.44 |   64.38 |  69.50 |  74.44 |  168.65 |  191.43 |
| ENSG00000000460 | C1orf112  |  95.86 |  56.26 |  72.03 |  68.59 |  67.87 |  44.56 |   59.62 |  50.50 |  91.56 |   88.35 |  163.57 |
| ENSG00000000938 | FGR       |  28.00 |  11.00 |  15.00 |  16.00 |   6.00 |  58.00 |   44.00 |  28.00 |  56.00 |   16.00 |  113.00 |
| ENSG00000000971 | CFH       | 103.00 |  41.00 |  90.00 |  70.00 |  41.00 | 198.00 |   59.00 | 151.00 | 100.00 |   97.00 |  105.00 |
| ENSG00000001036 | FUCA2     | 104.60 |  70.29 |  65.07 |  75.40 |  73.47 |  83.00 |   62.29 |  71.67 |  89.90 |   81.02 |  228.64 |
| ENSG00000001084 | GCLC      | 482.00 | 451.00 | 406.00 | 442.00 | 500.00 | 466.97 |  395.00 | 344.00 | 388.00 |  559.00 | 1042.00 |
| ENSG00000001167 | NFYA      | 533.00 | 325.00 | 224.00 | 333.97 | 307.00 | 340.00 |  508.89 | 344.00 | 450.00 |  260.00 |  491.00 |
| ENSG00000001460 | STPG1     | 173.02 | 227.16 | 157.14 | 155.07 | 183.09 | 140.37 |  131.45 | 146.10 | 182.06 |  205.00 |  214.00 |
| ENSG00000001461 | NIPAL3    | 809.98 | 963.84 | 656.86 | 670.93 | 757.91 | 776.63 | 1438.55 | 675.90 | 692.94 | 1072.00 |  992.00 |

A Knitr kable

### Meta data

``` r
meta = read.csv('/Users/skander/Library/CloudStorage/OneDrive-Personnel/Documents/UDEM/Master/BulkData/data/bulkExtra/Totalbulk/meta_data_3batch.csv', header = T,sep = ';', row.names = 1)
names(meta)[names(meta) == "type"] <- "lesion_type"
meta$sex = 'M'
meta[which(meta$sample=="AB203"),]$sex = 'F'
meta[which(meta$sample=="AB253"),]$sex = 'F'

meta$disease="MS"
meta[which(meta$sample %in% c("AB257","AB251","AB252")),]$disease="HC"
#meta$region="WM"
meta$group=paste(meta$disease,meta$tissue,sep=".")
meta$reads=colSums(dat[-1])[rownames(meta)]

meta$condition = 'Normal'
meta[which(meta$lesion_type %in% c("WML","CL","PVML")),'condition']="Lesion"
meta[grepl('^NA',meta$lesion_type),'condition']='Normal_Appearing'

knitr::kable(meta, caption ="A Knitr kable")
```

|          | id       | batch  | lesion_type | tissue   | condition        | sample | typeL        | sex | disease | group       |    reads |
|:---------|:---------|:-------|:------------|:---------|:-----------------|:-------|:-------------|:----|:--------|:------------|---------:|
| S13_b1   | S13_b1   | batch1 | WML         | WM       | Lesion           | AB187  |              | M   | MS      | MS.WM       | 10675040 |
| S14_b1   | S14_b1   | batch1 | NAGM        | Cortical | Normal_Appearing | AB187  |              | M   | MS      | MS.Cortical | 11601256 |
| S23      | S23      | batch1 | CL          | Cortical | Lesion           | AB203  |              | F   | MS      | MS.Cortical | 10388672 |
| S26      | S26      | batch1 | WML         | WM       | Lesion           | AB187  |              | M   | MS      | MS.WM       | 11312959 |
| S29      | S29      | batch1 | CL          | Cortical | Lesion           | AB187  |              | M   | MS      | MS.Cortical | 10849903 |
| S33      | S33      | batch1 | CL          | Cortical | Lesion           | AB200  |              | M   | MS      | MS.Cortical | 15027141 |
| S37      | S37      | batch1 | WML         | WM       | Lesion           | AB200  |              | M   | MS      | MS.WM       | 12031415 |
| S38      | S38      | batch1 | NAGM        | Cortical | Normal_Appearing | AB200  |              | M   | MS      | MS.Cortical | 12856284 |
| S3_b1    | S3_b1    | batch1 | PVML        | PV       | Lesion           | AB200  |              | M   | MS      | MS.PV       | 10431353 |
| S43      | S43      | batch1 | CL          | Cortical | Lesion           | AB203  |              | F   | MS      | MS.Cortical | 12121099 |
| S45      | S45      | batch1 | PVML        | PV       | Lesion           | AB203  |              | F   | MS      | MS.PV       | 10093676 |
| S53      | S53      | batch1 | NAWM        | WM       | Normal_Appearing | AB187  |              | M   | MS      | MS.WM       | 10907244 |
| S55      | S55      | batch1 | NAGM        | Cortical | Normal_Appearing | AB200  |              | M   | MS      | MS.Cortical | 20126720 |
| S62      | S62      | batch1 | WML         | WM       | Lesion           | AB200  |              | M   | MS      | MS.WM       | 12859996 |
| S64      | S64      | batch1 | NAGM        | Cortical | Normal_Appearing | AB203  |              | F   | MS      | MS.Cortical | 13502425 |
| S67      | S67      | batch1 | WML         | WM       | Lesion           | AB203  |              | F   | MS      | MS.WM       | 11484555 |
| S70      | S70      | batch1 | WML         | WM       | Lesion           | AB203  |              | F   | MS      | MS.WM       | 10809418 |
| S71      | S71      | batch1 | NAWM        | WM       | Normal_Appearing | AB187  |              | M   | MS      | MS.WM       | 13315868 |
| S72      | S72      | batch1 | CL          | Cortical | Lesion           | AB187  |              | M   | MS      | MS.Cortical | 12938331 |
| S76      | S76      | batch1 | CL          | Cortical | Lesion           | AB200  |              | M   | MS      | MS.Cortical | 14486529 |
| S77      | S77      | batch1 | WML         | WM       | Lesion           | AB200  |              | M   | MS      | MS.WM       | 14583455 |
| S79      | S79      | batch1 | NAWM        | WM       | Normal_Appearing | AB203  |              | F   | MS      | MS.WM       | 14931094 |
| S9_b1    | S9_b1    | batch1 | NAGM        | Cortical | Normal_Appearing | AB187  |              | M   | MS      | MS.Cortical | 14719787 |
| Sample1  | Sample1  | batch2 | PVML        | PV       | Lesion           | AB200  |              | M   | MS      | MS.PV       | 20735000 |
| Sample10 | Sample10 | batch2 | CL          | Cortical | Lesion           | AB187  |              | M   | MS      | MS.Cortical |  9925960 |
| Sample11 | Sample11 | batch2 | WML         | WM       | Lesion           | AB187  |              | M   | MS      | MS.WM       | 19113957 |
| Sample12 | Sample12 | batch2 | CL          | Cortical | Lesion           | AB187  |              | M   | MS      | MS.Cortical | 13267587 |
| Sample15 | Sample15 | batch2 | WML         | WM       | Lesion           | AB187  |              | M   | MS      | MS.WM       | 14887133 |
| Sample16 | Sample16 | batch2 | WML         | WM       | Lesion           | AB187  |              | M   | MS      | MS.WM       | 12598401 |
| Sample17 | Sample17 | batch2 | PVML        | PV       | Lesion           | AB203  |              | F   | MS      | MS.PV       | 10162274 |
| Sample18 | Sample18 | batch2 | CL          | Cortical | Lesion           | AB203  |              | F   | MS      | MS.Cortical | 11609221 |
| Sample21 | Sample21 | batch2 | NAWM        | WM       | Normal_Appearing | AB203  |              | F   | MS      | MS.WM       | 11016594 |
| Sample24 | Sample24 | batch2 | NAWM        | WM       | Normal_Appearing | AB203  |              | F   | MS      | MS.WM       | 14238903 |
| Sample25 | Sample25 | batch2 | NAGM        | Cortical | Normal_Appearing | AB187  |              | M   | MS      | MS.Cortical | 15168403 |
| Sample34 | Sample34 | batch2 | WML         | WM       | Lesion           | AB200  |              | M   | MS      | MS.WM       | 18291181 |
| Sample46 | Sample46 | batch2 | NAGM        | Cortical | Normal_Appearing | AB203  |              | F   | MS      | MS.Cortical | 16625491 |
| Sample52 | Sample52 | batch2 | NAWM        | WM       | Normal_Appearing | AB187  |              | M   | MS      | MS.WM       | 12701779 |
| Sample56 | Sample56 | batch2 | PVML        | PV       | Lesion           | AB200  |              | M   | MS      | MS.PV       | 13165303 |
| Sample57 | Sample57 | batch2 | CL          | Cortical | Lesion           | AB200  |              | M   | MS      | MS.Cortical | 11305279 |
| Sample6  | Sample6  | batch2 | NAGM        | Cortical | Normal_Appearing | AB200  |              | M   | MS      | MS.Cortical | 22662880 |
| Sample60 | Sample60 | batch2 | NAWM        | WM       | Normal_Appearing | AB200  |              | M   | MS      | MS.WM       | 14547748 |
| Sample65 | Sample65 | batch2 | PVML        | PV       | Lesion           | AB203  |              | F   | MS      | MS.PV       | 13327909 |
| Sample66 | Sample66 | batch2 | WML         | WM       | Lesion           | AB203  |              | F   | MS      | MS.WM       |  8626915 |
| Sample7  | Sample7  | batch2 | WML         | WM       | Lesion           | AB200  |              | M   | MS      | MS.WM       | 17559199 |
| Sample74 | Sample74 | batch2 | NAGM        | Cortical | Normal_Appearing | AB187  |              | M   | MS      | MS.Cortical | 17633780 |
| Sample75 | Sample75 | batch2 | NAWM        | WM       | Normal_Appearing | AB200  |              | M   | MS      | MS.WM       | 20237852 |
| S1       | S1       | batch3 | GM          | Cortical | Normal           | AB257  | GM           | M   | HC      | HC.Cortical | 14746041 |
| S2       | S2       | batch3 | GM          | Cortical | Normal           | AB257  | GM           | M   | HC      | HC.Cortical | 16135512 |
| S3_b3    | S3_b3    | batch3 | WM          | WM       | Normal           | AB257  | WM           | M   | HC      | HC.WM       | 17280964 |
| S4       | S4       | batch3 | WM          | WM       | Normal           | AB257  | WM           | M   | HC      | HC.WM       | 16128150 |
| S5       | S5       | batch3 | PVM         | PV       | Normal           | AB257  | PVM          | M   | HC      | HC.PV       | 15969401 |
| S6       | S6       | batch3 | PVM         | PV       | Normal           | AB257  | PVM          | M   | HC      | HC.PV       | 14119977 |
| S7       | S7       | batch3 | GM          | Cortical | Normal           | AB252  | GM           | M   | HC      | HC.Cortical | 16993959 |
| S8       | S8       | batch3 | GM          | Cortical | Normal           | AB252  | GM           | M   | HC      | HC.Cortical | 16902227 |
| S9_b3    | S9_b3    | batch3 | WM          | WM       | Normal           | AB252  | WM           | M   | HC      | HC.WM       | 14200304 |
| S10      | S10      | batch3 | WM          | WM       | Normal           | AB252  | WM           | M   | HC      | HC.WM       | 14993001 |
| S11      | S11      | batch3 | PVM         | PV       | Normal           | AB252  | PVM          | M   | HC      | HC.PV       | 13322475 |
| S12      | S12      | batch3 | PVM         | PV       | Normal           | AB252  | PVM          | M   | HC      | HC.PV       | 11873861 |
| S13_b3   | S13_b3   | batch3 | GM          | Cortical | Normal           | AB251  | GM           | M   | HC      | HC.Cortical | 14540503 |
| S14_b3   | S14_b3   | batch3 | GM          | Cortical | Normal           | AB251  | GM           | M   | HC      | HC.Cortical | 15068466 |
| S15      | S15      | batch3 | WM          | WM       | Normal           | AB251  | WM           | M   | HC      | HC.WM       | 13301374 |
| S16      | S16      | batch3 | WM          | WM       | Normal           | AB251  | WM           | M   | HC      | HC.WM       | 15610086 |
| S17      | S17      | batch3 | PVM         | PV       | Normal           | AB251  | PVM          | M   | HC      | HC.PV       | 12503687 |
| S18      | S18      | batch3 | PVM         | PV       | Normal           | AB251  | PVM          | M   | HC      | HC.PV       | 15465598 |
| S102f    | S102f    | batch3 | CL          | Cortical | Lesion           | AB253  | CL           | F   | MS      | MS.Cortical | 11101420 |
| S104f    | S104f    | batch3 | CL          | Cortical | Lesion           | AB253  | CL           | F   | MS      | MS.Cortical | 16624705 |
| S106f    | S106f    | batch3 | CL          | Cortical | Lesion           | AB253  | CL           | F   | MS      | MS.Cortical | 18467343 |
| S108f    | S108f    | batch3 | NAGM        | Cortical | Normal_Appearing | AB253  | NAGM         | F   | MS      | MS.Cortical | 12319237 |
| S109f    | S109f    | batch3 | NAGM        | Cortical | Normal_Appearing | AB253  | NAGM         | F   | MS      | MS.Cortical | 14185427 |
| S110f    | S110f    | batch3 | CL          | Cortical | Lesion           | AB253  | CL           | F   | MS      | MS.Cortical | 12999955 |
| S112f    | S112f    | batch3 | NAWM        | WM       | Normal_Appearing | AB253  | NAWM         | F   | MS      | MS.WM       | 13363398 |
| S113f    | S113f    | batch3 | NAWM        | WM       | Normal_Appearing | AB253  | NAWM         | F   | MS      | MS.WM       | 16157577 |
| S114f    | S114f    | batch3 | CL          | Cortical | Lesion           | AB253  | CL           | F   | MS      | MS.Cortical | 16266665 |
| S115f    | S115f    | batch3 | NAGM        | Cortical | Normal_Appearing | AB253  | NAGM         | F   | MS      | MS.Cortical | 17885078 |
| S118f    | S118f    | batch3 | PVML        | PV       | Lesion           | AB253  | Mixed.PD     | F   | MS      | MS.PV       | 10443427 |
| S119f    | S119f    | batch3 | PVML        | PV       | Lesion           | AB253  | Active.PD    | F   | MS      | MS.PV       | 13248535 |
| S120f    | S120f    | batch3 | PVML        | PV       | Lesion           | AB253  | Active.PD    | F   | MS      | MS.PV       | 12139042 |
| S121f    | S121f    | batch3 | PVML        | PV       | Lesion           | AB253  | Mixed.PD     | F   | MS      | MS.PV       | 10019098 |
| S19a     | S19a     | batch3 | WML         | WM       | Lesion           | AB253  | Active       | F   | MS      | MS.WM       | 17274976 |
| S22a     | S22a     | batch3 | PVML        | PV       | Lesion           | AB253  | PVML         | F   | MS      | MS.PV       | 11231127 |
| S28b     | S28b     | batch3 | WML         | WM       | Lesion           | AB187  | Active.PD    | M   | MS      | MS.WM       | 15599429 |
| S30b     | S30b     | batch3 | NAWM        | WM       | Normal_Appearing | AB187  | NAWM         | M   | MS      | MS.WM       | 16356084 |
| S31b     | S31b     | batch3 | WML         | WM       | Lesion           | AB187  | Active.PD    | M   | MS      | MS.WM       | 13453907 |
| S51c     | S51c     | batch3 | CL          | Cortical | Lesion           | AB187  | CL           | M   | MS      | MS.Cortical | 16757776 |
| S69c     | S69c     | batch3 | PVML        | PV       | Lesion           | AB203  | PVML         | F   | MS      | MS.PV       | 13814668 |
| S73d     | S73d     | batch3 | WML         | WM       | Lesion           | AB187  | Active.PD    | M   | MS      | MS.WM       | 13243821 |
| S78d     | S78d     | batch3 | NAWM        | WM       | Normal_Appearing | AB200  | NAWM         | M   | MS      | MS.WM       | 14215553 |
| S81d     | S81d     | batch3 | PVML        | PV       | Lesion           | AB203  | PVML         | F   | MS      | MS.PV       | 15303502 |
| S89f     | S89f     | batch3 | WML         | WM       | Lesion           | AB253  | Early.Active | F   | MS      | MS.WM       | 10392586 |
| S90f     | S90f     | batch3 | WML         | WM       | Lesion           | AB253  | Early.Active | F   | MS      | MS.WM       | 18177867 |
| S91f     | S91f     | batch3 | WML         | WM       | Lesion           | AB253  | Early.Active | F   | MS      | MS.WM       | 11856582 |
| S92f     | S92f     | batch3 | NAGM        | Cortical | Normal_Appearing | AB253  | NAGM         | F   | MS      | MS.Cortical | 12412586 |
| S94f     | S94f     | batch3 | WML         | WM       | Lesion           | AB253  | Mixed.PD     | F   | MS      | MS.WM       | 13343462 |
| S96f     | S96f     | batch3 | WML         | WM       | Lesion           | AB253  | Active.AD    | F   | MS      | MS.WM       |  8964589 |

A Knitr kable

``` r
go = read.gmt('/Users/skander/Library/CloudStorage/OneDrive-Personnel/Documents/UDEM/Master/BulkData/GO/c5.go.bp.v7.5.1.symbols.gmt')
mf = read.gmt('/Users/skander/Library/CloudStorage/OneDrive-Personnel/Documents/UDEM/Master/BulkData/GO/c5.go.mf.v7.5.1.symbols.gmt')
kegg = read.gmt('/Users/skander/Library/CloudStorage/OneDrive-Personnel/Documents/UDEM/Master/BulkData/GO/c2.cp.kegg.v7.5.1.symbols.gmt')
tft = read.gmt('/Users/skander/Library/CloudStorage/OneDrive-Personnel/Documents/UDEM/Master/BulkData/GO/c3.tft.v7.5.1.symbols.gmt')

sign_cols=c(DN='deepskyblue',NS='grey80',UP='tomato')

plot.gene <- function(x){
  gene_id = rownames(subset(resA,gene_name==x))[1]
  meta$gene = y$E[gene_id,rownames(meta)]
  ggplot(meta,aes(lesion_type,gene))+geom_boxplot()+geom_point(aes(col=batch))+theme_bw()+ ggtitle(x)
}

fisher_enrichment<-function(cluster_markers,pathway,universe,pathway_name){
enrichment_tables=list() 
  for(cluster in 1:length(cluster_markers)){
    cluster_name=names(cluster_markers)[cluster]
    output <- lapply(pathway, function(x) {
    Freq.table <- table(factor(universe %in% as.character(cluster_markers[[cluster]]), 
                          levels = c(TRUE,FALSE)), 
                          factor(universe %in% x, 
                          levels = c(TRUE,FALSE)))

      fit <- fisher.test(Freq.table, alternative = "greater")
      interSection <- intersect(cluster_markers[[cluster]], x)
      interSection <- paste(interSection, collapse = ",")
      return(value = c(NOM_pval = fit$p.value, INTERSECT = interSection,"OR"=fit$estimate))})
    
    term_names=character()
    for (pathway.term in 1:length(output)){
      term_names[pathway.term]=pathway[[pathway.term]][1]
    }
  
    results=data.frame(do.call("rbind",output))
    results$fdr=p.adjust(as.numeric(as.character(results$NOM_pval)),method = "BH")
    results=results[order(results$fdr),]
    enrichment_tables[[cluster_name]]=results
  }

  return(enrichment_tables)
}

find_enrichments<-function(resA,n=5,show='both'){
  resA.t = resA$t
  names(resA.t)=resA$gene_name 
  if(show=='up'){
    resA.bp = subset(fgsea(go,resA.t),NES>0)
    resA.mf = subset(fgsea(mf,resA.t),NES>0)
    resA.kegg= subset(fgsea(kegg,resA.t),NES>0)
    resA.tft= subset(fgsea(tft,resA.t),NES>0)
  }else if(show=='down'){
    resA.bp = subset(fgsea(go,resA.t),NES<0)
    resA.mf = subset(fgsea(mf,resA.t),NES<0)
    resA.kegg= subset(fgsea(kegg,resA.t),NES<0)   
    resA.tft= subset(fgsea(tft,resA.t),NES<0)
    
  }else{
    resA.bp = fgsea(go,resA.t)
    resA.mf = fgsea(mf,resA.t)
    resA.kegg= fgsea(kegg,resA.t)   
    resA.tft= subset(fgsea(tft,resA.t))
  }
  
  resA.bp$GO = 'BP'
  resA.mf$GO = 'MF'
  resA.kegg$GO = 'KEGG'
  out=list(enr=list(bp=resA.bp,mf=resA.mf,kegg=resA.kegg,tft=resA.tft))
  
  enr = rbind(
    head(resA.bp[order(resA.bp$padj),],n),
    head(resA.mf[order(resA.mf$padj),],n),
    head(resA.kegg[order(resA.kegg$padj),],n)
  )
  enr=data.frame(enr)
  enr$pathway=gsub('_',' ',gsub('GO_|KEGG_','',enr$pathway))
  enr$pathway=factor(enr$pathway,levels=as.character(enr[order(enr$GO,enr$NES),'pathway']))
  p=ggplot(enr,aes(pathway,NES,fill=GO))+geom_col(width=.6,alpha=.6)+geom_col(width=.6,fill=NA,aes(col=GO))+coord_flip()+theme_bw()+geom_vline(xintercept=c(n,n*2)+.5,linetype='dashed')+geom_hline(yintercept=0)
  out$plot=p
  return(out)
}

##top pathway
topPath<- function(fish ,specific, p ,n){
  l=n*20
  result <- matrix(NA, l, l)
  
  colnames(result) <- rownames(result) <- rownames(getElement(fish, specific))[1:l]
  for(i in 1:l){
    for(j in 1:l){
      minlist = 0
      if(length(unlist(strsplit(fish[[specific]][["INTERSECT"]][i], split = ",")))>=length(unlist(strsplit(fish[[specific]][["INTERSECT"]][j], split = ",")))) minlist=length(unlist(strsplit(fish[[specific]][["INTERSECT"]][j], split = ","))) 
      else minlist=length(unlist(strsplit(fish[[specific]][["INTERSECT"]][i], split = ",")))
      result[i, j] <- length(intersect(unlist(strsplit(fish[[specific]][["INTERSECT"]][i], split = ",")), unlist(strsplit(fish[[specific]][["INTERSECT"]][j], split = ","))))/(minlist)      
      result[j, i] <- result[i, j]     
    }
  }
  listBPspecific = list(c(rownames(result)[1]))
  i=1
  while(i < l){
    j=i
    while(j < l){
      if(result[i, j] < p){ 
        i=j
        listBPspecific <- c(listBPspecific, rownames(result)[j])
        j=l
      }else {
        j=j+1
      }
    }
    if(length(listBPspecific)>n-1) 
      i=l
  }
  return(listBPspecific)
}
```

# ALL batchs

``` r
d0 <- DGEList(dat[-1])
d0 <- calcNormFactors(d0)
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

mm <- model.matrix(~ 0 + lesion_type+batch+sample+sex, data=meta[colnames(d),])

# Without filtering low expressed genes
y <- voom(d0[,rownames(meta)], mm, plot = T)
```

    ## Coefficients not estimable: sampleAB257 sexM

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# Filtering low expressed genes
y <- voom(d[,rownames(meta)], mm, plot = T)
```

    ## Coefficients not estimable: sampleAB257 sexM

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
##gene name reference(ref_G)
conv_B=AnnotationDbi::select(org.Hs.eg.db,keys=rownames(y$E),keytype='ENSEMBL',columns=c('GENENAME','SYMBOL'))
ref_G=conv_B[!duplicated(conv_B[,1]),]
rownames(ref_G)=ref_G[,1]
##
```

## PCA All DATA

Perform a *PCA* on the normalized samples and identify main factors of
variation.

### Data including batch effect

``` r
BD<-pca(y$E, metadata = meta)
BD <- biplot(BD, colby = 'batch',hline = 0, vline = 0,legendPosition = 'right', lab="",title = "Data including batch effect")
BD
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
BDCorr <- removeBatchEffect(y$E,batch = meta$batch)
```

### Batch effect removed from data

``` r
BDCorr <- removeBatchEffect(y$E,batch = meta$batch)
BD1 <- pca(BDCorr, metadata = meta)
BDb <- biplot(BD1, colby = 'batch',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Batch effect removed from data")
BDb
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Effect of samples on data (sex, type of ms,…) (col by Sample)

``` r
BDsample <- biplot(BD1, colby = 'sample',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Effect of samples on data (sex, type of ms,...)")
BDsample
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Effect of samples on data (sex, type of ms,…) (col by Sex)

``` r
BDsex <- biplot(BD1, colby = 'sex',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Effect of samples on data (sex, type of ms,...)")
BDsex
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Batch+Patient effect removed from data (col by Sex)

``` r
BDCorrBS <- removeBatchEffect(y$E,batch = meta$batch, batch2 = meta$sample)
BD2 <- pca(BDCorrBS, metadata = meta)
BDsexCorr <- biplot(BD2, colby = 'sex',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Batch+Patient effect removed from data")
BDsexCorr
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Batch+Patient effect removed from data (col by Sample)

``` r
BDsampleCorr <- biplot(BD2, colby = 'sample',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Batch+Patient effect removed from data")
BDsampleCorr
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
pc=prcomp(t(BDCorrBS))
autoplot(pc)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
meta$pc1=pc$x[rownames(meta),1]
meta$pc2=pc$x[rownames(meta),2]
ggplot(meta,aes(pc1,pc2))+
  geom_point(aes(col=tissue))+
  geom_point(data=subset(meta,grepl("^Active",typeL)),col="black",shape="*",size=3)+
  geom_point(data=subset(meta,grepl("^Early",typeL)),col="red",shape="*",size=3)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

### Batch+Patient effect removed from data (col by Lesion type and Sex)

``` r
BDg <- biplot(BD2, colby = 'lesion_type',shape = 'sex',hline = 0, vline = 0,legendPosition = 'right', lab="",title="PCA")
BDg
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Clustering

``` r
hc=hclust(dist(t(BDCorrBS)), method='ward.D2')
hc_clus=data.frame(cutree(hc, k=2))
hc_clus[,1]=paste0('C',hc_clus[,1])
meta$clus_hm=hc_clus[rownames(meta),1]
t=data.frame(table(meta$clus_hm))

plot(hc)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggplot(t,aes(Var1,Freq))+geom_col(fill=c('darkolivegreen3','brown3'))+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
autoplot(prcomp(t(BDCorrBS)),data=meta, col='clus_hm')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

``` r
summarize_col<-function(col){
  t=table(meta[,c('clus_hm' ,col)])
  for(i in 1:ncol(t)){t[,i]=t[,i]/sum(t[,i])}
  for(i in 1:nrow(t)){t[i, ]=t[i,]/sum(t[i,])}
  t=data.frame(t)
  p=ggplot(t,aes(clus_hm,Freq,fill=t[,col]))+geom_col(col='white')+labs(fill=col)+theme_bw()
return(p)
}

p1=summarize_col('sample')
p2=summarize_col('sex')+scale_fill_manual(values=c('deeppink','deepskyblue'))
p3=summarize_col('disease')
p4=summarize_col('batch')+scale_fill_manual(values=c('grey30','grey60','grey90'))
p5=summarize_col('condition')+scale_fill_manual(values=c('brown1','darkolivegreen1','darkorange'))
p6=summarize_col('lesion_type')
p7=summarize_col('tissue')

cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,ncol=3)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Comparisons between cluster

``` r
mmG <- model.matrix(~ 0 + clus_hm+batch+sample+sex, data=meta[colnames(d),])

# Filtering low expressed genes
yG <- voom(d[,rownames(meta)], mmG, plot = F)
```

    ## Coefficients not estimable: sexM

``` r
fitG <- lmFit(yG, mmG)
```

    ## Coefficients not estimable: sexM

``` r
contrG <- makeContrasts(
  A = clus_hmC1 - clus_hmC2,
  levels = colnames(coef(fitG)))
tmpG <- contrasts.fit(fitG, contrG)
tmpG <- eBayes(tmpG) 

res <- topTable(tmpG, sort.by = "P", n = Inf,coef='A')
res$gene_id =  rownames(res)
res$comp = 'clus_hmC1_vs_clus_hmC2'
res$gene_name = ref_G[rownames(res),3]
res=res[order(-abs(res$logFC)),]
res$show = F
res[1:15,'show']=T
res$sign = 'NS'
res[which(res$adj.P.Val<0.05&res$logFC>1),'sign']='UP'
res[which(res$adj.P.Val<0.05&res$logFC<(-1)),'sign']='DN'
res$gene_name=ref_G[rownames(res),3]

ggplot(subset(res,-log10(adj.P.Val)<6),aes(logFC,-log10(adj.P.Val),col=sign)) +geom_point_interactive(aes(x=logFC,y=-log10(adj.P.Val),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(res,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('clus_hmC1_vs_clus_hmC2')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Gene set enrichment analysis(PCA)

``` r
pca = prcomp(t(BDCorrBS))$x[,1:2]
pvals=c()
est=c()
for(i in 1:nrow(BDCorrBS)){
  pvals = c(pvals,cor.test(BDCorrBS[i,],pca[,1])$p.value)
  est = c(est,cor.test(BDCorrBS[i,],pca[,1])$estimate)
}
names(pvals)=rownames (BDCorrBS)
names(pvals)=ref_G[names(pvals),3]
pvals=pvals[!is.na(names(pvals))]

names(est)=rownames(BDCorrBS)
names(est)=ref_G[names(est),3]
est=est[!is.na(names(est))]

df = data.frame(cbind(pvals,est[names(pvals)]))
df$gene = rownames(df)
pvals=c()
est=c()
for(i in 1:nrow(BDCorrBS)){
  pvals = c(pvals,cor.test(BDCorrBS[i,],pca[,2])$p.value)
  est = c(est,cor.test(BDCorrBS[i, ],pca[,2])$estimate)
}
names(pvals)=rownames(BDCorrBS)
names(pvals)=ref_G[names(pvals),3]
pvals=pvals[!is.na(names(pvals))]

names(est)=rownames(BDCorrBS)
names(est)=ref_G[names(est),3]
est=est[!is.na(names(est))]

df2=data.frame(cbind(pvals,est[names(pvals)]))

df$pval_pc2=df2[rownames(df),1]
df$est_pc2=df2[rownames(df),2]

vec1=df[,2]
vec2=df[,4]

names(vec1)=rownames(df)
names(vec2)=rownames(df)

gsea_res1=fgsea(pathways=go,stats=vec1)
gsea_res2=fgsea(pathways=go,stats=vec2)

head(gsea_res1[order(gsea_res1$padj),],30)
```

    ##                                                                            pathway
    ##  1:                                                           GOBP_IMMUNE_RESPONSE
    ##  2:                                                        GOBP_SYNAPTIC_SIGNALING
    ##  3:                                    GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING
    ##  4:                                                       GOBP_CELL_CELL_SIGNALING
    ##  5: GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERSPECIES_INTERACTION_BETWEEN_ORGANISMS
    ##  6:                                                          GOBP_DEFENSE_RESPONSE
    ##  7:                                        GOBP_DEFENSE_RESPONSE_TO_OTHER_ORGANISM
    ##  8:                                             GOBP_REGULATION_OF_IMMUNE_RESPONSE
    ##  9:                                       GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS
    ## 10:                                                    GOBP_INNATE_IMMUNE_RESPONSE
    ## 11:                              GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS
    ## 12:                                               GOBP_ION_TRANSMEMBRANE_TRANSPORT
    ## 13:                                     GOBP_INORGANIC_ION_TRANSMEMBRANE_TRANSPORT
    ## 14:                                                           GOBP_CELL_ACTIVATION
    ## 15:                                            GOBP_CATION_TRANSMEMBRANE_TRANSPORT
    ## 16:                                                   GOBP_TRANSMEMBRANE_TRANSPORT
    ## 17:                                                   GOBP_IMMUNE_EFFECTOR_PROCESS
    ## 18:                                                  GOBP_ADAPTIVE_IMMUNE_RESPONSE
    ## 19:                                          GOBP_REGULATION_OF_MEMBRANE_POTENTIAL
    ## 20:                                                     GOBP_LYMPHOCYTE_ACTIVATION
    ## 21:                                                          GOBP_CATION_TRANSPORT
    ## 22:                                                     GOBP_INFLAMMATORY_RESPONSE
    ## 23:                                    GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE
    ## 24:                                                    GOBP_NERVOUS_SYSTEM_PROCESS
    ## 25:                                                             GOBP_ION_TRANSPORT
    ## 26:                                                        GOBP_NEURON_DEVELOPMENT
    ## 27:                                                       GOBP_CYTOKINE_PRODUCTION
    ## 28:                              GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY
    ## 29:                                                                  GOBP_BEHAVIOR
    ## 30:                                                     GOBP_RESPONSE_TO_BACTERIUM
    ##                                                                            pathway
    ##             pval         padj  log2err         ES       NES size
    ##  1: 1.000000e-50 3.824500e-47       NA -0.3673220 -2.941422 1082
    ##  2: 1.000000e-50 3.824500e-47       NA  0.5289155  3.428126  612
    ##  3: 6.441139e-49 1.642276e-45 1.825919  0.5380582  3.338252  367
    ##  4: 5.698934e-48 1.089779e-44 1.803094  0.3726329  2.506656 1249
    ##  5: 1.259889e-43 1.927378e-40 1.720824 -0.3278873 -2.627554 1055
    ##  6: 1.355556e-42 1.728108e-39 1.702678 -0.3195581 -2.558941 1082
    ##  7: 2.965783e-40 3.240753e-37 1.653315 -0.3703934 -2.923288  673
    ##  8: 2.021935e-38 1.933223e-35 1.615302 -0.3772780 -2.953648  611
    ##  9: 5.778670e-38 4.911228e-35 1.602431 -0.3137426 -2.533956 1007
    ## 10: 1.051382e-35 8.042023e-33 1.556544 -0.3872205 -2.965707  529
    ## 11: 3.631323e-34 2.525090e-31 1.522922 -0.3605529 -2.822709  611
    ## 12: 1.456722e-33 9.285388e-31 1.509263  0.3735897  2.474652  849
    ## 13: 3.357742e-33 1.975643e-30 1.495479  0.4047524  2.633272  643
    ## 14: 3.804524e-32 2.078629e-29 1.474563 -0.3257178 -2.627535  762
    ## 15: 7.791261e-32 3.973024e-29 1.467524  0.3937101  2.571578  678
    ## 16: 9.826653e-32 4.697754e-29 1.467524  0.3400290  2.278790 1132
    ## 17: 7.720329e-30 3.473694e-27 1.417276 -0.3977285 -3.015989  411
    ## 18: 8.183022e-29 3.477330e-26 1.395187 -0.4510710 -3.239243  291
    ## 19: 1.240910e-28 4.995643e-26 1.387746  0.4646893  2.866429  348
    ## 20: 1.003666e-27 3.838520e-25 1.365180 -0.3495934 -2.689599  536
    ## 21: 1.919838e-27 6.992782e-25 1.357574  0.3465337  2.302087  898
    ## 22: 2.079479e-26 7.229970e-24 1.334497 -0.3417271 -2.622647  532
    ## 23: 3.000329e-26 9.978050e-24 1.326716 -0.3963151 -2.942017  368
    ## 24: 4.328352e-26 1.379482e-23 1.326716  0.3577825  2.353565  763
    ## 25: 7.155987e-26 2.189446e-23 1.318889  0.3214756  2.156398 1158
    ## 26: 2.184421e-25 6.426399e-23 1.303093  0.3365603  2.245914  950
    ## 27: 3.835010e-25 1.086444e-22 1.303093 -0.3375009 -2.598704  544
    ## 28: 7.747654e-25 2.116493e-22 1.295123 -0.4001712 -2.923401  333
    ## 29: 9.407937e-25 2.481425e-22 1.287104  0.4052365  2.558314  469
    ## 30: 1.118099e-24 2.850781e-22 1.287104 -0.3712641 -2.791966  407
    ##             pval         padj  log2err         ES       NES size
    ##                                      leadingEdge
    ##  1:          GSN,NFATC2,TRIM56,CSF1,ELF1,SPN,...
    ##  2:  STXBP1,KCNB1,GABRB3,DLGAP1,CLSTN3,CPEB3,...
    ##  3:   STXBP1,KCNB1,DLGAP1,CLSTN3,CPEB3,PRKCE,...
    ##  4:  STXBP1,KCNB1,GABRB3,DLGAP1,CLSTN3,CPEB3,...
    ##  5:          GSN,LITAF,TRIM56,CSF1,NFKB1,SP1,...
    ##  6:          GSN,TRIM56,PARP4,CSF1,NFKB1,SPN,...
    ##  7:          GSN,TRIM56,CSF1,SPN,IFI16,CIITA,...
    ##  8:       NFATC2,TRIM56,ELF1,SPN,IFI16,HMOX1,...
    ##  9:        NFATC2,EVI2B,TRIM56,CSF1,CD9,ELF1,...
    ## 10:     GSN,TRIM56,CSF1,IFI16,CIITA,CALCOCO2,...
    ## 11:        NFATC2,EVI2B,TRIM56,CSF1,ELF1,SPN,...
    ## 12:      KCNB1,GABRB3,GNB5,SCN2A,SCN8A,FGF14,...
    ## 13:      KCNB1,GABRB3,GNB5,SCN2A,SCN8A,FGF14,...
    ## 14:           GSN,NFATC2,CSF1,TLN1,CD9,IKZF1,...
    ## 15:       KCNB1,GNB5,SCN2A,SCN8A,FGF14,PRKCE,...
    ## 16:      KCNB1,GABRB3,GNB5,SCN2A,SCN8A,FGF14,...
    ## 17:             SPN,HMOX1,C3,CD84,LAMP1,CD74,...
    ## 18:            SPN,KLHL6,C3,CD84,CD74,UNC13D,...
    ## 19:    KCNB1,GABRB3,SCN2A,SCN8A,GABRB2,PTK2B,...
    ## 20:           GSN,NFATC2,IKZF1,SPN,CBFB,CD44,...
    ## 21:         KCNB1,NSF,GNB5,SCN2A,SCN8A,FGF14,...
    ## 22:       PARP4,CSF1,NFKB1,IFI16,CIITA,HMOX1,...
    ## 23:     NFATC2,TRIM56,ELF1,IFI16,KLHL6,IKBKB,...
    ## 24:   GABRB3,DLGAP1,RBFOX2,SCN2A,SCN8A,CPEB3,...
    ## 25:       STXBP1,KCNB1,GABRB3,NSF,GNB5,SCN2A,...
    ## 26: CHN1,STXBP1,KALRN,ARHGAP44,RBFOX2,CAMK1D,...
    ## 27:       LITAF,NFATC2,TRIM56,NFKB1,ELF1,SPN,...
    ## 28:       NFATC2,ELF1,KLHL6,ALPK1,IKBKB,FYB1,...
    ## 29:         SCN2A,CPEB3,PRKCE,PLCB1,GPI,GAD1,...
    ## 30:           LITAF,NFKB1,SPN,KLHL6,CARD8,C3,...
    ##                                      leadingEdge

``` r
head(gsea_res2[order(gsea_res2$padj),],30)
```

    ##                                                                             pathway
    ##  1:                                                           GOBP_AXONEME_ASSEMBLY
    ##  2:                                      GOBP_GLANDULAR_EPITHELIAL_CELL_DEVELOPMENT
    ##  3:                                         GOBP_TYPE_B_PANCREATIC_CELL_DEVELOPMENT
    ##  4:                             GOBP_AMYLOID_PRECURSOR_PROTEIN_BIOSYNTHETIC_PROCESS
    ##  5:                                                GOBP_DNA_METHYLATION_ON_CYTOSINE
    ##  6:                             GOBP_NEPHRON_TUBULE_EPITHELIAL_CELL_DIFFERENTIATION
    ##  7:                                                          GOBP_PROTEIN_SULFATION
    ##  8:                            GOBP_REGULATION_OF_CARDIAC_MUSCLE_TISSUE_DEVELOPMENT
    ##  9:                                                                  GOBP_SECRETION
    ## 10:                                        GOBP_REGULATION_OF_VASCULAR_PERMEABILITY
    ## 11:                               GOBP_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY
    ## 12: GOBP_MESENCHYMAL_TO_EPITHELIAL_TRANSITION_INVOLVED_IN_METANEPHROS_MORPHOGENESIS
    ## 13:                                                   GOBP_CARBON_DIOXIDE_TRANSPORT
    ## 14:                               GOBP_NEGATIVE_REGULATION_OF_CENTRIOLE_REPLICATION
    ## 15:                                            GOBP_CARBOHYDRATE_MEDIATED_SIGNALING
    ## 16:                               GOBP_PROTEIN_LOCALIZATION_TO_EXTRACELLULAR_REGION
    ## 17:                            GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION
    ## 18:                                            GOBP_REGULATION_OF_PROTEIN_SECRETION
    ## 19:                                                    GOBP_REGULATION_OF_SECRETION
    ## 20:                                              GOBP_ACTIVATION_OF_IMMUNE_RESPONSE
    ## 21:                                               GOBP_ADHESION_OF_SYMBIONT_TO_HOST
    ## 22:                                                         GOBP_AGGRESOME_ASSEMBLY
    ## 23:                                                            GOBP_AMIDE_TRANSPORT
    ## 24:                                                 GOBP_AMINE_BIOSYNTHETIC_PROCESS
    ## 25:                                           GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY
    ## 26:                                             GOBP_BRAIN_RENIN_ANGIOTENSIN_SYSTEM
    ## 27:                                           GOBP_CARDIOLIPIN_BIOSYNTHETIC_PROCESS
    ## 28:                                                           GOBP_CELL_AGGREGATION
    ## 29:                   GOBP_CELL_DIFFERENTIATION_INVOLVED_IN_METANEPHROS_DEVELOPMENT
    ## 30:                                 GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS
    ##                                                                             pathway
    ##             pval        padj   log2err         ES       NES size
    ##  1: 1.463132e-06 0.008088763 0.6435518  0.7832505  1.270402   68
    ##  2: 2.115262e-06 0.008088763 0.6272567  0.9207636  1.439050   25
    ##  3: 2.527220e-05 0.064427271 0.5756103  0.9214855  1.415747   19
    ##  4: 8.186655e-05 0.156528847 0.5384341  0.9447592  1.422740   14
    ##  5: 4.376469e-04 0.443285231 0.4984931 -0.9607868 -1.592956    3
    ##  6: 5.216484e-04 0.443285231 0.4772708  0.9646596  1.424486   10
    ##  7: 4.131947e-04 0.443285231 0.4984931  0.9874814  1.414605    7
    ##  8: 3.039188e-04 0.443285231 0.4984931  0.9983698  1.365602    3
    ##  9: 5.083207e-04 0.443285231 0.4772708  0.6323508  1.060588  716
    ## 10: 5.815903e-04 0.444800267 0.4772708  0.7816942  1.245243   39
    ## 11: 6.555150e-04 0.455761672 0.4772708  0.9041625  1.361605   14
    ## 12: 7.663241e-04 0.488403888 0.4772708  0.9731266  1.405826    8
    ## 13: 9.079975e-04 0.496026037 0.4772708  0.9959366  1.362274    3
    ## 14: 8.586214e-04 0.496026037 0.4772708  0.9810650  1.405413    7
    ## 15: 1.148167e-03 0.585412130 0.4550599  0.9677172  1.398011    8
    ## 16: 1.647531e-03 0.707632721 0.4550599  0.6491163  1.080996  272
    ## 17: 1.721163e-03 0.707632721 0.4550599  0.7532518  1.204842   43
    ## 18: 1.757979e-03 0.707632721 0.4550599  0.6643081  1.101583  187
    ## 19: 1.610715e-03 0.707632721 0.4550599  0.6347526  1.061791  468
    ## 20: 9.938984e-03 0.710924602 0.3807304  0.6451826  1.072694  241
    ## 21: 7.421221e-03 0.710924602 0.4070179  0.8888837  1.322323   11
    ## 22: 1.065246e-02 0.710924602 0.3807304  0.9617611  1.349136    5
    ## 23: 1.023790e-02 0.710924602 0.3807304  0.6480384  1.075441  210
    ## 24: 4.945072e-03 0.710924602 0.4070179  0.7955733  1.244981   26
    ## 25: 4.201468e-03 0.710924602 0.4070179  0.7706258  1.213900   30
    ## 26: 9.249524e-03 0.710924602 0.3807304  0.9835537  1.345336    3
    ## 27: 7.428514e-03 0.710924602 0.4070179  0.9270477  1.339258    8
    ## 28: 3.092436e-03 0.710924602 0.4317077  0.8747772  1.322330   15
    ## 29: 2.090411e-03 0.710924602 0.4317077  0.8664976  1.319743   17
    ## 30: 1.098517e-02 0.710924602 0.3807304  0.7419375  1.176073   34
    ##             pval        padj   log2err         ES       NES size
    ##                                      leadingEdge
    ##  1:          DNAAF8,BBS2,GAS8,TTLL3,TTLL1,DNAAF2
    ##  2:                            GSK3B,PDPK1,CLOCK
    ##  3:                            GSK3B,PDPK1,CLOCK
    ##  4:              NECAB3,NECAB2,ABCA7,NCSTN,ITM2A
    ##  5:                           EHMT2,DNMT1,DNMT3A
    ##  6:                                  STAT1,MEF2C
    ##  7:                               CHST5,HS3ST3B1
    ##  8:                                   CREB1,JPH2
    ##  9:      NECAB3,PFKFB2,CREB1,SYT12,CD84,AQP1,...
    ## 10:      OCLN,TACR1,CEACAM1,TRPV4,SLIT2,TJP2,...
    ## 11:                  OCLN,TACR1,TRPV4,TJP2,TGFB1
    ## 12:                                   STAT1,SIX2
    ## 13:                                     AQP1,CA2
    ## 14:                                 TRIM37,RBM14
    ## 15:                                CLEC7A,MLXIPL
    ## 16:   NECAB3,PFKFB2,C2CD2L,CAVIN1,UCP2,JAGN1,...
    ## 17: COLGALT1,RUNX1,TGFB1,CFLAR,TNFRSF1B,HAS2,...
    ## 18:      PFKFB2,C2CD2L,UCP2,JAGN1,CD33,UQCC2,...
    ## 19:        PFKFB2,CREB1,SYT12,CD84,AQP1,SYT2,...
    ## 20:              C2,C3AR1,C5,IFI16,PTPRC,BTK,...
    ## 21:                                HSP90AB1,CD81
    ## 22:                                  TRIM37,PRKN
    ## 23:  PFKFB2,SLC16A10,C2CD2L,UCP2,JAGN1,UQCC2,...
    ## 24:                       NR4A2,KL,SAT2,OAZ1,SRM
    ## 25:                    DNAAF8,DNAAF2,TTC12,DNAI4
    ## 26:                                    TACR1,AGT
    ## 27:                                  PLSCR3,PGS1
    ## 28:                                COL11A1,WNT7A
    ## 29:                                   STAT1,WNT4
    ## 30:   EPHB6,NR4A2,NIN,MYCBP2,HSP90AB1,PLXNA4,...
    ##                                      leadingEdge

## Lesion only

``` r
lesions=subset(meta,condition=='Lesion')
autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,col='tissue')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
lesions$cd3e = BDCorrBS['ENSG00000198851', rownames(lesions)]
lesions$vcam1 = BDCorrBS['ENSG00000162692', rownames(lesions)]
lesions$icam1 = BDCorrBS['ENSG00000090339', rownames(lesions)]
lesions$mcam = BDCorrBS['ENSG00000076706', rownames(lesions)] 
lesions$ighg2 = BDCorrBS['ENSG00000211893', rownames(lesions)]
lesions$ms4a1 = BDCorrBS['ENSG00000156738', rownames(lesions)]
lesions$itga8 = BDCorrBS['ENSG00000077943', rownames(lesions)]
lesions$irf1 = BDCorrBS['ENSG00000125347', rownames(lesions)] 
lesions$cd24 = BDCorrBS['ENSG00000272398', rownames(lesions)] 
lesions$ccl4 = BDCorrBS['ENSG00000275302', rownames(lesions)] 

a1=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='vcam1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a2=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='icam1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a3=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='mcam',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a4=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='cd3e',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a5=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='ighg2',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a6=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='ms4a1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a7=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='itga8',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a8=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='irf1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
#a9=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='cd24',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a9=autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=lesions,fill='ccl4',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()

cowplot::plot_grid(a1,a2,a3,a4,a5,a6,a7,a8,a9,ncol=3)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Clustering only Lesion

``` r
hc=hclust(dist(t(BDCorrBS[,rownames(lesions)])))
subtypeLesion=data.frame(cutree(hc, k=2))
subtypeLesion[,1]=paste0('C',subtypeLesion[,1])
colnames(subtypeLesion)='subtypeLesion'

meta$clus_hm=subtypeLesion[rownames(meta),1]
t=data.frame(table(meta$clus_hm))

plot(hc)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggplot(t,aes(Var1,Freq))+geom_col(fill=c('darkolivegreen3','brown3'))+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
autoplot(prcomp(t(BDCorrBS[,rownames(lesions)])),data=subtypeLesion,col='subtypeLesion')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

``` r
summarize_col<-function(col){
  t=table(lesions[,c('lesion_type' ,col)])
  for(i in 1:ncol(t)){t[,i]=t[,i]/sum(t[,i])}
  for(i in 1:nrow(t)){t[i, ]=t[i,]/sum(t[i,])}
  t=data.frame(t)
  p=ggplot(t,aes(lesion_type,Freq,fill=t[,col]))+geom_col(col='white')+labs(fill=col)+theme_bw()
  return(p)
}

lesions$lesion_type=subtypeLesion[rownames(lesions),1]

p1=summarize_col('sample')
p2=summarize_col('sex')+scale_fill_manual(values=c('deeppink','deepskyblue'))
p3=summarize_col('batch')+scale_fill_manual(values=c('grey30','grey60','grey90'))
p4=summarize_col('tissue')

cowplot::plot_grid(p1,p2,p3,p4,ncol=2)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
#subset(lesions,lesion_type=='C2'&tissue=='PV')
#subset(lesions,lesion_type=='C1'&tissue=='Cortical')
meta[c('S72','S81d'),]
```

    ##        id  batch lesion_type   tissue condition sample typeL sex disease
    ## S72   S72 batch1          CL Cortical    Lesion  AB187         M      MS
    ## S81d S81d batch3        PVML       PV    Lesion  AB203  PVML   F      MS
    ##            group    reads       pc1       pc2 clus_hm
    ## S72  MS.Cortical 12938331 -67.06396 10.947404      C1
    ## S81d       MS.PV 15303502  66.81026 -3.480749      C2

``` r
mmG <- model.matrix(~ 0 + lesion_type+batch, data=lesions)

# Filtering low expressed genes
yG <- voom(d[,rownames(lesions)], mmG, plot = F)

fitG <- lmFit(yG, mmG)
contrG <- makeContrasts(
  A = lesion_typeC1 - lesion_typeC2,
  levels = colnames(coef(fitG)))
tmpG <- contrasts.fit(fitG, contrG)
tmpG <- eBayes(tmpG) 

res <- list()
res[['C1_vs_C2']] <- topTable(tmpG, sort.by = "P", n = Inf,coef='A')

res=do.call(rbind,res)
res=res[order(-abs(res$logFC)),]
res$show = F
res[1:15,'show']=T
res$sign = 'NS'
res[which(res$adj.P.Val<0.01&res$logFC>0.5),'sign']='UP'
res[which(res$adj.P.Val<0.01&res$logFC<(-0.5)),'sign']='DN'
res$gene_name=ref_G[gsub(".+\\.","",rownames(res)),3]
res$comp="C1_vs_C2"

ggplot(res,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(res,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+facet_grid(~comp)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
meta$subtype = 'NoLesion'
meta[rownames(lesions),]$subtype = subtypeLesion[,1]
meta[grepl('^NA' ,meta$lesion_type), 'subtype']="NormalAppear"

autoplot(prcomp(t(BDCorrBS)), data = meta, col ='subtype' ,shape='tissue')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
autoplot(prcomp(t(BDCorr[,rownames(subset(meta,tissue%in%c('WM')))]),scale=T),data=subset(meta, tissue%in%c('WM')),col='lesion_type')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

``` r
autoplot(prcomp(t(BDCorr[,rownames(subset(meta,tissue%in%c('Cortical')))])),data=subset(meta,tissue%in%c('Cortical')),col='lesion_type' )+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-22-3.png)<!-- -->

``` r
autoplot(prcomp(t(BDCorr[,rownames(subset(meta,tissue%in%c('PV')))])),data=subset(meta,tissue%in%c('PV')),col='lesion_type')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-22-4.png)<!-- -->

``` r
mmG <- model.matrix(~ 0 + subtype+batch, data=meta[colnames(d),])

# Filtering low expressed genes
yG <- voom(d[,rownames(meta)], mmG, plot = F)

fitG <- lmFit(yG, mmG)
contrG <- makeContrasts(
  A = subtypeC1 - subtypeC2,
  B = subtypeC1 - (subtypeNoLesion+subtypeNormalAppear)/2,
  C = subtypeC2 - (subtypeNoLesion+subtypeNormalAppear)/2,
  levels = colnames(coef(fitG)))
tmpG <- contrasts.fit(fitG, contrG)
tmpG <- eBayes(tmpG) 

res <- list()
res[['C1_vs_C2']] <- topTable(tmpG, sort.by = "P", n = Inf,coef='A')
res[['C1_vs_NoLesion+NormalAppear']] <- topTable(tmpG, sort.by = "P", n = Inf,coef='B')
res[['C2_vs_NoLesion+NormalAppear']] <- topTable(tmpG, sort.by = "P", n = Inf,coef='C')

res[['C1_vs_C2']]$comp='C1_vs_C2'
res[['C1_vs_NoLesion+NormalAppear']]$comp='C1_vs_NoLesion+NormalAppear'
res[['C2_vs_NoLesion+NormalAppear']]$comp='C2_vs_NoLesion+NormalAppear'
res=do.call(rbind,res)
res$sign = 'NS'
res[which(res$adj.P.Val<0.01&res$logFC>0.5),'sign']='UP'
res[which(res$adj.P.Val<0.01&res$logFC<(-0.5)),'sign']='DN'
res$gene_name=ref_G[gsub(".+\\.","",rownames(res)),3]

ggplot(res,aes(logFC,-log10(P.Value),col=sign))+geom_point()+facet_grid(~comp)+theme_bw()+scale_color_manual(values=sign_cols)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
plot_gene_tissue<-function(gene) {
  id = rownames(ref_G[which(ref_G[,3]==gene),])
  datap = meta 
  meta$gene = BDCorrBS[id, rownames(meta)]
  ggplot(meta,aes(condition,gene))+geom_boxplot(outlier.shape=NA)+facet_grid(~tissue,space='free_x',scales='free_x')+theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))+geom_point(aes(col=subtype),position=position_jitter(width=.1))+ggtitle(gene)
}
plot_gene_tissue('IGHG2')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
plot_gene_tissue('GBP1')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

``` r
plot_gene_tissue('VCAM1')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->

``` r
plot_gene_tissue('MCAM')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-4.png)<!-- -->

``` r
plot_gene_tissue('CD3E')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-5.png)<!-- -->

``` r
plot_gene_tissue('VCAM1')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-6.png)<!-- -->

``` r
plot_gene_tissue('PLP1')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-7.png)<!-- -->

``` r
plot_gene_tissue('CNP')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-8.png)<!-- -->

``` r
plot_gene_tissue('IRF1')
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-24-9.png)<!-- -->

``` r
noLesions=subset(meta,condition=='Normal'|condition=='Normal_Appearing')
autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,col='condition')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
noLesions$cd3e = BDCorrBS['ENSG00000198851', rownames(noLesions)]
noLesions$vcam1 = BDCorrBS['ENSG00000162692', rownames(noLesions)]
noLesions$icam1 = BDCorrBS['ENSG00000090339', rownames(noLesions)]
noLesions$mcam = BDCorrBS['ENSG00000076706', rownames(noLesions)] 
noLesions$ighg2 = BDCorrBS['ENSG00000211893', rownames(noLesions)]
noLesions$ms4a1 = BDCorrBS['ENSG00000156738', rownames(noLesions)]
noLesions$itga8 = BDCorrBS['ENSG00000077943', rownames(noLesions)]
noLesions$irf1 = BDCorrBS['ENSG00000125347', rownames(noLesions)] 
noLesions$cd24 = BDCorrBS['ENSG00000272398', rownames(noLesions)] 
noLesions$ccl4 = BDCorrBS['ENSG00000275302', rownames(noLesions)] 

a1=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='vcam1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a2=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='icam1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a3=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='mcam',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a4=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='cd3e',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a5=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='ighg2',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a6=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='ms4a1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a7=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='itga8',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a8=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='irf1',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
#a9=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='cd24',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()
a9=autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=noLesions,fill='ccl4',shape=21)+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+theme_bw()

cowplot::plot_grid(a1,a2,a3,a4,a5,a6,a7,a8,a9,ncol=3)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
hc=hclust(dist(t(BDCorrBS[,rownames(noLesions)])))
subtypeLesion=data.frame(cutree(hc, k=2))
subtypeLesion[,1]=paste0('C',subtypeLesion[,1])
colnames(subtypeLesion)='subtypeLesion'

plot(hc)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
autoplot(prcomp(t(BDCorrBS[,rownames(noLesions)])),data=subtypeLesion,col='subtypeLesion')+theme_bw()
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

``` r
summarize_col<-function(col){
  t=table(noLesions[,c('sample_type' ,col)])
  for(i in 1:ncol(t)){t[,i]=t[,i]/sum(t[,i])}
  for(i in 1:nrow(t)){t[i,]=t[i,]/sum(t[i,])}
  t=data.frame(t)
  p=ggplot(t,aes(sample_type,Freq,fill=t[,col]))+geom_col(col='white')+labs(fill=col)+theme_bw()
  return(p)

}

noLesions$sample_type=subtypeLesion[rownames(noLesions),1]
noLesions$lesion_type=factor(noLesions$lesion_type,levels = c('NAGM','NAWM','GM','PVM','WM'))

p1=summarize_col('sample')
p2=summarize_col('sex')+scale_fill_manual(values=c('deeppink','deepskyblue'))
p3=summarize_col('batch')+scale_fill_manual(values=c('grey30','grey60','grey90'))
p4=summarize_col('disease')+scale_fill_manual(values=c('brown1','darkolivegreen1','darkorange'))
p5=summarize_col('lesion_type')
p6=summarize_col('tissue')

cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
mmG <- model.matrix(~ 0 + sample_type+batch, data=noLesions)

# Filtering low expressed genes
yG <- voom(d[,rownames(noLesions)], mmG, plot = F)

fitG <- lmFit(yG, mmG)
contrG <- makeContrasts(
  A = sample_typeC1 - sample_typeC2,
  levels = colnames(coef(fitG)))
tmpG <- contrasts.fit(fitG, contrG)
tmpG <- eBayes(tmpG) 

resN <- list()
resN[['C1_vs_C2']] <- topTable(tmpG, sort.by = "P", n = Inf,coef='A')

resN[['C1_vs_C2']]$comp='C1_vs_C2'
resN=do.call(rbind,resN)
resN$sign = 'NS'
resN[which(resN$adj.P.Val<0.01&resN$logFC>0.5),'sign']='UP'
resN[which(resN$adj.P.Val<0.01&resN$logFC<(-0.5)),'sign']='DN'
resN$gene_name=ref_G[gsub(".+\\.","",rownames(resN)),3]

ggplot(resN,aes(logFC,-log10(P.Value),col=sign))+geom_point()+facet_grid(~comp)+theme_bw()+scale_color_manual(values=sign_cols)
```

![](BulkAnalysisUnsupervisedMSLesions_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->
