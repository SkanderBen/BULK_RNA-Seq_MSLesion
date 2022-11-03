Bulk analysis of MS lesions(Supervised)
================

This is an analysis performed on the bulk RNA-Seq dataset

The workflow is as follows:  
1. Inspect variance in the dataset (Filtering, Normalization, PCA)  
2. Perform differential expression analysis (Using limma)  
3. Annotate biological signatures found by DEA (GSEA)

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

knitr::kable(meta, caption ="A Knitr kable")
```

|          | id       | batch  | lesion_type | tissue   | condition | sample | typeL        | sex | disease | group       |    reads |
|:---------|:---------|:-------|:------------|:---------|:----------|:-------|:-------------|:----|:--------|:------------|---------:|
| S13_b1   | S13_b1   | batch1 | WML         | WM       | Lesion    | AB187  |              | M   | MS      | MS.WM       | 10675040 |
| S14_b1   | S14_b1   | batch1 | NAGM        | Cortical | Normal    | AB187  |              | M   | MS      | MS.Cortical | 11601256 |
| S23      | S23      | batch1 | CL          | Cortical | Lesion    | AB203  |              | F   | MS      | MS.Cortical | 10388672 |
| S26      | S26      | batch1 | WML         | WM       | Lesion    | AB187  |              | M   | MS      | MS.WM       | 11312959 |
| S29      | S29      | batch1 | CL          | Cortical | Lesion    | AB187  |              | M   | MS      | MS.Cortical | 10849903 |
| S33      | S33      | batch1 | CL          | Cortical | Lesion    | AB200  |              | M   | MS      | MS.Cortical | 15027141 |
| S37      | S37      | batch1 | WML         | WM       | Lesion    | AB200  |              | M   | MS      | MS.WM       | 12031415 |
| S38      | S38      | batch1 | NAGM        | Cortical | Normal    | AB200  |              | M   | MS      | MS.Cortical | 12856284 |
| S3_b1    | S3_b1    | batch1 | PVML        | PV       | Lesion    | AB200  |              | M   | MS      | MS.PV       | 10431353 |
| S43      | S43      | batch1 | CL          | Cortical | Lesion    | AB203  |              | F   | MS      | MS.Cortical | 12121099 |
| S45      | S45      | batch1 | PVML        | PV       | Lesion    | AB203  |              | F   | MS      | MS.PV       | 10093676 |
| S53      | S53      | batch1 | NAWM        | WM       | Normal    | AB187  |              | M   | MS      | MS.WM       | 10907244 |
| S55      | S55      | batch1 | NAGM        | Cortical | Normal    | AB200  |              | M   | MS      | MS.Cortical | 20126720 |
| S62      | S62      | batch1 | WML         | WM       | Lesion    | AB200  |              | M   | MS      | MS.WM       | 12859996 |
| S64      | S64      | batch1 | NAGM        | Cortical | Normal    | AB203  |              | F   | MS      | MS.Cortical | 13502425 |
| S67      | S67      | batch1 | WML         | WM       | Lesion    | AB203  |              | F   | MS      | MS.WM       | 11484555 |
| S70      | S70      | batch1 | WML         | WM       | Lesion    | AB203  |              | F   | MS      | MS.WM       | 10809418 |
| S71      | S71      | batch1 | NAWM        | WM       | Normal    | AB187  |              | M   | MS      | MS.WM       | 13315868 |
| S72      | S72      | batch1 | CL          | Cortical | Lesion    | AB187  |              | M   | MS      | MS.Cortical | 12938331 |
| S76      | S76      | batch1 | CL          | Cortical | Lesion    | AB200  |              | M   | MS      | MS.Cortical | 14486529 |
| S77      | S77      | batch1 | WML         | WM       | Lesion    | AB200  |              | M   | MS      | MS.WM       | 14583455 |
| S79      | S79      | batch1 | NAWM        | WM       | Normal    | AB203  |              | F   | MS      | MS.WM       | 14931094 |
| S9_b1    | S9_b1    | batch1 | NAGM        | Cortical | Normal    | AB187  |              | M   | MS      | MS.Cortical | 14719787 |
| Sample1  | Sample1  | batch2 | PVML        | PV       | Lesion    | AB200  |              | M   | MS      | MS.PV       | 20735000 |
| Sample10 | Sample10 | batch2 | CL          | Cortical | Lesion    | AB187  |              | M   | MS      | MS.Cortical |  9925960 |
| Sample11 | Sample11 | batch2 | WML         | WM       | Lesion    | AB187  |              | M   | MS      | MS.WM       | 19113957 |
| Sample12 | Sample12 | batch2 | CL          | Cortical | Lesion    | AB187  |              | M   | MS      | MS.Cortical | 13267587 |
| Sample15 | Sample15 | batch2 | WML         | WM       | Lesion    | AB187  |              | M   | MS      | MS.WM       | 14887133 |
| Sample16 | Sample16 | batch2 | WML         | WM       | Lesion    | AB187  |              | M   | MS      | MS.WM       | 12598401 |
| Sample17 | Sample17 | batch2 | PVML        | PV       | Lesion    | AB203  |              | F   | MS      | MS.PV       | 10162274 |
| Sample18 | Sample18 | batch2 | CL          | Cortical | Lesion    | AB203  |              | F   | MS      | MS.Cortical | 11609221 |
| Sample21 | Sample21 | batch2 | NAWM        | WM       | Normal    | AB203  |              | F   | MS      | MS.WM       | 11016594 |
| Sample24 | Sample24 | batch2 | NAWM        | WM       | Normal    | AB203  |              | F   | MS      | MS.WM       | 14238903 |
| Sample25 | Sample25 | batch2 | NAGM        | Cortical | Normal    | AB187  |              | M   | MS      | MS.Cortical | 15168403 |
| Sample34 | Sample34 | batch2 | WML         | WM       | Lesion    | AB200  |              | M   | MS      | MS.WM       | 18291181 |
| Sample46 | Sample46 | batch2 | NAGM        | Cortical | Normal    | AB203  |              | F   | MS      | MS.Cortical | 16625491 |
| Sample52 | Sample52 | batch2 | NAWM        | WM       | Normal    | AB187  |              | M   | MS      | MS.WM       | 12701779 |
| Sample56 | Sample56 | batch2 | PVML        | PV       | Lesion    | AB200  |              | M   | MS      | MS.PV       | 13165303 |
| Sample57 | Sample57 | batch2 | CL          | Cortical | Lesion    | AB200  |              | M   | MS      | MS.Cortical | 11305279 |
| Sample6  | Sample6  | batch2 | NAGM        | Cortical | Normal    | AB200  |              | M   | MS      | MS.Cortical | 22662880 |
| Sample60 | Sample60 | batch2 | NAWM        | WM       | Normal    | AB200  |              | M   | MS      | MS.WM       | 14547748 |
| Sample65 | Sample65 | batch2 | PVML        | PV       | Lesion    | AB203  |              | F   | MS      | MS.PV       | 13327909 |
| Sample66 | Sample66 | batch2 | WML         | WM       | Lesion    | AB203  |              | F   | MS      | MS.WM       |  8626915 |
| Sample7  | Sample7  | batch2 | WML         | WM       | Lesion    | AB200  |              | M   | MS      | MS.WM       | 17559199 |
| Sample74 | Sample74 | batch2 | NAGM        | Cortical | Normal    | AB187  |              | M   | MS      | MS.Cortical | 17633780 |
| Sample75 | Sample75 | batch2 | NAWM        | WM       | Normal    | AB200  |              | M   | MS      | MS.WM       | 20237852 |
| S1       | S1       | batch3 | GM          | Cortical |           | AB257  | GM           | M   | HC      | HC.Cortical | 14746041 |
| S2       | S2       | batch3 | GM          | Cortical |           | AB257  | GM           | M   | HC      | HC.Cortical | 16135512 |
| S3_b3    | S3_b3    | batch3 | WM          | WM       |           | AB257  | WM           | M   | HC      | HC.WM       | 17280964 |
| S4       | S4       | batch3 | WM          | WM       |           | AB257  | WM           | M   | HC      | HC.WM       | 16128150 |
| S5       | S5       | batch3 | PVM         | PV       |           | AB257  | PVM          | M   | HC      | HC.PV       | 15969401 |
| S6       | S6       | batch3 | PVM         | PV       |           | AB257  | PVM          | M   | HC      | HC.PV       | 14119977 |
| S7       | S7       | batch3 | GM          | Cortical |           | AB252  | GM           | M   | HC      | HC.Cortical | 16993959 |
| S8       | S8       | batch3 | GM          | Cortical |           | AB252  | GM           | M   | HC      | HC.Cortical | 16902227 |
| S9_b3    | S9_b3    | batch3 | WM          | WM       |           | AB252  | WM           | M   | HC      | HC.WM       | 14200304 |
| S10      | S10      | batch3 | WM          | WM       |           | AB252  | WM           | M   | HC      | HC.WM       | 14993001 |
| S11      | S11      | batch3 | PVM         | PV       |           | AB252  | PVM          | M   | HC      | HC.PV       | 13322475 |
| S12      | S12      | batch3 | PVM         | PV       |           | AB252  | PVM          | M   | HC      | HC.PV       | 11873861 |
| S13_b3   | S13_b3   | batch3 | GM          | Cortical |           | AB251  | GM           | M   | HC      | HC.Cortical | 14540503 |
| S14_b3   | S14_b3   | batch3 | GM          | Cortical |           | AB251  | GM           | M   | HC      | HC.Cortical | 15068466 |
| S15      | S15      | batch3 | WM          | WM       |           | AB251  | WM           | M   | HC      | HC.WM       | 13301374 |
| S16      | S16      | batch3 | WM          | WM       |           | AB251  | WM           | M   | HC      | HC.WM       | 15610086 |
| S17      | S17      | batch3 | PVM         | PV       |           | AB251  | PVM          | M   | HC      | HC.PV       | 12503687 |
| S18      | S18      | batch3 | PVM         | PV       |           | AB251  | PVM          | M   | HC      | HC.PV       | 15465598 |
| S102f    | S102f    | batch3 | CL          | Cortical |           | AB253  | CL           | F   | MS      | MS.Cortical | 11101420 |
| S104f    | S104f    | batch3 | CL          | Cortical |           | AB253  | CL           | F   | MS      | MS.Cortical | 16624705 |
| S106f    | S106f    | batch3 | CL          | Cortical |           | AB253  | CL           | F   | MS      | MS.Cortical | 18467343 |
| S108f    | S108f    | batch3 | NAGM        | Cortical |           | AB253  | NAGM         | F   | MS      | MS.Cortical | 12319237 |
| S109f    | S109f    | batch3 | NAGM        | Cortical |           | AB253  | NAGM         | F   | MS      | MS.Cortical | 14185427 |
| S110f    | S110f    | batch3 | CL          | Cortical |           | AB253  | CL           | F   | MS      | MS.Cortical | 12999955 |
| S112f    | S112f    | batch3 | NAWM        | WM       |           | AB253  | NAWM         | F   | MS      | MS.WM       | 13363398 |
| S113f    | S113f    | batch3 | NAWM        | WM       |           | AB253  | NAWM         | F   | MS      | MS.WM       | 16157577 |
| S114f    | S114f    | batch3 | CL          | Cortical |           | AB253  | CL           | F   | MS      | MS.Cortical | 16266665 |
| S115f    | S115f    | batch3 | NAGM        | Cortical |           | AB253  | NAGM         | F   | MS      | MS.Cortical | 17885078 |
| S118f    | S118f    | batch3 | PVML        | PV       |           | AB253  | Mixed.PD     | F   | MS      | MS.PV       | 10443427 |
| S119f    | S119f    | batch3 | PVML        | PV       |           | AB253  | Active.PD    | F   | MS      | MS.PV       | 13248535 |
| S120f    | S120f    | batch3 | PVML        | PV       |           | AB253  | Active.PD    | F   | MS      | MS.PV       | 12139042 |
| S121f    | S121f    | batch3 | PVML        | PV       |           | AB253  | Mixed.PD     | F   | MS      | MS.PV       | 10019098 |
| S19a     | S19a     | batch3 | WML         | WM       |           | AB253  | Active       | F   | MS      | MS.WM       | 17274976 |
| S22a     | S22a     | batch3 | PVML        | PV       |           | AB253  | PVML         | F   | MS      | MS.PV       | 11231127 |
| S28b     | S28b     | batch3 | WML         | WM       |           | AB187  | Active.PD    | M   | MS      | MS.WM       | 15599429 |
| S30b     | S30b     | batch3 | NAWM        | WM       |           | AB187  | NAWM         | M   | MS      | MS.WM       | 16356084 |
| S31b     | S31b     | batch3 | WML         | WM       |           | AB187  | Active.PD    | M   | MS      | MS.WM       | 13453907 |
| S51c     | S51c     | batch3 | CL          | Cortical |           | AB187  | CL           | M   | MS      | MS.Cortical | 16757776 |
| S69c     | S69c     | batch3 | PVML        | PV       |           | AB203  | PVML         | F   | MS      | MS.PV       | 13814668 |
| S73d     | S73d     | batch3 | WML         | WM       |           | AB187  | Active.PD    | M   | MS      | MS.WM       | 13243821 |
| S78d     | S78d     | batch3 | NAWM        | WM       |           | AB200  | NAWM         | M   | MS      | MS.WM       | 14215553 |
| S81d     | S81d     | batch3 | PVML        | PV       |           | AB203  | PVML         | F   | MS      | MS.PV       | 15303502 |
| S89f     | S89f     | batch3 | WML         | WM       |           | AB253  | Early.Active | F   | MS      | MS.WM       | 10392586 |
| S90f     | S90f     | batch3 | WML         | WM       |           | AB253  | Early.Active | F   | MS      | MS.WM       | 18177867 |
| S91f     | S91f     | batch3 | WML         | WM       |           | AB253  | Early.Active | F   | MS      | MS.WM       | 11856582 |
| S92f     | S92f     | batch3 | NAGM        | Cortical |           | AB253  | NAGM         | F   | MS      | MS.Cortical | 12412586 |
| S94f     | S94f     | batch3 | WML         | WM       |           | AB253  | Mixed.PD     | F   | MS      | MS.WM       | 13343462 |
| S96f     | S96f     | batch3 | WML         | WM       |           | AB253  | Active.AD    | F   | MS      | MS.WM       |  8964589 |

A Knitr kable

Define functions for the analysis: - find_enrichments(): For a given
comparison (DEA) extract enrichments for GO(BP), GO(MF), KEGG and TFT  
- fisher_enrichment(): Perform an enrichment of genesets using a fisher
test approach - plot.gene(): Plot gene expression for a given gene
across conditions.

The gmt files contain the genesets for *GSEA*

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
    freq.table <- table(factor(universe %in% as.character(cluster_markers[[cluster]]), 
                          levels = c(TRUE,FALSE)), 
                          factor(universe %in% x, 
                          levels = c(TRUE,FALSE)))

      fit <- fisher.test(freq.table, alternative = "greater")
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

### Counts list et metadata for batch3

``` r
datB3 <- dat[,48:95]
metaB3 <- meta[47:94,]
```

## Differential expression analysis with limma-voom

### Filtering to remove lowly expressed genes + Normalization for composition bias

use the *edgeR/limma* workflow to filter out lowly expressed genes.
*Voom* is used to normalize the raw counts. A linear mixed model is used
to fit the expression values. Defining the model is done with the
model.matrix function of edgeR. In the model are included factors of
interest (group). If a batch/patient effect is detected, it can be added
to the model to account for its effect.

## Only Batch3

### Without filtering low expressed genes

``` r
d0 <- DGEList(datB3)
d0 <- calcNormFactors(d0)
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

mmB3 <- model.matrix(~ 0 + lesion_type+sex+sample+reads, data=metaB3[colnames(d),])

y <- voom(d0[,rownames(metaB3)], mmB3, plot = T)
```

    ## Coefficients not estimable: sampleAB253 sampleAB257

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Filtering low expressed genes

``` r
# Filtering low expressed genes
yB3 <- voom(d[,rownames(metaB3)], mmB3, plot = T)
```

    ## Coefficients not estimable: sampleAB253 sampleAB257

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## ALL batchs

### Without filtering low expressed genes

``` r
d0 <- DGEList(dat[-1])
d0 <- calcNormFactors(d0)
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

mm <- model.matrix(~ 0 + lesion_type+batch+sample+sex, data=meta[colnames(d),])

y <- voom(d0[,rownames(meta)], mm, plot = T)
```

    ## Coefficients not estimable: sampleAB257 sexM

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Filtering low expressed genes

``` r
# Filtering low expressed genes
y <- voom(d[,rownames(meta)], mm, plot = T)
```

    ## Coefficients not estimable: sampleAB257 sexM

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Reads By Batch

``` r
meta$id = factor(meta$id, levels=meta[order(meta$sample),'id'])
ggplot(meta,aes(id, reads, fill=sample))+theme_bw()+geom_col(width=0.6)+ theme(axis.text.x=element_text(angle=60,hjust=1,size=5))+facet_grid(~batch, space = "free_x",scales = "free_x")
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## PCA Batch3

Perform a *PCA* on the normalized samples and identify main factors of
variation.

### col by Sample

``` r
pcaB3<-pca(yB3$E, metadata = metaB3)
pcaB3S <- biplot(pcaB3, colby = 'sample',hline = 0, vline = 0,legendPosition = 'right', lab="",title = "")
pcaB3S
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
\### col by lesion type

``` r
pcaB3<-pca(yB3$E, metadata = metaB3)
pcaB3L <- biplot(pcaB3, colby = 'lesion_type',hline = 0, vline = 0,legendPosition = 'right', lab="",title = "")
pcaB3L
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
\### col by tissue

``` r
pcaB3<-pca(yB3$E, metadata = metaB3)
pcaB3t <- biplot(pcaB3, colby = 'tissue',hline = 0, vline = 0,legendPosition = 'right', lab="",title = "")
pcaB3t
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
\### col by type

``` r
pcaB3<-pca(yB3$E, metadata = metaB3)
pcaB3tl <- biplot(pcaB3, colby = 'typeL',shape='tissue',hline = 0, vline = 0,legendPosition = 'right', lab="",title = "")
pcaB3tl
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## PCA All DATA

Perform a *PCA* on the normalized samples and identify main factors of
variation.

### Data including batch effect

``` r
BD<-pca(y$E, metadata = meta)
BD <- biplot(BD, colby = 'batch',hline = 0, vline = 0,legendPosition = 'right', lab="",title = "Data including batch effect")
BD
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

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

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### Effect of samples on data (sex, type of ms,…) (col by Sample)

``` r
BDsample <- biplot(BD1, colby = 'sample',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Effect of samples on data (sex, type of ms,...)")
BDsample
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### Effect of samples on data (sex, type of ms,…) (col by Sex)

``` r
BDsex <- biplot(BD1, colby = 'sex',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Effect of samples on data (sex, type of ms,...)")
BDsex
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Batch+Patient effect removed from data (col by Sex)

``` r
BDCorrBS <- removeBatchEffect(y$E,batch = meta$batch, batch2 = meta$sample)
BD2 <- pca(BDCorrBS, metadata = meta)
BDsexCorr <- biplot(BD2, colby = 'sex',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Batch+Patient effect removed from data")
BDsexCorr
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Batch+Patient effect removed from data (col by Sample)

``` r
BDsampleCorr <- biplot(BD2, colby = 'sample',hline = 0, vline = 0,legendPosition = 'right', lab="",title="Batch+Patient effect removed from data")
BDsampleCorr
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### Batch+Patient effect removed from data (col by Lesion type and Sex)

``` r
BDg <- biplot(BD2, colby = 'lesion_type',shape = 'sex',hline = 0, vline = 0,legendPosition = 'right', lab="",title="PCA")
## PCA
BDg
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

## gene name reference(ref_G)

``` r
conv_B=AnnotationDbi::select(org.Hs.eg.db,keys=rownames(y$E),keytype='ENSEMBL',columns=c('GENENAME','SYMBOL'))
ref_G=conv_B[!duplicated(conv_B[,1]),]
rownames(ref_G)=ref_G[,1]
```

## Specify Contrast(s) of interest

Use a linear model to fit voomed counts (logCPM). Instanciate the
comparisons of interest by using the contrast approach from *limma*.
Extract the results, annotate the tables and merge them into a global
result dataframe.

## Comparisons between conditions Batch3

``` r
fitB3 <- lmFit(yB3, mmB3)
```

    ## Coefficients not estimable: sampleAB253 sampleAB257

``` r
contr <- makeContrasts(
  A = lesion_typeWML - lesion_typeNAWM,
  B = lesion_typeWML - lesion_typeCL,
  C = lesion_typeCL - lesion_typeNAGM,
  D = lesion_typeNAWM - lesion_typeNAGM,
  E = lesion_typePVML - lesion_typeNAWM,
  F = lesion_typePVML - lesion_typeWML,
  G = lesion_typePVML - lesion_typeCL,
  H = lesion_typeNAGM - lesion_typeGM,
  I = lesion_typePVML - lesion_typePVM,
  J = lesion_typeWML - lesion_typeWM,
  levels = colnames(coef(fitB3)))
tmp <- contrasts.fit(fitB3, contr)
tmp <- eBayes(tmp)

resA_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='A')
resB_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='B')
resC_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='C')
resD_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='D')
resE_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='E')
resF_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='F')
resG_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='G')
resH_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='H')
resI_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='I')
resJ_B3 <- topTable(tmp, sort.by = "P", n = Inf,coef='J')
```

## Comparisons between conditions All data

``` r
fit <- lmFit(y, mm)
```

    ## Coefficients not estimable: sampleAB257 sexM

``` r
contr <- makeContrasts(
  A = lesion_typeWML - lesion_typeNAWM,
  B = lesion_typeWML - lesion_typeCL,
  C = lesion_typeCL - lesion_typeNAGM,
  D = lesion_typeNAWM - lesion_typeNAGM,
  E = lesion_typePVML - lesion_typeNAWM,
  F = lesion_typePVML - lesion_typeWML,
  G = lesion_typePVML - lesion_typeCL,
  H = lesion_typeNAGM - lesion_typeGM,
  I = lesion_typePVML - lesion_typePVM,
  J = lesion_typeWML - lesion_typeWM,
  levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

resA <- topTable(tmp, sort.by = "P", n = Inf,coef='A')
resB <- topTable(tmp, sort.by = "P", n = Inf,coef='B')
resC <- topTable(tmp, sort.by = "P", n = Inf,coef='C')
resD <- topTable(tmp, sort.by = "P", n = Inf,coef='D')
resE <- topTable(tmp, sort.by = "P", n = Inf,coef='E')
resF <- topTable(tmp, sort.by = "P", n = Inf,coef='F')
resG <- topTable(tmp, sort.by = "P", n = Inf,coef='G')
resH <- topTable(tmp, sort.by = "P", n = Inf,coef='H')
resI <- topTable(tmp, sort.by = "P", n = Inf,coef='I')
resJ <- topTable(tmp, sort.by = "P", n = Inf,coef='J')
```

## Comparisons between group MS & HC

``` r
mmG <- model.matrix(~ 0 + group+batch+sample+sex, data=meta[colnames(d),])

# Filtering low expressed genes
yG <- voom(d[,rownames(meta)], mm, plot = F)
```

    ## Coefficients not estimable: sampleAB257 sexM

``` r
fitG <- lmFit(yG, mmG)
```

    ## Coefficients not estimable: sampleAB257 sexM

``` r
contrG <- makeContrasts(
  A = groupMS.WM - groupHC.WM,
  B = groupMS.Cortical - groupHC.Cortical,
  C = groupMS.PV - groupHC.PV,
  levels = colnames(coef(fitG)))
tmpG <- contrasts.fit(fitG, contrG)
tmpG <- eBayes(tmpG)

resA_group <- topTable(tmpG, sort.by = "P", n = Inf,coef='A')
resB_group <- topTable(tmpG, sort.by = "P", n = Inf,coef='B')
resC_group <- topTable(tmpG, sort.by = "P", n = Inf,coef='C')
```

# Batch3

## Comparisons between conditions Batch3 (Volcano plot)

Plot results in the form of volcano plots. Thresholds used are
adj.P.Val\<0.05 & \|logFC\|\>1. Red genes are upregulated in the first
group as found in the subtitle of each graph.

### PVML_vs_NAWM

``` r
resE_B3$gene_id =  rownames(resE_B3)
resE_B3$comp = 'PVML_vs_NAWM'
resE_B3$gene_name = ref_G[rownames(resE_B3),3]
resE_B3=resE_B3[order(-abs(resE_B3$logFC)),]
resE_B3$show = F
resE_B3[1:15,'show']=T
resE_B3$sign = 'NS'
resE_B3[which(resE_B3$adj.P.Val<0.05&resE_B3$logFC>1),'sign']='UP'
resE_B3[which(resE_B3$adj.P.Val<0.05&resE_B3$logFC<(-1)),'sign']='DN'
ggplot(resE_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resE_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_NAWM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### WML_vs_NAWM

``` r
resA_B3$gene_id =  rownames(resA_B3)
resA_B3$comp = 'WML_vs_NAWM'
resA_B3$gene_name = ref_G[rownames(resA_B3),3]
resA_B3=resA_B3[order(-abs(resA_B3$logFC)),]
resA_B3$show = F
resA_B3[1:15,'show']=T
resA_B3$sign = 'NS'
resA_B3[which(resA_B3$adj.P.Val<0.05&resA_B3$logFC>1),'sign']='UP'
resA_B3[which(resA_B3$adj.P.Val<0.05&resA_B3$logFC<(-1)),'sign']='DN'
ggplot(resA_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resA_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('WML_vs_NAWM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

### WML_vs_CL

``` r
resB_B3$gene_id =  rownames(resB_B3)
resB_B3$comp = 'WML_vs_CL'
resB_B3$gene_name = ref_G[rownames(resB_B3),3]
resB_B3=resB_B3[order(-abs(resB_B3$logFC)),]
resB_B3$show = F
resB_B3[1:15,'show']=T
resB_B3$sign = 'NS'
resB_B3[which(resB_B3$adj.P.Val<0.05&resB_B3$logFC>1),'sign']='UP'
resB_B3[which(resB_B3$adj.P.Val<0.05&resB_B3$logFC<(-1)),'sign']='DN'
ggplot(resB_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resB_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('WML_vs_CL')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

### CL_vs_NAGM

``` r
resC_B3$gene_id =  rownames(resC_B3)
resC_B3$comp = 'CL_vs_NAGM'
resC_B3$gene_name = ref_G[rownames(resC_B3),3]
resC_B3=resC_B3[order(-abs(resC_B3$logFC)),]
resC_B3$show = F
resC_B3[1:15,'show']=T
resC_B3$sign = 'NS'
resC_B3[which(resC_B3$adj.P.Val<0.05&resC_B3$logFC>1),'sign']='UP'
resC_B3[which(resC_B3$adj.P.Val<0.05&resC_B3$logFC<(-1)),'sign']='DN'
ggplot(resC_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resC_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('CL_vs_NAGM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

### NAWM_vs_NAGM

``` r
resD_B3$gene_id =  rownames(resD_B3)
resD_B3$comp = 'NAWM_vs_NAGM'
resD_B3$gene_name = ref_G[rownames(resD_B3),3]
resD_B3=resD_B3[order(-abs(resD_B3$logFC)),]
resD_B3$show = F
resD_B3[1:15,'show']=T
resD_B3$sign = 'NS'
resD_B3[which(resD_B3$adj.P.Val<0.05&resD_B3$logFC>1),'sign']='UP'
resD_B3[which(resD_B3$adj.P.Val<0.05&resD_B3$logFC<(-1)),'sign']='DN'
ggplot(resD_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resD_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('NAWM_vs_NAGM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

### PVML_vs_WML

``` r
resF_B3$gene_id =  rownames(resF_B3)
resF_B3$comp = 'PVML_vs_WML'
resF_B3$gene_name = ref_G[rownames(resF_B3),3]
resF_B3=resF_B3[order(-abs(resF_B3$logFC)),]
resF_B3$show = F
resF_B3[1:15,'show']=T
resF_B3$sign = 'NS'
resF_B3[which(resF_B3$adj.P.Val<0.05&resF_B3$logFC>1),'sign']='UP'
resF_B3[which(resF_B3$adj.P.Val<0.05&resF_B3$logFC<(-1)),'sign']='DN'
ggplot(resF_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resF_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_WML')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

### PVML_vs_CL

``` r
resG_B3$gene_id =  rownames(resG_B3)
resG_B3$comp = 'PVML_vs_CL'
resG_B3$gene_name = ref_G[rownames(resG_B3),3]
resG_B3=resG_B3[order(-abs(resG_B3$logFC)),]
resG_B3$show = F
resG_B3[1:15,'show']=T
resG_B3$sign = 'NS'
resG_B3[which(resG_B3$adj.P.Val<0.05&resG_B3$logFC>1),'sign']='UP'
resG_B3[which(resG_B3$adj.P.Val<0.05&resG_B3$logFC<(-1)),'sign']='DN'
ggplot(resG_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resG_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_CL')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

### NAGM_vs_GM

``` r
resH_B3$gene_id =  rownames(resH_B3)
resH_B3$comp = 'NAGM_vs_GM'
resH_B3$gene_name = ref_G[rownames(resH_B3),3]
resH_B3=resH_B3[order(-abs(resH_B3$logFC)),]
resH_B3$show = F
resH_B3[1:15,'show']=T
resH_B3$sign = 'NS'
resH_B3[which(resH_B3$adj.P.Val<0.05&resH_B3$logFC>1),'sign']='UP'
resH_B3[which(resH_B3$adj.P.Val<0.05&resH_B3$logFC<(-1)),'sign']='DN'
ggplot(resH_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resH_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('NAGM_vs_GM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

### PVML_vs_PVM

``` r
resI_B3$gene_id =  rownames(resI_B3)
resI_B3$comp = 'PVML_vs_PVM'
resI_B3$gene_name = ref_G[rownames(resI_B3),3]
resI_B3=resI_B3[order(-abs(resI_B3$logFC)),]
resI_B3$show = F
resI_B3[1:15,'show']=T
resI_B3$sign = 'NS'
resI_B3[which(resI_B3$adj.P.Val<0.05&resI_B3$logFC>1),'sign']='UP'
resI_B3[which(resI_B3$adj.P.Val<0.05&resI_B3$logFC<(-1)),'sign']='DN'
ggplot(resI_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resI_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_PVM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

### WML_vs_WM

``` r
resJ_B3$gene_id =  rownames(resJ_B3)
resJ_B3$comp = 'WML_vs_WM'
resJ_B3$gene_name = ref_G[rownames(resJ_B3),3]
resJ_B3=resJ_B3[order(-abs(resJ_B3$logFC)),]
resJ_B3$show = F
resJ_B3[1:15,'show']=T
resJ_B3$sign = 'NS'
resJ_B3[which(resJ_B3$adj.P.Val<0.05&resJ_B3$logFC>1),'sign']='UP'
resJ_B3[which(resJ_B3$adj.P.Val<0.05&resJ_B3$logFC<(-1)),'sign']='DN'
ggplot(resJ_B3,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resJ_B3,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('WML_vs_WM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

# ALL DATA

## Comparisons between conditions (Volcano plot)

Plot results in the form of volcano plots. Thresholds used are
adj.P.Val\<0.05 & \|logFC\|\>1. Red genes are upregulated in the first
group as found in the subtitle of each graph.

### PVML_vs_NAWM

``` r
resE$gene_id =  rownames(resE)
resE$comp = 'PVML_vs_NAWM'
resE$gene_name = ref_G[rownames(resE),3]
resE=resE[order(-abs(resE$logFC)),]
resE$show = F
resE[1:15,'show']=T
resE$sign = 'NS'
resE[which(resE$adj.P.Val<0.05&resE$logFC>1),'sign']='UP'
resE[which(resE$adj.P.Val<0.05&resE$logFC<(-1)),'sign']='DN'
ggplot(resE,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resE,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_NAWM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

### WML_vs_NAWM

``` r
resA$gene_id =  rownames(resA)
resA$comp = 'WML_vs_NAWM'
resA$gene_name = ref_G[rownames(resA),3]
resA=resA[order(-abs(resA$logFC)),]
resA$show = F
resA[1:15,'show']=T
resA$sign = 'NS'
resA[which(resA$adj.P.Val<0.05&resA$logFC>1),'sign']='UP'
resA[which(resA$adj.P.Val<0.05&resA$logFC<(-1)),'sign']='DN'
ggplot(resA,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resA,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('WML_vs_NAWM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

### WML_vs_CL

``` r
resB$gene_id =  rownames(resB)
resB$comp = 'WML_vs_CL'
resB$gene_name = ref_G[rownames(resB),3]
resB=resB[order(-abs(resB$logFC)),]
resB$show = F
resB[1:15,'show']=T
resB$sign = 'NS'
resB[which(resB$adj.P.Val<0.05&resB$logFC>1),'sign']='UP'
resB[which(resB$adj.P.Val<0.05&resB$logFC<(-1)),'sign']='DN'
ggplot(resB,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resB,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('WML_vs_CL')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

### CL_vs_NAGM

``` r
resC$gene_id =  rownames(resC)
resC$comp = 'CL_vs_NAGM'
resC$gene_name = ref_G[rownames(resC),3]
resC=resC[order(-abs(resC$logFC)),]
resC$show = F
resC[1:15,'show']=T
resC$sign = 'NS'
resC[which(resC$adj.P.Val<0.05&resC$logFC>1),'sign']='UP'
resC[which(resC$adj.P.Val<0.05&resC$logFC<(-1)),'sign']='DN'
ggplot(resC,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resC,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('CL_vs_NAGM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

### NAWM_vs_NAGM

``` r
resD$gene_id =  rownames(resD)
resD$comp = 'NAWM_vs_NAGM'
resD$gene_name = ref_G[rownames(resD),3]
resD=resD[order(-abs(resD$logFC)),]
resD$show = F
resD[1:15,'show']=T
resD$sign = 'NS'
resD[which(resD$adj.P.Val<0.05&resD$logFC>1),'sign']='UP'
resD[which(resD$adj.P.Val<0.05&resD$logFC<(-1)),'sign']='DN'
ggplot(resD,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resD,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('NAWM_vs_NAGM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

### PVML_vs_WML

``` r
resF$gene_id =  rownames(resF)
resF$comp = 'PVML_vs_WML'
resF$gene_name = ref_G[rownames(resF),3]
resF=resF[order(-abs(resF$logFC)),]
resF$show = F
resF[1:15,'show']=T
resF$sign = 'NS'
resF[which(resF$adj.P.Val<0.05&resF$logFC>1),'sign']='UP'
resF[which(resF$adj.P.Val<0.05&resF$logFC<(-1)),'sign']='DN'
ggplot(resF,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resF,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_WML')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

### PVML_vs_CL

``` r
resG$gene_id =  rownames(resG)
resG$comp = 'PVML_vs_CL'
resG$gene_name = ref_G[rownames(resG),3]
resG=resG[order(-abs(resG$logFC)),]
resG$show = F
resG[1:15,'show']=T
resG$sign = 'NS'
resG[which(resG$adj.P.Val<0.05&resG$logFC>1),'sign']='UP'
resG[which(resG$adj.P.Val<0.05&resG$logFC<(-1)),'sign']='DN'
ggplot(resG,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resG,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_CL')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

### NAGM_vs_GM

``` r
resH$gene_id =  rownames(resH)
resH$comp = 'NAGM_vs_GM'
resH$gene_name = ref_G[rownames(resH),3]
resH=resH[order(-abs(resH$logFC)),]
resH$show = F
resH[1:15,'show']=T
resH$sign = 'NS'
resH[which(resH$adj.P.Val<0.05&resH$logFC>1),'sign']='UP'
resH[which(resH$adj.P.Val<0.05&resH$logFC<(-1)),'sign']='DN'
ggplot(resH,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resH,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('NAGM_vs_GM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

### PVML_vs_PVM

``` r
resI$gene_id =  rownames(resI)
resI$comp = 'PVML_vs_PVM'
resI$gene_name = ref_G[rownames(resI),3]
resI=resI[order(-abs(resI$logFC)),]
resI$show = F
resI[1:15,'show']=T
resI$sign = 'NS'
resI[which(resI$adj.P.Val<0.05&resI$logFC>1),'sign']='UP'
resI[which(resI$adj.P.Val<0.05&resI$logFC<(-1)),'sign']='DN'
ggplot(resI,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resI,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('PVML_vs_PVM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

### WML_vs_WM

``` r
resJ$gene_id =  rownames(resJ)
resJ$comp = 'WML_vs_WM'
resJ$gene_name = ref_G[rownames(resJ),3]
resJ=resJ[order(-abs(resJ$logFC)),]
resJ$show = F
resJ[1:15,'show']=T
resJ$sign = 'NS'
resJ[which(resJ$adj.P.Val<0.05&resJ$logFC>1),'sign']='UP'
resJ[which(resJ$adj.P.Val<0.05&resJ$logFC<(-1)),'sign']='DN'
ggplot(resJ,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resJ,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('WML_vs_WM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

## Comparisons between group MS & HC (Volcano plot)

Plot results in the form of volcano plots. Thresholds used are
adj.P.Val\<0.05 & \|logFC\|\>1. Red genes are upregulated in the first
group as found in the subtitle of each graph.

### groupMS.WM_vs_groupHC.WM

``` r
resA_group$gene_id =  rownames(resA_group)
resA_group$comp = 'groupMS.WM_vs_groupHC.WM'
resA_group$gene_name = ref_G[rownames(resA_group),3]
resA_group=resA_group[order(-abs(resA_group$logFC)),]
resA_group$show = F
resA_group[1:15,'show']=T
resA_group$sign = 'NS'
resA_group[which(resA_group$adj.P.Val<0.05&resA_group$logFC>1),'sign']='UP'
resA_group[which(resA_group$adj.P.Val<0.05&resA_group$logFC<(-1)),'sign']='DN'
ggplot(resA_group,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resA_group,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('groupMS.WM_vs_groupHC.WM')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

### groupMS.Cortical_vs_groupHC.Cortical

``` r
resB_group$gene_id =  rownames(resB_group)
resB_group$comp = 'groupMS.Cortical_vs_groupHC.Cortical'
resB_group$gene_name = ref_G[rownames(resB_group),3]
resB_group=resB_group[order(-abs(resB_group$logFC)),]
resB_group$show = F
resB_group[1:15,'show']=T
resB_group$sign = 'NS'
resB_group[which(resB_group$adj.P.Val<0.05&resB_group$logFC>1),'sign']='UP'
resB_group[which(resB_group$adj.P.Val<0.05&resB_group$logFC<(-1)),'sign']='DN'
ggplot(resB_group,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resB_group,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('groupMS.Cortical_vs_groupHC.Cortical')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

### groupMS.PV_vs_groupHC.PV

``` r
resC_group$gene_id =  rownames(resC_group)
resC_group$comp = 'groupMS.PV_vs_groupHC.PV'
resC_group$gene_name = ref_G[rownames(resC_group),3]
resC_group=resC_group[order(-abs(resC_group$logFC)),]
resC_group$show = F
resC_group[1:15,'show']=T
resC_group$sign = 'NS'
resC_group[which(resC_group$adj.P.Val<0.05&resC_group$logFC>1),'sign']='UP'
resC_group[which(resC_group$adj.P.Val<0.05&resC_group$logFC<(-1)),'sign']='DN'
ggplot(resC_group,aes(logFC,-log10(P.Value),col=sign))+geom_point_interactive(aes(x=logFC,y=-log10(P.Value),color=sign,tooltip = gene_name,data_id = gene_name))+theme_bw()+scale_color_manual(values=sign_cols)+ geom_label_repel(data=subset(resC_group,show==T),aes(label=gene_name),col='black',size=3,force=30,fill = alpha(c("white"),0.5))+ggtitle('groupMS.PV_vs_groupHC.PV')
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

``` r
resE100=head(resE,100)
pheatmap::pheatmap(y$E[rownames(resE100),],scale = "row")
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

``` r
pheatmap::pheatmap(cor(y$E[rownames(resE100),]))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-49-2.png)<!-- -->

## Gene ontology term enrichment analysis of comparisons between conditions

This extracts biological signatures from the transcriptomic changes
found in the differential expression analysis

``` r
enr.a=find_enrichments(resA)
enr.b=find_enrichments(resB)
enr.c=find_enrichments(resC)
enr.d=find_enrichments(resD)
enr.e=find_enrichments(resE)
enr.f=find_enrichments(resF)
enr.g=find_enrichments(resG)
enr.h=find_enrichments(resH)
enr.i=find_enrichments(resI)
enr.j=find_enrichments(resJ)
```

### Biological Processes

``` r
enrA = enr.a$enr$bp
enrA=enrA[order(enrA$padj),]
enrA$comp = unique(resA$comp)

enrB = enr.b$enr$bp
enrB=enrB[order(enrB$padj),]
enrB$comp = unique(resB$comp)

enrC = enr.c$enr$bp
enrC=enrC[order(enrC$padj),]
enrC$comp = unique(resC$comp)

enrD = enr.d$enr$bp
enrD=enrD[order(enrD$padj),]
enrD$comp = unique(resD$comp)

enrE = enr.e$enr$bp
enrE=enrE[order(enrE$padj),]
enrE$comp = unique(resE$comp)

enrF = enr.f$enr$bp
enrF=enrF[order(enrF$padj),]
enrF$comp = unique(resF$comp)

enrG = enr.g$enr$bp
enrG=enrG[order(enrG$padj),]
enrG$comp = unique(resG$comp)

enrH = enr.h$enr$bp
enrH=enrH[order(enrH$padj),]
enrH$comp = unique(resH$comp)

enrI = enr.i$enr$bp
enrI=enrI[order(enrI$padj),]
enrI$comp = unique(resI$comp)

enrJ = enr.j$enr$bp
enrJ=enrJ[order(enrJ$padj),]
enrJ$comp = unique(resJ$comp)

paths = c(head(enrA,10)$pathway,
          head(enrB,10)$pathway,
          head(enrC,10)$pathway,
          head(enrD,10)$pathway,
          head(enrE,10)$pathway,
          head(enrF,10)$pathway,
          head(enrG,10)$pathway,
          head(enrH,10)$pathway,
          head(enrI,10)$pathway,
          head(enrJ,10)$pathway)

paths=unique(paths)
rb=rbind(enrA[which(enrA$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrB[which(enrB$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrC[which(enrC$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrD[which(enrD$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrE[which(enrE$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrF[ which(enrF$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrG[which(enrG$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrH[which(enrH$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrI[ which(enrI$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrJ[which(enrJ$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')])
rb[which(rb$padj>=0.1),'NES']=NA
rb$pathway <- gsub("GOBP", "",rb$pathway)
rb$pathway <- gsub("_", " ",rb$pathway)
ss  = strsplit(as.character(rb$leadingEdge),',')
rb$n_genes = unlist(lapply(ss,length))
rb$pval <- as.numeric(rb$pval)
rb[which(rb$pval>=0.1),'pval']=NA
ggplot(na.omit(rb),aes(pathway,comp,fill=NES))+geom_point(aes(size=pval),shape=21)+scale_x_discrete(label = function(x) stringr::str_trunc(x,50))+theme_bw()+coord_flip()+theme(axis.text.x=element_text( angle=60,hjust=1))+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+ scale_size(trans = 'reverse')+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

### Molecular Function

``` r
enrA = enr.a$enr$mf
enrA=enrA[order(enrA$padj),]
enrA$comp = unique(resA$comp)

enrB = enr.b$enr$mf
enrB=enrB[order(enrB$padj),]
enrB$comp = unique(resB$comp)

enrC = enr.c$enr$mf
enrC=enrC[order(enrC$padj),]
enrC$comp = unique(resC$comp)

enrD = enr.d$enr$mf
enrD=enrD[order(enrD$padj),]
enrD$comp = unique(resD$comp)

enrE = enr.e$enr$mf
enrE=enrE[order(enrE$padj),]
enrE$comp = unique(resE$comp)

enrF = enr.f$enr$mf
enrF=enrF[order(enrF$padj),]
enrF$comp = unique(resF$comp)

enrG = enr.g$enr$mf
enrG=enrG[order(enrG$padj),]
enrG$comp = unique(resG$comp)

enrH = enr.h$enr$mf
enrH=enrH[order(enrH$padj),]
enrH$comp = unique(resH$comp)

enrI = enr.i$enr$mf
enrI=enrI[order(enrI$padj),]
enrI$comp = unique(resI$comp)

enrJ = enr.j$enr$mf
enrJ=enrJ[order(enrJ$padj),]
enrJ$comp = unique(resJ$comp)

paths = c(head(enrA,10)$pathway,
          head(enrB,10)$pathway,
          head(enrC,10)$pathway,
          head(enrD,10)$pathway,
          head(enrE,10)$pathway,
          head(enrF,10)$pathway,
          head(enrG,10)$pathway,
          head(enrH,10)$pathway,
          head(enrI,10)$pathway,
          head(enrJ,10)$pathway)

paths=unique(paths)
rb=rbind(enrA[which(enrA$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrB[which(enrB$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrC[which(enrC$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrD[which(enrD$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrE[which(enrE$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrF[ which(enrF$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrG[which(enrG$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrH[which(enrH$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrI[ which(enrI$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrJ[which(enrJ$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')])

rb[which(rb$padj>=0.1),'NES']=NA
rb$pathway <- gsub("GOMF", "",rb$pathway)
rb$pathway <- gsub("_", " ",rb$pathway)
ss  = strsplit(as.character(rb$leadingEdge),',')
rb$n_genes = unlist(lapply(ss,length))
rb$pval <- as.numeric(rb$pval)
rb[which(rb$pval>=0.1),'pval']=NA
ggplot(na.omit(rb),aes(pathway,comp,fill=NES))+geom_point(aes(size=pval),shape=21)+scale_x_discrete(label = function(x) stringr::str_trunc(x,50))+theme_bw()+coord_flip()+theme(axis.text.x=element_text( angle=60,hjust=1))+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+ scale_size(trans = 'reverse')+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

### KEGG pathway database

``` r
enrA = enr.a$enr$kegg
enrA=enrA[order(enrA$padj),]
enrA$comp = unique(resA$comp)

enrB = enr.b$enr$kegg
enrB=enrB[order(enrB$padj),]
enrB$comp = unique(resB$comp)

enrC = enr.c$enr$kegg
enrC=enrC[order(enrC$padj),]
enrC$comp = unique(resC$comp)

enrD = enr.d$enr$kegg
enrD=enrD[order(enrD$padj),]
enrD$comp = unique(resD$comp)

enrE = enr.e$enr$kegg
enrE=enrE[order(enrE$padj),]
enrE$comp = unique(resE$comp)

enrF = enr.f$enr$kegg
enrF=enrF[order(enrF$padj),]
enrF$comp = unique(resF$comp)

enrG = enr.g$enr$kegg
enrG=enrG[order(enrG$padj),]
enrG$comp = unique(resG$comp)

enrH = enr.h$enr$kegg
enrH=enrH[order(enrH$padj),]
enrH$comp = unique(resH$comp)

enrI = enr.i$enr$kegg
enrI=enrI[order(enrI$padj),]
enrI$comp = unique(resI$comp)

enrJ = enr.j$enr$kegg
enrJ=enrJ[order(enrJ$padj),]
enrJ$comp = unique(resJ$comp)

paths = c(head(enrA,10)$pathway,
          head(enrB,10)$pathway,
          head(enrC,10)$pathway,
          head(enrD,10)$pathway,
          head(enrE,10)$pathway,
          head(enrF,10)$pathway,
          head(enrG,10)$pathway,
          head(enrH,10)$pathway,
          head(enrI,10)$pathway,
          head(enrJ,10)$pathway)

paths=unique(paths)
rb=rbind(enrA[which(enrA$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrB[which(enrB$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrC[which(enrC$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrD[which(enrD$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrE[which(enrE$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrF[ which(enrF$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrG[which(enrG$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrH[which(enrH$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrI[ which(enrI$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')],
         enrJ[which(enrJ$pathway%in%paths),c('pval','leadingEdge','NES','comp','pathway')])

rb[which(rb$padj>=0.1),'NES']=NA
rb$pathway <- gsub("KEGG", "",rb$pathway)
rb$pathway <- gsub("_", " ",rb$pathway)
ss  = strsplit(as.character(rb$leadingEdge),',')
rb$n_genes = unlist(lapply(ss,length))
rb$pval <- as.numeric(rb$pval)
rb[which(rb$pval>=0.1),'pval']=NA
ggplot(na.omit(rb),aes(pathway,comp,fill=NES))+geom_point(aes(size=pval),shape=21)+scale_x_discrete(label = function(x) stringr::str_trunc(x,50))+theme_bw()+coord_flip()+theme(axis.text.x=element_text( angle=60,hjust=1))+scale_fill_gradientn(colors=rev(brewer.pal(11,'RdYlBu')))+ scale_size(trans = 'reverse')+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

# Plot the overlap between DEGs with an upset plot/Venn Diagram

## Overlap of DEGs in comparaisons

``` r
res=rbind(resA,resB,resC,resD,resE,resF,resG,resH,resI,resJ)
res= na.omit(res)

degsListvs = list(
    WML_vs_NAWM=subset(res,comp=='WML_vs_NAWM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    PVML_vs_WML=subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    PVML_vs_NAWM=subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    WML_vs_CL=subset(res,comp=='WML_vs_CL' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    #NAWM_vs_NAGM=subset(res,comp=='NAWM_vs_NAGM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    NAGM_vs_GM=subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    PVML_vs_PVM=subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    WML_vs_WM=subset(res,comp=='WML_vs_WM' & P.Value<0.05 & abs(logFC)>1)$gene_name)
    

degsListUP = list(
    WML_vs_NAWM=subset(res,comp=='WML_vs_NAWM' & P.Value<0.05 & (logFC)>1)$gene_name,
    PVML_vs_WML=subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & (logFC)>1)$gene_name,
    PVML_vs_NAWM=subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & (logFC)>1)$gene_name,
    WML_vs_CL=subset(res,comp=='WML_vs_CL' & P.Value<0.05 & (logFC)>1)$gene_name,
    #NAWM_vs_NAGM=subset(res,comp=='NAWM_vs_NAGM' & P.Value<0.05 & (logFC)>1)$gene_name,
    NAGM_vs_GM=subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & (logFC)>1)$gene_name,
    PVML_vs_PVM=subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & (logFC)>1)$gene_name,
    WML_vs_WM=subset(res,comp=='WML_vs_WM' & P.Value<0.05 & (logFC)>1)$gene_name)

degsListDOWN = list(
    WML_vs_NAWM=subset(res,comp=='WML_vs_NAWM' & P.Value<0.05 & (logFC)<1)$gene_name,
    PVML_vs_WML=subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & (logFC)<1)$gene_name,
    PVML_vs_NAWM=subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & (logFC)<1)$gene_name,
    WML_vs_CL=subset(res,comp=='WML_vs_CL' & P.Value<0.05 & (logFC)<1)$gene_name,
    #NAWM_vs_NAGM=subset(res,comp=='NAWM_vs_NAGM' & P.Value<0.05 & (logFC)<1)$gene_name,
    NAGM_vs_GM=subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & (logFC)<1)$gene_name,
    PVML_vs_PVM=subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & (logFC)<1)$gene_name,
    WML_vs_WM=subset(res,comp=='WML_vs_WM' & P.Value<0.05 & (logFC)<1)$gene_name)

degsListvs_UpSet = list(
    WML_vs_NAWM=subset(res,comp=='WML_vs_NAWM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    PVML_vs_WML=subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    PVML_vs_NAWM=subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    WML_vs_CL=subset(res,comp=='WML_vs_CL' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    NAWM_vs_NAGM=subset(res,comp=='NAWM_vs_NAGM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    NAGM_vs_GM=subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    PVML_vs_PVM=subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & abs(logFC)>1)$gene_name,
    WML_vs_WM=subset(res,comp=='WML_vs_WM' & P.Value<0.05 & abs(logFC)>1)$gene_name)
    

degsListUP_UpSet = list(
    WML_vs_NAWM=subset(res,comp=='WML_vs_NAWM' & P.Value<0.05 & (logFC)>1)$gene_name,
    PVML_vs_WML=subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & (logFC)>1)$gene_name,
    PVML_vs_NAWM=subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & (logFC)>1)$gene_name,
    WML_vs_CL=subset(res,comp=='WML_vs_CL' & P.Value<0.05 & (logFC)>1)$gene_name,
    NAWM_vs_NAGM=subset(res,comp=='NAWM_vs_NAGM' & P.Value<0.05 & (logFC)>1)$gene_name,
    NAGM_vs_GM=subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & (logFC)>1)$gene_name,
    PVML_vs_PVM=subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & (logFC)>1)$gene_name,
    WML_vs_WM=subset(res,comp=='WML_vs_WM' & P.Value<0.05 & (logFC)>1)$gene_name)

degsListDOWN_UpSet = list(
    WML_vs_NAWM=subset(res,comp=='WML_vs_NAWM' & P.Value<0.05 & (logFC)<1)$gene_name,
    PVML_vs_WML=subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & (logFC)<1)$gene_name,
    PVML_vs_NAWM=subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & (logFC)<1)$gene_name,
    WML_vs_CL=subset(res,comp=='WML_vs_CL' & P.Value<0.05 & (logFC)<1)$gene_name,
    NAWM_vs_NAGM=subset(res,comp=='NAWM_vs_NAGM' & P.Value<0.05 & (logFC)<1)$gene_name,
    NAGM_vs_GM=subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & (logFC)<1)$gene_name,
    PVML_vs_PVM=subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & (logFC)<1)$gene_name,
    WML_vs_WM=subset(res,comp=='WML_vs_WM' & P.Value<0.05 & (logFC)<1)$gene_name)
```

### Venn Diagram

#### Overlap between DEGs (All)

``` r
#label = "count"
ggVennDiagram(degsListvs,label = "percent", label_alpha = 0)+ scale_fill_distiller(palette = "Reds", direction = 1)
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

#### Overlap between DEGs (UP)

``` r
#label = "count"
ggVennDiagram(degsListUP,label = "percent", label_alpha = 0)+ scale_fill_distiller(palette = "Reds", direction = 1)
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

#### Overlap between DEGs (DOWN)

``` r
#label = "count"
ggVennDiagram(degsListDOWN,label = "percent", label_alpha = 0)+ scale_fill_distiller(palette = "Reds", direction = 1)
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

### UpSet plot

#### Overlap between DEGs (All)

``` r
upset(fromList(degsListvs_UpSet),sets = c('WML_vs_NAWM','PVML_vs_WML','PVML_vs_NAWM',"WML_vs_CL","NAWM_vs_NAGM","NAGM_vs_GM", "PVML_vs_PVM", "WML_vs_WM"),order.by='freq',decreasing=c(F), text.scale = c(1.3, 2, 1.5, 1.5, 2, 2))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

#### Overlap between DEGs (UP)

``` r
upset(fromList(degsListUP_UpSet),sets = c('WML_vs_NAWM','PVML_vs_WML','PVML_vs_NAWM',"WML_vs_CL","NAWM_vs_NAGM","NAGM_vs_GM", "PVML_vs_PVM", "WML_vs_WM"),order.by='freq',decreasing = c(F),text.scale = c(1.3, 2, 1.5, 1.5, 2, 2))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

#### Overlap between DEGs (DOWN)

``` r
upset(fromList(degsListDOWN_UpSet),sets = c('WML_vs_NAWM','PVML_vs_WML','PVML_vs_NAWM',"WML_vs_CL","NAWM_vs_NAGM","NAGM_vs_GM", "PVML_vs_PVM", "WML_vs_WM"),order.by='freq',decreasing = c(F),text.scale = c(1.3, 2, 1.5, 1.5, 2, 2))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

# Gene Ontology enrichment analysis of overlapping genes in comparisons between conditions(The top 5 specific Genesets from the upsets plot)

``` r
WML_vs_CL=na.omit(subset(res,comp=='WML_vs_CL' & P.Value<0.05 & logFC>1)$gene_name)
PVML_vs_PVM=na.omit(subset(res,comp=='PVML_vs_PVM' & P.Value<0.05 & logFC>1)$gene_name)
PVML_vs_NAWM=na.omit(subset(res,comp=='PVML_vs_NAWM' & P.Value<0.05 & logFC>1)$gene_name)
NAGM_vs_GM=na.omit(subset(res,comp=='NAGM_vs_GM' & P.Value<0.05 & logFC>1)$gene_name)
PVML_vs_WML=na.omit(subset(res,comp=='PVML_vs_WML' & P.Value<0.05 & logFC>1)$gene_name)

specificWML_vs_CL = as.character(WML_vs_CL[!WML_vs_CL%in%PVML_vs_PVM & !WML_vs_CL%in%PVML_vs_NAWM & !WML_vs_CL%in%NAGM_vs_GM & !WML_vs_CL%in%PVML_vs_WML])

specificPVML_vs_PVM = as.character(PVML_vs_PVM[!PVML_vs_PVM%in%WML_vs_CL & !PVML_vs_PVM%in%PVML_vs_NAWM & !PVML_vs_PVM%in%NAGM_vs_GM & !PVML_vs_PVM%in%PVML_vs_WML])

specificPVML_vs_NAWM = as.character(PVML_vs_NAWM[!PVML_vs_NAWM%in%WML_vs_CL & !PVML_vs_NAWM%in%PVML_vs_PVM & !PVML_vs_NAWM%in%NAGM_vs_GM & !PVML_vs_NAWM%in%PVML_vs_WML])

specificNAGM_vs_GM = as.character(NAGM_vs_GM[!NAGM_vs_GM%in%WML_vs_CL & !NAGM_vs_GM%in%PVML_vs_PVM & !NAGM_vs_GM%in%PVML_vs_NAWM & !NAGM_vs_GM%in%PVML_vs_WML])
  
commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = as.character(PVML_vs_NAWM[!PVML_vs_NAWM%in%WML_vs_CL & PVML_vs_NAWM%in%PVML_vs_WML & PVML_vs_NAWM%in%PVML_vs_PVM & !PVML_vs_NAWM%in%NAGM_vs_GM])
  

list = list(commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM=commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM, specificWML_vs_CL=specificWML_vs_CL,specificPVML_vs_PVM=specificPVML_vs_PVM, specificPVML_vs_NAWM=specificPVML_vs_NAWM,specificNAGM_vs_GM=specificNAGM_vs_GM)

fishBP = fisher_enrichment(list,go,unique(ref_G$SYMBOL))
fishMF = fisher_enrichment(list,mf,unique(ref_G$SYMBOL))
#fishKEGG = fisher_enrichment(list,kegg,unique(ref_G$SYMBOL))
```

## biological processes

``` r
  enr.specificWML_vs_CL = fishBP$specificWML_vs_CL
  enr.specificWML_vs_CL$pathway = rownames(enr.specificWML_vs_CL)
  enr.specificWML_vs_CL$fdr <- as.numeric(enr.specificWML_vs_CL$fdr) 
  enr.specificWML_vs_CL$OR.odds.ratio <- as.numeric(enr.specificWML_vs_CL$OR.odds.ratio) 
  enr.specificWML_vs_CL$NOM_pval <- as.numeric(enr.specificWML_vs_CL$NOM_pval)
  topspecificWML_vs_CL <- topPath(fishBP ,"specificWML_vs_CL", 0.5, 20)
  enr.specificWML_vs_CL <- enr.specificWML_vs_CL[match(topspecificWML_vs_CL, rownames(enr.specificWML_vs_CL)),]
  rownames(enr.specificWML_vs_CL) <- gsub("GOBP_", "",rownames(enr.specificWML_vs_CL)) 
  enr.specificWML_vs_CL=enr.specificWML_vs_CL[order(enr.specificWML_vs_CL$NOM_pval),]
  enr.specificWML_vs_CL$comp = "specificWML_vs_CL"
  
  enr.specificPVML_vs_PVM = fishBP$specificPVML_vs_PVM
  enr.specificPVML_vs_PVM$pathway = rownames(enr.specificPVML_vs_PVM)
  enr.specificPVML_vs_PVM$fdr <- as.numeric(enr.specificPVML_vs_PVM$fdr) 
  enr.specificPVML_vs_PVM$OR.odds.ratio <- as.numeric(enr.specificPVML_vs_PVM$OR.odds.ratio) 
  enr.specificPVML_vs_PVM$NOM_pval <- as.numeric(enr.specificPVML_vs_PVM$NOM_pval)
  topspecificPVML_vs_PVM <- topPath(fishBP ,"specificPVML_vs_PVM", 0.7, 20)
  enr.specificPVML_vs_PVM <- enr.specificPVML_vs_PVM[match(topspecificPVML_vs_PVM, rownames(enr.specificPVML_vs_PVM)),]
  rownames(enr.specificPVML_vs_PVM) <- gsub("GOBP_", "",rownames(enr.specificPVML_vs_PVM)) 
  enr.specificPVML_vs_PVM=enr.specificPVML_vs_PVM[order(enr.specificPVML_vs_PVM$NOM_pval),]
  enr.specificPVML_vs_PVM$comp = "specificPVML_vs_PVM"
  
  enr.specificPVML_vs_NAWM = fishBP$specificPVML_vs_NAWM
  enr.specificPVML_vs_NAWM$pathway = rownames(enr.specificPVML_vs_NAWM)
  enr.specificPVML_vs_NAWM$fdr <- as.numeric(enr.specificPVML_vs_NAWM$fdr) 
  enr.specificPVML_vs_NAWM$OR.odds.ratio <- as.numeric(enr.specificPVML_vs_NAWM$OR.odds.ratio) 
  enr.specificPVML_vs_NAWM$NOM_pval <- as.numeric(enr.specificPVML_vs_NAWM$NOM_pval)
  enr.specificPVML_vs_NAWM=enr.specificPVML_vs_NAWM[order(enr.specificPVML_vs_NAWM$NOM_pval),]
  topspecificPVML_vs_NAWM = subset(enr.specificPVML_vs_NAWM,NOM_pval<=0.1)
  topspecificPVML_vs_NAWM = topspecificPVML_vs_NAWM$pathway
  enr.specificPVML_vs_NAWM <- enr.specificPVML_vs_NAWM[match(topspecificPVML_vs_NAWM, rownames(enr.specificPVML_vs_NAWM)),]
  rownames(enr.specificPVML_vs_NAWM) <- gsub("GOBP_", "",rownames(enr.specificPVML_vs_NAWM)) 
  enr.specificPVML_vs_NAWM=enr.specificPVML_vs_NAWM[order(enr.specificPVML_vs_NAWM$NOM_pval),]
  enr.specificPVML_vs_NAWM$comp = "specificPVML_vs_NAWM"
  enr.specificPVML_vs_NAWM <- head(enr.specificPVML_vs_NAWM,20)
  
  enr.specificNAGM_vs_GM = fishBP$specificNAGM_vs_GM
  enr.specificNAGM_vs_GM$pathway = rownames(enr.specificNAGM_vs_GM)
  enr.specificNAGM_vs_GM$fdr <- as.numeric(enr.specificNAGM_vs_GM$fdr) 
  enr.specificNAGM_vs_GM$OR.odds.ratio <- as.numeric(enr.specificNAGM_vs_GM$OR.odds.ratio) 
  enr.specificNAGM_vs_GM$NOM_pval <- as.numeric(enr.specificNAGM_vs_GM$NOM_pval)
  enr.specificNAGM_vs_GM=enr.specificNAGM_vs_GM[order(enr.specificNAGM_vs_GM$NOM_pval),]
  topspecificNAGM_vs_GM = subset(enr.specificNAGM_vs_GM,NOM_pval<=0.1)
  topspecificNAGM_vs_GM = topspecificNAGM_vs_GM$pathway
  enr.specificNAGM_vs_GM <- enr.specificNAGM_vs_GM[match(topspecificNAGM_vs_GM, rownames(enr.specificNAGM_vs_GM)),]
  rownames(enr.specificNAGM_vs_GM) <- gsub("GOBP_", "",rownames(enr.specificNAGM_vs_GM)) 
  enr.specificNAGM_vs_GM=enr.specificNAGM_vs_GM[order(enr.specificNAGM_vs_GM$NOM_pval),]
  enr.specificNAGM_vs_GM$comp = "specificNAGM_vs_GM"
  enr.specificNAGM_vs_GM <- head(enr.specificNAGM_vs_GM,20)
  
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = fishBP$commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$pathway = rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM)
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$fdr <- as.numeric(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$fdr) 
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$OR.odds.ratio <- as.numeric(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$OR.odds.ratio) 
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval <- as.numeric(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval)
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM=enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM[order(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval),]
  topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = subset(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM,NOM_pval<=0.1)
  topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$pathway
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM <- enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM[match(topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM, rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM)),]
  rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM) <- gsub("GOBP_", "",rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM)) 
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM=enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM[order(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval),]
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$comp = "common PVML_vs_NAWM & PVML_vs_WML & PVML_vs_PVM"  
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM <- head(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM,20)
```

### Specific WML_vs_CL

``` r
ggplot(enr.specificWML_vs_CL, aes(x=reorder(rownames(enr.specificWML_vs_CL),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "WML_vs_CL") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

### Specific PVML_vs_PVM

``` r
ggplot(enr.specificPVML_vs_PVM, aes(x=reorder(rownames(enr.specificPVML_vs_PVM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "PVML_vs_PVM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

### Specific PVML_vs_NAWM

``` r
ggplot(enr.specificPVML_vs_NAWM, aes(x=reorder(rownames(enr.specificPVML_vs_NAWM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "PVML_vs_NAWM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

### Specific NAGM_vs_GM

``` r
ggplot(enr.specificNAGM_vs_GM, aes(x=reorder(rownames(enr.specificNAGM_vs_GM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "NAGM_vs_GM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

### common PVML_vs_NAWM & PVML_vs_WML & PVML_vs_PVM

``` r
ggplot(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM, aes(x=reorder(rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "common PVML_vs_NAWM & PVML_vs_WML & PVML_vs_PVM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

## Molecular Function

``` r
  enr.specificWML_vs_CL = fishMF$specificWML_vs_CL
  enr.specificWML_vs_CL$pathway = rownames(enr.specificWML_vs_CL)
  enr.specificWML_vs_CL$fdr <- as.numeric(enr.specificWML_vs_CL$fdr) 
  enr.specificWML_vs_CL$OR.odds.ratio <- as.numeric(enr.specificWML_vs_CL$OR.odds.ratio) 
  enr.specificWML_vs_CL$NOM_pval <- as.numeric(enr.specificWML_vs_CL$NOM_pval)
  enr.specificWML_vs_CL=enr.specificWML_vs_CL[order(enr.specificWML_vs_CL$NOM_pval),]
  topspecificWML_vs_CL = subset(enr.specificWML_vs_CL,NOM_pval<=0.1)
  topspecificWML_vs_CL = topspecificWML_vs_CL$pathway
  enr.specificWML_vs_CL <- enr.specificWML_vs_CL[match(topspecificWML_vs_CL, rownames(enr.specificWML_vs_CL)),]
  rownames(enr.specificWML_vs_CL) <- gsub("GOMF_", "",rownames(enr.specificWML_vs_CL)) 
  enr.specificWML_vs_CL=enr.specificWML_vs_CL[order(enr.specificWML_vs_CL$NOM_pval),]
  enr.specificWML_vs_CL$comp = "specificWML_vs_CL"
  enr.specificWML_vs_CL <- head(enr.specificWML_vs_CL,20)
  
  enr.specificPVML_vs_PVM = fishMF$specificPVML_vs_PVM
  enr.specificPVML_vs_PVM$pathway = rownames(enr.specificPVML_vs_PVM)
  enr.specificPVML_vs_PVM$fdr <- as.numeric(enr.specificPVML_vs_PVM$fdr) 
  enr.specificPVML_vs_PVM$OR.odds.ratio <- as.numeric(enr.specificPVML_vs_PVM$OR.odds.ratio) 
  enr.specificPVML_vs_PVM$NOM_pval <- as.numeric(enr.specificPVML_vs_PVM$NOM_pval)
  enr.specificPVML_vs_PVM=enr.specificPVML_vs_PVM[order(enr.specificPVML_vs_PVM$NOM_pval),]
  topspecificPVML_vs_PVM = subset(enr.specificPVML_vs_PVM,NOM_pval<=0.1)
  topspecificPVML_vs_PVM = topspecificPVML_vs_PVM$pathway
  enr.specificPVML_vs_PVM <- enr.specificPVML_vs_PVM[match(topspecificPVML_vs_PVM, rownames(enr.specificPVML_vs_PVM)),]
  rownames(enr.specificPVML_vs_PVM) <- gsub("GOMF_", "",rownames(enr.specificPVML_vs_PVM)) 
  enr.specificPVML_vs_PVM=enr.specificPVML_vs_PVM[order(enr.specificPVML_vs_PVM$NOM_pval),]
  enr.specificPVML_vs_PVM$comp = "specificPVML_vs_PVM"
  enr.specificPVML_vs_PVM <- head(enr.specificPVML_vs_PVM,20)
  
  enr.specificPVML_vs_NAWM = fishMF$specificPVML_vs_NAWM
  enr.specificPVML_vs_NAWM$pathway = rownames(enr.specificPVML_vs_NAWM)
  enr.specificPVML_vs_NAWM$fdr <- as.numeric(enr.specificPVML_vs_NAWM$fdr) 
  enr.specificPVML_vs_NAWM$OR.odds.ratio <- as.numeric(enr.specificPVML_vs_NAWM$OR.odds.ratio) 
  enr.specificPVML_vs_NAWM$NOM_pval <- as.numeric(enr.specificPVML_vs_NAWM$NOM_pval)
  enr.specificPVML_vs_NAWM=enr.specificPVML_vs_NAWM[order(enr.specificPVML_vs_NAWM$NOM_pval),]
  topspecificPVML_vs_NAWM = subset(enr.specificPVML_vs_NAWM,NOM_pval<=0.1)
  topspecificPVML_vs_NAWM = topspecificPVML_vs_NAWM$pathway
  enr.specificPVML_vs_NAWM <- enr.specificPVML_vs_NAWM[match(topspecificPVML_vs_NAWM, rownames(enr.specificPVML_vs_NAWM)),]
  rownames(enr.specificPVML_vs_NAWM) <- gsub("GOMF_", "",rownames(enr.specificPVML_vs_NAWM)) 
  enr.specificPVML_vs_NAWM=enr.specificPVML_vs_NAWM[order(enr.specificPVML_vs_NAWM$NOM_pval),]
  enr.specificPVML_vs_NAWM$comp = "specificPVML_vs_NAWM"
  enr.specificPVML_vs_NAWM <- head(enr.specificPVML_vs_NAWM,20)
  
  enr.specificNAGM_vs_GM = fishMF$specificNAGM_vs_GM
  enr.specificNAGM_vs_GM$pathway = rownames(enr.specificNAGM_vs_GM)
  enr.specificNAGM_vs_GM$fdr <- as.numeric(enr.specificNAGM_vs_GM$fdr) 
  enr.specificNAGM_vs_GM$OR.odds.ratio <- as.numeric(enr.specificNAGM_vs_GM$OR.odds.ratio) 
  enr.specificNAGM_vs_GM$NOM_pval <- as.numeric(enr.specificNAGM_vs_GM$NOM_pval)
  enr.specificNAGM_vs_GM=enr.specificNAGM_vs_GM[order(enr.specificNAGM_vs_GM$NOM_pval),]
  topspecificNAGM_vs_GM = subset(enr.specificNAGM_vs_GM,NOM_pval<=0.1)
  topspecificNAGM_vs_GM = topspecificNAGM_vs_GM$pathway
  enr.specificNAGM_vs_GM <- enr.specificNAGM_vs_GM[match(topspecificNAGM_vs_GM, rownames(enr.specificNAGM_vs_GM)),]
  rownames(enr.specificNAGM_vs_GM) <- gsub("GOMF_", "",rownames(enr.specificNAGM_vs_GM)) 
  enr.specificNAGM_vs_GM=enr.specificNAGM_vs_GM[order(enr.specificNAGM_vs_GM$NOM_pval),]
  enr.specificNAGM_vs_GM$comp = "specificNAGM_vs_GM"
  enr.specificNAGM_vs_GM <- head(enr.specificNAGM_vs_GM,20)

  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = fishMF$commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$pathway = rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM)
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$fdr <- as.numeric(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$fdr) 
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$OR.odds.ratio <- as.numeric(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$OR.odds.ratio) 
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval <- as.numeric(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval)
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM=enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM[order(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval),]
  topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = subset(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM,NOM_pval<=0.1)
  topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM = topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$pathway
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM <- enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM[match(topcommonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM, rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM)),]
  rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM) <- gsub("GOMF_", "",rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM)) 
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM=enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM[order(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$NOM_pval),]
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM$comp = "common PVML_vs_NAWM & PVML_vs_WML & PVML_vs_PVM"
  enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM <- head(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM,20)
```

### Specific WML_vs_CL

``` r
ggplot(enr.specificWML_vs_CL, aes(x=reorder(rownames(enr.specificWML_vs_CL),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "WML_vs_CL") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

### Specific PVML_vs_PVM

``` r
ggplot(enr.specificPVML_vs_PVM, aes(x=reorder(rownames(enr.specificPVML_vs_PVM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "PVML_vs_PVM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

### Specific PVML_vs_NAWM

``` r
ggplot(enr.specificPVML_vs_NAWM, aes(x=reorder(rownames(enr.specificPVML_vs_NAWM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "PVML_vs_NAWM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

### Specific NAGM_vs_GM

``` r
ggplot(enr.specificNAGM_vs_GM, aes(x=reorder(rownames(enr.specificNAGM_vs_GM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "NAGM_vs_GM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

### common PVML_vs_NAWM & PVML_vs_WML & PVML_vs_PVM

``` r
ggplot(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM, aes(x=reorder(rownames(enr.commonPVML_vs_NAWM.PVML_vs_WML.PVML_vs_PVM),-rank(NOM_pval)), y=-log10(NOM_pval),fill=-log10(NOM_pval))) + geom_col(width=.65)+ scale_fill_gradient(low="#F6BDC0",high="#DC1C13") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) + labs(title = "common PVML_vs_NAWM & PVML_vs_WML & PVML_vs_PVM") + xlab("") + coord_flip() + theme_bw()+ theme(axis.text = element_text(face="bold"))
```

![](BulkAnalysisSupervisedMSLesions_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->
