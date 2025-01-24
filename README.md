# RankBat
**RankBat** a computational approach based on ranked genes for batch effect correction (RankBat) which integrates heterogeneous transcriptomic data for identification of MB subtypes across diverse cohorts

## Workflow

## Run steps
1. use RCA() function in step1.R file:
```bash
all_rank_t_genes<- RCA(GSE85217, sampAnnote_GSE85217)
```
2. use RRA() function in step2.R file:
```bash
all_reversed_gp_genes<-RRA(GSE85217, sampAnnote_GSE85217, all_rank_t_genes)
```
3. use LaSelect() function in step3.R file:
```bash
MB_RANK_GP<-LaSelect(GSE85217, sampAnnote_GSE85217, all_rank_t_genes,all_reversed_gp_genes)
```
4. use MBS() function in step4.R file:
```bash
myMat<- (GSE21140, MB_RANK_GP, sampAnnote_GSE21140)
```
