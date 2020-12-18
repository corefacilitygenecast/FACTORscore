timestart<-Sys.time()
cat("###################################################",sep = "\n")
cat("# Author     : Fang Lu                            #",sep = "\n")
cat("# E-mail     : fang.lu@genecast.com.cn            #",sep = "\n")
cat("# Department : Bioinformatics Core Facility       #",sep = "\n")
cat("# Date       : Tur Oct 20 2020                    #",sep = "\n")
cat("# Version    : V1.0                               #",sep = "\n")
cat("# You are using the program scripted by Fang Lu   #",sep = "\n")
cat("# Please let me know if you have any suggestions  #",sep = "\n")
cat("# Thank you !                                     #",sep = "\n")
cat("# Best wishes !                                   #",sep = "\n")
cat("###################################################",sep = "\n")

cat("R packages:",sep = "\n")
cat("   R                    : 3.5.1", sep = "\n")
cat("   getopt               : 1.20.2", sep = "\n")
cat("   ConsensusClusterPlus : 1.46.0", sep = "\n")
cat("   ggsci                : 2.9", sep = "\n")
cat("   ggplot2              : 3.3.2",sep = "\n")
cat("   ComplexHeatmap       : 1.18.1", sep = "\n")
cat("   circlize             : 0.4.9", sep = "\n")
cat("   ggpubr               : 0.2", sep = "\n")
cat("   tidyverse            : 1.3.0", sep = "\n")
cat("   FactoMineR           : 1.41", sep = "\n")
cat("   factoextra           : 1.0.5", sep = "\n")
cat("   survminer            : 0.4.4.999", sep = "\n")
cat("   survival             : 3.1.8", sep = "\n")
cat("   broom                : 0.5.6", sep = "\n")
cat("   stringr              : 1.4.0", sep = "\n")


cat("--------------------------------------",sep = "\n")
cat(paste("start at :",timestart,sep = ""),sep = "\n")

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(ConsensusClusterPlus))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(stringr))


spec = matrix(c(
    'help',    'h', 0, 'logical',   'this help',
    'input',   'i', 1, 'character', 'input matrix',
    'time',    't', 1, 'character', 'time and status file',
    'fpkm',    'f', 1, 'character', 'T: RNA-seq count by gene; TR: fpkm/tpm by transcript; C: other data',
    'outdir',  'o', 1, 'character', 'pathway of output files',
    'prefix',  'p', 1, 'character', 'prefix of the result'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
if(!is.null(opt$h) || is.null(opt$input) || is.null(opt$time) || is.null(opt$outdir) || is.null(opt$prefix)){
        cat(paste(getopt(spec,usage = T)))
        cat("输入文件格式：
        左上角必填：sample
        CNV： 行名为samples, 列名为genes；
        RNA： （1）以基因为单位做定量：行名为samples,列名为genes;
              （2）以转录本为单位做定量：行名为genes，列名为samples
        生存时间格式：sample    OS  re", sep = "\n")
        q(status=1)
}

dir.create(opt$outdir)
mat <- read.csv(opt$i, header = TRUE, sep="\t", check.names=F, encoding="UTF-8")

if(opt$f == "TR"){
    file.rename(opt$i,paste(opt$i,".old",sep=""))
    mat_out <- mat %>% group_by(sample) %>% summarise_all(sum)
    mat_out <- as.data.frame(mat_out)
    rownames(mat_out)<-mat_out[,1]
    mat <- mat_out[,-1]
    mat_tmp <- as.data.frame(t(mat))
    mat_tmp <- mat_tmp[which(rowSums(mat_tmp) > 0),  , drop=FALSE]
    mat_tmp <- mat_tmp[, which(colSums(mat_tmp) > 0), drop=FALSE]
    sample <- data.frame(sample = rownames(mat_tmp))
    mat_out <- cbind(sample,mat_tmp)    
    write.table(mat_out, opt$i, quote=FALSE, sep = "\t",row.names=FALSE)
    data <- as.data.frame(mat)  
}else{
    file.rename(opt$i,paste(opt$i,".old",sep=""))
    rownames(mat)<-mat[,1]
    mat <- mat[,-1]
    mat <- mat[which(rowSums(mat) > 0),  , drop=FALSE]
    mat <- mat[, which(colSums(mat) > 0), drop=FALSE]
    sample <- data.frame(sample = rownames(mat))
    mat_out <- cbind(sample,mat)    
    write.table(mat_out, opt$i, quote=FALSE, sep = "\t",row.names=FALSE)
    data <- as.data.frame(t(mat)) 
}

fpkm <- data
colname <- colnames(fpkm)
fpkm <- as.matrix(fpkm)
fpkm <- apply(fpkm,1,scale)
fpkm <- as.data.frame(t(fpkm))
colnames(fpkm) <- colname
data <- as.matrix(fpkm)
print(head(data))
write.table(data,file="Normalized.data.xls",sep="\t",quote=F,col.names=NA,row.names=T)

setwd(opt$o)
ConsensusClusterPlus(data,maxK=6,reps=50,pItem=0.8,pFeature=1,title="euclidean_pam",clusterAlg="pam",distance="euclidean",plot="pdf",writeTable=TRUE,verbose=F)

get_group <- function(k=2){
    group_new <- c()
    group_info <- read.table(paste(opt$o, "/euclidean_pam", "/euclidean_pam.k=", k, ".consensusClass.csv", sep = ""), header = FALSE, sep = ",", check.names=F, stringsAsFactors = FALSE, encoding="UTF-8")
    for(i in 1:nrow(group_info)){
        group_new <- rbind(group_new,c(gsub(pattern='"',replacement='',group_info[i,1]),group_info[i,2]))
    }
    group_new <- as.matrix(group_new)
    colnames(group_new) <- c("sample","groups")
    write.table(group_new,paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""),row.names = FALSE,quote=FALSE,sep = "\t")
}

complexheatmap <- function(input, group, k){
    diff=read.csv(input,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",encoding='UTF-8',check.names=F,row.names=1)
    diff<-as.data.frame(t(diff))
    
    if( opt$f == "T" | opt$f == "TR" ){
        diff <- log10(diff + 1)
    }else{
        diff <- diff
    }
    
    group_data<-as.data.frame(read.table(group,header=1,sep="\t",encoding='UTF-8',check.names=F))
    sample<-as.character(group_data$sample) 
    newdata<-(diff[,sample])
    #data<-as.matrix(newdata[,c(3:dim(newdata)[2])])
    data<-newdata
    group_data<-group_data[order(group_data$groups),]
    fpkm<-data
    compares <- unique(group_data$groups)
    #####################################

    fpkm <- subset(fpkm,select=as.character(group_data$sample))
    #fpkm <- fpkm[rownames(geneinfo),]
    #fpkm <- fpkm[rowSums(fpkm)>0,]
    fpkm <- na.omit(fpkm)
    #fpkm <- log10(fpkm+1)
    colname <- colnames(fpkm)
    fpkm <- as.matrix(fpkm)
    fpkm <- apply(fpkm,1,scale)
    fpkm <- as.data.frame(t(fpkm))
    #fpkm<-as.data.frame(fpkm)
    colnames(fpkm) <- colname
    fpkm <- fpkm[apply(fpkm,1,function(y) !sum(is.na(y))),]
    fpkm <- fpkm[,apply(fpkm,2,function(y) !sum(is.na(y)))]
    write.table(fpkm,file=paste(opt$o,"/",opt$p,'_data_nrom.txt',sep=''),sep='\t',quote=F,row.names=T)
    #####################################
    data<-as.data.frame(t(fpkm))
    number<-dim(data)[2]
    matchinfo<-match(rownames(data),group_data$sample,nomatch=0)
    data$Group<-group_data[matchinfo,"groups"]
    data<-data[order(data$Group),]
    horder<-"0"
    for(i in 1:length(compares)){
        data1<-data[which(data$Group==compares[i]),]
        data1<-as.matrix(data1[,1:number,drop=FALSE])
        h1<-Heatmap(data1)
        order1<-as.matrix(row_order(h1)[[1]])
        rownames(order1)<-rownames(data1)
        if(max(horder)==0){
            horder<-order1
        }else{
            order2<-apply(order1,2,function(x) (x+max(horder)))
            if(dim(data1)[1]==1){
                order2<-as.matrix(order2)
                rownames(order2)<-rownames(data1)
            }
            horder<-c(horder,order2)
        }
    }

    fpkm<-as.data.frame(t(as.matrix(data[,1:number,drop=FALSE])))

    #####################################

    color_n = length(compares)
    if(color_n==1){
        col_nw = "darkorange"
    }else if (color_n == 2){
        col_nw = c("#003399","#CC0033")
    }else{
        col_nw = pal_jama("default",alpha=0.8)(color_n)
    }
    group.assign<-setNames(col_nw,compares)
    match.id<-match(colnames(fpkm),group_data$sample,nomatch=0)
    df<-group_data[match.id,"groups",drop=FALSE]
    ha<-HeatmapAnnotation(df=df,col=list(groups=group.assign),annotation_height = unit(0.5, "cm"))

    col_fun = colorRamp2(c(-3,0,3), c("green","white", "red"))
    p<-Heatmap(fpkm, column_title = "Heatmap", column_title_gp = gpar(fontsize = 25, fontface = "bold", col="black"),row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 1),top_annotation=ha,cluster_rows=T,cluster_columns=F,show_column_names = F,show_row_names = F,column_order=horder, heatmap_legend_param = list( title = "Factor", at = c(-3, -1 ,1, 3)),col=col_fun)

    out_pdf<-paste(opt$o,"/euclidean_pam/",opt$prefix,".k=",k,".heatmap.pdf",sep="")
    pdf(out_pdf,width=15,height=30)
    draw(p)
    dev.off()
}

wilcox_test <- function(df,group){
    merged_data <-  merge(df, group,by="sample")
    merged_data[is.na(merged_data)] <- 0
    for_wilcox <- subset(merged_data, select = -sample)
    sig_gene <- c()
    for(gene in 1:(ncol(for_wilcox)-1))
    {
        gene_wilcox <- for_wilcox[,c(colnames(for_wilcox)[gene],"groups")]
        gene_wilcox <- rename(gene_wilcox, c("gene" = colnames(gene_wilcox)[1]))
        out <- compare_means(gene~groups, gene_wilcox, method = "wilcox.test", exact=FALSE ,correct=FALSE)
        if(out$p < 0.05){
            sig_gene <- c(sig_gene, colnames(for_wilcox)[gene])
        }
    }
    return(sig_gene)
}
KW_test <- function(df,group){
    merged_data <-  merge(df, group, by="sample")
    merged_data[is.na(merged_data)] <- 0
    for_KW <- subset(merged_data, select = -sample)
    sig_gene <- c()
    for(gene in 1:(ncol(for_KW)-1))
    {
        gene_KW <- for_KW[,c(colnames(for_KW)[gene],"groups")]
        gene_KW <- rename(gene_KW, c("gene" = colnames(gene_KW)[1]))
        out <- compare_means(gene~groups, gene_KW, method = "kruskal.test", exact=FALSE ,correct=FALSE)
        if(out$p < 0.05){
            sig_gene <- c(sig_gene, colnames(for_KW)[gene])
        }
    }
    return(sig_gene)
}

get_GeneGroup <- function(k=2,g=2){
    GeneGroup_new <- c()
    GeneGroup_info <- read.table(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/", "/euclidean_pam.k=", g, ".consensusClass.csv", sep = ""), header = FALSE, sep = ",", check.names=F, stringsAsFactors = FALSE, encoding="UTF-8")
    for(i in 1:nrow(GeneGroup_info)){
        GeneGroup_new <- rbind(GeneGroup_new,c(gsub(pattern='"',replacement='',GeneGroup_info[i,1]),GeneGroup_info[i,2]))
    }
    GeneGroup_new <- as.matrix(GeneGroup_new)
    colnames(GeneGroup_new) <- c("sample","groups")
    write.table(GeneGroup_new,paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,".group",sep = ""),row.names = FALSE,quote=FALSE,sep = "\t")
}

complexheatmap_v2 <- function(input, sample_group, gene_group, k, g_k){
    diff = read.csv(input, fill=T, header=TRUE, stringsAsFactors=FALSE, sep="\t", encoding='UTF-8', check.names=F, row.names=1)
    
    if( opt$f == "T" | opt$f == "TR"){
        diff <- log10(diff + 1)
    }else{
        diff <- diff
    }

    group_data <- as.data.frame(read.table(sample_group, header=1, sep="\t", encoding='UTF-8', check.names=F))
    group_data <- group_data[order(group_data$groups),]
    gene_order <- read.csv(gene_group, header = TRUE, stringsAsFactors = FALSE, sep = "\t", encoding = 'UTF-8', check.names = F)
    gene_order <- gene_order[order(gene_order$groups),]
    gene_rank <- as.character(gene_order$sample)

    diff <- as.data.frame(t(diff))
    sample <- as.character(group_data$sample)
    newdata <- (diff[,sample])
    data <- newdata
    fpkm <- data

    #####################################

    fpkm <- subset(fpkm, select=as.character(group_data$sample))
    
    fpkm <- na.omit(fpkm)
    #fpkm <- log10(fpkm+1)
    colname <- colnames(fpkm)
    fpkm <- as.matrix(fpkm)
    fpkm <- apply(fpkm,1,scale)
    fpkm <- as.data.frame(t(fpkm))
    colnames(fpkm) <- colname
    fpkm <- fpkm[apply(fpkm,1,function(y) !sum(is.na(y))),]
    fpkm <- fpkm[,apply(fpkm,2,function(y) !sum(is.na(y)))]

    match_genegroup_to_matrix <- match(rownames(data), gene_order$sample, nomatch=0)
    tmp <- data[match_genegroup_to_matrix, ,drop = FALSE]

    data<-as.data.frame(t(fpkm))
    number<-dim(data)[2]
    matchinfo<-match(rownames(data), group_data$sample, nomatch=0)
    data$Group<-group_data[matchinfo,"groups"]
    data<-data[order(data$Group),]
    horder<-"0"
    compares <- unique(group_data$groups)
    for(i in 1:length(compares)){
        data1<-data[which(data$Group==compares[i]),]
        data1<-as.matrix(data1[,1:number,drop=FALSE])
        h1<-Heatmap(data1)
        order1<-as.matrix(row_order(h1)[[1]])
        rownames(order1)<-rownames(data1)
        if(max(horder)==0){
            horder<-order1
        }else{
            order2<-apply(order1,2,function(x) (x+max(horder)))
            if(dim(data1)[1]==1){
                order2<-as.matrix(order2)
                rownames(order2)<-rownames(data1)
            }
            horder<-c(horder,order2)
        }
    }

    fpkm<-as.data.frame(t(as.matrix(data[,1:number,drop=FALSE])))

    color_n = length(compares)
    if(color_n==1){
        col_nw = "darkorange"
    }else if (color_n == 2){
        col_nw = c("#003399","#CC0033")
    }else{
        col_nw = pal_jama("default",alpha=0.8)(color_n)
    }
    group.assign<-setNames(col_nw,compares)
    match.id<-match(colnames(fpkm),group_data$sample,nomatch=0)
    df<-group_data[match.id,"groups",drop=FALSE]
    ha<-HeatmapAnnotation(df=df,col=list(groups=group.assign),annotation_height = unit(0.5, "cm"))

    compares_gene <- unique(gene_order$groups)
    color_gene = length(compares_gene)
    if(color_gene==1){
        col_gene = "darkorange"
    }else if (color_gene == 2){
        col_gene = c("#3CB371","#FFA500")
    }else{
        col_gene = pal_jama("default",alpha=0.8)(color_gene)
    }
    gene_group.assign <- setNames(col_gene,compares_gene)
    match.id2 <- match(rownames(fpkm), gene_order$sample, nomatch = 0)
    df_gene <- gene_order[match.id2,"groups", drop=FALSE]
    h_gene <- rowAnnotation(df = df_gene, col = list(groups = gene_group.assign), annotation_height = unit(0.5, "cm",),annotation_legend_param = list(groups = list(title = "GeneGroups")))
  
    col_fun = colorRamp2(c(-3,0,3), c("green","white", "red"))
    p<-Heatmap(fpkm, column_title = "Heatmap", column_title_gp = gpar(fontsize = 25, fontface = "bold", col="black"),row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 1),top_annotation=ha, cluster_row = F,row_order = gene_rank,cluster_columns=F,show_column_names = F,show_row_names = F,column_order=horder, heatmap_legend_param = list( title = "Factor", at = c(-3, -1 ,1, 3)),col=col_fun)
    
    out_pdf<-paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","k=",k,"geneK=",g_k,".heatmap_new.pdf",sep="")
    pdf(out_pdf,width=15,height=30)
    draw(h_gene+p)
    dev.off()
}

cox_regression <- function(timeinfo, groupinfo, pcainfo, dir){
    timeinfo  <- read.table(timeinfo, sep="\t", header=T, encoding='UTF-8', check.names=F)
    rownames(timeinfo) <- timeinfo$sample
    colnames(timeinfo) <- c("sample", "OS", "re")
    groupinfo <- read.table(groupinfo, sep="\t", header=T, encoding='UTF-8', check.names=F)
    rownames(groupinfo) <- groupinfo$sample
    pca_info<-read.table(pcainfo, sep = "\t", header = T, encoding = 'UTF-8', check.names = F)
    rownames(pca_info) <- pca_info$sample

    groupinfo$groups <- factor(groupinfo$groups)
    
    tmpdata <- merge(timeinfo, groupinfo, by="sample")
    newdata <- merge(tmpdata, pca_info, by="sample")

    data <- na.omit(newdata)
    #data <- data[order(data$groups, decreasing= T),]
    write.table(data, file=(paste(dir,"/ConsensusGroup_PC1_svv.txt", sep="")), sep="\t", quote=F, row.names=F)

    data <- select(data,-sample)
    fit=eval(parse(text=paste("survfit(Surv(OS,re) ~ ", "groups", ", data = data)", sep="")))
    fit3=eval(parse(text=paste("coxph(Surv(OS,re) ~ .", ", data = data)", sep=""))) 
    write.table(summary(fit)$table, file=(paste(dir,"/ConsensusGroup_PC1_surv_summary.txt", sep="")), sep="\t",quote=F,row.names=T)

    #tidy the coxph analysis result and plot forest
    coef<-as.data.frame(tidy(fit3))
    write.table(coef,file=(paste(dir,"/ConsensusGroup_PC1_summary.txt", sep="")), sep="\t",quote=F,row.names=F)
    pdf(paste(dir,"/ConsensusGroup_PC1_forest.pdf",sep=""),onefile=FALSE,width=30,height=10)
    pforest<-ggforest(fit3, fontsize = 2, cpositions = c(0.0005, 0.15, 0.28),data=data)
    print(pforest)
    dev.off()

    #extract coxph analysis result
    fit4<-summary(fit3)
    sctest<-as.matrix(fit4$sctest)
    colnames(sctest)<-"sctest"
    sctest_p <- sctest[3,1]
    waldtest<-as.matrix(fit4$waldtest)
    colnames(waldtest)<-"waldtest"
    wald_p <- waldtest[3,1]
    logtest<-as.matrix(fit4$logtest)
    colnames(logtest)<-"logtest"
    logrank_p <- logtest[3,1]

    head<-c(rep(0,8))
    head<-c("Item","HR[exp(coef)]","HR[exp(-coef)]","95% CI lower","95% CI upper","Score_logrank_test_pval","Wald_test_pval","Likelihood_ratio_test_pval")

    item_num <- dim(fit4$conf.int)[1]
    variable_name <- rownames(fit4$conf.int)
    out<-c()
    for(n in 1:item_num){
        TT <- c(variable_name[n],fit4$conf.int[n,1], fit4$conf.int[n,2], fit4$conf.int[n,3], fit4$conf.int[n,4], sctest_p, wald_p, logrank_p)
        out <- rbind(out,TT)
    }
    result <- as.data.frame(out)
    colnames(result) <- head
    outfile<-paste(dir, "/ConsensusGroup_PC1_survival.HR.pval.xls", sep="")
    write.table(result, file=outfile, sep="\t", row.names=F, quote=F)

    # plot survival graph
    num_type=dim(coef)[1]
    #pvalue<-coef$p.value[num_type]
    pvalue <- as.character(sctest_p)
    nn<-round(as.numeric(pvalue),4)
    HR<-round(max(as.data.frame(fit4$conf.int)[,"exp(coef)"]),4)
    #HR<-as.data.frame(fit4$conf.int)[,"exp(coef)"][num_type]
    #hr <- as.character(fit4$conf.int[1])
    #HR <-round(as.numeric(hr),4)
    pv<-paste("P = ",nn,"\nHR = ",HR,sep="")
    #pv <- paste("P = ",nn)
    name <- paste("groups", "_Survival_Curves", sep="")

    str = gsub("sig=", "", rownames(as.matrix(fit$strata)), fixed = TRUE)[1]
    str2 = gsub("sig=", "", rownames(as.matrix(fit$strata)), fixed = TRUE)[2]
    if(is.null(opt$quantile)){
        legend_lable=c(paste(str," (n=",as.matrix(fit$n)[1],")\n",sep=""),paste(str2," (n=",as.matrix(fit$n)[2],")\n",sep=""))
        if(length(rownames(as.matrix(fit$strata)))>=3){
            for(j in 3:length(rownames(as.matrix(fit$strata)))){
                name<-gsub("sig=", "", rownames(as.matrix(fit$strata)), fixed = TRUE)[j]
                legend_temp=paste(name," (n=",as.matrix(fit$n)[j],")\n",sep="")
                legend_lable=c(legend_lable,legend_temp)
            }	
        }
        if(length(rownames(as.matrix(fit$strata))) == 2){
        k2<-ggsurvplot(fit,data = data, censor.shape="|", censor.size = 4,pval.size = 8, font.x = c(25, "bold.italic", "black"), font.y = c(25, "bold.italic", "black"), font.tickslab = c(10, "plain", "black"), legend.title = "groups",font.legend = c(10,"bold.italic", "black"),pval=pv,legend = "right",palette = c("#003399","#CC0033"), legend.labs = legend_lable, risk.table=T, surv.median.line = "hv")
        }else{
        k2<-ggsurvplot(fit,data = data, censor.shape="|", censor.size = 4,pval.size = 8,font.x = c(25, "bold.italic", "black"),font.y = c(25, "bold.italic", "black"),font.tickslab = c(10, "plain", "black"),legend.title = "groups",font.legend = c(10,"bold.italic", "black"),pval=pv,legend = "right",palette = "jama",legend.labs = legend_lable,risk.table=T,surv.median.line = "hv")
        }
    }else{
        if(str_detect(str,">")){
            legend_lable=c(paste(opt$quantile,":",str," (n=",as.matrix(fit$n)[1],", median PFS=",me_pfsl,"m)\n",sep=""),paste(opt$quantile,":",str2," (n=",as.matrix(fit$n)[2],", median PFS=",me_pfss,"m)\n",sep=""))
            k2<-ggsurvplot(fit, data = data, censor.shape="|", censor.size = 4, pval.size = 8, font.x = c(25, "bold.italic", "black"), font.y = c(25, "bold.italic", "black"), font.tickslab = c(25, "plain", "black"), legend.title = "groups", font.legend = c(23,"bold.italic", "black"), palette = c("#003399","#CC0033"), legend.labs = legend_lable, pval=pv, legend = "right",, surv.median.line = "hv")
        }else{
            legend_lable=c(paste(opt$quantile,":",str," (n=",as.matrix(fit$n)[1],", median PFS=",me_pfss,"m)\n",sep=""),paste(opt$quantile,":",str2," (n=",as.matrix(fit$n)[2],", median PFS=",me_pfsl,"m)\n",sep=""))
            k2<-ggsurvplot(fit,data = data, censor.shape="|", censor.size = 4,pval.size = 8,font.x = c(25, "bold.italic", "black"),font.y = c(25, "bold.italic", "black"),font.tickslab = c(25, "plain", "black"),legend.title = "groups",font.legend = c(23,"bold.italic", "black"),palette = c("#003399","#CC0033"),legend.labs = legend_lable,pval=pv,legend = "right",, surv.median.line = "hv")
        }
    }

    if(length(rownames(as.matrix(fit$n)))<=6 && length(rownames(as.matrix(fit$n)))!=1){
        pdf(paste(dir,"/ConsensusGroup_PC1_svv.pdf",sep=""),onefile=FALSE,width=20)
        print(k2)
        dev.off()
    }
}
SVV_maxstat <- function(timeinfo, scoreinfo, outdir, factor){
    timeinfo<-read.table(timeinfo,sep="\t",header=T,encoding='UTF-8',check.names=F)
    diff<-read.table(factor,sep="\t",header=F,encoding='UTF-8',check.names=F)
    exp<-as.data.frame(t(read.table(scoreinfo,sep="\t",header=T,encoding='UTF-8',check.names=F,row.names=1)))
    #exp<-log(exp1+1)
    #rownames(exp)<-exp$Patient

    len=dim(diff)[1]
    for(i in 1:len){
        id=as.character(diff[i,])
        if(str_detect(id,"[+]")){
            allsig<-unlist(strsplit(id, "+", fixed = TRUE))
            gene<-as.data.frame(t(exp[allsig,,drop=FALSE]))
        }else{
            gene<-as.data.frame(t(exp[which(rownames(exp)==id),]))
        }
        if(str_detect(id,"Adjuvant_Therapy")){
            gene1<-gene[which(gene$Adjuvant_Therapy!=""),]
            gene1$Adjuvant_Therapy<-as.character(gene1$Adjuvant_Therapy)
            gene<-gene1
        }
        gene<-gene[which(gene!=""),,drop=F]
        #gene<-na.omit(gene)
        commonsample<-intersect(rownames(gene),timeinfo$sample)
        gene_tmp<-as.data.frame(gene[commonsample,,drop=FALSE])
        #commonsample1<-intersect(rownames(gene_tmp),groupinfo$sample)
        #gene_new<-as.data.frame(gene_tmp[commonsample1,,drop=FALSE])
        #match.id1<-match(rownames(gene_new),groupinfo$sample,nomatch=0)
        #gene_new$Group=groupinfo[match.id1,"groups"]
        gene_new <- gene_tmp
        match.id2<-match(rownames(gene_new),timeinfo$sample,nomatch=0)
        gene_new$OS<-timeinfo[match.id2,"OS"]
        gene_new$re<-timeinfo[match.id2,"re"]
        if(str_detect(id,"[+]")){
            info<-strsplit(id,"+",fixed=T)[[1]]
            data<-gene_new
            data$patient<-rownames(gene_new)
            data<-na.omit(data[c("patient",info,"Group","OS","re")])
            fit=eval(parse(text=paste("survfit(Surv(OS,re) ~ ",id,",data = data)",sep="")))
            fit3=eval(parse(text=paste("coxph(Surv(OS,re) ~ ",id,",data = data)",sep="")))
        }else{
            #colnames(gene_new)<-c("gene","Group","OS","re")
            colnames(gene_new)<-c("gene","OS","re")
            gene_new<-gene_new[order(gene_new$gene,decreasing= T),]
            if(is.null(opt$quantile)){
                gene_new$tmpsig=as.numeric(as.character(gene_new$gene))
                res.cut <- tryCatch(surv_cutpoint(gene_new, time = "OS", event = "re",variables = "tmpsig"),error=function(e){cat("Can not get best cutoff,because of error:",conditionMessage(e),"\n\n")})
                if(is.null(res.cut)){
                    next
                }
                cutoff=summary(res.cut)[1,1]
                gene_new$sig=gene_new$tmpsig>cutoff
                me_pfsl<-median(gene_new[which(gene_new$sig=="TRUE"),]$OS)
                n_pfsl<-dim(gene_new[which(gene_new$sig=="TRUE"),])[1]
                me_pfss<-median(gene_new[which(gene_new$sig=="FALSE"),]$OS)
                n_pfss<-dim(gene_new[which(gene_new$sig=="FALSE"),])[1]
                gene_new[which(gene_new$sig=="TRUE"),]$sig<-paste("group>",signif(cutoff,3),sep="")
                gene_new[which(gene_new$sig=="FALSE"),]$sig<-paste("group<=",signif(cutoff,3),sep="")
                data<-gene_new
            }else{
                medienfpkm=quantile(gene_new$gene,opt$quantile)
                gene_new$sig=gene_new$gene>medienfpkm
                me_pfsl<-median(gene_new[which(gene_new$sig=="TRUE"),]$OS)
                n_pfsl<-dim(gene_new[which(gene_new$sig=="TRUE"),])[1]
                me_pfss<-median(gene_new[which(gene_new$sig=="FALSE"),]$OS)
                n_pfss<-dim(gene_new[which(gene_new$sig=="FALSE"),])[1]
                gene_new[which(gene_new$sig=="TRUE"),]$sig<-paste("group>",signif(medienfpkm,3),sep="")
                gene_new[which(gene_new$sig=="FALSE"),]$sig<-paste("group<=",signif(medienfpkm,3),sep="")
                data<-gene_new
            }
            data$patient<-rownames(gene_new)
            data<-na.omit(data)
            fit <- survfit(Surv(OS,re) ~sig,data = data)
            fit3<-coxph(Surv(OS,re) ~sig,data = data)
        }
        write.table(data,file=(paste(outdir,"/Score","_svv.txt",sep="")),sep="\t",quote=F,row.names=F)
        coef<-as.data.frame(tidy(fit3))
        write.table(coef,file=(paste(outdir,"/Score","_summary.txt",sep="")),sep="\t",quote=F,row.names=F)
        fit4<-summary(fit3)
        pdf(paste(outdir,"/Score","_forest.pdf",sep=""),onefile=FALSE,width=32,height=12)
        #pforest<-ggforest(fit3, fontsize = 0.7, data=data)
        pforest<-ggforest(fit3, fontsize = 2, cpositions = c(0.00, 0.05, 0.28),data=data)
        print(pforest)
        dev.off()
        #fit4<-summary(fit3)
        sctest<-as.matrix(fit4$sctest)
        colnames(sctest)<-"sctest"
        waldtest<-as.matrix(fit4$waldtest)
        colnames(waldtest)<-"waldtest"
        logtest<-as.matrix(fit4$logtest)
        colnames(logtest)<-"logtest"
        TT<-c(rep(0,8))
        if(dim(waldtest)[1]!=3){
        TT<-c(id,fit4$conf.int[1],fit4$conf.int[2],fit4$conf.int[3],fit4$conf.int[4],sctest[3,1],"NA",logtest[3,1])
        }else{
        TT<-c(id,fit4$conf.int[1],fit4$conf.int[2],fit4$conf.int[3],fit4$conf.int[4],sctest[3,1],waldtest[3,1],logtest[3,1])
        }
        head<-c(rep(0,8))
        head<-c("Item","HR[exp(coef)]","HR[exp(-coef)]","95% CI lower","95% CI upper","Score_logrank_test_pval","Wald_test_pval","Likelihood_ratio_test_pval")
        result<-data.frame(value=rep(0,8))
        result$value<-TT
        
        result1<-t(result)
        colnames(result1)<-head
        result1<-as.data.frame(result1)
            outfile<-paste(outdir,"/Score","_survival.HR.pval.xls",sep="")
            write.table(result1,file=outfile,sep="\t",row.names=F,quote=F)
        #pdf(paste(opt$outdir,"/",opt$p,"_",id,"_svv.pdf",sep=""),onefile=FALSE,width=15)
        num_type=dim(coef)[1]
        #pvalue<-coef$p.value[num_type]
        pvalue <- as.character(result1$Score_logrank_test_pval)
            nn<-round(as.numeric(pvalue),4)
        HR<-max(as.data.frame(fit4$conf.int)[,"exp(coef)"])
        #HR<-as.data.frame(fit4$conf.int)[,"exp(coef)"][num_type]
        #hr <- as.character(fit4$conf.int[1])

        #HR <-round(as.numeric(hr),4)
        pv<-paste("P = ",nn,"\nHR = ",HR,sep="")
        name <-paste(id,"_Survival_Curves",sep="")
        if(dim(waldtest)[1]==3){
        str=gsub("sig=", "", rownames(as.matrix(fit$strata)), fixed = TRUE)[1]
        str2=gsub("sig=", "", rownames(as.matrix(fit$strata)), fixed = TRUE)[2]
        if(is.null(opt$quantile)){
            legend_lable=c(paste(str," (n=",as.matrix(fit$n)[1],")\n",sep=""),paste(str2," (n=",as.matrix(fit$n)[2],")\n",sep=""))
            if(length(rownames(as.matrix(fit$strata)))>=3){
                for(j in 3:length(rownames(as.matrix(fit$strata)))){
                    name<-gsub("sig=", "", rownames(as.matrix(fit$strata)), fixed = TRUE)[j]
                    legend_temp=paste(name," (n=",as.matrix(fit$n)[j],")\n",sep="")
                    legend_lable=c(legend_lable,legend_temp)
                }	
            }
        if(length(rownames(as.matrix(fit$strata)))<=9){
            k2<-ggsurvplot(fit,data = data, censor.shape="|", censor.size = 4,pval.size = 10,font.x = c(25, "bold.italic", "black"),font.y = c(25, "bold.italic", "black"),font.tickslab = c(10, "plain", "black"),legend.title = id,font.legend = c(10,"bold.italic", "black"),pval=pv,legend = "right",legend.labs = legend_lable,risk.table=T,palette = c("#003399","#CC0033"),surv.median.line = "hv")
        } #palette = "NEJM",palette = "jama",
        }else{
            if(str_detect(str,">")){
                legend_lable=c(paste(opt$quantile,":",str," (n=",as.matrix(fit$n)[1],", median PFS=",me_pfsl,"m)\n",sep=""),paste(opt$quantile,":",str2," (n=",as.matrix(fit$n)[2],", median PFS=",me_pfss,"m)\n",sep=""))
                k2<-ggsurvplot(fit,data = data, censor.shape="|", censor.size = 4,pval.size = 10,font.x = c(25, "bold.italic", "black"),font.y = c(25, "bold.italic", "black"),font.tickslab = c(25, "plain", "black"),legend.title = id,font.legend = c(23,"bold.italic", "black"),palette = c("#003399","#CC0033"),legend.labs = legend_lable,pval=pv,legend = "right",surv.median.line = "hv")
            }else{
                legend_lable=c(paste(opt$quantile,":",str," (n=",as.matrix(fit$n)[1],", median PFS=",me_pfss,"m)\n",sep=""),paste(opt$quantile,":",str2," (n=",as.matrix(fit$n)[2],", median PFS=",me_pfsl,"m)\n",sep=""))
                k2<-ggsurvplot(fit,data = data, censor.shape="|", censor.size = 4,pval.size = 10,font.x = c(25, "bold.italic", "black"),font.y = c(25, "bold.italic", "black"),font.tickslab = c(25, "plain", "black"),legend.title = id,font.legend = c(23,"bold.italic", "black"),palette = c("#003399","#CC0033"),legend.labs = legend_lable,pval=pv,legend = "right",surv.median.line = "hv")
            }
        }
        if(length(rownames(as.matrix(fit$strata)))<=9){
            pdf(paste(outdir,"/Score","_svv.pdf",sep=""),onefile=FALSE,width=15)
            print(k2)
            dev.off()
        }
        }
    }
}

for(k in 2:6){
    get_group(k)
    complexheatmap(opt$i, paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""), k)
    dir.create(paste(opt$outdir,"/euclidean_pam/k=",k,sep=""))
    mat_k <- as.data.frame(read.csv(opt$i, header = TRUE, sep="\t", check.names=F, encoding="UTF-8"))
    group_k <- as.data.frame(read.table(paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""),header=1,sep="\t",encoding='UTF-8',check.names=F))
    if(k == 2){
        sig_gene <- wilcox_test(mat_k,group_k)
        mat_sigGene <- subset(mat_k,select=c("sample",sig_gene))
        mat_sigGene_raw <- mat_sigGene
        write.table(sig_gene,paste(opt$o,"/euclidean_pam/k=",k,"/wilcox_sig_gene.txt",sep = ""),col.names = FALSE, row.names = FALSE,quote=FALSE,sep = "\t")
        write.table(mat_sigGene,paste(opt$o,"/euclidean_pam/k=",k,"/wilcox_sig_gene.mat",sep = ""),row.names = FALSE,quote=FALSE,sep = "\t")
        
        rownames(mat_sigGene)<-mat_sigGene[,1]
        mat_sigGene<-mat_sigGene[,-1]
        data_k <-as.data.frame(t(mat_sigGene))
        fpkm_k <- data_k
        colname <- colnames(fpkm_k)
        fpkm_k <- as.matrix(fpkm_k)
        fpkm_k <- apply(fpkm_k,1,scale)
        fpkm_k <- as.data.frame(t(fpkm_k))
        colnames(fpkm_k) <- colname
        data_k <-as.matrix(t(fpkm_k))
        write.table(data_k,file=paste(opt$o,"/euclidean_pam/k=",k,"/Normalized.data.xls",sep=""),sep="\t",quote=F,col.names=NA,row.names=T)
        path_k <- paste(opt$o,"/euclidean_pam/k=",k,sep="")
        setwd(path_k)
        ConsensusClusterPlus(data_k,maxK=6,reps=50,pItem=0.8,pFeature=1,title="euclidean_pam",clusterAlg="pam",distance="euclidean",plot="pdf",writeTable=TRUE,verbose=F)
        for(g in 2:6){
            get_GeneGroup(k,g)
            complexheatmap_v2(paste(opt$o,"/euclidean_pam/k=",k,"/wilcox_sig_gene.mat",sep = ""),paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,".group",sep = ""), k, g)
            dir.create(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,sep = ""))
            GeneGroup_info <- read.table(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,".group",sep = ""), header = TRUE, sep = "\t")
            PC_group <- c()
            for(gg in 1:g){
                Gene_gg <- as.character(GeneGroup_info[which(GeneGroup_info$groups == gg),]$sample)
                mat_gg <- subset(mat_sigGene_raw,select=c("sample",Gene_gg))
                write.table(mat_gg,file=paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".mat",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE)
                rownames(mat_gg)<-mat_gg[,1]
                mat_gg<-mat_gg[,-1]
                gg_pca <- PCA(mat_gg, scale.unit = TRUE, ncp = 2, graph = FALSE)
                gg_PC1 <- get_pca_ind(gg_pca)$coord[,1, drop=FALSE]
                colnames(gg_PC1) <- paste("group",gg,"_PC1",sep ="")
                PC_group <- cbind(PC_group,gg_PC1)
            }
            PC_for_Score <- as.data.frame(PC_group) %>% rownames_to_column('sample')
            write.table(as.data.frame(PC_group) %>% rownames_to_column('sample'), file = paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".PC1",sep = ""), row.names = FALSE, sep = "\t", quote=FALSE)
            cox_regression(opt$t, paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".PC1",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g, sep = ""))
            HR_result <- read.table(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/ConsensusGroup_PC1_survival.HR.pval.xls", sep = ""),header = TRUE, sep = "\t")
            
            fomula<-""
            for(hr in 2:nrow(HR_result)){
                if( fomula == "" ){
                    if( HR_result[hr,2] > 1){
                        fomula <- HR_result[hr,1]
                    }else{
                        fomula <- paste("-", HR_result[hr,1], sep = "")
                    }
                }else{
                    if( HR_result[hr,2] > 1){
                        fomula <- paste(fomula,"+",HR_result[hr,1],sep = "")
                    }else{
                        fomula <- paste(fomula,"-",HR_result[hr,1], sep = "")
                    }
                }  
            }
            eval(parse(text=paste("Score <- PC_for_Score %>% mutate(Score = ",fomula,")",sep="")))
            write.table(Score, file = paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".Score",sep = ""), row.names = FALSE, sep = "\t", quote=FALSE)
            factor <- "Score"
            write.table(factor, file = paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/factor.list",sep = ""), col.names = FALSE, row.names = FALSE,quote = FALSE)
            SVV_maxstat(opt$t, paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".Score",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/factor.list",sep = ""))
        }
    }else{
        sig_gene <- KW_test(mat_k,group_k)
        mat_sigGene <- subset(mat_k,select=c("sample",sig_gene))
        mat_sigGene_raw <- mat_sigGene
        write.table(sig_gene,paste(opt$o,"/euclidean_pam/k=",k,"/KW_sig_gene.txt",sep = ""),col.names = FALSE, row.names = FALSE,quote=FALSE,sep = "\t")
        write.table(mat_sigGene,paste(opt$o,"/euclidean_pam/k=",k,"/KW_sig_gene.mat",sep = ""),row.names = FALSE,quote=FALSE,sep = "\t")

        rownames(mat_sigGene)<-mat_sigGene[,1]
        mat_sigGene<-mat_sigGene[,-1]
        data_k <-as.data.frame(t(mat_sigGene))
        fpkm_k <- data_k
        colname <- colnames(fpkm_k)
        fpkm_k <- as.matrix(fpkm_k)
        fpkm_k <- apply(fpkm_k,1,scale)
        fpkm_k <- as.data.frame(t(fpkm_k))
        colnames(fpkm_k) <- colname
        data_k <-as.matrix(t(fpkm_k))
        write.table(data_k,file=paste(opt$o,"/euclidean_pam/k=",k,"/Normalized.data.xls",sep=""),sep="\t",quote=F,col.names=NA,row.names=T)
        path_k <- paste(opt$o,"/euclidean_pam/k=",k,sep="")
        setwd(path_k)
        ConsensusClusterPlus(data_k,maxK=6,reps=50,pItem=0.8,pFeature=1,title="euclidean_pam",clusterAlg="pam",distance="euclidean",plot="pdf",writeTable=TRUE)
        for(g in 2:6){
            get_GeneGroup(k,g)
            complexheatmap_v2(paste(opt$o,"/euclidean_pam/k=",k,"/KW_sig_gene.mat",sep = ""),paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,".group",sep = ""), k, g)
            dir.create(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,sep = ""))
            GeneGroup_info <- read.table(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,".group",sep = ""), header = TRUE, sep = "\t")
            PC_group <- c()
            for(gg in 1:g){
                Gene_gg <- as.character(GeneGroup_info[which(GeneGroup_info$groups == gg),]$sample)
                mat_gg <- subset(mat_sigGene_raw,select=c("sample",Gene_gg))
                write.table(mat_gg,file=paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".mat",sep = ""),sep = "\t", quote = FALSE, row.names = FALSE)
                rownames(mat_gg)<-mat_gg[,1]
                mat_gg<-mat_gg[,-1]
                gg_pca <- PCA(mat_gg, scale.unit = TRUE, ncp = 2, graph = FALSE)
                gg_PC1 <- get_pca_ind(gg_pca)$coord[,1, drop=FALSE]
                colnames(gg_PC1) <- paste("group",gg,"_PC1",sep ="")
                PC_group <- cbind(PC_group,gg_PC1)
            }
            PC_for_Score <- as.data.frame(PC_group) %>% rownames_to_column('sample')
            write.table(as.data.frame(PC_group) %>% rownames_to_column('sample'), file = paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".PC1",sep = ""), row.names = FALSE, sep = "\t", quote=FALSE)
            cox_regression(opt$t, paste(opt$o,"/euclidean_pam/",opt$p,".","k=",k,".group",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".PC1",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g, sep = ""))
            HR_result <- read.table(paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/ConsensusGroup_PC1_survival.HR.pval.xls", sep = ""),header = TRUE, sep = "\t")
            
            fomula<-""
            for(hr in k:nrow(HR_result)){
                if( fomula == "" ){
                    if( HR_result[hr,2] > 1){
                        fomula <- HR_result[hr,1]
                    }else{
                        fomula <- paste("-", HR_result[hr,1], sep = "")
                    }
                }else{
                    if( HR_result[hr,2] > 1){
                        fomula <- paste(fomula,"+",HR_result[hr,1],sep = "")
                    }else{
                        fomula <- paste(fomula,"-",HR_result[hr,1], sep = "")
                    }
                }  
            }
            eval(parse(text=paste("Score <- PC_for_Score %>% mutate(Score = ",fomula,")",sep="")))
            write.table(Score, file = paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".Score",sep = ""), row.names = FALSE, sep = "\t", quote=FALSE)
            factor <- "Score"
            write.table(factor, file = paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/factor.list",sep = ""), col.names = FALSE, row.names = FALSE,quote = FALSE)
            SVV_maxstat(opt$t, paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/GeneGroup=",gg,".Score",sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,sep = ""), paste(opt$o,"/euclidean_pam/k=",k,"/euclidean_pam/","gene_k=",g,"/factor.list",sep = ""))
        }
    }
}

if(opt$f == "TR"){
    file.rename(opt$i,paste(opt$i,".transRowCol",sep=""))
    file.rename(paste(opt$i,".old",sep=""),opt$i)
}else{
    file.rename(opt$i,paste(opt$i,".delALL0",sep=""))
    file.rename(paste(opt$i,".old",sep=""),opt$i)
}

timeend<-Sys.time()
cat(paste("End at :",timeend,sep =""),sep = "\n")
cat("--------------------------------------",sep = "\n")
q()
