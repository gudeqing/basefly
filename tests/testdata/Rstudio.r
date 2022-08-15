library(ggplot2)

# 读入数据
data = read.csv('all.data.group.csv', header=T)

# 将数据转换成数据框
data = as.data.frame(data)

# 筛选出肿瘤样本
tumor = data[data['sample_type']=='tumor', ]

# 按照样本id排序
tumor = tumor[order(tumor$sample), ]

# 对tumor查看化疗效果分布
# 查看诊断分布
count = as.data.frame(table(tumor$RECIST))
count= count[order(count[, 2], decreasing=TRUE), ]
p = ggplot(data=count, aes(x=Var1, y=Freq, fill=Var1)) + geom_bar(stat='identity', width=0.5) + theme_minimal()

# 查看分组metatsic
unique(tumor$metastatic)
table(tumor$metastatic)

# 比较primary liver 中 CD274 ('PD-L1')的表达差异
r = t.test(tumor$CD274~tumor$metastatic, paired=T)
r2 = t.test(tumor$CD274~tumor$metastatic, paired=T)

# 查看检验结果
r
# 画boxplot, 直观比较PD-L1在原发灶和转移灶中的区别
p<-ggplot(tumor, aes(x=metastatic, y=CD274, color=metastatic)) + geom_boxplot(outlier.color='blue')
p

# 画boxplot, 直观比较原发灶中和转移灶中，有化疗组和无化疗组中，CD274的表达分布
p<-ggplot(tumor, aes(x=metastatic, y=CD274, fill=chemotherapy)) + geom_boxplot(outlier.color='blue')
p

# 聚类热图展示HLA基因的表达
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
mat = as.matrix(tumor[, grep("HLA\\.", colnames(tumor))])
# mat = as.matrix(tumor[, 21:85])
row.names(mat) = tumor$sample
mat = t(mat)
mat_scale = t(apply(mat, 1, scale))
ht_lst = Heatmap(
	mat_scale, name = "expr", row_km = 5, 
    col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
    top_annotation =  HeatmapAnnotation(
		df=tumor[, c('RECIST', 'metastatic', 'chemotherapy')],
		col=list(
			RECIST = c("CR"=brewer.pal(7,"Set2")[1],
						"PR"=brewer.pal(7,"Set2")[2],
						"SD"=brewer.pal(7,"Set2")[3],
						"PD"=brewer.pal(7,"Set2")[4],
						"unknown"=brewer.pal(7,"Set2")[7]
						),
			metastatic = c("liver"=brewer.pal(7,"Set2")[5],
						"primary"=brewer.pal(7,"Set2")[6]),
			chemotherapy = c("chemotherapy"=brewer.pal(7,"Set2")[5],
						"noChemotherapy"=brewer.pal(7,"Set2")[6])
		)
	),
    show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE,
	clustering_method_columns = "complete",
	clustering_distance_columns = "pearson"
)
draw(ht_lst)
