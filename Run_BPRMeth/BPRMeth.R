library("BPRMeth")

args <- commandArgs(TRUE);
f_annotation <- args[1];
f_methylLevel <- args[2];
outdir <- args[3];

anno_dt <- read_anno(f_annotation)
met_dt <-read_met(file=f_methylLevel,type="bulk_seq")
met_region <- create_region_object(met_dt = met_dt, anno_dt = anno_dt)


for (i in 4:8){
	K <- i;
	# Create basis object
	basis_obj <- create_rbf_object(M = 3)
	# Perform clustering
	cl_obj <- cluster_profiles_vb(X = met_region$met, K = K, model = "binomial", alpha_0 = .5, beta_0 = .1, basis = basis_obj, vb_max_iter = 20)
	names(cl_obj$labels)<- names(met_region$met)
	cluster_plot <- plot_cluster_profiles(cluster_obj = cl_obj)

	f_clusterLabel <- paste(outdir,"/cluster_",K,".txt",sep="")
	f_clusterPlot <- paste(outdir,"/cluster_",K,".pdf",sep="")

	pdf(f_clusterPlot);
	print(cluster_plot);
	dev.off();
	write.table(cl_obj$labels,f_clusterLabel,row.names=TRUE, sep="\t",quote = FALSE)
}
