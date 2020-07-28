library(umap)
library(colorspace)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(viridis)

compute_pca <- function(met_mat, vmr_count = 2000, max_ratio_of_na_cells = 0.75)
{
  count_na_by_region = rowSums(is.na(met_mat))
  head(count_na_by_region)
  ratio_na_by_region = count_na_by_region / ncol(met_mat)
  head(ratio_na_by_region)
  sum(ratio_na_by_region > max_ratio_of_na_cells)
  use_regions = ratio_na_by_region <= max_ratio_of_na_cells
  met_mat_filtered = met_mat[use_regions, ]
  dim(met_mat_filtered)
  met_mat_filtered[1:5, 1:5]

  met_mat_imputed = met_mat_filtered
  na_index <- which(is.na(met_mat_imputed), arr.ind=TRUE)
  met_mat_imputed[na_index] <- rowMeans(met_mat_imputed, na.rm=TRUE)[na_index[,1]]

  met_mat_imputed[1:5, 1:5]
  var_region = apply(X = met_mat_imputed, FUN = var, MARGIN = 1)
  var_region_sorted = sort(var_region, decreasing = T)
  vmr = names(head(var_region_sorted, vmr_count))
  head(vmr)
  length(vmr)

  met_mat_vmr = met_mat_imputed[vmr, ]
  dim(met_mat_vmr)

  pca <- prcomp(t(met_mat_vmr))
  class(pca)
  #plot(pca$x[,1:2])

  return(pca$x)
}

compute_umap <- function(pca, num_dims = 10)
{
  pca_input = pca[, 1:num_dims]
  set.seed(123)
  umap_object = umap(pca_input)

  return(umap_object$layout)
}



plot_dim_red <- function(dim_red, groups = NULL, title = '',
                         reduc_type = 'UMAP')
{

  group_flag = !is.null(groups)

  if(!group_flag)
  {
    groups = rep(1, nrow(dim_red))
  }

  par(mar=c(4.1, 4.1, 6.1, 10.1), xpd=TRUE)



  #color_set = get_divergent_color_set()
  #cell_colors = color_set[groups]

  #plot(dim_red[,1:2], col = group_colors,
  #     pch = 20, cex = 0.75,
  #     xlab = "X", ylab = "Y", main = title, cex.main = 0.9 )

  point_size = 500 / nrow(dim_red)

  data = data.frame(dim_red)
  colnames(data) = c('X1', 'X2')

  if(group_flag)
  {
    gg <- ggplot(data, aes(x = X1, y = X2, color = groups))
  }else
  {
    gg <- ggplot(data, aes(x = X1, y = X2))
  }

  gg1 = gg + geom_point(size=point_size) +
        theme_bw() +
        labs(color = '') +
        ggtitle(title) +
        theme(plot.title = element_text( hjust = 0.5)) +
        xlab(paste0(reduc_type, '-1') ) +
        ylab(paste0(reduc_type, '-2') )

  return(gg1)

}#plot_dim_red




plot_feature <- function(dim_red, feature_matrix = NULL,
                         feature = 'NULL', title = '',
                         legend_title = '', reduc_type = 'UMAP',
                         min_cells_with_info = 1){


    idx = which(rownames(feature_matrix) == feature)

    if(length(idx) == 0 ){
      print(paste0(feature, ' not found in activity matrix, skipping. '))
      return()
    }
    non_na_cell_count = sum(!is.na(feature_matrix[feature, ]))

    if(non_na_cell_count < min_cells_with_info ){
      msg = paste0(feature,  ': Only ', non_na_cell_count,
                   ' cells with information. Less than threshold: ',
                   min_cells_with_info, ' . Skipping.')
      print(msg)
      return()
    }

    num_color_slices = 100
    rbPal <- colorRampPalette(c('navy',  'yellow', 'red'))
    viridis_hcl <- colorspace::sequential_hcl(num_color_slices,
                                              h = c(300, 75),
                                              c = c(35, 95),
                                              l = c(15, 90),
                                              power = c(0.8, 1.2))

    feature_vector = feature_matrix[feature, ]
    quant = quantile(x = feature_vector, probs = c(0.05, 0.95), na.rm = T)
    roof = quant[2]
    feature_vector[feature_vector < quant[1]] = quant[1]
    feature_vector[feature_vector > quant[2]] = quant[2]

    #print('---1')

    data = data.frame(dim_red)
    #data$Col <- rbPal(5 + 1)[as.numeric(cut(gene_activity_vector,breaks = 5 + 1))]
    #data$Col <- rbPal(num_color_slices + 1)[as.numeric(cut(gene_activity_vector,breaks = num_color_slices + 1))]
    data$Col <- viridis_hcl[as.numeric(cut(feature_vector, breaks = num_color_slices + 1))]

    head(data)



    #data$Col[gene_expression == 0] = 'gray'
    data$Col[is.na(feature_vector)] = 'gray'

    #opacity = 0.85
    #trans_colors = alpha(data$Col, opacity)
    #na_opacity = 0.15
    #trans_colors[trans_colors == "gray"] = alpha(data$Col, na_opacity)


    point_size = 500 / length(feature_vector)

    gg <- ggplot(data, aes(X1, X2))
    gg1 = gg + geom_point(size=point_size, aes(colour = feature_vector) ) +
       scale_color_viridis(discrete=F) +
       theme_bw() +
       theme(plot.title = element_text( hjust = 0.5)) +
       ggtitle(paste0(title, '\n', feature)) +
       labs(colour = legend_title) +
       xlab(paste0(reduc_type, '-1') ) +
       ylab(paste0(reduc_type, '-2') )



    return(gg1)


}#plot_dim_red



reduce_dims_for_sample <- function(met_mat_for_dim_red, met_mat_for_features,
                                   name_for_dim_red, name_for_features,
                                   methylation_calls_dir, plot_dir,
                                   methylation_type = 'CpG',
                                   min_call_count_threshold = 10,
                                   features = NULL)
{

  pca <- compute_pca(met_mat_for_dim_red, vmr_count = 2000, max_ratio_of_na_cells = 0.75)
  dim(pca)
  umap <- compute_umap(pca, num_dims = 10)
  plot(umap)

  clusters = dens_clus(umap)
  names(clusters) = rownames(umap)
  head(clusters)

  title = paste0(sample_name, ' - ' , methylation_type ,
                 '\nDR region: ', name_for_dim_red,
                 '\nFeature region: ', name_for_features)


  plot_file = paste0(plot_dir, '/UMAP.DR_',name_for_dim_red,'.clusters.eps')
  gg1 = plot_dim_red(umap, groups = clusters, title = title )
  ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')


  feature_plot_dir = paste0(plot_dir, '/regions_',name_for_dim_red,'/')
  dir.create(feature_plot_dir, recursive = T, showWarnings = F )

  #plot_file = paste0(plot_dir, '/genes.eps')
  #cairo_ps(plot_file, fallback_resolution = 300, onefile = T)
  #setEPS()

  if(!is.null(features))
  {
    for(feature in features)
    {
      #print(feature)
      plot_file = paste0(feature_plot_dir, feature,'.eps')
      #cairo_ps(plot_file, fallback_resolution = 2400)
      #postscript(plot_file, onefile = F, width = 7, height = 6)
      gg1 = plot_feature(dim_red = umap,
                         feature_matrix = met_mat_for_features,
                         feature = feature, title = title,
                         reduc_type = 'UMAP', legend_title = 'Methylation')
      print(gg1)
      ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')
      #dev.off()
    }#for(gene in marker_genes)

  }#if(!is.na(features))

  dim_red_list = list(met_mat_for_dim_red = met_mat_for_dim_red,
                met_mat_for_features = met_mat_for_features,
                name_for_dim_red = name_for_dim_red,
                name_for_features  = name_for_features,
                methylation_type = methylation_type,
                pca = pca,
                umap = umap,
                clusters = clusters)

  return(dim_red_list)

}#reduce_dims_for_sample

