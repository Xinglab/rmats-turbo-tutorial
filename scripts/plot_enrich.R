# /**
#  * @author [Yuanyuan Wang]
#  * @email [wyynju1993@gmail.com]
#  * @create date 2021-06-03 21:39:48
#  * @modify date 2021-06-03 21:39:48
#  * @desc [description]
#  */

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


run_hypergeometric <- function(observed_gene_list, bg_gene_list, dir_database, termFrom = 'GO') {
      # result = data.frame(Term = NA, Overlap = NA, P.value = NA, Adjusted.P.value = NA, Odds.Ratio = NA, Genes = NA, db = NA, Overlap_number = NA, Term_gene_number = NA)[numeric(0), ]
      result = c()
      if (termFrom == 'GO') {
            db_vec = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
            db_vec_text = c("Biological_Process", "Cellular_Component", "Molecular_Function"); names(db_vec_text) = db_vec
      } else if (termFrom == 'pathway') {
            db_vec = c("WikiPathway_2021_Human", "BioPlanet_2019", "Elsevier_Pathway_Collection")
            db_vec_text = c("WikiPathway", "BioPlanet", "Elsevier"); names(db_vec_text) = db_vec
      }
      
      for (db in db_vec) {
      database = file(paste0(dir_database, '/', db, '.txt'), open = 'r')
      while(length(oneLine <- readLines(database, n = 1, warn = FALSE)) > 0) {
            x = strsplit(oneLine, '\t')[[1]]
            term = x[1]
            genes = x[-c(1)]; genes = genes[!genes == '']
            total = toupper(bg_gene_list)
            group1 = toupper(observed_gene_list)
            group2 = intersect(genes, total)
            overlap = intersect(group1, group2)
            if (length(overlap) == 0) {
                  next
            }
            pval = phyper(length(overlap) - 1, length(group2), length(total) - length(group2), length(group1), lower.tail = FALSE)
            odds.ratio = ( length(overlap) / length(group1) ) / ( length(group2) / length(total) )
            Genes = paste(overlap, collapse = ';')
            Overlap = paste(c(length(overlap), length(group2), length(genes)), collapse = '/')
            result = rbind(result, c(Term = term, Overlap = Overlap, P.value = pval, Adjusted.P.value = NA, Odds.Ratio = odds.ratio, Genes = Genes, db = as.character(db_vec_text[db]), Overlap_number = length(overlap), Term_gene_number = length(group2)))
      }
      close(database)
      }
      # c(term, Overlap, pval, NA, odds.ratio, Genes, db, length(overlap), length(group2))
      result = data.frame(result, stringsAsFactors = FALSE)
      result$Adjusted.P.value = p.adjust(as.numeric(as.character(result$P.value)), method = 'fdr')
      result[, c(3,4,5,8,9)] = apply(result[, c(3,4,5,8,9)], 2, as.numeric)
      return(result)
}


plot_enrich <- function(df_enrich, cutoff = 0.05) {
      
      cutoff = 0.05
      aspect_ratio = 0.7
      bar_width = 0.6
      x_axis.tilt = 0
      
      {
      if (nrow(df_enrich) == 0) {
      print("STOP: There is no enriched GO term.")
      quit("no")}
      }
      
      db_vec = unique(df_enrich$db)
      tmp = lapply(db_vec, FUN = function(x) {
            df_subset = df_enrich[df_enrich$db == x, ]
            df_subset$score = df_subset$Adjusted.P.value
            # if ('updated_fdr' %in% names(df_subset)) {
            #     df_subset$score = df_subset$updated_fdr
            # } else {
            # df_subset$score = df_subset$Adjusted.P.value
            # }
            df_subset <- df_subset[order(df_subset$score), ]
            if (df_subset$score[1] <= cutoff) {
            sig_subset <- subset(df_subset, score <= cutoff)
            retain_index = min(10, max(min(5, nrow(df_subset)), nrow(sig_subset)))
            } else {
            retain_index = min(5, nrow(df_subset))
            }
            retain_subset <- df_subset[1:retain_index, ]
            retain_subset
      })
      
      top_significant_hits = do.call(rbind, tmp)
      
      # shorten name
      index = nchar(top_significant_hits$Term) > 100
      top_significant_hits$Term[index] = paste(substr(top_significant_hits$Term[index], 1, 100), substrRight(top_significant_hits$Term[index], 12), sep = '... ')
      
      library(ggplot2)
      library(scales)
      
      p <- ggplot(top_significant_hits, aes(x = reorder(Term, nrow(top_significant_hits):1), y = -log10(score), fill = db, alpha = Odds.Ratio)) + 
      geom_col(position = position_dodge2(width = bar_width, preserve = "single"))+
      facet_grid(rows = vars(db), scales = 'free_y', space = 'free_y', switch="both") +
      scale_fill_Publication(guide = "none") +
      labs(x = "Terms", y = "-log10(adjusted p value)", alpha = 'Odds.Ratio') + 
      # theme(aspect.ratio = aspect_ratio)+ # somehow affects the facet_grid scale and space
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
      scale_x_discrete(position = "bottom", expand = c(0, 0)) +
      scale_alpha_continuous(range = c(0.3, 1.5)) +
      # scale_fill_gradient2(midpoint = 0, low = muted("blue"), mid = "white", high = muted("red"),
      #                         breaks = seq(0, ceiling(max(top_significant_hits$Odds.Ratio)), length.out=5), 
      #                         labels = seq(0, ceiling(max(top_significant_hits$Odds.Ratio)), length.out=5), 
      #                         limits = c(0, NA)) +
      geom_text(aes(label = paste(Overlap_number, Term_gene_number, sep = '/'), y = -log10(score) + 0.1), position = position_dodge(0.9), vjust = 0.5, hjust = 0, col = 'darkgrey', show.legend = FALSE) +
      scale_y_continuous(limits = c(0, max(-log10(top_significant_hits$score)) * 1.2), expand = c(0,0)) +
      # theme_minimal(base_family = "Roboto Condensed") +
      theme(axis.text.x = element_text(hjust=0)) +
      theme(panel.background = element_blank())+
      theme(axis.text.x=element_text(angle = x_axis.tilt, hjust=0))+
      theme(legend.position = 'bottom', 
            legend.direction = 'horizontal',
            legend.key.width = unit(1.5,"line"),
            legend.key.size = unit(0.7,"line")) +
      coord_flip()
      return(list(plot = p, termNum = nrow(top_significant_hits)))
}


save_enrich <- function(df_enrich, p, suffix, outFolder = './', termFrom = 'GO'){
      if (p$termNum <= 20) {
            height = 4
      } else if (p$termNum <= 25) {
            height = 4.5
      } else {
            height = 5
      }
      pdf(file = file.path(outFolder, paste0('plot_enrich_', termFrom, '_', suffix, '.pdf')), width = 10, height = height); print(p$plot); dev.off()
      write.table(df_enrich, file = file.path(outFolder, paste0('output_enrich_', termFrom, '_', suffix, '.txt')), row.names = T, col.names = T, sep = '\t', quote = FALSE)
}



#################### MAIN ####################
# Rscript {params.sc_enrich} {output.fn_enrich_input} {input.bg} {params.out_suffix} {params.out_path}
# Rscript /mnt/isilon/xing_lab/wangy14/snakemake/RNAseq_snakemake/scripts/plot/plot_enrich.R /mnt/isilon/xing_lab/wangy14/project/PRMT9_mouse_new/07_enrich/WT_KO/input_enrich_GO_SE.txt /mnt/isilon/xing_lab/wangy14/project/PRMT9_mouse_new/02_deseq2Filtering/WT_KO/filtered_deseq2.txt SE /mnt/isilon/xing_lab/wangy14/project/PRMT9_mouse_new/07_enrich/WT_KO 

args <- commandArgs(trailingOnly = TRUE)
fn_enrich_input <- args[1]
fn_enrich_input_bg <- args[2]
out_suffix = args[3]
dir_database = args[4]
if (length(args) > 4) {termFrom = args[5]} else {termFrom = 'GO'} # termFrom = 'GO' or 'pathway'

observed_gene_list = unique(tryCatch(read.table(fn_enrich_input, header = F, stringsAsFactors = FALSE), error=function(e) NULL)$V1)
bg_gene_list = unique(tryCatch(read.table(fn_enrich_input_bg, header = F, stringsAsFactors = FALSE), error=function(e) NULL)$V1)


{
      if (length(observed_gene_list) < 5) { 
            # opt <- options(show.error.messages=FALSE)
            # on.exit(options(opt))
            # stop("STOP: Foreground gene list is too short.") 
            # print("Attention: script did NOT end!")
            quit()
      }
}

df_enrich = run_hypergeometric(observed_gene_list, bg_gene_list, dir_database, termFrom)
p = plot_enrich(df_enrich, cutoff = 0.05)
save_enrich(df_enrich, p, out_suffix, './' , termFrom)
