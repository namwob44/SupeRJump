#' Compare TFs across different lineage Fate pairs
#'
#' @param res_new 
#' @param folder_name 
#'
#' @return
#' @export
#'
#' @examples
PlotTheCrossFates<-function(res_new,folder_name){
  res_new_wide<-res_new%>%tidyr::pivot_wider(id_cols = source,values_from = sign_score,names_from = condition)
  
  for(pair_iter in 1:(dim(res_new_wide)[2]-2)){
    for(other_iter2 in (pair_iter+1):(dim(res_new_wide)[2]-1)){
      
      color_values<-as.numeric(as.matrix(res_new_wide[,1+pair_iter])*as.matrix(res_new_wide[,1+other_iter2]))
      
      res_wide_check<-res_new_wide%>%
        mutate(color_col=case_when(
          (res_new_wide[,1+pair_iter]>quantile(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))[4]+1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+pair_iter]))))&(res_new_wide[,1+other_iter2]>quantile(as.numeric(as.matrix(res_new_wide[,1+other_iter2])))[4]+1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+other_iter2]))))~"top_right",
          
          (res_new_wide[,1+pair_iter]< quantile(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))[2]-1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+pair_iter]))))&(res_new_wide[,1+other_iter2]>quantile(as.numeric(as.matrix(res_new_wide[,1+other_iter2])))[4]+1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+other_iter2]))))~"top_left",
          
          (res_new_wide[,1+pair_iter]>quantile(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))[4]+1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+pair_iter]))))&
            (res_new_wide[,1+other_iter2]<quantile(as.numeric(as.matrix(res_new_wide[,1+other_iter2])))[2]-1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+other_iter2]))))~"bottom_right",
          
          (res_new_wide[,1+pair_iter]<quantile(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))[2]-1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+pair_iter]))))&
            (res_new_wide[,1+other_iter2]< quantile(as.numeric(as.matrix(res_new_wide[,1+other_iter2])))[2]-1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+other_iter2]))))~"bottom_left",
          
          (res_new_wide[,1+other_iter2]>quantile(as.numeric(as.matrix(res_new_wide[,1+other_iter2])))[4]+1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+other_iter2]))))~"top_middle",
          (res_new_wide[,1+other_iter2]<quantile(as.numeric(as.matrix(res_new_wide[,1+other_iter2])))[2]-1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+other_iter2]))))~"bottom_middle",
          (res_new_wide[,1+pair_iter]<quantile(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))[2]-1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+pair_iter]))))~"left_middle",
          (res_new_wide[,1+pair_iter]>quantile(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))[4]+1.5*IQR(as.numeric(as.matrix(res_new_wide[,1+pair_iter])))) ~"right_middle",
          # ~"middle_middle"
          TRUE~"PANIC"))
      res_wide_check$color_col<-factor(res_wide_check$color_col,levels=c("top_right","bottom_left","right_middle","left_middle","top_middle","bottom_middle","top_left","bottom_right","PANIC"))
      ggplot(res_wide_check%>%as.data.frame,
             aes(x=!!sym(colnames(res_wide_check)[1+pair_iter]),
                 y=!!sym(colnames(res_wide_check)[1+other_iter2]),
                 color=color_col,label=source,fontface="bold"))+
        #color=as.numeric(as.matrix(res_new_wide[,1+pair_iter])*as.matrix(res_new_wide[,1+other_iter2]))))+
        geom_point(alpha=0.6)+
        geom_label_repel(data=res_wide_check%>%
                           dplyr::select(source,!!sym(colnames(res_wide_check)[1+pair_iter]),!!sym(colnames(res_wide_check)[1+other_iter2]),color_col),
                         aes(x=!!sym(colnames(res_wide_check)[1+pair_iter]),
                             y=!!sym(colnames(res_wide_check)[1+other_iter2]),color=color_col),
                         #label=source,
                         #fontface="bold"),
                         max.overlaps = 20,
                         label.padding=0.1,
                         size=4,
                         force=1,show.legend = F)+
        scale_color_manual(values=c("top_right"=brewer.pal(n=5,name="Reds")[5],
                                    "top_left"="yellowgreen",
                                    "top_middle"="orange",
                                    "bottom_left"=brewer.pal(n=5,name="Blues")[5],
                                    "bottom_right"="purple",
                                    "bottom_middle"="violet",
                                    "left_middle"="darkgreen",
                                    "right_middle"="palevioletred1",
                                    "PANIC"="grey80"),
                           labels=c("top_right"="Synergistic Up",
                                    "top_left"="Antagonist 1",
                                    "top_middle"=paste0(as.character(colnames(res_wide_check)[1+other_iter2])," Up Only"),
                                    "bottom_left"="Synergistic Down",
                                    "bottom_right"="Antagonist 2",
                                    "bottom_middle"=paste0(as.character(colnames(res_wide_check)[1+other_iter2])," Down Only"), 
                                    "left_middle"=paste0(as.character(colnames(res_wide_check)[1+pair_iter])," Down Only"), 
                                    "right_middle"=paste0(as.character(colnames(res_wide_check)[1+pair_iter])," Up Only"),
                                    "PANIC"="nonsignificant"),name="")+
        theme_bw(base_size = 20)
      ggsave(filename = paste0(folder_name,as.character(colnames(res_wide_check)[1+pair_iter]),"_vs_",as.character(colnames(res_wide_check)[1+other_iter2]),".pdf"),width=9,height=6,plot=last_plot())
      
    }
  }
  
}