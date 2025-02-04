
#' Compare the TF activites in a cross plot with the axes the different conditions
#'
#' @param res_new_combo 
#' @param model_list_check 
#' @param condition_colname 
#' @param foldername 
#'
#' @return
#' @export
#'
#' @examples
PlotTheTFByCondition<-function(res_new_combo,model_list_check,condition_colname,foldername){
  
  for_level_check_labels<-colnames(res_new_combo)[grep(condition_colname,colnames(res_new_combo))]
  
  
  for(iter in 1:(length(model_list_check)-1)){
    res_wide_check<-res_new_combo%>%dplyr::filter(condition==as.character(model_list_check[iter]))%>%
      mutate(color_col=case_when(
        (`sign_score.x`>quantile(as.numeric(as.matrix(`sign_score.x`)))[4]+1.5*IQR(as.numeric(as.matrix(`sign_score.x`))))&(`sign_score.y`>quantile(as.numeric(as.matrix(`sign_score.y`)))[4]+1.5*IQR(as.numeric(as.matrix(`sign_score.y`))))~"top_right",
        
        (`sign_score.x`< quantile(as.numeric(as.matrix(`sign_score.x`)))[2]-1.5*IQR(as.numeric(as.matrix(`sign_score.x`))))&(`sign_score.y`>quantile(as.numeric(as.matrix(`sign_score.y`)))[4]+1.5*IQR(as.numeric(as.matrix(`sign_score.y`))))~"top_left",
        
        (`sign_score.x`>quantile(as.numeric(as.matrix(`sign_score.x`)))[4]+1.5*IQR(as.numeric(as.matrix(`sign_score.x`))))&
          (`sign_score.y`<quantile(as.numeric(as.matrix(`sign_score.y`)))[2]-1.5*IQR(as.numeric(as.matrix(`sign_score.y`))))~"bottom_right",
        
        (`sign_score.x`<quantile(as.numeric(as.matrix(`sign_score.x`)))[2]-1.5*IQR(as.numeric(as.matrix(`sign_score.x`))))&
          (`sign_score.y`< quantile(as.numeric(as.matrix(`sign_score.y`)))[2]-1.5*IQR(as.numeric(as.matrix(`sign_score.y`))))~"bottom_left",
        
        (`sign_score.y`>quantile(as.numeric(as.matrix(`sign_score.y`)))[4]+1.5*IQR(as.numeric(as.matrix(`sign_score.y`))))~"top_middle",
        (`sign_score.y`<quantile(as.numeric(as.matrix(`sign_score.y`)))[2]-1.5*IQR(as.numeric(as.matrix(`sign_score.y`))))~"bottom_middle",
        (`sign_score.x`<quantile(as.numeric(as.matrix(`sign_score.x`)))[2]-1.5*IQR(as.numeric(as.matrix(`sign_score.x`))))~"left_middle",
        (`sign_score.x`>quantile(as.numeric(as.matrix(`sign_score.x`)))[4]+1.5*IQR(as.numeric(as.matrix(`sign_score.x`)))) ~"right_middle",
        # ~"middle_middle"
        TRUE~"PANIC"))
    res_wide_check$color_col<-factor(res_wide_check$color_col,levels=c("top_right","bottom_left","right_middle","left_middle","top_middle","bottom_middle","top_left","bottom_right","PANIC"))
    
    ggplot(res_wide_check%>%as.data.frame,
           aes(x=sign_score.x,y=sign_score.y,color=color_col,label=source,fontface="bold"))+
      geom_point()+
      xlab(paste0(as.character(res_wide_check[1,for_level_check_labels[1]])," Score"))+
      ylab(paste0(as.character(res_wide_check[1,for_level_check_labels[2]])," Score"))+
      geom_label_repel(data=res_wide_check%>%as.data.frame%>%
                         dplyr::select(source,sign_score.x,sign_score.y,color_col),
                       aes(x=sign_score.x,
                           y=sign_score.y,
                           label=source),
                       #fontface="bold"),
                       max.overlaps = 20,
                       label.padding=0.1,
                       size=3,show.legend = F)+
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
                                  "top_middle"=paste0(as.character(res_wide_check[1,for_level_check_labels[2]])," Up Only"),
                                  "bottom_left"="Synergistic Down",
                                  "bottom_right"="Antagonist 2",
                                  "bottom_middle"=paste0(as.character(res_wide_check[1,for_level_check_labels[2]])," Down Only"), 
                                  "left_middle"=paste0(as.character(res_wide_check[1,for_level_check_labels[1]])," Down Only"), 
                                  "right_middle"=paste0(as.character(res_wide_check[1,for_level_check_labels[1]])," Up Only"),
                                  "PANIC"="nonsignificant"),name="")+
      theme_bw()
    ggsave(filename = paste0(foldername,as.character(model_list_check[iter]),".pdf"),width=9,height=6,plot=last_plot())
    
  }
  
  
}
