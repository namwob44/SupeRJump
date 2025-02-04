#' Compare biased cells for the same population and the same lineage for different conditions. It Saves figures and produces a table
#'
#' @param reduced_data_frame 
#' @param lineages_to_compare 
#' @param folder_name 
#'
#' @return
#' @export
#'
#' @examples
GetDiffBetweenConditionalPoisedCells<-function(reduced_data_frame,lineages_to_compare,folder_name){
  
  if(length(as.character(unique(reduced_data_frame$Group)))<2){
    return(NULL)
  }
  if(!dir.exists(folder_name)){dir.create(folder_name)}
  
  p_values_TF_table<-reduced_data_frame%>%
    tidyr::pivot_longer(cols = !c(Poised,Group,Combined_Ordering,Cluster.Assignment),
                        names_to = "TF_names",
                        values_to = "Score")%>%
    group_by(TF_names)%>%
    reframe(ttest =list(pairwise.wilcox.test(x=Score,g=Group,p.adjust.method="bonf",pool.sd = F))) %>%
    mutate(ttest = purrr::map(ttest, broom::tidy)) %>%
    unnest(cols = c(ttest))
  
  mean_groups_TF_table<-reduced_data_frame%>%
    tidyr::pivot_longer(cols = !c(Poised,Group,Combined_Ordering,Cluster.Assignment),
                        names_to = "TF_names",
                        values_to = "Score")%>%
    group_by(TF_names,Group)%>%
    summarise(across(Score,mean))%>%
    mutate(log_score = log2(Score))
  
  Log2FC_pvalue_df<-data.frame(TF_names=mean_groups_TF_table%>%
                                 dplyr::filter(Group==as.character(p_values_TF_table$group1[1]))%>%
                                 pull(TF_names),
                               log2FC=(mean_groups_TF_table%>%dplyr::filter(Group==as.character(p_values_TF_table$group1[1]))%>%pull(Score))-(mean_groups_TF_table%>%dplyr::filter(Group==as.character(p_values_TF_table$group2[1]))%>%pull(Score)))%>%
    inner_join(p_values_TF_table,by="TF_names")%>%
    dplyr::select(TF_names,group1,group2,log2FC,p.value)%>%
    mutate(significant_point=(-1)*log10(p.value)*sign(log2FC))%>%
    mutate(SigCall=case_when(
      (log2FC>quantile(log2FC)[4])&(p.value<0.01)~"Up",
      (log2FC<quantile(log2FC)[2])&(p.value<0.01)~"Down",
      TRUE~"NS"
    ))%>%
    mutate(SigCall=factor(SigCall,levels=unique(SigCall)))
  
  
  ggplot(Log2FC_pvalue_df,aes(x=log2FC,y=(-1)*log10(p.value),color=SigCall,label=TF_names,fontface="bold"))+
    geom_point(alpha=0.6)+
    xlab(paste0("Differential TF expression (",as.character(Log2FC_pvalue_df$group1[1]),"-",as.character(Log2FC_pvalue_df$group2[1]),")"))+
    geom_label_repel(data=Log2FC_pvalue_df,
                     aes(x=log2FC,
                         y=(-1)*log10(p.value),color=SigCall),
                     max.overlaps = 20,
                     label.padding=0.1,
                     size=3,
                     force=2,show.legend = F)+
    scale_color_manual(values=c("Up"=brewer.pal(n=5,name="Reds")[5],
                                "Down"=brewer.pal(n=5,name="Blues")[5],
                                "NS"="grey80"),
                       name="")+
    theme_bw()
  ggsave(filename = paste0(folder_name,lineages_to_compare,"_poised_",as.character(Log2FC_pvalue_df$group1[1]),"_vs_",as.character(Log2FC_pvalue_df$group2[1]),".pdf"),width=9,height=6,plot=last_plot())
  return(Log2FC_pvalue_df)
  
  
}