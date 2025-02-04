#' PRoduce a plot that shows where jumps are present on PCs and cell types
#'
#' @param down_Z 
#' @param Y 
#' @param model_fits_list 
#' @param foldername 
#' @param top_n_PCs 
#'
#' @return
#' @export
#'
#' @examples
PlotJumpDensityAlongPCLoading<-function(down_Z,Y,model_fits_list,foldername,top_n_PCs=10){
  
  
  list_of_lambdas<-lapply(model_fits_list[["model_fits"]],function(x){
    return(x[["solution"]][1])})%>%as.matrix
  idx_for_jump_models<-(which(list_of_lambdas<max(down_Z$PT_score)&(list_of_lambdas>0)))
  
  Y$type<-"synth_sink"
  Y<-Y%>%as.data.frame
  rownames(Y)<-LETTERS[1:(dim(Y)[1])]
  
  wide_plot_df<-rbind(down_Z,Y[,colnames(down_Z)])
  wide_plot_df<-wide_plot_df%>%
    mutate(color_value=ifelse(type!="synth_sink","skip",as.character(group)))
  
  long_loading_pc_df_norm<-wide_plot_df%>%
    tibble::rownames_to_column(var="Cell")%>%
    as.data.frame%>%
    pivot_longer(cols=!c(Cell,type,group,PT_score,color_value),
                 names_to = "PC_name",
                 values_to="score")
  
  long_loading_pc_df_norm$color_value<-factor(long_loading_pc_df_norm$color_value,levels=c("skip",levels(Y$group)))
  
  long_loading_pc_df_norm$score<-as.numeric(long_loading_pc_df_norm$score)
  long_loading_pc_df_norm$PC_name<-factor(long_loading_pc_df_norm$PC_name,levels=rev(unique(long_loading_pc_df_norm$PC_name)))
  
  sink_levels<-as.character(Y$group)
  jump_table<-matrix(FALSE,nrow=length(unique(long_loading_pc_df_norm$PC_name)),ncol=length(unique(long_loading_pc_df_norm$group)))
  rownames(jump_table)<-unique(long_loading_pc_df_norm$PC_name)
  colnames(jump_table)<-levels(long_loading_pc_df_norm$group)
  for(iter in 1:length(idx_for_jump_models)){
    jump_table[idx_for_jump_models[iter]%%nrow(jump_table),sink_levels[ceiling(idx_for_jump_models[iter]/nrow(jump_table))]]<-TRUE
  }
  
  jump_df<-jump_table%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="PC_names")%>%
    pivot_longer(cols=!PC_names,names_to = "Sink_Type",values_to="Jump")%>%
    mutate(color_value=ifelse(Jump==FALSE,"skip",as.character(Sink_Type)))
  
  jump_df$color_value<-factor(jump_df$color_value,levels=c("skip",levels(Y$group)))
  jump_df$PC_names<-factor(jump_df$PC_names,levels=rev(unique(jump_df$PC_names)))
  jump_df$Sink_Type<-factor(jump_df$Sink_Type,levels=levels(Y$group))
  
  
  phase_density_plot<-ggplot(long_loading_pc_df_norm,
                             aes(x=score,y=PC_name,color=color_value,group=color_value,alpha=ifelse(type!="synth_sink",0.1,1)))+#,size=ifelse(type!="synth_sink",1,2)))+
    geom_point(pch = "|",na.rm = F,cex=4)+
    theme_bw(base_size=8)+
    scale_color_manual(values=c("grey80",pals::tol()[1:(length(unique(long_loading_pc_df_norm$color_value))-1)]))+
    xlab("PC_loading")+
    guides(alpha="none",size="none",color="none")
  
  Jump_heatmap<-ggplot(jump_df,aes(y=PC_names,x=Sink_Type,fill=color_value))+
    geom_tile(color="grey80")+scale_fill_manual(values=c("white",pals::tol()[1:(length(unique(jump_df$color_value))-1)]),name="Sink with Jumps")+
    theme_bw(base_size=10)+
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1),
          axis.text.y=element_blank())+xlab("")+ylab("") 
  
  
  cowplot::plot_grid(phase_density_plot, Jump_heatmap,nrow = 1, align = "h", 
                     axis = "l", rel_widths = c(0.5, 0.5))
  
  #ggsave(filename =paste0(foldername,"Jumps_to_sinks_per_gene_profile.pdf"),width=9,height=6,plot=last_plot())
  
}