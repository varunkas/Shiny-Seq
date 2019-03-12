#####download ppt
source("functions.R")



powerpoint_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    fluidRow(
      column(2,
             textInput(ns('download_all_tables_path'), 'Insert download path for all tables'))),
    
    fluidRow(
      column(2,
             selectInput(ns('datachoiceX')  ,label = h5("Select Data Type"), 
                         choices = list("Excel" = 1, "CSV" = 2),
                         selected = 1)),
      column(1, 
             br(),
             br(),
             actionButton(ns('download_all_tables'), 'Download all tables'))),
    fluidRow(column(5,bsAlert("message_download_path"))),
    fluidRow(column(5, bsAlert("progress"))),
    br(),
    hr(),
    br(),
    br(),
    fluidRow(
      column(3,
             actionButton(ns('generateppt'),'Generate PowerPoint Presentation')),
      column(1,
             uiOutput(ns('generate_download_button')))),
    fluidRow(column(5, bsAlert("progress2")))
  )
}

powerpoint_module<-function(input,output,session,
                            normal,dds.fc,conchoice,top_unorm,top_norm,pca_input_norm,pca_input_raw,
                            top_b,batch_corrected,
                            batch_choice,combination,ok3,DE_genes,DE_TF,
                            input_scale,p_values,input_scale_volx,input_scale_voly,hypothesis_choice,
                            input_Distance_de,input_Linkage_de,input_Distance_tf,input_Linkage_tf,
                            dds,input_Distance_anova,input_Linkage_anova,
                            Enriched_Kegg_table,Enriched_Kegg_obj,
                            Enriched_BP_table,Enriched_BP_obj,
                            Enriched_hall_table,Enriched_hall_obj,
                            anova_table){
  q<-reactiveValues(plot_names=c(),dat=NULL,head=list(),top_pca_batch=1)
  
  # for(run in 1:2){ #execute twice
    observeEvent(input$generateppt,
                 {
                   createAlert(session,"progress2", alertId = "ppt",
                               title="Please Wait", 
                               content = "PowerPoint Presentation is generated. Please wait till this message disappears and klick the Download ppt button")
                   # Create a Progress object
                   progress <- shiny::Progress$new()
                   # Make sure it closes when we exit this reactive, even if there's an error
                   on.exit(progress$close())
                   progress$set(message = "Processing Data", value = 0)
                   combo<-combination()
                   n<-7+length(combo())
                   
                   # Increment the progress bar, and update the detail text.
                   progress$inc(1/(n), detail = paste("Doing part", 1,"/",(n)))
                   # Pause for 0.1 seconds to simulate a long computation.
                   Sys.sleep(0.1)
                   
                   
                   # Slide 2 : Add boxplot
                   #+++++++++++++++++++++++
                   
                   png(filename="./plots/Boxplot_for_raw_data.png")
                   boxplot_output(assay(dds.fc()),colData(dds.fc()),as.numeric(conchoice()))
                   dev.off()
                   
                   q$plot_names=c(q$plot_names,"./plots/Boxplot_for_raw_data.png")
                   
                   #plot_names=c(q$plot_names,"./plots/Boxplot for raw data.png")
                   
                   png(filename="./plots/Boxplot_for_normalized_data.png")
                   boxplot_output(normal(),colData(dds.fc()),as.numeric(conchoice()))
                   dev.off()
                   
                   q$plot_names=c(q$plot_names,"./plots/Boxplot_for_normalized_data.png")
                   
                   # Increment the progress bar, and update the detail text.
                   progress$inc(1/(n), detail = paste("Doing part", 2,"/",(n)))
                   # Pause for 0.1 seconds to simulate a long computation.
                   Sys.sleep(0.1)
                   
                   #slide 3: Add PCA
                   #+++++++++++++++++++++++
                   
                   point_size = 4
                   #2d raw
                   raw_pca<-pcaplot(pca_input_raw(),conchoice(),top_unorm(),"2D",NULL, point_size)[[2]]
                   ggsave(filename="./plots/un_normalized_2D_pca.png", plot=raw_pca,dpi = 600,  height = 7.5, width = 11, units = "cm", device = "png")
                   
                   q$plot_names=c(q$plot_names,"./plots/un_normalized_2D_pca.png")
                   
                   #2d normalized
                   norm_pca<-pcaplot(pca_input_norm(),conchoice(),top_norm(),"2D",NULL, point_size)[[2]]
                   ggsave(filename="./plots/normalized_2D_pca.png", plot=norm_pca,dpi = 600,  height = 7.5, width = 11, units = "cm", device = "png")
                   
                   q$plot_names=c(q$plot_names,"./plots/normalized_2D_pca.png")
                   
                   # Increment the progress bar, and update the detail text.
                   progress$inc(1/(n), detail = paste("Doing part", 3,"/",(n)))
                   # Pause for 0.1 seconds to simulate a long computation.
                   Sys.sleep(0.1)
                   
                   # #  # Slide 4 & 5: Add batch and normalized(2D and 3D)
                   #  #+++++++++++++++++++++++
                   
                   if(!is.null(batch_choice()) && as.numeric(batch_choice()!=1))
                   {
                     pca_input_batch<-list(assay(dds.fc()),colData(dds.fc()),batch_corrected())
                     #print(head(pca_input_batch))
                     print(top_b())
                     top_pca_batch<-1
                     top_batch<-top_b()
                     
                     if(!is.null(top_batch))
                     {
                       top_pca_batch<-top_batch()
                     }
                     q$top_pca_batch<-top_pca_batch
                     batch_corrected_pca<-pcaplot(pca_input_batch,conchoice(),top_pca_batch,"2D",NULL, point_size)[[2]]
                     ggsave(filename="./plots/batch_corrected_2D_pca.png", plot=batch_corrected_pca,dpi = 600, height = 7.5, width = 11, units = "cm")
                     
                     q$plot_names=c(q$plot_names,"./plots/batch_corrected_2D_pca.png")
                   }
                   # Increment the progress bar, and update the detail text.
                   progress$inc(1/(n), detail = paste("Doing part", 4,"/",(n)))
                   # Pause for 0.1 seconds to simulate a long computation.
                   Sys.sleep(0.1)
                   
                   if(!is.null(ok3()) && ok3()>0)
                   {
                     # Slide 7: Add heatmap of 1000 most variable genes
                     #+++++++++++++++++++++++
                     rld<-NULL
                     if(!is.null(batch_choice()) && as.numeric(batch_choice()==1))
                     {
                       rld<-pca_input_norm()[[3]]
                     }
                     else rld<-batch_corrected()
                     #   
                     #   comboX<-combination()
                     #  numX<-length(comboX())
                     #  
                     
                     anova_heatmap_1000 <- heatmap_genes("ANOVA",dds(),dds.fc(),rld,
                                                         1,
                                                         #NULL,
                                                         DE_genes,NULL,NULL,
                                                         #numX,
                                                         NULL,
                                                         "Heatmap of 1000 most variable genes",
                                                         input_Distance_anova(),input_Linkage_anova())
                     
                     png("./plots/Heatmap_of_1000_most_variable_genes.png",width = 23, height = 12.5, 
                         units = "cm", res = 600)
                     print(anova_heatmap_1000)
                     dev.off()
                     
                     q$plot_names=c(q$plot_names,"./plots/Heatmap_of_1000_most_variable_genes.png")
                     #vector graphic to make image editable in ppt
                     # Increment the progress bar, and update the detail text.
                     progress$inc(1/(n), detail = paste("Doing part", 5,"/",(n)))
                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     
                     # Slide 8: summary slide
                     #+++++++++++++++++++++++
                     # Slide 9: comparison
                     #+++++++++++++++++++++++
                     combo<-combination()
                     num<- length(combo())
                     #enrichment kegg data
                     res<-Enriched_Kegg_obj()
                     result<-Enriched_Kegg_table()
                     res_bp<-Enriched_BP_obj()
                     result_bp<-Enriched_BP_table()
                     res_hall<-Enriched_hall_obj()
                     result_hall<-Enriched_hall_table()
                     result_de<-DE_genes()
                     result_tf<-DE_TF()
                     
                     #summary
                     #DE genes
                     
                     q$dat<-summary_analysis(result_de,result_tf,result,result_bp,result_hall,combo())
                     
                     #TF
                     
                     #dat[[2]]
                     
                     #Kegg
                     
                     #dat[[3]]
                     
                     #GO Terms
                     
                     #dat[[4]]
                     
                     #Hallmark
                     
                     #dat[[5]]
                     
                     
                     # Increment the progress bar, and update the detail text.
                     progress$inc(1/(n), detail = paste("Doing part", 7,"/",(n)))
                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     
                     #
                     for (i in 1:num)
                     {
                       name<-str_replace_all(combo()[[i]],"[^[:alnum:]]"," ")
                       #p-value plot
                       png(filename = str_replace_all(paste("./plots/p value plot ",name,".png",sep = "")," ","_"), height = 10, width = 20, units = "cm", res = 600)
                       p_value_all(result_de,i,combo())
                       dev.off()
                       q$plot_names=c(q$plot_names,str_replace_all(paste("./plots/p value plot ",name,".png",sep = "")," ","_"))
                       #ma plot
                       png(filename=str_replace_all(paste("./plots/MA plot ",name,".png",sep = "")," ","_"), height = 10, width = 20, units = "cm", res = 600)
                       ma_plot(i,combo(),result_de,input_scale(),p_values())
                       dev.off()
                       
                       q$plot_names=c(q$plot_names,str_replace_all(paste("./plots/MA plot ",name,".png",sep = "")," ","_"))
                       
                       #volcano plot
                       #nam<-paste("condition",str_replace_all(combination()[[i]][2],"[^[:alnum:]]","."),sep = "")
                       
                       nam<- str_replace_all(paste("./plots/volcano plot ",name,".png",sep = "")," ","_")
                       
                       lim_y<-input_scale_voly()
                       if(is.null(lim_y))
                       {
                         num<- length(combo())
                         #Get the deseq2 dataset
                         res_vol<-result_de[[as.numeric(i)]][5] [[1]]
                         num<-(-log10(na.omit(res_vol$padj)))
                         idx<-which(is.finite(num))
                         max<-max(round(num[idx],1))
                         #max<-(-log10(min(na.omit(res$padj))))
                         #max<-round(max,1)
                         lim_y<-max+5
                       }
                       lim_x<-input_scale_volx()
                       if(is.null(lim_x))
                       {
                         num<- length(combo())
                         #Get the deseq2 dataset
                         res_vol<-result_de[[as.numeric(i)]][5] [[1]]
                         max<-max(na.omit(abs(res_vol$log2FoldChange)))
                         max<-round(max,1)
                         print(max)
                         lim_x<-max+5
                       }
                       
                       p<-volcano_plot(i,combo(),result_de,
                                       lim_x,lim_y,
                                       p_values(),hypothesis_choice())[[1]]
                       ggsave(filename=nam, plot=p,
                              device="png", height = 10, width = 20, units = "cm"
                              
                       )
                       q$plot_names=c(q$plot_names,nam)
                       
                       #comparison of heatmap of de genes and heatmap of tf
                       nam1<- str_replace_all(paste("./plots/Heatmap of DE genes for ",name,".png",sep = "")," ","_")
                       nam<- str_replace_all(paste("./plots/Heatmap of TF for ",name,".png",sep = "")," ","_")
                       rld<-NULL
                       print(batch_choice())
                       if(!is.null(batch_choice()) && as.numeric(batch_choice()==1))
                       {
                         print("ppt 416")
                         rld<-pca_input_norm()[[3]]
                       }
                       else rld<-batch_corrected()
                       print("ppt 419")
                       print(head(rld))
                       
                       h_de_genes<-heatmap_genes("DE",NULL,dds.fc(),rld,i,
                                        DE_genes,NULL,NULL,
                                        num,paste("Heatmap of DE genes for ",name),
                                        input_Distance_de(),input_Linkage_de())
                       #print(h)
                       # if(is.null(h_de_genes)) {
                       #   png(nam1,width = 980, height = 680)
                       #   dev.off()
                       #   
                       # }
                       # else
                       # {
                         png(nam1,width = 12, height = 12, units = "cm", res = 600)
                         print(h_de_genes)
                         dev.off()
                         
                         q$plot_names=c(q$plot_names,nam1)
                       # }
                       
                       # tf heatmap
                       
                       if(!is.null(result_tf))
                       {
                         h_de_tf<-heatmap_genes("TF",NULL,dds.fc(),rld,i,
                                          DE_TF,NULL,NULL,
                                          num,paste("Heatmap of TF genes for ",name),
                                          input_Distance_tf(),input_Linkage_tf())
                         print("ppt line 434")
                         
                         # if(is.null(h_de_tf))
                         # {
                         #   png(nam, width = 980, height = 680)
                         #   dev.off()
                         # }
                         
                         # else
                         # {
                           png(nam,width = 12, height = 12, units = "cm", res = 600)
                           print(h_de_tf)
                           dev.off()
                           
                           q$plot_names=c(q$plot_names,nam)
                           
                         # }
                       }
                       
                       #
                       #up regulated genes
                       
                       #barplot of kegg pathways
                       #barplot of go terms
                       #barplot of hallmark
                       df<-as.data.frame(result_de[[i]][1])
                       print(head(df))
                       print(head(df$FoldChange))
                       genes<-df[,1]
                       df<-df[-1]
                       rownames(df)<-genes
                       print("temp")
                       print(head(order(unlist(df$FoldChange),decreasing=TRUE)))
                       df_final<- df[order(unlist(df$FoldChange),decreasing=TRUE),]
                       colnames(df_final)[1]<-"Overall mean"
                       print(colnames(df_final))
                       for (j in 8:length(colnames(df_final)))
                       {
                         colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                       }
                       a_tab<-as.data.frame(df_final[,-4])
                       print(colnames(a_tab))
                       c<-colnames(a_tab)
                       print(length(c))
                       temp<-as.vector(c[7:length(c)])
                       print(length(temp))
                       print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
                       print('howdy')
                       print(c)
                       df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
                       data<-df_de[,1:(3+length(temp))]
                       print(head(data))
                       
                       if(nrow(data)>0)
                       {
                         q$head[[length(q$head)+1]]<-head(data,10)
                       }
                       ######
                       
                       #up regulated tf
                       
                       df<-as.data.frame(result_tf[[i]][[1]])
                       
                       ######
                       if(nrow(df)>0)
                       {
                         df_final<- df[order(unlist(df$FoldChange),decreasing = TRUE),]
                         colnames(df_final)[1]<-"Overall mean"
                         for (j in 8:length(colnames(df_final)))
                         {
                           colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                         }
                         data<-as.data.frame(df_final[,-c(1,3,4,5)])
                         print(head(data))
                         
                         title<-NULL
                         if(nrow(data)>=10) title<-paste("Top 10 Up regulated Transcription Factors",name)
                         else title<-paste("All Up regulated Transcription Factors",name)
                         print(title)
                         
                         q$head[[length(q$head)+1]]<-head(data,10)
                         
                       }
                       else q$head[[length(q$head)+1]]<-0
                       
                       #up regulated kegg
                       r<-result[[i]][[1]]#as.data.frame(result[[i]][[1]])
                       print('kegg')
                       print(typeof(r))
                       print(head(r))
                       if(nrow(r)==0){
                         title<-paste("No up regulated kegg pathways",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = {
                                  plot(c(-5,5),c(-5,5))
                                  text("No upregulated kegg pathways", x = 0, y = 0)
                                },
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       
                       else{
                         title<-NULL
                         if(nrow(r)>=10) title<-paste("Top 10 up regulated kegg pathways",name)
                         else title<-paste("All up regulated kegg pathways",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         keggplot <- enrichment_plot("kegg",res,result,i,1,10,"")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = keggplot,
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       
                       #biological process
                       #up regulated
                       r<-as.data.frame(result_bp[[i]][[1]])
                       
                       if(nrow(r)==0){
                         title<-paste("No up regulated biological processes",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = {
                                  plot(c(-5,5),c(-5,5),type = "n")
                                  text("No up regulated biological processes", x = 0, y = 0)
                                },
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       else
                       {
                         title<-NULL
                         if(nrow(r)>=10) title<-paste("Top 10 up regulated biological processes",name)
                         else title<-paste("All up regulated biological processes",name)
                         print(title)
                         bioprocplot <- enrichment_plot("biological process",res_bp,result_bp,i,1,10,"")
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = bioprocplot,
                                dpi = 1000,height = 6,width = 12)
                         
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                         
                       }
                       #up regulated hallmark
                       r<-result_hall[[i]][[1]]#as.data.frame(result[[i]][[1]])
                       print('hallmark')
                       print(typeof(r))
                       print(head(r))
                       
                       if(nrow(r)==0){
                         title<-paste("No up regulated Hallmark gene sets",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = {
                                  plot(c(-5,5),c(-5,5),type = "n")
                                  text("No up regulated Hallmark gene sets", x = 0, y = 0)
                                },
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       else
                       {
                         title<-NULL
                         if(nrow(r)>=10) title<-paste("Top 10 up regulated Hallmark gene sets",name)
                         else title<-paste("All up regulated Hallmark gene sets",name)
                         print(title)
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         hallplot <- enrichment_plot("hallmark",res_hall,result_hall,i,1,10,"")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = hallplot,
                                dpi = 1000,height = 6,width = 12)
                         
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                         
                       }
                       #Down regulated genes
                       #pepare table
                       df<-as.data.frame(result_de[[i]][2])
                       genes<-df[,1]
                       df<-df[-1]
                       rownames(df)<-genes
                       df_final<- df[order(unlist(df$FoldChange)),]
                       colnames(df_final)[1]<-"Overall mean"
                       print(colnames(df_final))
                       for (j in 8:length(colnames(df_final)))
                       {
                         colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                       }
                       a_tab<-as.data.frame(df_final[,-4])
                       print(colnames(a_tab))
                       c<-colnames(a_tab)
                       print(length(c))
                       temp<-as.vector(c[7:length(c)])
                       print(length(temp))
                       print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
                       print('howdy')
                       print(c)
                       df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
                       data<-df_de[,1:(3+length(temp))]
                       print(head(data))
                       ######
                       #down regulated genes
                       if(nrow(data)>0)
                       {
                         
                         q$head[[length(q$head)+1]] <- head(data,10)
                         
                       }
                       #down regulated tf
                       
                       df<-as.data.frame(result_tf[[i]][[2]])
                       # genes<-df[,1]
                       # df<-df[-1]
                       # rownames(df)<-genes
                       
                       ######
                       if(nrow(df)>0)
                       {
                         df_final<- df[order(unlist(df$FoldChange)),]
                         colnames(df_final)[1]<-"Overall mean"
                         for (j in 8:length(colnames(df_final)))
                         {
                           colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                         }
                         data<-as.data.frame(df_final[,-c(1,3,4,5)])
                         print("ppt line 642")
                         print(head(data))
                         
                         title<-NULL
                         if(nrow(data)>=10) title<-paste("Top 10 Down regulated Transcription Factors",name)
                         else title<-paste("All Down regulated Transcription Factors",name)
                         print(title)
                         
                         q$head[[length(q$head)+1]] <- head(data,10)
                       }
                       else q$head[[length(q$head)+1]]<-0
                       #down regulated kegg
                       r<-as.data.frame(result[[i]][[2]])
                       #up reg col
                       r2<-as.data.frame(result[[i]][[1]])
                       
                       if(nrow(r)==0){
                         title<-paste("No down regulated kegg pathways",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = {
                                  plot(c(-5,5),c(-5,5),type = "n")
                                  text("No down regulated kegg pathways", x = 0, y = 0)
                                },
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       else
                       {
                         #kegg<-function(){keggplot_all(i,2,res)}
                         
                         title<-NULL
                         if(nrow(r)>=10) title<-paste("Top 10 down regulated kegg pathways",name)
                         else title<-paste("All down regulated kegg pathways",name)
                         print(title)
                         #png(paste(title,'.png'))
                         num<-2
                         if(nrow(r2)==0) num<-1
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                enrichment_plot("kegg",res,result,i,num,10,""),
                                dpi = 1000,height = 6,width = 12)
                         
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                         
                       }
                       #down regulated bp
                       r<-as.data.frame(result_bp[[i]][[2]])
                       #up reg col
                       r2<-as.data.frame(result_bp[[i]][[1]])
                       
                       if(nrow(r)==0){
                         title<-paste("No down regulated biological pathways",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = {
                                  plot(c(-5,5),c(-5,5),type = "n")
                                  text("No down regulated biological pathways", x = 0, y = 0)
                                },
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       else
                       {
                         title<-NULL
                         if(nrow(r)>=10) title<-paste("Top 10 down regulated biological processes",name)
                         else title<-paste("All down regulated biological processes",name)
                         print(title)
                         num<-2
                         if(nrow(r2)==0) num<-1
                         p<-function(){enrichment_plot("biological process",res_bp,result_bp,i,num,10,"")}
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),p(),dpi = 1000,height = 6,width = 12)
                         
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                         
                       }
                       #down regulated hall
                       r<-as.data.frame(result_hall[[i]][[2]])
                       #up reg col
                       r2<-as.data.frame(result_hall[[i]][[1]])
                       
                       if(nrow(r)==0){
                         title<-paste("No down regulated Hallmark gene sets",name)
                         print(title)
                         #png(paste(title,'.png'))
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                plot = {
                                  plot(c(-5,5),c(-5,5),type = "n")
                                  text("No down regulated Hallmark gene sets", x = 0, y = 0)
                                },
                                dpi = 1000,height = 6,width = 12)
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                       }
                       else
                       {
                         
                         title<-NULL
                         if(nrow(r)>=10) title<-paste("Top 10 down regulated Hallmark gene sets",name)
                         else title<-paste("All down regulated Hallmark gene sets",name)
                         print(title)
                         #png(paste(title,'.png'))
                         num<-2
                         if(nrow(r2)==0) num<-1
                         nam<-str_replace_all(paste('./plots/',title,sep = "")," ","_")
                         ggsave(paste(nam,'.png',sep = ""),
                                enrichment_plot("hallmark",res_hall,result_hall,i,num,10,NULL),
                                dpi = 1000,height = 6,width = 12)
                         #dev.off()
                         q$plot_names=c(q$plot_names,paste(nam,'.png',sep = ""))
                         
                       }
                       # Increment the progress bar, and update the detail text.
                       progress$inc(1/(n), detail = paste("Doing part", 7+i,"/",(n)))
                       # Pause for 0.1 seconds to simulate a long computation.
                       Sys.sleep(0.1)
                     }
                     print("line 545")
                     print(q$head)
                     print(q$dat)
                   }
                   closeAlert(session, "ppt")
                   output$generate_download_button<-renderUI({
                     downloadButton(session$ns('downloadppt'), 'Download ppt')
                   })   
                 })
  # }
  ###########
  output$downloadppt <- downloadHandler(
    filename = "result.pptx",
    content = function(file) {
      
      library(officer)
      library(magrittr)
      
      my_pres <- read_pptx()
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Processing Data", value = 0)
      combo<-combination()
      n<-7+length(combo())
      
      my_pres <- my_pres %>% 
        add_slide(layout = "Title Slide", master = "Office Theme") %>%
        ph_with_text(type = "ctrTitle", str = "Data Analysis results of project") %>%
        ph_with_text(type = "subTitle", str = "Presented By:")
     
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 1,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      # Slide 2 : Add boxplot
      #+++++++++++++++++++++++
      print(q$plot_names)
      my_pres<-my_pres %>%
        add_slide(layout = "Comparison", master = "Office Theme") %>%
        ph_with_text(type = "body", index = 4, str = "Raw Data") %>%
        ph_with_img(type = "body", index = 3, src = q$plot_names[1]) %>%
        ph_with_text(type = "body", str = "Normalized Data", index = 1) %>%
        ph_with_img(type = "body", index = 2, src = q$plot_names[2]) %>%
        ph_with_text(type = "title", str = "Box Plot")
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 2,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      # Slide 3 : Add PCA
      #++++++++++++++++++++++
      title<-""
      subtitle1<-"Raw data"
      subtitle2<-"Normalized data"
      if(top_norm()==2 && top_unorm()==2) title<-"2D PCA of All genes"
      else if(top_norm()==1 && top_unorm()==1) title<-"2D PCA of top 500 most variable genes"
      else
      {
        title<-"2D PCA"
        if(top_unorm()==2)
        {
          subtitle1<-"Raw data of all genes"
        }
        else
        {
          subtitle1<-"Raw data of top 500 most variable genes"
        }
        if(top_norm()==2)
        {
          subtitle2<-"Raw data of all genes"
        }
        else
        {
          subtitle2<-"Raw data of top 500 most variable genes"
        }
      }
      
      my_pres <- my_pres %>%
        add_slide(layout = "Comparison", master = "Office Theme") %>%
        ph_with_text(type = "body", index = 4, str = subtitle1) %>%
        ph_with_img(type = "body", index = 3, src = q$plot_names[3], height = 2.95, width = 4.33) %>%
        ph_with_text(type = "body", str = subtitle2, index = 1) %>%
        ph_with_img(type = "body", index = 2, src = q$plot_names[4], height = 2.95, width = 4.33) %>%
        ph_with_text(type = "title", str = title)
      
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 2,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      
      title<-""
      subtitle1<-"Normalized data"
      subtitle2<-"Batch corrected data"
      print(q$top_pca_batch)
      if(top_norm()==2 && q$top_pca_batch==2) title<-"2D PCA of All genes"
      else if(top_norm()==1 && q$top_pca_batch==1) title<-"2D PCA of top 500 most variable genes"
      else
      {
        title<-"2D PCA"
        if(top_norm()==2)
        {
          subtitle1<-"Normalized data of all genes"
        }
        else
        {
          subtitle1<-"Normalized data of top 500 most variable genes"
        }
        if(q$top_pca_batch==2)
        {
          subtitle2<-"Batch corrected data of all genes"
        }
        else
        {
          subtitle2<-"Batch corrected data of top 500 most variable genes"
        }
      }
      src=5
      if(q$plot_names[5]=="./plots/batch_corrected_2D_pca.png")
      {
        my_pres <- my_pres %>%
          add_slide(layout = "Comparison", master = "Office Theme") %>%
          ph_with_text(type = "body", index = 1, str = subtitle1) %>%
          ph_with_img(type = "body", index = 2, src = q$plot_names[4], height = 2.95, width = 4.33) %>%
          ph_with_text(type = "body", str = subtitle2, index = 4) %>%
          ph_with_img(type = "body", index = 3, src = q$plot_names[5], height = 2.95, width = 4.33) %>%
          ph_with_text(type = "title", str = title)
        
      }
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 3,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      if(q$plot_names[6]=="./plots/Heatmap_of_1000_most_variable_genes.png") src=6
      my_pres <- my_pres %>% 
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        ph_with_text(type = "title", str = "Heatmap of 1000 most variable genes") %>%
        ph_with_img(type = "body", src = q$plot_names[src])
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 4,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      name<-unlist(lapply(combo(),function(x) str_replace_all(x,"[^[:alnum:]]","")))
      #str_replace_all(combo()[[i]],"[^[:alnum:]]"," ")
      #DE genes
      dat<-q$dat[[1]]
      print(dat)
      print(colnames(dat))
      dat<-cbind.data.frame(Comparison=name,dat)
      colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
      colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))
      print(dat)
      # ft <- flextable(data = dat) %>% 
      #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
      #   autofit() 
      
      my_pres <- my_pres %>%
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        ph_with_text(type = "title", str = "summary: DE genes") %>%
        #  ph_with_flextable(ft,type = "body")
        #ph_with_table(type = "body", value = q$dat[[1]],first_column=T)
        ph_with_table_at(value=dat, first_column = T, width=9, height=1.5, left=0.5, top = 3)
      
      #TF
      dat<-q$dat[[2]]
      dat<-cbind.data.frame(Comparison=name,dat)
      colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
      colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))      # ft <- flextable(data = dat) %>% 
      #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
      #   autofit() 
      my_pres <- my_pres %>%
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        ph_with_text(type = "title", str = "summary: Transcription factors") %>%
        #  ph_with_flextable(ft,type = "body")
        #ph_with_table(type = "body", value = q$dat[[2]],first_column=T)
        ph_with_table_at(value=dat, first_column = T, width=9, height=1.5, left=0.5, top = 3)
      
      #Kegg
      dat<-q$dat[[3]]
      dat<-cbind.data.frame(Comparison=name,dat)
      colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
      colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))      # ft <- flextable(data = dat) %>% 
      #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
      #   autofit() 
      
      my_pres <- my_pres %>%
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        ph_with_text(type = "title", str = "summary: Enriched Kegg pathways") %>%
        # ph_with_flextable(ft,type = "body")
        #ph_with_table(type = "body", value = q$dat[[3]],first_column=T)
        ph_with_table_at(value=dat, first_column = T, width=9, height=1.5, left=0.5, top = 3)
      
      #GO(Terms)
      dat<-q$dat[[4]]
      dat<-cbind.data.frame(Comparison=name,dat)
      colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
      colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))      # ft <- flextable(data = dat) %>% 
      #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
      #   autofit()
      
      my_pres <- my_pres %>%
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        ph_with_text(type = "title", str = "summary: Enriched GOTerms(Biological processes)") %>%
        # ph_with_flextable(ft,type = "body")
        #ph_with_table(type = "body", value = q$dat[[4]],first_column=T)
        ph_with_table_at(value=dat, first_column = T, width=9, height=1.5, left=0.5, top = 3)
      
      #Hallmark
      dat<-q$dat[[5]]
      dat<-cbind.data.frame(Comparison=name,dat)
      colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
      colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))      # ft <- flextable(data = dat) %>% 
      #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
      #   autofit()
      
      my_pres <- my_pres %>%
        add_slide(layout = "Title and Content", master = "Office Theme") %>%
        ph_with_text(type = "title", str = "summary: Enriched Hallmark") %>%
        #ph_with_table(type = "body", value = q$dat[[5]],first_column=T)
        #ph_with_flextable(ft,type = "body")
        ph_with_table_at(value=dat, first_column = T, width=9, height=1.5, left=0.5, top = 3)
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 5,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      print(src)
      head_idx<-0
      num<-length(combo())
      
      for (i in 1:num)
      {
        name<-str_replace_all(combo()[[i]],"[^[:alnum:]]"," ")
        
        print(name)
        #slide with comparison
        my_pres <- my_pres %>% 
          add_slide(layout = "Title Slide", master = "Office Theme") %>%
          ph_with_text(type = "ctrTitle", str = name)
       
        #p value plot
        print(q$plot_names[src+1])
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = "p-value plot") %>%
          ph_with_img(type = "body", src = q$plot_names[src+1])
        #ma plot
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = "MA plot") %>%
          ph_with_img(type = "body", src = q$plot_names[src+2])
        #volcano plot
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = "Volcano plot") %>%
          ph_with_img(type = "body", src = q$plot_names[src+3])
        result_tf<-DE_TF()
        if(startsWith(q$plot_names[src+5], "./plots/Heatmap_of_TF_for") == TRUE)
        {
          my_pres <- my_pres %>%
            add_slide(layout = "Comparison", master = "Office Theme") %>%
            ph_with_text(type = "body", index = 4, str = "DE genes") %>%
            ph_with_img(type = "body", index = 3, src = q$plot_names[src+4]) %>%
            ph_with_text(type = "body", index = 1, str = "Transcription factor") %>%
            ph_with_img(type = "body", index = 2, src = q$plot_names[src+5]) %>%
            ph_with_text(type = "title", str = "Heatmaps")
          src<-src+5
          
        }
        else
        {
          my_pres <- my_pres %>% 
            add_slide(layout = "Title and Content", master = "Office Theme") %>%
            ph_with_text(type = "title", str = "Heat map of DE genes") %>%
            ph_with_img(type = "body", src = q$plot_names[src+4])
          src=src+4
          
        }
        print(q$plot_names[src+5])
        print(q$head)
        print(head_idx)
        print(q$head[[head_idx+1]])
        
        #up regulated DE genes top 10
        dat<-q$head[[head_idx+1]]
        dat[,1] <- as.numeric(dat[,1])
        dat[,-c(2,3)] <- round(dat[,-c(2,3)], 2)
        dat[,c(2,3)] <- signif(dat[,c(2,3)], 3)
        dat<-cbind.data.frame(Gene=rownames(dat),dat)
        colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
        colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))       
        # ft <- flextable(data = dat) %>% 
        #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
        #   autofit()
        print("hey2")
        
        my_pres <- my_pres %>%
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = "Top 10 most up reglated genes") %>%
          #  ph_with_table(type = "body", value = dat,first_column=T)
          ph_with_table_at(value = dat, height = 4, width = 9, left = 0.5, top = 2, first_column = T)
        #ph_with_flextable(ft,type = "body")
        print("hey3")
        
        #up regulated TF top 10
        if(q$head[[head_idx+2]]!=0)
        {
          dat<-q$head[[head_idx+2]]
          dat[,1] <- as.numeric(dat[,1])
          dat[,-c(2,3)] <- round(dat[,-c(2,3)], 2)
          dat[,c(2,3)] <- signif(dat[,c(2,3)], 3)
          dat<-cbind.data.frame(TF=rownames(dat),dat)
          colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
          colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))        
          # ft <- flextable(data = dat) %>% 
          #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
          #   autofit()
          print("hey5")
          
          my_pres <- my_pres %>%
            add_slide(layout = "Title and Content", master = "Office Theme") %>%
            ph_with_text(type = "title", str = "Up regulated transcription factors") %>%
            # ph_with_table(type = "body", value = dat,first_column=T)
            ph_with_table_at(value = dat, height = 4, width = 9, left = 0.5, top = 2, first_column = T)
          #ph_with_flextable(ft,type = "body")
          print("hey6")
        }
        #up regulated kegg
        str_name<-"Top 10"
        temp<-substr(strsplit(q$plot_names[src+1],"plots/")[[1]][2],1,3)
        if(temp=="All") str_name<-"All"
        if(temp=="No_") str_name <- "No"
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = paste(str_name,"up regulated kegg pathways for",name,sep=" ")) %>%
          ph_with_img(type = "body", src = q$plot_names[src+1])
        
        #up regulated GO Terms
        str_name<-"Top 10"
        temp<-substr(strsplit(q$plot_names[src+2],"plots/")[[1]][2],1,3)
        if(temp=="All") str_name<-"All"
        if(temp=="No_") str_name <- "No"
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = paste(str_name,"up regulated GO Terms for",name,sep=" ")) %>%
          ph_with_img(type = "body", src = q$plot_names[src+2])
        
        #up regulated hallmark
        str_name<-"Top 10"
        temp<-substr(strsplit(q$plot_names[src+3],"plots/")[[1]][2],1,3)
        if(temp=="All") str_name<-"All"
        if(temp=="No_") str_name <- "No"
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = paste(str_name,"up regulated molecular signatures for",name,sep=" ")) %>%
          ph_with_img(type = "body", src = q$plot_names[src+3])
        src=src+3
        
        #down regulated DE genes top 10
        print(q$plot_names[src+1])
        dat<-q$head[[head_idx+3]]
        dat[,1] <- as.numeric(dat[,1])
        dat[,-c(2,3)] <- round(dat[,-c(2,3)], 2)
        dat[,c(2,3)] <- signif(dat[,c(2,3)], 3)
        dat<-cbind.data.frame(Gene=rownames(dat),dat)
        colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
        colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))        # ft <- flextable(data = dat) %>% 
        #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
        #   autofit()
        my_pres <- my_pres %>%
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = "Top 10 down regulated genes") %>%
          #  ph_with_table(type = "body", value = dat,first_column=T)
          ph_with_table_at(value = dat, height = 4, width = 9, left = 0.5, top = 2, first_column = T)
        #ph_with_flextable(ft,type = "body")
        
        #down regulated TF top 10
        if(q$head[[head_idx+4]]!=0)
        {
          print("hey8")
          dat<-q$head[[head_idx+4]]
          dat[,1] <- as.numeric(dat[,1])
          dat[,-c(2,3)] <- round(dat[,-c(2,3)], 2)
          dat[,c(2,3)] <- signif(dat[,c(2,3)], 3)
          dat<-cbind.data.frame(TF=rownames(dat),dat)
          colnames(dat)<- gsub(pattern = "-", replacement = " ",x = colnames(dat))
          colnames(dat)<- gsub(pattern = "_", replacement = " ",x = colnames(dat))          # ft <- flextable(data = dat) %>% 
          #   theme_zebra(odd_header = "transparent", odd_body = "lightblue") %>% 
          #   autofit()
          my_pres <- my_pres %>%
            add_slide(layout = "Title and Content", master = "Office Theme") %>%
            ph_with_text(type = "title", str = "Down regulated transcription factors") %>%
            #  ph_with_table(type = "body", value = dat,first_column=T)
            ph_with_table_at(value = dat, height = 4, width = 9, left = 0.5, top = 2, first_column = T)
          #ph_with_flextable(ft,type = "body")
        }
        head_idx<-head_idx+4
        
        #down regulated kegg
        str_name<-"Top 10"
        temp<-substr(strsplit(q$plot_names[src+1],"plots/")[[1]][2],1,3)
        if(temp=="All") str_name<-"All"
        if(temp=="No_") str_name <- "No"
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = paste(str_name,"down regulated kegg pathways for",name,sep=" ")) %>%
          ph_with_img(type = "body", src = q$plot_names[src+1])
        
        #down regulated GO Terms
        str_name<-"Top 10"
        temp<-substr(strsplit(q$plot_names[src+2],"plots/")[[1]][2],1,3)
        if(temp=="All") str_name<-"All"
        if(temp=="No_") str_name <- "No"
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = paste(str_name,"down regulated GO Terms for",name,sep=" ")) %>%
          ph_with_img(type = "body", src = q$plot_names[src+2])
        
        #down regulated hallmark
        str_name<-"Top 10"
        temp<-substr(strsplit(q$plot_names[src+3],"plots/")[[1]][2],1,3)
        if(temp=="All") str_name<-"All"
        if(temp=="No_") str_name <- "No"
        my_pres <- my_pres %>% 
          add_slide(layout = "Title and Content", master = "Office Theme") %>%
          ph_with_text(type = "title", str = paste(str_name,"down regulated molecular signatures for",name,sep=" ")) %>%
          ph_with_img(type = "body", src = q$plot_names[src+3])
        src=src+3
        
      }
      #print(filename)
      print(file)
      print(my_pres,
            target = file) %>% 
        invisible()
      #}
      
    })
  
  
  #download all tables
  #incase this throws error
  #download and install rtools and set system variable 'PATH' in control panel
  #http://stackoverflow.com/questions/29129681/create-zip-file-error-running-command-had-status-127
  observeEvent(input$download_all_tables,
               {
                 if(input$download_all_tables_path!=""){
                   
                   closeAlert(session, "message1")
                   
                   createAlert(session,"progress", alertId = "message2",
                               title="Please Wait!", content = "Download in progress. Please wait until this message disappears")
                   
                   if(as.numeric(input$datachoiceX==1)){
                     wb <- createWorkbook()
                     addWorksheet(wb, sheetName = "Normalize table")
                     writeData(wb = wb, sheet = 1, x = normal(), colNames = T, rowNames = T)
                     saveWorkbook(wb, file = paste(input$download_all_tables_path,'/normalized_table.xlsx',sep=""), overwrite = T)
                     
                     # write.xlsx2(normal(), file=paste(input$download_all_tables_path,'/normalized_table.xlsx',sep=""), sheetName = "Normalized table",
                     #             col.names = TRUE, row.names = TRUE, append = FALSE)
                   }
                   else {
                     write.csv(normal(), file=paste(input$download_all_tables_path,'/normalized_table.csv',sep=""))
                   }
                   
                   
                   #batch corrected table
                   
                   if(!is.null(batch_choice()) && as.numeric(batch_choice()!=1))
                   {
                     if(as.numeric(input$datachoiceX==1)){
                       wb <- createWorkbook()
                       addWorksheet(wb, sheetName = "Batch corrected table")
                       writeData(wb = wb, sheet = 1, x = batch_corrected(), colNames = T, rowNames = T)
                       saveWorkbook(wb, file = paste(input$download_all_tables_path,'/batch_corrected.xlsx',sep=""), overwrite = T)
                       
                       # write.xlsx2(batch_corrected(),
                       #             #[[1]],
                       #   file=paste(input$download_all_tables_path,'/batch_corrected.xlsx',sep=""), sheetName = "Batch corrected table",
                       #             col.names = TRUE, row.names = TRUE, append = FALSE)
                     }
                     else {
                       write.csv(batch_corrected(),
                                 #[[1]], 
                                 file=paste(input$download_all_tables_path,'/batch_corrected.csv',sep=""))
                     }
                   }
                   
                   #total number of comparisons
                   combo<-combination()
                   num<- length(combo())
                   #de genes table
                   result_de<-DE_genes()
                   #tf table
                   result_tf<-DE_TF()
                   #kegg pathways
                   res<-Enriched_Kegg_obj()
                   result<-Enriched_Kegg_table()
                   #go terms(biological process)
                   res_bp<-Enriched_BP_obj()
                   result_bp<-Enriched_BP_table()
                   #hallmark
                   res_hall<-Enriched_hall_obj()
                   result_hall<-Enriched_hall_table()
                   anova<-anova_table()
                   a_tab<-anova_table()[,-c(2,3)]
                   
                   
                   
                   #excel table for the comparisons selected
                   
                   for (i in 1:num)
                   {
                     
                     #up regulated 
                     #de genes
                     df<-as.data.frame(result_de[[i]][1])
                     genes<-df[,1]
                     df<-df[-1]
                     rownames(df)<-genes
                     df_final<- df[order(unlist(df$FoldChange),decreasing=TRUE),]
                     colnames(df_final)[1]<-"Overall mean"
                     print(colnames(df_final))
                     for (j in 8:length(colnames(df_final)))
                     {
                       colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                     }
                     a_tab<-as.data.frame(df_final[,-4])
                     print(colnames(a_tab))
                     c<-colnames(a_tab)
                     print(length(c))
                     temp<-as.vector(c[7:length(c)])
                     print(length(temp))
                     print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
                     print('howdy')
                     print(c)
                     print("test1")
                     df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
                     data<-df_de[,1:(3+length(temp))]
                     print(head(data))
                     sheet_numb = 0
                     ######
                     if(nrow(data)>0)
                     {
                       
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_name<-paste('Up regulated DE Genes')
                         sheet_numb <- sheet_numb +1
                         # Write the first data set in a new workbook
                         M <- as.matrix(data)
                         wb1 <- createWorkbook()
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         
                         
                         # write.xlsx(data, file = paste(input$download_all_tables_path,"/",condition,".xlsx",sep=""),
                         #            sheetName = sheet_name, append = FALSE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Up regulated DE Genes for', condition, '.csv',sep=""))
                       }
                     }
                     
                     
                     #up regulated tf
                     
                     df<-as.data.frame(result_tf[[i]][[1]])
                     df_final<- df[order(unlist(df$FoldChange),decreasing = TRUE),]
                     colnames(df_final)[1]<-"Overall mean"
                     for (j in 8:length(colnames(df_final)))
                     {
                       colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                     }
                     data<-as.data.frame(df_final[,-c(1,3,4,5)])
                     print(head(data))
                     print("test3")
                     ######
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Up regulated TFs')
                         #write.csv(data, file)
                         sheet_numb <- sheet_numb + 1
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Up regulated TF for', condition, '.csv',sep=""))
                       }
                     }
                     
                     
                     #up regulated kegg
                     data<-result[[i]][[1]]#as.data.frame(result[[i]][[1]])
                     print('kegg')
                     print(typeof(data))
                     print(head(data))
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Up regulated kegg pathways')
                         #write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a third data in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Up regulated kegg pathways for', condition, '.csv',sep=""))
                       }
                       
                     }
                     
                     #biological process
                     #up regulated
                     data<-as.data.frame(result_bp[[i]][[1]])
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Up regulated BPs')
                         # write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Up regulated BP for', condition, '.csv',sep=""))
                       }
                       
                     }
                     #up regulated hallmark
                     data<-result_hall[[i]][[1]]#as.data.frame(result[[i]][[1]])
                     print('hallmark')
                     print(typeof(data))
                     print(head(data))
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Up regulated Hallmarks')
                         # write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a fourth data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Up regulated Hallmark gene sets for', condition, '.csv',sep=""))
                       }
                     }
                     #Down regulated
                     #de genes
                     df<-as.data.frame(result_de[[i]][2])
                     genes<-df[,1]
                     df<-df[-1]
                     rownames(df)<-genes
                     df_final<- df[order(unlist(df$FoldChange),decreasing=TRUE),]
                     colnames(df_final)[1]<-"Overall mean"
                     print(colnames(df_final))
                     for (j in 8:length(colnames(df_final)))
                     {
                       colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                     }
                     a_tab<-as.data.frame(df_final[,-4])
                     print(colnames(a_tab))
                     c<-colnames(a_tab)
                     print(length(c))
                     temp<-as.vector(c[7:length(c)])
                     print(length(temp))
                     print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
                     print('howdy')
                     print(c)
                     print("test5")
                     df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
                     data<-df_de[,1:(3+length(temp))]
                     print(head(data))
                     ######
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste("Down regulated genes")
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Down regulated genes for', condition, '.csv',sep=""))
                       }
                     }
                     #Down regulated tf
                     
                     df<-as.data.frame(result_tf[[i]][[2]])
                     df_final<- df[order(unlist(df$FoldChange),decreasing = TRUE),]
                     colnames(df_final)[1]<-"Overall mean"
                     for (j in 8:length(colnames(df_final)))
                     {
                       colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
                     }
                     data<-as.data.frame(df_final[,-c(1,3,4,5)])
                     print(head(data))
                     ######
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(str_replace_all(combo()[[i]],"[^[:alnum:]]",".")," ","_")
                         sheet_name<-paste('Down regulated TFs')
                         # write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         # write.xlsx(data, file = file, append = TRUE, sheetName = sheet_name)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Down regulated TF for', condition, '.csv',sep=""))
                       }
                     }
                     
                     #Down regulated kegg
                     data<-result[[i]][[2]]#as.data.frame(result[[i]][[1]])
                     print('kegg')
                     print(typeof(data))
                     print(head(data))
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Down regulated kegg pathways')
                         # write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Down regulated kegg pathways for', condition, '.csv',sep=""))
                       }
                     }
                     
                     #biological process
                     #Down regulated
                     data<-as.data.frame(result_bp[[i]][[2]])
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Down regulated BPs')
                         # write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Down regulated BP for', condition, '.csv',sep=""))
                       }
                     }
                     #up regulated hallmark
                     data<-result_hall[[i]][[2]]#as.data.frame(result[[i]][[1]])
                     print('hallmark')
                     print(typeof(data))
                     print(head(data))
                     if(nrow(data)>0)
                     {
                       if (as.numeric(input$datachoiceX==1)){
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         sheet_name<-paste('Down regulated Hallmarks')
                         # write.csv(data, file)
                         file<-paste(input$download_all_tables_path,"/",condition,".xlsx",sep="")
                         sheet_numb <- sheet_numb + 1
                         # Add a second data set in a new worksheet
                         M <- as.matrix(data)
                         addWorksheet(wb1, sheetName = sheet_name)
                         writeData(wb = wb1, sheet = sheet_numb, x = M)
                         # write.xlsx(data, file = file, 
                         #            sheetName=sheet_name, append=TRUE)
                         
                         saveWorkbook(wb1, file = paste(input$download_all_tables_path,"/",condition,".xlsx",sep=""), overwrite = T )
                       }
                       else{
                         condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
                         write.csv(data, file=paste(input$download_all_tables_path,'/Down regulated Hallmark gene sets for', condition, '.csv',sep=""))
                       }
                     }
                   }
                   
                   #anova table
                   
                   ########reordering columns in anova table for display#######
                   cond<-unique(colData(dds())[,as.numeric(conchoice())])
                   print(cond)
                   c<-colnames(a_tab)
                   print(length(c))
                   temp<-as.vector(c[4:(3+length(cond))])
                   temp2<-as.vector(c[(length(cond)+4):length(c)])
                   print(temp2)
                   print(c(temp,c[2],c[3],temp2,c[1]))
                   print('howdy')
                   print(head(a_tab[,c(temp,c[2],c[3],temp2,c[1])]))
                   
                   ######################
                   anova=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
                   print(head(anova))
                   # Add a second data set in a new worksheet
                   if(as.numeric(input$datachoiceX==1)){
                     M <- as.matrix(anova)
                     wb <- createWorkbook()
                     addWorksheet(wb, sheetName = "ANOVA table")
                     writeData(wb = wb, sheet = 1, x = M, colNames = T, rowNames = T)
                     saveWorkbook(wb, file=paste(input$download_all_tables_path,'/ANOVA_table.xlsx',sep=""), overwrite = T)
                     # write.xlsx(anova, file = paste(input$download_all_tables_path,'/ANOVA_table.xlsx',sep=""), 
                     #           sheetName="ANOVA table", append=FALSE)
                   }
                   else {
                     write.csv(x=anova, file = paste(input$download_all_tables_path,'/ANOVA_table.csv',sep=""))
                   }
                   
                   closeAlert(session, "message2")
                   
                   
                   ###########################################################################
                 }
                 else{
                   createAlert(session,"message_download_path", alertId = "message1",
                               title="Error", 
                               content = "Please enter path.")
                 }
               }
  )}
