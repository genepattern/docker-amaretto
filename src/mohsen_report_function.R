##################################  HTML Report Functions

amaretto_html_report <- function(AMARETTOinit,AMARETTOresults,CNV_matrix,MET_matrix,hyper_geo_test_bool=TRUE)
{

    suppressMessages(suppressWarnings(library("AMARETTO")))
    #file_wd=dirname(rstudioapi::getSourceEditorContext()$path)
    #setwd(file_wd)
    file_wd='./'
    setwd(file_wd)
    ########################################################
    # Evaluate AMARETTO Results
    ########################################################    
    # AMARETTOtestReport<-AMARETTO_EvaluateTestSet(AMARETTOresults,
    #                                              AMARETTOinit$MA_matrix_Var,AMARETTOinit$RegulatorData)
    ######################################################################################################################################################################################
    ######################################################################################################################################################################################
    NrModules<-AMARETTOresults$NrModules
    if (hyper_geo_test_bool)
      {
      saveRDS(AMARETTOresults, file = "hyper_geo_test/AMARETTOtestReport.RData")
      ########################################################
      # Save AMARETTO results in different formats including .gmt 
      ######################################################## 
      suppressMessages(suppressWarnings(source("/usr/local/bin/amaretto/hyper_geo_test/ProcessTCGA_modules.R")))
      rileen("/usr/local/bin/amaretto/hyper_geo_test/AMARETTOtestReport.RData", AMARETTOinit, AMARETTOresults)
      }
    ######################################################################################################################################################################################
    ######################################################################################################################################################################################
    ######################################################################################################################################################################################
    ######################################################################################################################################################################################
    
    ##################################################################################################################################################################
    #REPORT
    ##################################################################################################################################################################
    unlink("report_htm/*")
    unlink("report_html/htmls/*")
    unlink("report_html/htmls/images/*")
    unlink("report_html/htmls/data/*")
    unlink("report_html/htmls/tables/*")
    unlink("report_html/htmls/tables/module_hyper_geo_test/*")

    dir.create("report_html")
    dir.create("report_html/htmls")
    dir.create("report_html/htmls/images")
    dir.create("report_html/htmls/data")
    dir.create("report_html/htmls/tables")
    dir.create("report_html/htmls/tables/module_hyper_geo_test")
    ########################################################
    # Save images of all the modules
    ########################################################   
    address1=paste("./","report_html",sep="")
    address2=paste("./","htmls",sep="")
    address3=paste("./","htmls/images",sep="")
    for (ModuleNr in 1:NrModules )
    {
      html_address=paste("report_html","/htmls/images","/module",as.character(ModuleNr),".jpeg",sep="")
      jpeg(file =html_address )
      AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults=AMARETTOresults, CNV_matrix, MET_matrix, ModuleNr=ModuleNr) 
      dev.off()
    }
    
    ##############################################################################
    # Create HTMLs for each module
    ##############################################################################   

    if (hyper_geo_test_bool)
    {
        ###################################################
        library("GSEABase")
        library("rstudioapi")
        suppressMessages(suppressWarnings(source("/usr/local/bin/amaretto/hyper_geo_test/HyperGTestGeneEnrichment.R")))
        suppressMessages(suppressWarnings(source("/usr/local/bin/amaretto/hyper_geo_test/word_Cloud.R")))

        b<- HyperGTestGeneEnrichment("/usr/local/bin/amaretto/hyper_geo_test/H.C2CP.genesets_forRileen.gmt", "/usr/local/bin/amaretto/hyper_geo_test/TCGA_modules_target_only.gmt", "hyper_geo_test/output.txt",show.overlapping.genes=TRUE)
        df1<-read.table("hyper_geo_test/output.txt",sep="\t",header=TRUE, fill=TRUE)
        #df2<-df1[order(-df1$p.value),]
        df2=df1
        df3<-read.table("hyper_geo_test/output.genes.txt",sep="\t",header=TRUE, fill=TRUE)
        ###################################################
        print(head(df3))
    }
    

    library(R2HTML)
    
    number_of_significant_gene_overlappings<-c()
    for (ModuleNr in 1:NrModules )
    {
      
      module_name=paste("module",as.character(ModuleNr),sep="")
      
      ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
      currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
      RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
      module_regulators_weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
      module_regulators_weights<-data.frame(module_regulators_weights)
      positiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] > 0)]
      negetiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] < 0)]
      ModuleGenes=rownames(ModuleData)
      RegulatoryGenes=rownames(RegulatorData)
      module_all_genes_data <- rbind(ModuleData, RegulatorData)
      module_all_genes_data <-module_all_genes_data[order(rownames(module_all_genes_data)),]
      module_all_genes_data <- unique(module_all_genes_data)
      module_annotations<-create_gene_annotations(module_all_genes_data,ModuleGenes,module_regulators_weights)
      
      
      if (hyper_geo_test_bool)
      {
        ####################### Hyper Geometric Significance
        module_name2=paste("Module_",as.character(ModuleNr),sep="")
        print(module_name2)
        
        
        
        filter_indexes<-(df3$Testset==module_name2) & (df3$p.value<0.05)
        
        gene_descriptions<-df3$Description[filter_indexes]
        
        gene_names<-df3$Geneset[filter_indexes]
        print(length(gene_names))
        print('hassan')
        overlapping_gene_names<-df3$Overlapping.genes[filter_indexes]

        number_overlappings<-df3$n.Overlapping[filter_indexes]
        
        p_values<-df3$p.value[filter_indexes]
        q_values<-df3$q.value[filter_indexes]
        
        number_of_significant_gene_overlappings<-c(number_of_significant_gene_overlappings,length(gene_names))
        
        
        print(head(df3))
        print('hassan')
        
        mmm<-gene_descriptions
        mm<-as.vector(unique(mmm))
        
        descriptions=""
        for (var in mm) 
          {
            descriptions = paste(descriptions,var,sep=" ")
            descriptions =gsub(">",",",descriptions)
            # descriptions<-substring(descriptions, 1)
            # descriptions<-sub('.', '', descriptions)
          }
        if (nchar(descriptions)>0)
          {
            wordcloud_making(descriptions,module_name2)
          }
        #############################################
      }
      
      print(number_of_significant_gene_overlappings)
      
      address=address2
      fname=paste("module",as.character(ModuleNr),sep="")
      tite_page=paste("module",as.character(ModuleNr),sep="")
      graph1=paste("./images","/module",as.character(ModuleNr),".jpeg",sep = "")
      tmpfic<-HTMLInitFile("./report_html/htmls/",filename=fname,Title = tite_page,CSSFile="http://www.stat.ucl.ac.be/R2HTML/Pastel.css")
      ####### CSS ####
      bootstrap1='<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">'
      bootstrap2='<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>'
      bootstrap3='<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>'
      
      HTML(bootstrap1,file=tmpfic)
      HTML(bootstrap2,file=tmpfic)
      HTML(bootstrap3,file=tmpfic)
      
      HTML('<div class="container-fluid">',file=tmpfic)
      
      ################
      HTML("<h1 class='text-center text-primary'> Module Results  </h1>",file=tmpfic)
      HTML('<br /><br />')
      HTMLInsertGraph(graph1,file=tmpfic)
      ##################
            
        
            if (hyper_geo_test_bool)
            {
                
                HTML('<hr class="col-xs-12">')
                HTML("<h2 class='text-center text-primary'> Hyper Geometric Test </h2>",file=tmpfic)
                HTML("<p> Conditioned on P-value <0.05 </p>",file=tmpfic)
                HTML('<div class="row">',file=tmpfic)
                if (nchar(descriptions)==0){
                  HTML("<h4 class='text-center text-danger'> Not enough for wordcloud </h4>",file=tmpfic)
                } 
                if (nchar(descriptions)>0)
                {
                  graph2=paste("./images","/",module_name2,"_WordCloud.png",sep = "")
                  HTMLInsertGraph(graph2,file=tmpfic)
                }
                
                
                if (length(gene_names)>0)
                {
                
                  HTML('<div class="col-sm-1">',file=tmpfic)
                  HTML('</div>',file=tmpfic)
                  HTML('<div class="col-sm-10">',file=tmpfic)
                  
                      table_command2=
                        '
                      <table class="table table-hover .table-striped table-bordered">
                      <thead>
                      <tr>
                      <th scope="col">Gene Names</th>
                      <th scope="col">Gene Description</th>
                      <th scope="col">Number of Overlapping Genes</th>
                      <th scope="col">Overlapping Genes Names</th>
                      <th scope="col">p-value</th>
                      <th scope="col">q-value</th>
                      
                      </tr>
                      </thead>
                      <tbody>
                      '
                      
                      module_hypo_table_header<-c('Gene-Names','Gene-Description','Number-of-Overlapping-Genes','Overlapping-Genes-Names','p-value','q-value')
                      module_hypo_table<-c()
                      ##################
                      descriptions=strsplit(descriptions,",")[[1]]
                      
                      
                      HTML(table_command2,file=tmpfic)
                      
                          for (kk in 1:length(gene_names))
                            {
                            link_command=paste("<a href=http://software.broadinstitute.org/gsea/msigdb/cards/",as.character(gene_names[kk]),".html>",gene_names[kk],'</a>',sep="")
                              HTML(paste('<tr>',
                                         '<td valign="middle">',link_command,'</td>',
                                         '<td valign="middle">',gsub(">"," ",gene_descriptions[kk]),'</td>',
                                         '<td valign="middle">',number_overlappings[kk],'</td>',
                                         '<td valign="middle">',gsub("  ",' ',as.character(overlapping_gene_names[kk])),'</td>',
                                         '<td valign="middle">',round(p_values[kk],4),'</td>',
                                         '<td valign="middle">',round(q_values[kk],4),'</td>',
                                         '</tr>'),file=tmpfic)
                              rr<-c(as.character(gene_names[kk]),gsub(">"," ",gene_descriptions[kk]),number_overlappings[kk],gsub("  ",' ',as.character(overlapping_gene_names[kk])),round(p_values[kk],4),round(q_values[kk],4))
                              module_hypo_table<-rbind(module_hypo_table,rr)
                          }
                      
                      
                      HTML('</tbody></table>',file=tmpfic)
                      colnames(module_hypo_table)<-module_hypo_table_header
                      outfile=paste('./report_html/htmls/tables/module_hyper_geo_test/Module',ModuleNr,'_hypergeometric_test.tsv',sep='')
                      write.table(module_hypo_table,file=outfile,sep='\t',quote=F,col.names=T,row.names=F)
                  
                    HTML('</div>',file=tmpfic)
                  
                    HTML('<div class="col-sm-1">',file=tmpfic)
                    
                    HTML('</div>',file=tmpfic)
                ##################
                HTML('</div>',file=tmpfic)
                }
          }
      
          HTML('<hr class="col-xs-12">')
          
          ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
          currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
          positiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] > 0)]
          negetiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] < 0)]
          module_regulators_data=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
          #colnames(module_regulators_data) <- c("Expression data")
          HTML("<h2 class='text-center text-primary'>Module Genes Expression Data </h2>",file=tmpfic)
          
          all_gene_expression_file_name_save=paste(module_name,"_","data",".csv",sep="")
          all_gene_expression_file_address=paste("./report_html/htmls/data",'/',all_gene_expression_file_name_save,sep="")
          
          write.csv(module_all_genes_data, file =all_gene_expression_file_address)
          
          ModuleData<-round(ModuleData,2)
          
          HTML(paste('<a href=', paste('./data','/',all_gene_expression_file_name_save,sep=""),'  download>',' download all module gene data ','</a>',sep=""))
          
          
          HTML(ModuleData,file=tmpfic)
          
          HTML("<h2 class='text-center text-primary'>Regulators</h2>",file=tmpfic)
          #############
          annotations_file_name_save=paste(module_name,"_","annotations",".csv",sep="")
          annotations_file_address=paste('./report_html/htmls/data','/',annotations_file_name_save,sep="")
          write.csv(module_annotations, file =annotations_file_address)
          
          HTML(paste('<a href=', paste('./data','/',annotations_file_name_save,sep=""),'  download>',' download annotations data ','</a>',sep=""))
          #############
          
          # colnames(module_regulators_data)<-""
          
          HTML(paste('<p class="text-success">',paste(as.character(positiveRegulators)),'</p>'),file=tmpfic)
          HTML(paste('<p class="text-danger">',paste(as.character(negetiveRegulators)),'</p>'),file=tmpfic)
      
      HTML('</div>',file=tmpfic)
    }
    
    
    
    ##############################################################################
    #Create the landing page
    ##############################################################################   
    tmpfic<-HTMLInitFile(address1,filename="index",Title = "Amartto Report",CSSFile="http://www.stat.ucl.ac.be/R2HTML/Pastel.css")
    bootstrap1='<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">'
    bootstrap2='<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>'
    bootstrap3='<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>'
    
    HTML(bootstrap1,file=tmpfic)
    HTML(bootstrap2,file=tmpfic)
    HTML(bootstrap3,file=tmpfic)
    
    HTML('<div class="container-fluid">',file=tmpfic)
    
    
    #######################################  Create the TEXT  #########################
    HTML('<h1 class="text-primary text-center"> AMARETTO results  </h1>',file=tmpfic)
    HTML('<br /><br /><br /><br /><br />')
    #######################################  Create the table  #########################
    table_command5=
      '
    <table class="table table-hover.table-striped table-bordered">
    <thead>
    <tr>
    <th scope="col"># of samples</th>
    <th scope="col"># of modules</th>
    <th scope="col"> Var-Percentage</th>
    </tr>
    </thead>
    <tbody>
    '
    # HTML(table_command5,file=tmpfic)
    # HTML('<tr>',file=tmpfic)
    # HTML(paste('<td>',as.character(number_of_samples),'</td>'),file=tmpfic)
    # HTML(paste('<td>',as.character(NrModules),'</td>'),file=tmpfic)
    # HTML(paste('<td>',as.character(number_of_regulators),'</td>'),file=tmpfic)
    # #HTML(paste('<td>',as.character(number_of_samples),'</td>'),file=tmpfic)
    # HTML('</tr>',file=tmpfic)
    
    
    table_command1=
      '
    
    <table class="table table-hover ">
    <thead ">
    <tr>
    <th scope="col" class="align-middle">Module #</th>
    <th scope="col" class="align-middle"># of target genes</th>
    <th scope="col" class="align-middle"># of regulator genes</th>
    <th scope="col" class="align-middle"># of significant gene overlappings</th>
    </tr>
    </thead>
    <tbody>
    '
    
    ModuleNr<-1
    ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
    number_of_samples=length(colnames(ModuleData))
    
    
    
    
    HTML('<div class="col-sm-3">',file=tmpfic)
    HTML('</div>',file=tmpfic)
    HTML('<div class="col-sm-6">',file=tmpfic)
    
    HTML(paste('<p class=".text-success text-right"> # of samples = ',as.character(number_of_samples),'</p>'))
    HTML('<br /><br />')
    
    HTML(table_command1,file=tmpfic)
    
    
    amaretto_result_table_header<-c('Module_No','number_of_target_genes','number_of_regulator_genes','number_of_significant_gene_overlappings')
    amaretto_result_table<-c()
    
        for (ModuleNr in 1:NrModules )
        {
          module_name=paste("module",as.character(ModuleNr),sep="")
          ###################### find module info
          ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
          currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
          RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
          module_regulators_weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
          ModuleGenes=rownames(ModuleData)
          RegulatoryGenes=rownames(RegulatorData)
          
          
          
          number_of_genes=length(rownames(ModuleData))
          number_of_regulators=length(currentRegulators)
          number_of_samples=length(colnames(ModuleData))
          
          #########################  creating Link for each module ####################
          address='./htmls'
          htmladdress=paste("'",address,"/module",as.character(ModuleNr),".html","'",sep="")
          
  
          link_command=paste("<a href=",htmladdress,'>',module_name,'</a>',sep="")
          #HTML(link_command,file=tmpfic)
          ###########################################################################
          HTML('<tr>',file=tmpfic)
          HTML(paste('<td class="align-middle">',link_command,'</td>'),file=tmpfic)
          HTML(paste('<td class="align-middle">',as.character(number_of_genes),'</td>'),file=tmpfic)
          HTML(paste('<td class="align-middle">',as.character(number_of_regulators),'</td>'),file=tmpfic)
          HTML(paste('<td class="align-middle">',as.character( number_of_significant_gene_overlappings[ModuleNr]),'</td>'),file=tmpfic)
         
          #HTML(paste('<td>',as.character(number_of_samples),'</td>'),file=tmpfic)
          HTML('</tr>',file=tmpfic)
          rr<-c(module_name,as.character(number_of_genes),as.character(number_of_regulators),as.character( number_of_significant_gene_overlappings[ModuleNr]))
          amaretto_result_table<-rbind(amaretto_result_table,rr)
        }
    
    colnames(amaretto_result_table)<-amaretto_result_table_header
    outfile=paste('./report_html/htmls/tables/amaretto','.tsv',sep='')
    write.table(amaretto_result_table,file=outfile,sep='\t',quote=F,col.names=T,row.names=F)
    
    HTML('</tbody></table>',file=tmpfic)
    HTML('</div>',file=tmpfic)
    HTML('<div class="col-sm-3">',file=tmpfic)
    HTML('</div>',file=tmpfic)
    ####################################### ####################################### 
    
    HTML('</div>',file=tmpfic)
    # HTMLEndFile()
    #################################################################################[####
    zip(zipfile = 'reportZip', files = './report_html')
    ###############################
}



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
create_gene_annotations<-function(module_all_genes_data,Module_Genes_names,Module_regulators_weights)
{
  all_genes_names=rownames(module_all_genes_data)
 
  targets_bool<-c()
  regulators_bool<-c()
  regulators_weight<-c()

  for (i in 1:length(all_genes_names))
  {
    
    gene_name=all_genes_names[i]

    
    a=0
    b=0
    c=0
    
    
    if (is.element(gene_name, Module_Genes_names))
    {
      
      a<-1
      
    }
    
    
    
    if (is.element(gene_name,  rownames(Module_regulators_weights)))
    {
      b<-1
      c<-Module_regulators_weights$module_regulators_weights[rownames(Module_regulators_weights)==gene_name]
    }
    
    targets_bool<-c(targets_bool,a)
    regulators_bool<-c(regulators_bool,b)
    regulators_weight<-c(regulators_weight,c)
  }
  
  df=data.frame(all_genes_names, targets_bool, regulators_bool,regulators_weight)  
  return(df)
}

