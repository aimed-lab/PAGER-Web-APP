## server.R
source("lib/PAGER.R")
source("lib/beere.R")

#sudo apt-get install -y  libudunits2-dev libgdal-dev libgeos-dev libproj-dev
#sudo apt-get install libgdal-dev



# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
  # linux
  setwd(system("pwd", intern = T) )
  t <- list(family = 'Arial',
            size = 20, color = '#424B54')
  #User Interface functions.
  NumericCols <- NULL
  FactorCols <- NULL  
  
  ### event trigger functions ###
  observeEvent(input$loadbtn, {
    updateTabItems(session, "tabs", "filterTab")
  })
  observeEvent(input$uploadbtn, {
    updateTabItems(session, "tabs", "filterTab")
  })
  observeEvent(input$filterbtn, {
    updateTabItems(session, "tabs", "enrichTab")
  })
  observeEvent(input$filterbackbtn, {
    updateTabItems(session, "tabs", "loadTab")
  })
  observeEvent(input$enrichbtn, {
    updateTabItems(session, "tabs", "TanalysisTab")
  })
  observeEvent(input$enrichbackbtn, {
    updateTabItems(session, "tabs", "filterTab")
  })
  observeEvent(input$Tanalysisbackbtn, {
    updateTabItems(session, "tabs", "enrichTab")
  }) 
  
  ### end event trigger functions ###
  
  ##### data #####
  ### Step 1. load data ###
  input_dir = "dataset/"
  upload_df <- reactive({
    upload_df <- read.csv(
      input$DEtable$datapath, 
      header = input$fheader1,
      sep = "\t"
    )
    return(upload_df)
  })
  

  #server
  observe({ 
    toggle(id="uploadbtn", condition=!is.null(input$DEtable))
  })
  observeEvent(
    input$DEtable,
    {
      updateSelectInput(
        session,
        "symbolcolumn",
        choices = names(upload_df()),
        selected = names(upload_df())[1]
      )
      updateSelectInput(
        session,
        "log2FCcolumn",
        choices = names(upload_df()),
        selected = names(upload_df())[2]
      )
      updateSelectInput(
        session,
        "pvaluecolumn",
        choices = names(upload_df()),
        selected = names(upload_df())[3]
      )       
    }
  )
  
  rv <- reactiveValues(
    clicked = NULL
  )
  
  # If user moved slider say "slider" 
  observeEvent({
    input$loadbtn
    # add other sliders here
  }, {
    rv$clicked <- NULL  # reset if previous "loadbtn"
    rv$clicked <- "loadbtn"
    cat(paste0("individual: ",rv$clicked,"\n"))
  }, ignoreInit = TRUE) 
  
  # If user clicked radio button say "radio"
  observeEvent({
    input$uploadbtn
  }, 
  {
    rv$clicked <- NULL  # reset if previous "uploadbtn"
    rv$clicked <- "uploadbtn"
    cat(paste0("individual: ",rv$clicked,"\n"))
  }, ignoreInit = TRUE) 
  

  
  df_DEGs <- eventReactive({
    input$uploadbtn
    input$loadbtn
    1
  },{
    if(rv$clicked == 'loadbtn'){
      cat("loadbtn has been clicked\n")
    }else{
      cat("uploadbtn has been clicked\n")
    }
    reference = read.csv(paste0(input_dir,"HUMAN_PROTEIN_CODING_GENE.txt"), header = TRUE,sep="\t")
    reference = levels(reference$SYMBOL)[reference$SYMBOL]
    #cat(print(reference))
    if(rv$clicked == 'loadbtn'){
      if (input$diseaseSelect == 'Melanoma'){
        df <- read.csv(paste0(input_dir,input$diseaseSelect,".txt"), header = TRUE,sep="\t")
        res_DESeq2 = df
        reference = res_DESeq2$symbol
      }
      else{
        df <- read.csv(paste0(input_dir,input$diseaseSelect,".txt"), header = TRUE,sep="\t")
        colnames(df) = c('symbol','log2FoldChange','pvalue')
        res_DESeq2 = df
      }
    }
    
    else if(rv$clicked == 'uploadbtn'){
      upload_df = upload_df()
      names(upload_df)[names(upload_df) == input$symbolcolumn] <- 'symbol'
      names(upload_df)[names(upload_df) == input$log2FCcolumn] <- 'log2FoldChange'
      names(upload_df)[names(upload_df) == input$pvaluecolumn] <- 'pvalue'
      res_DESeq2 = upload_df      
    }
    res_DESeq2 = transform(res_DESeq2, pvalue=as.numeric(as.character(pvalue)))
    res_DESeq2%>%mutate(
      pvalue = ifelse(pvalue==0,10^-20,pvalue)
    ) %>% as.data.frame
    #res_DESeq2 = res_DESeq2%>%mutate(threshold = ifelse(log2FoldChange >= l2FC & pvalue<=p_cut,"UP", ifelse(log2FoldChange<=-l2FC & pvalue<=p_cut, "DN", "None"))) 
    
    #cat(print(reference))
    return(
      list(
        res_DESeq2 = res_DESeq2,
        reference = reference
      )
    )
  })

  ##########

  ### Step 2. filter data ###
  

  #v <- reactiveValues(data = NULL)
  df_DEGs_updated <- eventReactive({
    input$loadbtn
    input$uploadbtn
    input$sliderl2FC
    input$sliderp_cut
    1
  },{
    res_DESeq2 = df_DEGs()[["res_DESeq2"]]
    #cat(print(res_DESeq2))
    l2FC = input$sliderl2FC
    p_cut =2^(-input$sliderp_cut)
    res_DESeq2 = res_DESeq2%>%mutate(threshold = ifelse(log2FoldChange >= l2FC & pvalue<=p_cut,"UP", ifelse(log2FoldChange<=-l2FC & pvalue<=p_cut, "DN", "None")))
    res_DESeq2_filtered = res_DESeq2[res_DESeq2$threshold%in%c('UP','DN'),]
    
    #max_l2fc <- floor(max(abs(res_DESeq2$log2FoldChange[!is.na(res_DESeq2$log2FoldChange)])))
    return(
      list(
        res_DESeq2=res_DESeq2,
        res_DESeq2_filtered=res_DESeq2_filtered
      )
    )
  })    
  

  ### Step 3. enrich ###
  df_PAGER_GS <- eventReactive({
    input$apply_btn
    input$filterbtn
    1
    }
    ,{
      res_DESeq2_filtered = df_DEGs_updated()$res_DESeq2_filtered
      source = isolate({input$selectPAGERsource})
      olap = isolate({input$sliderPAGERovlp})
      sim = isolate({input$sliderPAGERsimilarity})
      
      nlogpval = 2^(-isolate({input$sliderPAGERnlogpval}))
      res_PAGER_GS = pathFun(markers = res_DESeq2_filtered$symbol,
                             source = source,similarity = sim,
                             overlap = olap,FDR = nlogpval)
      res_PAGER_GS = res_PAGER_GS%>%mutate(
        SIMILARITY_SCORE = round(as.double(SIMILARITY_SCORE), 4), 
        PAGER_LINK=paste0('http://discovery.informatics.uab.edu/PAGER/index.php/geneset/view/',GS_ID),
        pvalue = signif(res_PAGER_GS$pvalue, digits=3),
        FDR = signif(res_PAGER_GS$pFDR, digits=3)
      )
      res_PAGER_GS = res_PAGER_GS[,c('GS_ID','NAME','SOURCE','OLAP','SIMILARITY_SCORE','pvalue','FDR','PAGER_LINK','LINK','TYPE','GS_SIZE')]
      return(res_PAGER_GS)
  })
  
  df_mtype <- eventReactive({
    input$apply_btn
    input$filterbtn
    1
  }
  ,{
    res_path = df_PAGER_GS()
    m_type = pathInt(PAGER = res_path$GS_ID)$data
    m_type = m_type %>% mutate(
      SIMILARITY = signif(as.double(m_type$SIMILARITY),digits = 3),
      nlogPvalue = signif(as.double(m_type$PVALUE),digits = 3)
    )
    m_type = m_type%>%subset(select=c("GS_A_ID",'GS_B_ID','OLAP','SIMILARITY','nlogPvalue'))
    m_type = m_type[m_type$GS_A_ID!=m_type$GS_B_ID,]
    
    return(m_type)
  })
  
  graph_mtype <- eventReactive({
    input$apply_btn
    input$filterbtn
    1
  }
  ,{
    m_type = df_mtype()
    edgelist = cbind(m_type$GS_A_ID,m_type$GS_B_ID)
    graph <- graph_from_data_frame(edgelist)
    return(graph)
  })
  
  pd_mtype <- eventReactive({
    input$apply_btn
    input$filterbtn
    1
  }
  ,{
    graph = graph_mtype()
    clusterlouvain <- cluster_louvain(graph)
    res_path_mtype = res_path[v_orders,]
    res_path_mtype['cluster'] = clusterlouvain$membership
    return(res_path_mtype)
  })  
  
  df_enrichR_GS <- eventReactive({
    input$enrichR_apply_btn
    input$filterbtn
    1
    #input$enrichR
  }
  ,{
    res_DESeq2_filtered = df_DEGs_updated()$res_DESeq2_filtered
    src_enrichR = isolate({input$selectenrichRsource})
    sliderenrichRovlp = isolate({input$sliderenrichRovlp})
    sliderenrichRnlogpval = 2^(-isolate({input$sliderenrichRnlogpval}))
    res_enrich <- enrichr(genes = as.vector(res_DESeq2_filtered$symbol), databases = src_enrichR)
    res_enrichR = NULL
    for(idx in 1:length(res_enrich)){
      enRdata = res_enrich[[idx]]
      source = names(res_enrich)[idx]
      if(length(enRdata)>0){
        res_enrichR = rbind(res_enrichR, cbind(enRdata,source))
      }
    }
    res_enrichR = cbind(
      res_enrichR,
      ovlap=as.numeric(gsub("(\\/[0-9]+$)","",res_enrichR[,2]))
    )
    res_enrichR = as.data.frame(res_enrichR)
    res_enrichR_filtered = res_enrichR[res_enrichR$Adjusted.P.value<=sliderenrichRnlogpval & res_enrichR$ovlap>=sliderenrichRovlp ,] #
    res_enrichR_filtered = res_enrichR_filtered%>%mutate(
      P.value = signif(res_enrichR_filtered$P.value, digits=3),
      Adjusted.P.value = signif(res_enrichR_filtered$Adjusted.P.value, digits=3)
    )
    return(res_enrichR_filtered)
  })
  
  df_webgestaltR_GS <- eventReactive({
    input$webgestaltR_apply_btn
    input$filterbtn
  },{
    gestaltR_pathway_sources = input$selectwebgestaltRsource
    sliderwebgestaltRovlp = input$sliderwebgestaltRovlp
    sliderwebgestaltRnlogpval = 2^(-isolate({input$sliderwebgestaltRnlogpval}))
    
    geneList = df_DEGs_updated()$res_DESeq2_filtered$symbol
    referenceList = df_DEGs()[['reference']]
    
    # get the list of available genesets
    res_webgestaltR = WebGestaltR::WebGestaltR(
      fdrThr = sliderwebgestaltRnlogpval,
      enrichDatabase = gestaltR_pathway_sources,
      interestGene = as.vector(geneList),interestGeneType="genesymbol",
      referenceGene = as.vector(referenceList),referenceGeneType = "genesymbol",
      isOutput=FALSE#,
      #minNum=2,maxnum=3000
    )
    if(!is.null(res_webgestaltR)){
      res_webgestaltR_filtered = res_webgestaltR[
        res_webgestaltR$overlap > sliderwebgestaltRovlp &
        res_webgestaltR$FDR < sliderwebgestaltRnlogpval,
      ]
    }else{
      res_webgestaltR_filtered=data.frame(matrix(ncol = 12, nrow = 0))
      names(res_webgestaltR_filtered) = c('geneSet','description','link','size','overlap','expect','enrichmentRatio','pValue','FDR','overlapId','database','userId')
    }
    res_webgestaltR_filtered = res_webgestaltR_filtered%>%mutate(
      pValue = signif(res_webgestaltR_filtered$pValue, digits=3),
      FDR = signif(res_webgestaltR_filtered$FDR, digits=3)
    )
    return(res_webgestaltR_filtered)
  })
  
  ### step 4. similarity terms ###
  res_similarity_pair_filtered <- eventReactive({
    input$sliderSimilairtyScore
    input$enrichbtn
  },{
    # input #
    res_path = df_PAGER_GS()
    res_enrichR_filtered = df_enrichR_GS()
    res_webgestaltR_filtered = df_webgestaltR_GS()
    slidersimilarity = input$sliderSimilairtyScore
    ### distance-based similarity of names ###
    # filter common words or IDs

    IDs = c(' Homo sapiens',' R-HSA-[0-9]+$',' WP[0-9]+','hsa[0-9]+ ',' - Homo sapiens \\(human\\)',' $')
    IDs_str = paste0(IDs,collapse = "|")
    common_bow = c(
      #' signaling',' pathway',' in',
      #' by',' pathways',' metabolism',
      '\ $','\\','-',' h il1rPathway')#
    bow_str = paste0(common_bow,collapse = "|")
    
    res_path$NAME = gsub(IDs_str,"\\1",res_path$NAME,ignore.case = T)
    pager_names = tolower(gsub(paste0(bow_str),"\\1",res_path$NAME,ignore.case = T))
    res_enrichR_filtered$Term = gsub(IDs_str,"\\1",res_enrichR_filtered$Term,ignore.case = T)
    enirchR_Terms = tolower(gsub(paste0(bow_str),"\\1",res_enrichR_filtered$Term,ignore.case = T))
    webgestaltR_descs = tolower(gsub(paste0(bow_str),"\\1",res_webgestaltR_filtered$description,ignore.case = T))
    
    # compare enrichR and pager
    P_vs_E_pair = DataFrame(matrix(ncol=3,nrow=0))
    names(P_vs_E_pair) = c("pager","enrichR","P_vs_E")
    P_vs_E_mat_osa = stringdist::stringsimmatrix(pager_names,enirchR_Terms,method = 'osa')
    P_vs_E_mat = P_vs_E_mat_osa
    rownames(P_vs_E_mat) = pager_names
    colnames(P_vs_E_mat) = enirchR_Terms
    
    score = apply(P_vs_E_mat,1,function(x) x[which(x==max(x))][1])
    A = names(score)
    B = as.character(apply(P_vs_E_mat,1,function(x) names(x)[which(x==max(x))][1]))
    score_pair = data.frame(
      A,B,score
    )
    names(score_pair) = c("pager","enrichR","P_vs_E")
    P_vs_E_pair_filtered = unique(score_pair[score_pair$P_vs_E>=slidersimilarity,])
    
    
    # compare enrichR and webgestaltR
    W_vs_E_pair = DataFrame(matrix(ncol=3,nrow=0))
    names(W_vs_E_pair) = c("webgestaltR","enrichR","W_vs_E")
    W_vs_E_mat_osa = stringdist::stringsimmatrix(webgestaltR_descs,enirchR_Terms,method = 'osa')
    W_vs_E_mat = W_vs_E_mat_osa
    
    rownames(W_vs_E_mat) = webgestaltR_descs
    colnames(W_vs_E_mat) = enirchR_Terms
    score = apply(W_vs_E_mat,1,function(x) x[which(x==max(x))][1])
    A = names(score)
    B = as.character(apply(W_vs_E_mat,1,function(x) names(x)[which(x==max(x))][1]))
    score_pair = data.frame(
      A,B,score
    )
    names(score_pair) = c("webgestaltR","enrichR","W_vs_E")
    W_vs_E_pair_filtered = unique(score_pair[score_pair$W_vs_E>=slidersimilarity,])
    
    
    # compare pager and webgestaltR
    W_vs_P_pair = DataFrame(matrix(ncol=3,nrow=0))
    names(W_vs_P_pair) = c("webgestaltR","pager","W_vs_P")
    
    W_vs_P_mat_osa = stringdist::stringsimmatrix(webgestaltR_descs,pager_names,method = 'osa')
    W_vs_P_mat = W_vs_P_mat_osa
    
    rownames(W_vs_P_mat) = webgestaltR_descs
    colnames(W_vs_P_mat) = pager_names
    score = apply(W_vs_P_mat,1,function(x) x[which(x==max(x))][1])
    A = names(score)
    B = as.character(apply(W_vs_P_mat,1,function(x) names(x)[which(x==max(x))][1]))
    score_pair = data.frame(
      A,B,score
    )
    names(score_pair) = c("webgestaltR","pager","W_vs_P")
    W_vs_P_pair_filtered = unique(score_pair[score_pair$W_vs_P>=slidersimilarity,])
    
    ovlaps = merge(W_vs_P_pair_filtered,P_vs_E_pair_filtered,by.x = "pager", by.y = "pager",all=F)
    ovlaps = merge(ovlaps,W_vs_E_pair_filtered[c('enrichR','W_vs_E')],by.x = "enrichR", by.y = "enrichR",all=F)
    ovlaps = unique(ovlaps)
    ovlaps = ovlaps[,c('pager','enrichR','webgestaltR','W_vs_P','P_vs_E','W_vs_E')]
    
    if(nrow(ovlaps[c('webgestaltR','pager','W_vs_P')])>0){
      W_vs_P_uniq = anti_join(W_vs_P_pair_filtered, ovlaps[c('webgestaltR','pager','W_vs_P')])
      W_vs_E_uniq = anti_join(W_vs_E_pair_filtered, ovlaps[c('webgestaltR','enrichR','W_vs_E')])
      P_vs_E_uniq = anti_join(P_vs_E_pair_filtered, ovlaps[c('pager','enrichR','P_vs_E')])
    }else{
      W_vs_P_uniq = W_vs_P_pair_filtered
      W_vs_E_uniq = W_vs_E_pair_filtered
      P_vs_E_uniq = P_vs_E_pair_filtered
    }
    
    webgestaltR_uniq = setdiff(
      setdiff(
        setdiff(
          webgestaltR_descs,
          as.character(W_vs_P_uniq$webgestaltR)
        ),
        as.character(W_vs_E_uniq$webgestaltR)
      ),
      as.character(ovlaps$webgestaltR)
    )
    enrichR_uniq = setdiff(
      setdiff(
        setdiff(
          enirchR_Terms,
          as.character(W_vs_E_uniq$enrichR)
        ),
        as.character(P_vs_E_uniq$enrichR)
      ),
      as.character(ovlaps$enrichR)
    )
    PAGER_uniq = setdiff(
      setdiff(
        setdiff(
          pager_names,
          as.character(W_vs_P_uniq$pager)
        ),
        as.character(P_vs_E_uniq$pager)
      ),
      as.character(ovlaps$pager)
    )
    
    return(
      list(
        consensus=ovlaps,
        W_vs_P_uniq = W_vs_P_uniq,
        W_vs_E_uniq = W_vs_E_uniq,
        P_vs_E_uniq = P_vs_E_uniq,
        PAGER_uniq = data.frame(PAGER = PAGER_uniq),
        enrichR_uniq = data.frame(enrichR = enrichR_uniq),
        webgestaltR_uniq = data.frame(webgestaltR = webgestaltR_uniq)
      )
    )
  })
  
  ### end data ###
  
  ##################
  ##### output #####
  # step 1. load page #
  output$diseaseSelect <- renderText({
    paste0("You've selected ",input$diseaseSelect," dataset.")
  })
  
  output$tableCancer_list <- renderDT({
    DT::datatable(
      read.csv(paste0(input_dir,"ListOfAbbr.txt"), header = TRUE,sep="\t"),
      extensions = c('Responsive'), 
      class = 'cell-border stripe',
    )
  })
  
  # step 2. filter page #
  output$vulcanoPlot <- renderPlotly({
    #req(input$l2FC)
    #req(input$p_cut)
    res_DESeq2 = df_DEGs_updated()[['res_DESeq2']]
    volcano = ggplot(
      data = res_DESeq2, 
      aes(x=log2FoldChange, y=-log2(pvalue),
      text = paste('catergory:', threshold,
        '<br>symbol:', symbol,
        '<br>log2FC:', log2FoldChange,
        '<br>-log2(pvalue):', floor(-log2(pvalue))                                     
        )
      ),
      xlab="log2 fold change", ylab="-log2 p-value",
    )+
    #xlim(c(-10,10))+
    theme_bw() +
      geom_point(aes(colour = threshold), size=0.5) + 
      scale_colour_manual(values = c("UP"= "red", "DN"="blue",  "None"= "grey"),labels= c("UP"= "UP", "DN"="DN",  "None"= "None"))
    
    fig = plotly::ggplotly(volcano,tooltip = "text")
  })
  output$tableDEGs_updated <- renderDT({
    DT::datatable(
      df_DEGs_updated()$res_DESeq2_filtered,
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })
  ###
  
  # step 3. enrichment and histogram #
  ## PAGER ##
  output$tablePAGER <- renderDT({
    datatable(elementId = "tablePAGER",
      df_PAGER_GS(),
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })
  output$PAGERpie <- renderPlotly({
    res_path = df_PAGER_GS()
    res_source_count = data.frame(source = names(table(res_path$SOURCE)),frequency = as.vector(table(res_path$SOURCE)))
    fig <- plot_ly(
      textinfo = 'label+percent+value',
      res_source_count, labels = ~source, values = ~frequency, type = 'pie'
    ) 
    fig <- fig %>% layout(
      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
    )
  })
  output$PAGERhist <- renderPlotly({
    res_path = df_PAGER_GS()
    path_histo = ggplot(
      data = res_path, 
      aes(x=reorder(NAME, -log2(FDR)), 
          y=-log2(FDR),color = SOURCE,
      text = paste('PAG_ID', GS_ID,
      '<br>Name: ', reorder(NAME, -log2(FDR)),
      '<br>-log2(FDR):', floor(-log2(FDR))
    ))) +
      geom_bar(stat='identity') + theme_bw()+ ylab("-log2 FDR")+xlab('PAGs')+
      theme(axis.text.x = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
            axis.title.x = element_text(color = "grey20", size = 5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
      coord_flip()
    plotly::ggplotly(path_histo,tooltip = "text")
    
  })
  
  output$mtype <- renderPlotly({
    graph = graph_mtype()
    res_path = df_PAGER_GS()
    #cat(print(m_type))
    set.seed(7)
    L <- layout_nicely(graph)
    vs <- V(graph)
    es <- as.data.frame(get.edgelist(graph))
    Ne <- length(es[1]$V1)
    Xn <- L[,1]
    Yn <- L[,2]
    v_orders = unlist(lapply(names(vs), function(i) which(res_path$GS_ID == i)))
    network <- plot_ly(
      #width = -log(res_path[v_orders,]$pFDR)/log(10),
      #height = -log(res_path[v_orders,]$pFDR)/log(10),
      color = factor(res_path[v_orders,]$SOURCE),
      type = "scatter", 
      x = Xn, 
      y = Yn, 
      mode = "marker+text",
      marker = list(size = ~-log(res_path[v_orders,]$FDR)/log(10), sizeref = 4000, sizemode = 'area'),
      text = res_path[v_orders,]$NAME, 
      #label = res_path[v_orders,]$NAME, 
      hovertemplate = paste(
        "<b>%{text}</b><br><br>",
        "<extra></extra>"
      )
      
    )
    edge_shapes <- list()
    for(i in 1:Ne) {
      v0 <- es[i,]$V1
      v1 <- es[i,]$V2
      
      edge_shape = list(
        type = "line",
        line = list(color = "#030303", width = 0.3),
        x0 = Xn[match(v0,names(vs))],
        y0 = Yn[match(v0,names(vs))],
        x1 = Xn[match(v1,names(vs))],
        y1 = Yn[match(v1,names(vs))],
        opacity = 1
      )
      
      edge_shapes[[i]] <- edge_shape}
    
    network <- layout(
      network,
      shapes = edge_shapes,
      xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
      height = 1500
    )
    
    return(network)
  })
  
  
  
  output$mtype_table <- renderDT({
    datatable(elementId = "mtype_table",
              df_mtype(),
              extensions = c('Responsive','Buttons'), 
              class = 'cell-border stripe',
              options = list(
                pageLength = 5,
                buttons = c('copy', 'csv', 'excel'),
                dom = 'Blfrtip',
                responsive = TRUE
              )
    )  
  })
  
  ## enrichR ##
  output$tableenrichR_GS <- renderDT({
    datatable(
      df_enrichR_GS(),      
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })  
  output$enrichRhist <- renderPlotly({
    res_enrichR_filtered = df_enrichR_GS()
    enrichRpath_histo = ggplot(
      data = res_enrichR_filtered, 
       aes(x=reorder(Term, -log2(Adjusted.P.value)), 
           y=-log2(Adjusted.P.value),color = source,
           text = paste('Name: ', reorder(Term, -log2(Adjusted.P.value)),
           '<br>-log2(FDR):', floor(-log2(Adjusted.P.value))
      ))) +
      geom_bar(stat='identity') + theme_bw()+ ylab("-log2 FDR")+xlab('Genesets')+
      theme(
        axis.text.x = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain")
      )+
      coord_flip()
    plotly::ggplotly(enrichRpath_histo,tooltip = "text")
  })
  
  ## webgestaltR ##
  output$tablewebgestaltR <- renderDT({
    datatable(
      df_webgestaltR_GS(),
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })
  
  output$webgestaltRhist <- renderPlotly({
    
    res_webgestaltR_filtered = df_webgestaltR_GS()
    if(nrow(res_webgestaltR_filtered)>0){
    webgestaltR_histo = ggplot(
      data = res_webgestaltR_filtered, 
      aes(x=reorder(description, -log2(FDR)), 
          y=-log2(FDR),color = database,
          text = paste('Name: ', reorder(description, -log2(FDR)),'<br>-log2(FDR):', floor(-log2(FDR))
          ))) +
      geom_bar(stat='identity') + theme_bw()+ ylab("-log2 FDR")+xlab('Genesets')+
      theme(
        axis.text.x = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 5, angle = 90, hjust = .5, vjust = .5, face = "plain")
      )+
      coord_flip()
    plotly::ggplotly(webgestaltR_histo,tooltip = "text")
    }
  })


  
  
  # step 4. venn diagram #
  output$PAGER_uniq <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['PAGER_uniq']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  }) 
  output$enrichR_uniq <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['enrichR_uniq']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })   
  output$webgestaltR_uniq <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['webgestaltR_uniq']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })   
  output$W_vs_P_uniq <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['W_vs_P_uniq']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  }) 
  output$W_vs_E_uniq <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['W_vs_E_uniq']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })   
  output$P_vs_E_uniq <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['P_vs_E_uniq']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  }) 
  
  output$Consensus <- renderDT({
    datatable(
      res_similarity_pair_filtered()[['consensus']],
      extensions = c('Responsive','Buttons'), 
      class = 'cell-border stripe',
      options = list(
        pageLength = 5,
        buttons = c('copy', 'csv', 'excel'),
        dom = 'Blfrtip',
        responsive = TRUE
      )
    )
  })
  output$venn <- renderPlot({

    res_similarity_pair_filtered = res_similarity_pair_filtered()
    enrichR_num = nrow(res_similarity_pair_filtered$enrichR_uniq)
    PAGER_num = nrow(res_similarity_pair_filtered$PAGER_uniq)
    webgestaltR_num = nrow(res_similarity_pair_filtered$webgestaltR_uniq)
    ovlaps = res_similarity_pair_filtered[['consensus']]
    W_vs_P_uniq_num = nrow(res_similarity_pair_filtered[['W_vs_P_uniq']])
    W_vs_E_uniq_num = nrow(res_similarity_pair_filtered[['W_vs_E_uniq']])
    P_vs_E_uniq_num = nrow(res_similarity_pair_filtered[['P_vs_E_uniq']])
    cons_num = nrow(ovlaps)
    W_vs_P_num = W_vs_P_uniq_num+nrow(ovlaps)
    W_vs_E_num = W_vs_E_uniq_num+nrow(ovlaps)
    P_vs_E_num = P_vs_E_uniq_num+nrow(ovlaps)
    
    fig = draw.triple.venn(
      area1 = PAGER_num+W_vs_P_uniq_num+P_vs_E_uniq_num+cons_num,
      area2 = enrichR_num+P_vs_E_uniq_num+W_vs_E_uniq_num+cons_num,
      area3 = webgestaltR_num+W_vs_P_uniq_num+W_vs_E_uniq_num+cons_num,
      n12 = P_vs_E_num,
      n23 = W_vs_E_num,
      n13 = W_vs_P_num,
      n123 = cons_num,
      category = c("PAGER", "EnrichR", "WebGestaltR"),
      fill = c("blue2","red2", "green4"),
      lty = "blank",
      cex = 2,
      cat.cex = 2,
      print.mode = c("percent","raw"),
      sigdig=2, ind = T,scaled = T,
      #cat.col = c("blue", "red", "green")
    )
  })
}