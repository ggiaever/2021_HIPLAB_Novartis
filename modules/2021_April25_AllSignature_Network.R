library(shiny)
library(dplyr)
library(DT)
library(shinydashboard)
library(visNetwork)
library(igraph)
library(shinyjs)

source(file.path("modules","2021_April25_GOenrichment_functions.R"))
source(file.path("modules" ,"2021_April25_functions.R"))
source(file.path("modules/2021_April3_functions.R"))

dfsig   <-  read.delim("2021_April23_data/mkq5_slim.txt",check.names = F,stringsAsFactors = F)

fdat    <-  read.delim("2021_April23_data/december4_fdat_updated.txt",stringsAsFactors = F,check.names = F)
noness  <-  fdat  %>% filter(essential == "noness")
ess     <-  fdat  %>% filter(essential == "ess")


sigNetworkModuleUI <- function(id, label = "sigNetworkModule") {;
   ns <- NS(id)
   shiny_wdth = ns("shiny_wdth")
  
  tagList(
    fluidPage(
     fluidRow(
        column(width  = 6,
               ############### this is the query set ############
               box(title = "Select response signature:",
                   selectizeInput(ns("respMOD"),"",width = "400px", choices = sort(unique(dfsig$signature)), selected = "TOR signaling; rapamycin:18"), #,multiple = F,selected = "iron chelators"),
                   status = "primary", solidHeader = T,width="100%",height = 150)
              ),
        column(width  = 3,       
               box(title = "Reset response menu:",
                   br(),
                   fluidRow(
                     column(align ="center",
                            width = 12,
                            actionButton(ns("resetSig"), label = "respMenu"),
                            
                     )),
                     
                   solidHeader = TRUE, width = "100%",height = 150, status = "primary")#box
               ),
        column(width  = 3,
               box(title= "Download signatures:", align = "center",
                   br(),
                downloadButton(ns("downloadSignatures"), "signatures"),
                   solidHeader = T,status = "primary", width="100%", height = 150)
               ),
        ),
        
       fluidRow(
         column(width = 3,
               
               box(title = "Set fitness score threshold:",
                   sliderInput(ns("scorethreshSIGNAT"), label = "", min = 0, max = 5,
                               value = 1, step = 0.5),
               
                   status = "primary", solidHeader = T,width="100%", height = 150)
        ),
        column(width  = 3,
               box(title = "Set FDR threshold:",
                   sliderInput(ns("fdrSIGNAT"),label = "", min = 0, max = 0.5,
                               value = 0.2, step = 0.05), 
                   status = "primary", solidHeader = T,width="100%", height = 150)
               
        ),
        column(width = 3,  
               box(title= "Background genome:", align = "center",
                   prettyRadioButtons(ns("poolset"), label = "",
                                      choices = c("HIP" = "essential", "HOP" = "nonessential","HIPHOP" = "all"),
                                      status = "primary", shape = "square", animation = "tada",
                                      bigger = T,selected = "all", inline = T),
                   solidHeader = T,status = "primary",width = "100%",height = 150)
               
        ),
        column(width  = 3,
               box(title = "Compare signatures tab:",
                   br(),
                   column(align = "center",width  = 12,
                   actionButton(ns("compResp"), label = "Compare signatures"),
                   ),
                   status = "primary",solidHeader = T, width = "100%", height = 150)
               
               
        )
        ), 
      
      
      fluidRow(
        
        box(title = "GO enrichment analysis of the response signature:",
            
       
            visNetworkOutput(ns("network_proxy"), width = "100%",height = 855),
              
            width = 9, status = "primary", solidHeader = TRUE,height = 1000),
        
        box(title = "Enrichment details:",width = 3,
            uiOutput(ns("Gotable")),solidHeader = T,status = "primary",background = "navy", height = 270),
        
        
        box(title = "Top-contributing genes:",width = 3,
            uiOutput(ns("LeadingEdge")),solidHeader = T,status = "primary", height = 210),
        
        box(title = "Response signature genes:", width = 3,
            
            uiOutput(ns("sigOutput")), status = "primary", solidHeader = T, height = 480)
       
      ),
      
      
      fluidRow(
        box(title = "GO enrichment table:", DT::dataTableOutput(ns("enrichSig")),
            status = "primary", solidHeader = TRUE,width = 12, collapsible = TRUE)
      )
            ),
         
     )   
}
###########
sigNetworkModule = function(input,output,session, xRespDat, xRespInput = NULL, inputTab = NULL,#parent,
            message = "No GO enrichment, try relaxing the FDR or scorethreshold"){
  datSERV = reactive({xRespDat()})
  respSERV = reactive({xRespInput()})
  
  tabSERV = reactive({inputTab()})
  

#########################################
#############   NETWORK   ###############   
######################################### 

observeEvent(respSERV() ,{
    req(!is.null(respSERV()))
  
    drugMOD = input$respMOD
   
    sigRECVD =  respSERV()
    
    ns <- session$ns
    req(sigRECVD)
    req(length(sigRECVD > 0))
    
    
    ### if compound received from the server then update the compound in the module
    
    if(sigRECVD != drugMOD){
      
      updateSelectizeInput(session, 'respMOD', choices = sigRECVD)

      updateSliderInput(session,'scorethreshSIGNAT',label = "", min = 0, max = 5,
                        value = 1, step = 0.5)
      updateSliderInput(session,'fdrSIGNAT',label = "", min = 0, max = 0.5,
                        value = 0.2, step = 0.05)
      updatePrettyRadioButtons(session,"poolset", label = "",
                               choices = c("HIP" = "essential", "HOP" = "nonessential","HIPHOP" = "all"), 
                               selected = "all", inline = T)
     
    }
    
    
  },ignoreInit = F, ignoreNULL = T)


observeEvent(input$respMOD,{
    #req(datSERV())
    req(input$respMOD)
    #req(reactive({inputTab()}))
    
    
    ns <- session$ns 
    #req(tabSERV())
    updateSliderInput(session,'scorethreshSIGNAT',label = "", min = 0, max = 5,
                        value = 1, step = 0.5)
    updateSliderInput(session,'fdrSIGNAT',label = "", min = 0, max = 0.5,
                        value = 0.2, step = 0.05)
    updatePrettyRadioButtons(session,"poolset", label = "",
                               choices = c("HIP" = "essential", "HOP" = "nonessential","HIPHOP" = "all"), 
                               selected = "all", inline = T)
    
    #notSYNC = input$respMOD != respSERV()
    tabMOD = tabSERV() == "signature" 
    tabINP= tabSERV()
    ### if cmp was changed from the MOD side, then send the compound back
    
    
        
        returnRESP$sig = input$respMOD
        

  },ignoreInit = F, ignoreNULL = T)
#

#########################################
############# START NETWORK   ###############   
#########################################
  goEnrich <- reactive({ 
    req(datSERV())
    req(input$poolset)
    ns <- session$ns 
    x = xRespDat()
    
    we  <-  which(rownames(x)%in% ess$sgd_gene)
    wn  <-  which(rownames(x)%in% noness$sgd_gene)
    
    if(input$poolset == "nonessential") { xinp= x[wn,]}
    if  (input$poolset == "essential") {xinp = x[we,]} 
    
    if(input$poolset == "all") {xinp = x}
    
    w = which(colnames(xinp) %in% input$respMOD)
    
    validate(
      need(length(w) != 0, "Please select a compound")
              )
    
     thresh = input$scorethreshSIGNAT
    
    w2 = which(xinp[,w] >= thresh)
    
    validate(need(length(w2)!=0, message = "No scores above threshold"))

    req(length(w2)!=0)
    
    df = compSCORE(mat = xinp,coln = colnames(xinp)[w],sig = thresh)
   
    curr_exp = "network"
    
    FDR = input$fdrSIGNAT
    
    network = runGOENRICH(fdrThresh = FDR, curr_exp = "tst",score = df, 
                               bp_path = "2021_April5_GOBP_SGD.RDS",go_path = "2021_April5_GOID_GOBP_SGD.txt")
    
    validate(
      need(!is.null(network$enrichInfo), message = message)
    )
   
    enrichInfo = network$enrichInfo
    
    req(!(is.null(enrichInfo)))
    
    edgeMat = network$edgeMat
    
    return(network)
  })
#########################################
#############   NETWORK   ###############   
#########################################
net <- reactive({ 
  
    req(goEnrich()$enrichInfo)
    
    enrich = goEnrich()$enrichInfo
    
    edge = goEnrich()$edgeMat
    
    vis = visSetup(enrichInfo = enrich,edgeMat = edge, fontsize = 18, fontface = "Courier")
    
    vis
    
    })
 
   output$network_proxy <- renderVisNetwork({
    req(net()$nodes)
    vis = net()
    ns <- session$ns 
    n = net()$nodes
    
    req(n)
    
    w =  nrow(n)
    n <- n %>% arrange(term)
    
    names = n$id
    
    
    if(nrow(vis$edges)==0) {
      visNetwork(vis$nodes) %>% 
        visNodes(shadow=list(enabled=T,size=25),borderWidth=1) %>%
        visOptions(
          
          highlightNearest = list(enabled = T, degree = 5, hover = T),
          
          nodesIdSelection = list(enabled = TRUE, values = names,
            style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),
                                  
          
          selectedBy = list(variable="FDR",
            style = 'width: 200px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%
        
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    }
    
    
    else {
      ns <- session$ns
     
      wdth = as.numeric(input$shiny_wdth)/12*8
      visNetwork(vis$nodes, vis$edges) %>% 
        visNodes(shadow=list(enabled = T,size = 25)) %>%
        
        visOptions(
          
          highlightNearest = list(enabled = T, degree = 5, hover = T),
          
          
          nodesIdSelection = list(enabled = TRUE, values = names,
            style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),                         
          
          selectedBy = list(variable="FDR",
            style = 'width: 200px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%
        
        visIgraphLayout(layout = "layout_nicely", type = "full") %>%
        
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    }
  })
#########################################
#############   END NETWORK   ###########   
#########################################
   
 ########################################################################
 ################  START  GO NETWORK, PLOTS & DATATABLES  ###############
 ########################################################################
  output$Gotable = renderUI({
    
    req(input$network_proxy_selected)
    
    ns <- session$ns 
    
    DT::dataTableOutput(ns("detail"))
    
  })
 
#########################################
#############   NETWORK   ###############   
#########################################

  output$detail <- DT::renderDataTable({
    
    req(input$network_proxy_selected)
    
    vis = net()
    
    n = net()$nodes
    
    w = which(vis$nodes$id %in% c(input$network_proxy_selected))
    
    req(length(w) > 0)
    
    term = vis$nodes$label[w]
    
    nam = c("term","nGenes","geneSetFraction","FDR")
    
    m = match(nam,names(vis$nodes))
    
    n = vis$nodes[w,m]
    
    term = vis$nodes$label[w]
    
    nam = c("term","nGenes","geneSetFraction","FDR")
    
    m = match(nam,names(vis$nodes))
    
    names(vis$nodes)[m] = c("GO term","geneSet size","% of geneSet","FDR")
    
    n = vis$nodes[w,m]
    
    req(nrow(n)!=0)
    
    t = t(n[,2:4])
    
    datatable(t,width=220,selection = "single",caption = htmltools::tags$caption(term,
                                                            style = "caption-side: top; text-align: center; color:black;background:white;font-weight:bold;"),
              options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F,rowCallback = JS(
                "function(row, data) {",
                "for (i = 1; i < data.length; i++) {",
                "if (data[i]>1000 | data[i]<1){",
                "$('td:eq('+i+')', row).html(data[i].toExponential(1));",
                "}",
                "}",
                "}")),
              height = 400,colnames = "") %>%
      formatStyle( target = "row", color = "black",backgroundColor = "white",
                   columns = c(1,2),fontWeight = "bold")
  })  
##################################
##################################
##################################
#############   NETWORK   ###############   
######################################### 

   output$LeadingEdge = renderUI({
     req(input$network_proxy_selected)
     ns <- session$ns
     htmlOutput(ns("barPlot"))
   })

  
#########################################
#############   NETWORK   ###############   
#########################################
output$barPlot = renderText({

    req(input$network_proxy_selected)

    vis = net()

    #n = net()$nodes

    w = which(vis$nodes$id %in% c(input$network_proxy_selected))

    req(length(w) > 0)

    n = vis$nodes[w,]

    validate(need(nrow(n)!=0, message = "click node for detail"))

    leadInfo = genebarplot(n$overlapGenes)
    
    scoreRange = c(0,n$maxOverlapGeneScore)
    
    if(is.null(n$maxOverlapGeneScore)) scoreRange = c(0,max(leadInfo$score))
    
    hght = genOverlapGenePlot.gChart(oGeneStr  = n$overlapGenes,scoreRange = scoreRange)
    
    c('<img src="',hght[3],'">')
   
  }) 
#########################################
#############   NETWORK   ###############   
#########################################

###########################
##### START signature plot output 
###########################
   output$sigOutput = renderUI({
     
     ns <- session$ns 
     plotOutput(ns("barSigPlot"))
     
   })
   
#########################################
#############   NETWORK   ###############   
#########################################

   output$barSigPlot <- renderPlot({
     
     ns <- session$ns
     
     w = which(dfsig$signature %in% input$respMOD)
     
     validate(
       need(length(w)!= 0, "Please choose a signature")
     )
     
     d = dfsig[w,]
     d$score = as.numeric(d$score)
     d = d %>% arrange(desc(score))
     nrows = nrow(d)
     
     if(nrows > 10) d = d[1:10,]
     
     d = d %>% arrange(score)
     tit = stringWINDOW(paste(d$signature[1],d$site[1], sep=":"), width = 35)
     
     height = barplotHEIGHT(d)/96
     width = 150/96
     par(pin = c(width,height))
     
     barplot(d$score,names.arg = d$gene,las=1,horiz=T,col="navy",main = tit, xlab = "median fitness defect score")
     
   })
#########################################
#############   NETWORK   ###############   
#########################################
  hgtSigPlot = reactive({
    
    req(input$respMOD)

    w = which(dfsig$signature %in% input$respMOD)

    validate(
      need(length(w)!= 0, "Please choose a signature")
    )

    d = dfsig[w,]

    d = d %>% arrange(score)

    nrows = nrow(d)

    if(nrows > 10) d = d[1:10,]

    height = barplotHEIGHT(d)

    height = height*2

    })
###########################
##### END signature plot output 
########################### 

  observeEvent(input$network_proxy_selectedBy,
              
               
               {
                 req(input$network_proxy_selectedBy)
                 
                 n = net()$nodes
                 
                 w = which(n$FDR %in% as.numeric(input$network_proxy_selectedBy))
                 
                 id = n$id[w]
                 
                 ns <- session$ns 
                 visNetworkProxy(ns("network_proxy")) %>%
                   visSelectNodes(id = id)
                 
               }
  )
 ########################################################################
 ################  END  GO NETWORK, PLOTS & DATATABLES  #################
 ########################################################################
  
   
   
   
 #################################################
 #################################################
 ########## START GO ENCHICHMENT TABLE ###########
 #################################################
 ################################################# 
   
   enrid = reactive({
    
    req(net()$nodes)
    
    enrich = net()$nodes
    
    req(length(nrow(enrich)) > 0)
    
    row <- input$enrichSig_rows_selected
    
    out = outenrich()
    
    id = out$id[row]
    
    id
    
    
  })
  

############################  
outenrich = reactive({
    
    req(input$respMOD)
    req(goEnrich()$enrichInfo)
    
    enrich = goEnrich()$enrichInfo
    
    w = which(names(enrich) %in% c("querySetFraction", "geneSetFraction" ,
                                   "foldEnrichment" , "P" , "FDR" ))
    enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")] = 
      format(round(enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")],2),nsmall = 1,scientific = F)
    
    enrich[,c("P", "FDR")] = 
      format(signif(enrich[,c("P", "FDR")],2),nsmall = 1,scientific = T)
    
    enrich
  }
  )  
 ############################  
  observeEvent(enrid(), {
    req(enrid())
    id = enrid()
    
    ns <- session$ns
    visNetworkProxy(ns("network_proxy")) %>%
      visSelectNodes(id = id)
  })
############################   

  
############################ 
output$enrichSig = DT::renderDataTable({
    
    out = outenrich()
    w = which(names(out) %in% c("GOID","term","querySetFraction", "geneSetFraction" ,
                                "foldEnrichment" , "P" , "FDR","overlapGenes" ))
    
    out = out[,c("GOID","term","querySetFraction", "geneSetFraction" ,
                 "foldEnrichment" , "FDR","overlapGenes" )]
   
    
    opts = list(dom = 'Bfrtip',
    autoWidth = T,scrollX = T,columnDefs = list(
                  list(className = 'dt-left',targets = c(0,1,6)),
                  list(className = 'dt-center',targets = c(2:5)),
                  list(width = c('40px'),targets = c(2:4)),
                  list(width = c('300px'),targets = c(6)),
                  list(width = c('50px'),targets = c(0,5)),
                  list(width = c('400px'),targets = 1)
                    ))
                    
    df =  DT::datatable(out, options = opts,rownames = F, escape = F, selection = "single") %>% 
      formatStyle(c(1:7),fontWeight = 'bold', fontSize = '14px')
 })
      
 ############################ 
 ############################ END GO ENCHICHMENT TABLE and downloads
 ############################ 
 
output$downloadSignatures <- downloadHandler(
    filename = function() {
      paste0("enrich:",input$respMOD,"_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(dfsig, file, row.names = F,sep="\t")
    }
  )
 #################################################
 ########## END GO ENCHICHMENT TABLE #############
 #################################################
 
  
 ################################################
 ###### START UPDATES & ObserveEVENTS ###########
 ################################################
  
observeEvent(input$compResp,{
             
   req(input$compResp)
   tabMOD = tabSERV()
   if(tabMOD == "signature") {returnRESP$tab = "compsig"} 
   count = returnRESP$cnt
   count = count + 1
   
   returnRESP$cnt = count
   
    
  }, ignoreInit = F, ignoreNULL = T)    
############################ 

observeEvent(input$resetSig,{
    
    
    cmpSERV = colnames(datSERV())
    updateSelectizeInput(session, 'respMOD', choices = cmpSERV, selected = "iron homeostasis")
    
    updateSliderInput(session,'scorethreshSIGNAT',label = "", min = 0, max = 5,
                      value = 1, step = 0.5)
    updateSliderInput(session,'fdrSIGNAT',label = "", min = 0, max = 0.5,
                      value = 0.2, step = 0.05)
    
  },ignoreInit = F, ignoreNULL = T)
###############

  
  returnRESP = reactiveValues(
    sig = NULL,
    tab = NULL,
    cnt = 0
  )
   
 ################################################
 ###### END UPDATES & ObserveEVENTS #############
 ################################################ 
  
 

  return(
    
      list(
        sig = reactive({ returnRESP$sig }), 
        tab =  reactive({ returnRESP$tab }),
        cnt =  reactive({ returnRESP$cnt })
      )
      
    )
} 

            


            










 ############################ 
 ############################ 
 ############################ 
  
##################################
################################## 
##########   MAP   ##############   
#################################
#################################
#################################

  
  

  




