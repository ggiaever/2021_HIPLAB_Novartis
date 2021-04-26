###### 940 pixels by default 12 columns 


library(shiny)
library(dplyr)
library(DT)
library(shinydashboard)
library(visNetwork)
library(igraph)
library(shinyjs)


source(file.path("modules" ,"2021_April25_functions.R"))

source(file.path("modules/2021_April3_functions.R"))


pgo    <-  read.delim("2021_April5_GOID_GOBP_SGD.txt",stringsAsFactors = F,check.names = F)


geneSets    <-  readRDS("2021_April23_data/geneSets.RDS")
querySets    <-  readRDS("2021_April23_data/querySets.RDS")
matrixSets    <-  readRDS("2021_April23_data/matrixSets.RDS")
setnames = names(geneSets)

fdat    <-  read.delim("2021_April23_data/december4_fdat_updated.txt",stringsAsFactors = F,check.names = F)
bpref    <-  geneSets$`UBC:45`
dchallenge0  <-  querySets$`NIBR|42`

schallenge = split(dchallenge0$gene,dchallenge0$signature)
lnchallenge = sapply(schallenge,length)
lnbpref = sapply(bpref,length)
w = which(lnchallenge > 30)
lschallenge =  lapply(schallenge,function(x) paste(x,collapse = "|"))
swchallenge= lapply(schallenge[w], function(x) x = x[1:30])
lswchallenge=  lapply(swchallenge,function(x) paste(x,collapse = "|"))
lschallenge[w]=lswchallenge
lbpref =  lapply(bpref,function(x) paste(x,collapse = "|"))
w = which(lnbpref > 30)
swref= lapply(bpref[w], function(x) x = x[1:30])
lswref =  lapply(swref,function(x) paste(x,collapse = "|"))
lbpref[w]=lswref


CompareResponseModuleUI <- function(id, label = "CompareResponseModule") {

  ns <- NS(id)
  
 tagList(
    fluidPage(title = "Comparing Signatures",
              
          
      #########
      tabsetPanel(type = "tabs",
         id = ns("tabs"),
          tabPanel(
            "Comparing Signatures",
            
            
            fluidRow(
              column(
                width  = 3,
                
                box(
                  title  = "Compare response sets:",
                  fluidRow(
                    column(width = 2,
                           ),
                  
                    column(width = 10,
                      prettyRadioButtons(
                        ns("querySet"),
                        label = "Query signature set:",
                        choices = setnames,
                        outline = T,
                        fill = F,
                        status = "primary",
                        shape = "square",
                        bigger = T,
                        selected = setnames[2],
                        inline = T
                      ),
                      prettyRadioButtons(
                        ns("geneSet"),
                        label = "Reference signature set:",
                        choices = setnames,
                        outline = T,
                        fill = F,
                        status = "primary",
                        shape = "square",
                        bigger = T,
                        selected = setnames[1],
                        inline = T
                  ))),
                  solidHeader = T,
                  status = "primary",
                  width = "100%",
                  height = 200)#box
              ),
              
              
              
              column(
                width  = 6,
                box(
                  title = "Select individual query response signature:",
                 
                  selectizeInput(
                    ns("respMOD"), label = "",
                    choices = sort(unique(dchallenge0$signature)) ,
                    multiple = FALSE, width = "400px",
                    selected = "TOR signaling; rapamycin:18"
                  ),
                  status = "primary",
                  solidHeader = T,
                  width  =  "100%",
                  height = 200)#box
                
              ),#coln
              
              
              column(
                width = 3,
                box(
                  title = "All signatures tab:",
                  br(),
                  br(),
              fluidRow(
                
                
                column(
                  width = 12,align = "center",
                      
                      
                      actionButton(ns("allSigs"), label = "All signatures")
                      
                    )),
                    status = "primary", solidHeader = T, width = "100%", height = 200)#box
                )
                
               ),
              #############################
              #####################GOOD
              fluidRow(
                box(
                  title = "Query and reference response signatures (the first 20 genes are shown):",
                  htmlOutput(ns("alertsV"), placeholder = FALSE),
                  htmlOutput(ns("alertsIV"), placeholder = FALSE),
                  br(),
                  #uiOutput(ns("break1")),
                  DT::dataTableOutput(ns("silvertab")),
                  DT::dataTableOutput(ns("goldtab")),
                  
                  solidHeader = T, status = "primary", background = "navy", height = 300, width = 12)#box
                
                
              ),
              #####################GOOD
              fluidRow(
                box(
                  title = "Overlap analysis of the response signatures:",
                  
                  visNetworkOutput(ns("network_proxy"), width = "100%",height = 370),
                  
                  width = 6, status = "primary", solidHeader = TRUE, height = 500),#box
                
                
                box(
                  
                  title = "Comparison of query overlap with reference response signatures:",
                  width = 6,
                  br(),
                  br(),
                  
                  fluidRow(column(
                    width = 6,
                    offset = 1,
                    htmlOutput(ns("queryI"), placeholder = F)
                  )),
                  br(),
                  br(),
                  
                  fluidRow(
                    column(
                      width = 4,
                      align = "center",
                      
                      htmlOutput(ns("alertsI"), placeholder = FALSE),
                      
                      br(),
                      
                      uiOutput(ns("LeadingEdge"))
                    ),
                    column(
                      width = 3,
                      align = "center",
                      htmlOutput(ns("alertsII"), placeholder = TRUE),
                      
                      br(),
                      uiOutput(ns("sigOutput")),
                      
                    ),
                    
                    column(
                      width = 4,
                      align = "center",
                      
                      htmlOutput(ns("alertsIII"), placeholder = FALSE),
                      
                      br(),
                      
                      
                      uiOutput(ns("refOutput"))
                    )
                  ),
                  
                  
                  status = "primary", solidHeader = T, height = 500)#box
              ),
              
              #####################GOOD
              fluidRow(
                box(
                  title = "Query and reference (geneSet) intersections, foldEnrichment:",
                  DT::dataTableOutput(ns("enrichSig")),
                  
                  status = "primary", solidHeader = TRUE, width = 12, collapsible = TRUE)#box
              )
              
              ),
              
              #############################
              #############################
              #############################
              
              tabPanel(
                "GOenrichment",
                p("This page is for viewing GO enrichments of respoonse signatures"),
                
                fluidPage(
                  fluidRow(
                      box(
                        title = "Set fitness score threshold:",
                        sliderInput(
                          ns("scorethreshSIGNAT"),
                          label = "",
                          min = 0,
                          max = 5,
                          value = 1,
                          step = 0.5
                        ),
                        status = "primary",solidHeader = T, height = 200),#box
                      
                      box(
                        title = "Set FDR threshold:",
                        sliderInput(
                          ns("fdrSIGNAT"),
                          label = "",
                          min = 0,
                          max = 0.5,
                          value = 0.2,
                          step = 0.05
                        ),
                        status = "primary", solidHeader = T, height = 200)#box)
                    )
                 )
              ))))
      #########
      #########
      #########
      
   
}
###########
CompareResponseModule = function(input,output,session,inputTab = NULL,   
            message = "No GO enrichment, try relaxing the FDR or scorethreshold"){
  
  tabSERV = reactive({inputTab()})
  
  dfsig = reactiveValues(
    df = NULL
    
  )
    
 alerts = reactiveValues(
    go = NULL
    
  )
  
############################## FIRST UPDATE respMOD
# ### updates the respMOD depending on the choice of wuery and reference sets
# 
# ############################
#  
observeEvent(input$querySet,{
    
    req(input$respMOD)    
    req(input$querySet)  
  
   
if(input$querySet == setnames[1]) { 
    #queryset is input$respMOD
    dchallenge =  querySets[[1]]
      
      } else if(input$querySet == setnames[2]) { 
      
    dchallenge = querySets[[2]]
    
    
    } 
    
    ns <- session$ns
    
    choices = unique(sort(dchallenge$signature))
    
    dfsig$df = dchallenge
    
      updateSelectizeInput(session, 'respMOD', choices = choices, selected = choices[1])

      updateSliderInput(session,'scorethreshSIGNAT',label = "", min = 0, max = 5,
                        value = 3, step = 0.5)
      updateSliderInput(session,'fdrSIGNAT',label = "", min = 0, max = 0.5,
                        value = 0.2, step = 0.05)

  },ignoreInit = F, ignoreNULL = T)

############################ UPDATE FDR threshold and scorethresold whenever respMOD changes
### this output is hidden
### but is required to run the analysis for signature comparison
############################ 
observeEvent(input$querySet,{
req(input$respMOD)
req(input$querySet)
ns <- session$ns
hideTab(inputId = "tabs", target = "GOenrichment")  
})

############################   
observeEvent(input$respMOD,{
   
    req(input$respMOD)
    
    output$alertsIV =  renderText({paste("<h4><b>", "                                   ","</b></h4>")})
    ns <- session$ns 
   
    updateSliderInput(session,'scorethreshSIGNAT',label = "", min = 0, max = 5,
                        value = 1, step = 0.5)
    updateSliderInput(session,'fdrSIGNAT',label = "", min = 0, max = 0.5,
                        value = 0.2, step = 0.05)
    
    
    tabMOD = tabSERV() == "compsig" 
    tabINP= tabSERV()
   
        
        #returnRESP$sig = input$respMOD
        

  },ignoreInit = F, ignoreNULL = T)
 
 
 
 ########################################################################
 #####START NETWORK OF REFERENCE SIGNATURES THAT OVERLAP WITH QUERY #####
 ########################################################################
  
  goEnrich <- reactive({ 
    req(input$respMOD)
    
    ns <- session$ns 
   
    if(is.null(input$geneSet)) {geneSets = bpref } else if(input$geneSet == setnames[1]) { 
      geneSets = geneSets[[1]]} else if(input$geneSet == setnames[2]) { 
        geneSets = geneSets[[2]]} else if(input$geneSet == setnames[3]) { 
          geneSets = geneSets[[3]]}
    
   
  
    if(input$querySet == setnames[1]) { 
      
      xinp = matrixSets[[1]]
     
      
      } else if(input$querySet == setnames[2]) { 
      
      xinp = matrixSets[[2]]
      
    } 
    
    
    w = which(colnames(xinp) %in% input$respMOD)
    
    validate(
      need(length(w) != 0, "Please select a compound")
              )
    
    thresh = input$scorethreshSIGNAT
    
    w2 = which(xinp[,w] >= thresh)
    
    validate(need(length(w2)!=0, message = "No scores above threshold"))

    req(length(w2)!=0)
    
    curr_exp = "network"
    
    FDR = input$fdrSIGNAT
    
    network = runGORESP(fdrThresh = FDR, curr_exp = colnames(xinp)[w],coln = w,mat = xinp, bp_input = geneSets, minGeneSetSize = 2, sig = 1)
    
    df = compSCORE(mat = xinp,coln = colnames(xinp)[w],sig = thresh)
    
    
    gowork = runGORESP(fdrThresh = FDR, curr_exp = colnames(xinp)[w],coln = w, mat = xinp, bp_path = "2021_April5_GOBP_SGD.RDS",minGeneSetSize = 5, sig = 1, go_path = "2021_April5_GOID_GOBP_SGD.txt")
    
    
      output$alertsI = renderText({paste("<b>",  "  Query overlap with reference:","</b>")})
      output$alertsII = renderText({paste("<b>",  "  Query signature:","</b>")})
      output$alertsIII = renderText({paste("<b>",  "  Reference signature:","</b>")})
    
    if(is.null(network$enrichInfo) & !is.null(gowork$enrichInfo)) {
      
      network = gowork
      output$alertsI = renderText({"NO OVERLAP! GSEA used instead."})
      
    }  
      
      
      
      
   output$queryI = renderText({paste("<h4><b>",  "  QUERY = ",input$querySet,HTML('&emsp;'),"  REFERENCE = ",input$geneSet, "</b></h4>")})#"  SIGNIFICANT OVERLAP = ", overlapRESP$overlap,"</b>")})
   
    
    # validate(
    #   need(!is.null(network$enrichInfo), message = message)
    # )
   
    enrichInfo = network$enrichInfo
    
    if(is.null(network$enrichInfo)) {network$enrichInfo = 0
                            
                             return(network)}
    
    edgeMat = network$edgeMat
    return(network)
   
  })
 
######################################################################

net <- reactive({ 
  
    req(goEnrich()$enrichInfo!= 0)
    
    enrich = goEnrich()$enrichInfo
    
    edge = goEnrich()$edgeMat
    
    vis = visSetup(enrichInfo = enrich,edgeMat = edge, fontsize = 16, fontface = "Courier")
    
    vis
    
    })
############################ 
  
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
            style = 'width: 500px;  height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),               
                                  
          
          selectedBy = list(variable="FDR",
             style = 'width: 200px;  height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'))%>% 
        
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    }
    
    else {visNetwork(vis$nodes, vis$edges, height = 400) %>% 
        visNodes(shadow=list(enabled = T,size = 25)) %>%
    
        visOptions(
          
          highlightNearest = list(enabled = T, degree = 5, hover = T),
          
          
          nodesIdSelection = list(enabled = TRUE, values = names,
              style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),
          
          selectedBy = list(variable="FDR",
              style = 'width: 200px;  height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'))%>%
                            
        
        visIgraphLayout(type = "full",layout = "layout_nicely") %>%
        
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    }
  })

 ########################################################################
 #####END NETWORK OF REFERENCE SIGNATURES THAT OVERLAP WITH QUERY #######
 ########################################################################



############################
############################
 ##############################################################
 ##############  START SELECT NODE BY FDR #####################
 ##############################################################
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
 ##############################################################
 ##############  END SELECT NODE BY FDR #######################
 ##############################################################


 ##############################################################
 ##############  START OVERLAP SIGNATURE BARPLOT ##############
 ##############################################################
 output$LeadingEdge = renderUI({
    req(input$network_proxy_selected)
    ns <- session$ns 
    
    htmlOutput(ns("barPlot"))
  
    })

############################ 
############################ 
 output$barPlot = renderText({

    req(input$network_proxy_selected)

    vis = net()

    n = net()$nodes

    w = which(vis$nodes$id %in% c(input$network_proxy_selected))

    req(length(w) > 0)

    n = vis$nodes[w,]

    validate(need(nrow(n)!=0, message = "click node for detail"))

    leadInfo = geneBARPLOT(n$overlapGenes)
    
    scoreRange = c(0,n$maxOverlapGeneScore)
   
    if(is.null(n$maxOverlapGeneScore)) scoreRange = c(0,max(leadInfo$score))
    
    hght = genOverlapGenePlot.gChart(oGeneStr  = n$overlapGenes,scoreRange = scoreRange)
    
    hght[1:2]=paste0(hght[1:2],"px")
    h3 = paste0("<img src=",hght[3],">")
    
     h3
  })
 ##############################################################
 ##############  END OVERLAP SIGNATURE BARPLOT ################
 ##############################################################



 ##############################################################
 ##############  START REFERENCE SIGNATURE BARPLOT ############
 ##############################################################
output$refOutput = renderUI({
    
    
      ns <- session$ns 
      htmlOutput(ns("barRefPlot"))#, height = hgtSigPlot(),width = 150)
      #plotOutput(ns("barSigPlot"), height = hgtSigPlot())

  })

###########################
###########################

 


  output$barRefPlot  = renderText({
    
    
    q = querySets[[input$geneSet]]
  
    
    req(input$network_proxy_selected)

    vis = net()

    n = net()$nodes

    w = which(vis$nodes$id %in% c(input$network_proxy_selected))

    req(length(w) > 0)

    n = vis$nodes[w,]

    wq = which(q$signature %in% n$term)
    
    validate(
      need(length(wq)!= 0, "Not applicable.")
    )

    d = q[wq,]

    d = d %>% arrange(desc(score))

    nrows = nrow(d)

    if(nrows > 15) d = d[1:15,]
    
    overlapSet <- paste(d$gene, "(", d$score, ")", sep = "")
    overlapSet <- paste(overlapSet, collapse = "|")
    
    hght = genOverlapGenePlot.gChart(oGeneStr  = overlapSet,scoreRange = c(0,max(d$score)),maxGenes = 15)
    
    c('<img src="',hght[3],'">')
    
    

    })
 ##############################################################
 ##############  END REFERENCE SIGNATURE BARPLOT ##############
 ##############################################################

 ##############################################################
 ##############  START QUERY SIGNATURE BARPLOT ################
 ##############################################################
output$sigOutput = renderUI({
    
    
      ns <- session$ns 
      htmlOutput(ns("barSigPlot"))
      

  })

###########################
###########################

 
###########################
###########################
  output$barSigPlot  = renderText({
    
    if(!is.null(dfsig$df)) {dchallenge = dfsig$df} else { dchallenge = dchallenge0}

    w = which(dchallenge$signature %in% input$respMOD)

    validate(
      need(length(w)!= 0, "Please choose a signature")
    )

    d = dchallenge[w,]

    d = d %>% arrange(desc(score))

    nrows = nrow(d)

    if(nrows > 15) d = d[1:15,]
    
    
    overlapSet <- paste(d$gene, "(", d$score, ")", sep = "")
    overlapSet <- paste(overlapSet, collapse = "|")
    
    hght = genOverlapGenePlot.gChart(oGeneStr  = overlapSet,scoreRange = c(0,max(d$score)),maxGenes = 15)
    
    c('<img src="',hght[3],'">')
    
    })
 ##############################################################
 ##############  END QUERY SIGNATURE BARPLOT ##################
 ##############################################################

#############################################################
####  START GOLDTAB: TABLE OF OVERLAP BETWEEN SIGNATURES ####
#############################################################

output$goldtab = DT::renderDataTable({
    req(input$network_proxy_selected)
    req(input$geneSet)
    req(input$querySet)
    if(!is.null(dfsig$df)) {dchallenge = dfsig$df} else { dchallenge = dchallenge0}
   
    geneSet = geneSets[[input$geneSet]]
    querySets = querySets[[input$querySet]]
  
    schallenge = split(querySets$gene,querySets$signature)
    lschallenge2 = collapseSTR(schallenge, limit = max(sapply(schallenge,length)))
    lschallenge = collapseSTR(schallenge,limit = 15)
    lgeneSets2 = collapseSTR(geneSet, limit = max(sapply(geneSet,length)))  
    lgeneSets = collapseSTR(geneSet, limit = 15)
    #################
    #################
    vis = net()
    
    n = net()$nodes
    
    
    wnode = which(vis$nodes$id %in% c(isolate(input$network_proxy_selected)))
    
    
    req(length(wnode) > 0)
    
    term = vis$nodes$term[wnode]
    
    winput = which(names(lschallenge) %in% isolate(input$respMOD))
    
    wterm = which(names(lgeneSets) %in% term)
   
    
    if(length(wterm) >  0) {
      ngeneSets = lgeneSets[wterm]
      
      names(ngeneSets) = names(lgeneSets)[wterm]
      ngeneSets
      
      
    } else {
      
      leadInfo = geneBARPLOT(n$overlapGenes[wnode])
      
      overlapSet <- paste(leadInfo$gene, collapse = "|")
      
      ngeneSets = as.list(overlapSet)
      names(ngeneSets) = term[wnode]
      ngeneSets
      
    }
    
    goterm = vis$nodes$term[wnode]
    over = vis$nodes$overlapGenes[wnode]
    if(length(wterm) == 0){
    lenov = length(revSPLIT(over)[[1]])
    wover = which(pgo$term %in% goterm)
    gene = pgo$gene[wover]
    
  
    leng = length(revSPLIT(gene)[[1]])
    query = revSPLIT(lschallenge2)
    lnchallenge = sapply(query,length)
    lenquery = lnchallenge[winput]
    wgin = min(lenquery,leng)
    overGOCOEF = round(lenov/wgin*100,2)
    
 
    wover = which(pgo$term %in% goterm)
    
    if(pgo$geneSetSize[wover] > 20){
     ugene = revSPLIT(gene)
     gene1 = collapseSTR(ugene,limit = 15)
      } else {
      ugene = revSPLIT(gene) 
      gene1 = collapseSTR(ugene) 
     }
    
    oGenes <- unlist(strsplit(over, "\\|"))
    oGenes <- strSPLIT(oGenes,"\\(",1)
    over = collapseSTR(list(oGenes))
    
    } else { overGOCOEF = 0 }
    

    
    if(length(wterm) > 0){
    
    query = revSPLIT(lschallenge2)
    
    gterm = revSPLIT(lgeneSets2)
    lnchallenge = sapply(query,length)
    lnterm = sapply(gterm,length)
    
    lenquery = lnchallenge[winput]
    lengene = lnterm[wterm]
    lgterm = gterm[[wterm]]
    
    lsinput = query[[winput]]
    
   
    int = intersect(lsinput,lgterm)
   
    lenint = length(int)
    if(lenint > 0) lint = paste(sort(int),collapse = "|")
    if(lenint >= 20) lint = paste(int[1:20],collapse = "|")
    
   
    wmin = min(lenquery,lengene)
    overlapCOEF = round(lenint/wmin*100,2)
    
    } else {
      
    overlapCOEF = 0
    }
    
  observeEvent(input$network_proxy_selected,{
    
    
    if(overlapCOEF!= 0 & input$network_proxy_selected != "") {
      show("alertsIV")
      output$alertsIV = renderText({paste("<h4><b>",  "   overlapCOEF = ",overlapCOEF,"</b></h4>")})
      
    
      
    } else if(overlapCOEF== 0 & input$network_proxy_selected != ""& length(wterm) == 0) {
      show("alertsIV")
      output$alertsIV = renderText({paste("<h4><b>",  "   overlapCOEF = ",overGOCOEF," (No overlap! GO enrichment analysis (GSEA) used instead)","</b></h4>")})
      
      } else {hide("alertsIV")}
  }, ignoreInit = F, ignoreNULL = T)
  
  
  
 
   
     if(length(wterm) == 0 ){
    
        df = data.frame(SET = c("querySet","goSet","overLap"),RESPONSE= c(names(lschallenge[winput]),goterm, ""), 
                       SIGNATURE = c(unlist(lschallenge[winput],
                      use.names = F),unlist(gene1,use.names=F), unlist(over,use.names=F)),stringsAsFactors=F)
      
    } else if ( overlapCOEF != 0) {
      
      
      df = data.frame(SET = c("querySet","geneSet","overlap"),
                    RESPONSE = c(names(lschallenge[winput]), names(ngeneSets)," "),
                    SIGNATURE = c(unlist(lschallenge[winput],use.names = F),
                              unlist(ngeneSets,use.names = F),lint),
                      stringsAsFactors=F   )
    }
    
    
    
    opts = list(paging = F, searching = F, info = F,
      dom = 'Bfrtip',
                autowidth = T, scrollX = T, columnDefs = list(
                  list(className = 'dt-left',targets = c(0:2)),
                  list(width = c('500px'),targets = c(2)),
                  list(width = c('40px'),targets = c(0)),
                  list(width = c('400px'),targets  = c(1))))
    
    df =  DT::datatable(df, options = opts,rownames = F, selection = "single", escape = F) %>% 
      formatStyle(c(1:3),fontWeight = 'bold', fontSize = '14px')
    
    
})
 
 #######################################################
 #######################################################
 #### TABLE WHEN NO OVERLAP OR GO ENRICHMENT EXISTS ####
 #######################################################
 #######################################################
 

observeEvent(goEnrich(),{
 req(input$querySet)
 req(input$respMOD)
 
 enrich = isolate(goEnrich()$enrichInfo)
 
 
 if(length(isolate(goEnrich()$enrichInfo)) == 1) {
      
      output$alertsV = renderText({paste("<h4><b>",  "   overlapCOEF = ","No overlap OR GO enrichment!","</b></h4>")})
}
 
 if(length(isolate(goEnrich()$enrichInfo)) == 1) { 
   
       show("silvertab")} else {hide("silvertab")}
 
 if(length(isolate(goEnrich()$enrichInfo)) == 1){
       show("alertsV")} else {hide("alertsV")}
 
 if(length(isolate(goEnrich()$enrichInfo)) == 1) {    
        req(length(enrich) == 1)
        querySets = querySets[[input$querySet]]
       
        schallenge = split(querySets$gene,querySets$signature)
        
        
        lschallenge = collapseSTR(schallenge,limit = 15)
        winput = which(names(lschallenge) %in% input$respMOD)
        req(length(winput) > 0)    
        req(length(enrich) == 1)
        
        output$silvertab = DT::renderDataTable({
        df = data.frame(SET = "querySet",RESPONSE = names(lschallenge)[winput],
                       SIGNATURE = unlist(lschallenge[winput],use.names = F),stringsAsFactors=F   )
      opts = list(paging = F, searching = F, info = F, 
            #scrollY = TRUE,#dom = 'Bfrtip',
                      autowidth = T, scrollX = T, columnDefs = list(
                        list(className = 'dt-left',targets = c(0:2)),
                        list(width = c('10px'),targets = c(2)),
                        list(width = c('10px'),targets = c(0)),
                        list(width = c('700px'),targets  = c(1))))

          df =  DT::datatable(df, options = opts,rownames = F, selection = "single", escape = F) %>%
            formatStyle(c(1:3),fontWeight = 'bold', fontSize = '14px')
      })}
      
        
       }, ignoreInit = T, ignoreNULL = F) 
 
 
############################################
##############  END GOLDTAB ################
############################################



 
 
 #################################################
 #################################################
 ########## START GO ENCHICHMENT TABLE ###########
 #################################################
 ################################################# 
  enrID = reactive({
    
    req(net()$nodes)
    
    enrich = net()$nodes
    
    req(length(nrow(enrich)) > 0)
    
    row <- input$enrichSig_rows_selected
    
    out = outenrich()
    
    id = out$id[row]
    
    id
    
    
  })
  
############################ get entire enrichment table 
############################  
outenrich = reactive({
    
    req(input$respMOD)
    
    req(goEnrich()$enrichInfo!=0)
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
 #################################################
 #################################################
 #################################################
  observeEvent(enrID(), {
    req(enrID())
    id = enrID()
    
    ns <- session$ns
    visNetworkProxy(ns("network_proxy")) %>%
      visSelectNodes(id = id)
  })
############################   
############################
############################

############################ 
output$enrichSig = DT::renderDataTable({
    
    out = outenrich()
    
    
    #nam = c("term","interSect","filename","querySetFraction","nQuery","geneSetFraction","nGenes" ,"FDR", "overlapGenes"  )
    nam = c("term","interSect","filename","nQuery","nGenes" ,"FDR", "overlapGenes"  )
    
    w = which(names(out) %in% nam)
    
    out = out[,nam]
    
    opts = list(pageLength = 25,paging = F,searching = F,info = F,
    autoWidth = T,scrollX = T,columnDefs = list(
                  list(className = 'dt-left',targets = c(0,2,5,6)),
                  list(className = 'dt-center',targets = c(1,3:5)),
                  list(width = c('300px'),targets = c(0,2)),
                  
                  list(width = c('100px'),targets = c(5)),
                  list(width = c('20px'),targets = c(1,3:5)),
                  list(width = c('100px'),targets = c(6))
                
                  ))
    df =  DT::datatable(out, options = opts,rownames = F, selection = "single", escape = F) %>% 
      formatStyle(c(1:6),fontWeight = 'bold', fontSize = '14px')
    
})
 #################################################     
 #################################################
 ########## END GO ENCHICHMENT TABLE #############
 #################################################
 ################################################# 

    
###############
observeEvent({input$allSigs
              },{
   req(input$allSigs)
   tabMOD = tabSERV()
   if(tabMOD == "compsig") {returnRESP$tab = "signature"} 
   count = returnRESP$cnt
   count = count + 1
   
   returnRESP$cnt = count
   
    
  }, ignoreInit = F, ignoreNULL = T)  
############### 
############### 

############### 
############### 
############### 
############### 
  
  returnRESP = reactiveValues(
   #sig = NULL,
   tab = NULL,
   cnt = 0)
# ############################ 
  return(

      list(
       # sig = reactive({ returnRESP$sig }), 
        tab =  reactive({ returnRESP$tab }),
        cnt =  reactive({ returnRESP$cnt })
      )
      
    )
  
}
  
  
  



