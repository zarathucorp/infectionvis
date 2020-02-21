
#' @title epiODE equation
#' @description epiODE equation
#' @param time time
#' @param state initial states
#' @param modelparams parameters
#' @return epiODE output
#' @details epiODE
#' @examples 
#' init    <- c(S = 49656776, E = 0, I = 1, A = 0, Q = 0, R = 0, D = 0, Report = 1)
#' params  <- c(sigma1 = 1, sigma2 = 0, sigma3 = 0, sigma4 = 1, l_E = 0.1, 
#'              l_A = 0.5, l_Q = 0.1, gamma = 1/6, f = 0, p = 2/3, k = 1.9, alpha_Q = 1/6,
#'              beta = 0.33)
#' times   <- seq(0, 141, by = 1)
#' deSolve::ode(y = init, times = times, func = epiODE, parms = params)
#' @rdname epiODE
#' @export 

epiODE <- function(time, state, modelparams) {
  with(as.list(c(state, modelparams)), {
    N = S + E + I + A + Q + R + D
    lambda = beta*(I + l_A*A + l_E*E + l_Q*Q)/N
    dS <- -sigma1*lambda*S 
    dE <- sigma1*lambda*S - k*E
    dI <- sigma2*(1-p)*lambda*S + sigma1*(1-p)*k*E - f*I - sigma3*alpha_Q*I - sigma4*gamma*I
    dA <- sigma2*p*lambda*S + sigma1*p*k*E - gamma*A
    dQ <- sigma3*alpha_Q*I - gamma*Q - f*Q
    dR <- gamma*A + sigma4*gamma*I + gamma*Q
    dD <- f*(I + Q)
    dReport <- sigma2*(1-p)*lambda*S + sigma1*(1-p)*k*E
    return(list(c(dS, dE, dI, dA, dQ, dR, dD, dReport)))
  })
}




#' @title 'Shinydashboard' for Infection Simulation
#' @description RStudio addin- 'Shinydashboard' for Infection Simulation 
#' @param max.filesize Maximum file size to upload (MB), Default: 2048 (2 GB)
#' @return RStudio addin- 'Shinydashboard' for Infection Simulation
#' @details RStudio addin- 'Shinydashboard' for Infection Simulation
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[deSolve]{ode}}
#' @rdname infectionvis
#' @export 
#' @importFrom deSolve ode
#' @importFrom shinycustomloader withLoader
#' @importFrom shinydashboard box
#' @importFrom data.table fread
#' @importFrom pracma lsqnonlin
#' @import shiny 
#' @import shinydashboard
#' @import highcharter

infectionvis <- function(max.filesize = 2048){
  options(shiny.maxRequestSize = max.filesize * 1024^2)

  
  ui <- dashboardPage(
    dashboardHeader(
      title = "Infection Dashboard",
      titleWidth = 500,
      dropdownMenu(type = "messages",
                   messageItem(
                     from = "Sales Dept",
                     message = "Sales are steady this month."
                   ),
                   messageItem(
                     from = "New User",
                     message = "How do I register?",
                     icon = icon("question"),
                     time = "13:45"
                   ),
                   messageItem(
                     from = "Support",
                     message = "The new server is ready.",
                     icon = icon("life-ring"),
                     time = "2014-12-01"
                   )
      )),
    
    dashboardSidebar(
      width = 500,
      withMathJax(),
      sidebarMenu(
        fluidRow(
          column(8, fileInput("file", "Upload zip file")),
          column(4, helpText(a(h4("Example data", class = "btn btn-default action-button" , 
                                  style = "fontweight:600"), href="https://zarathu.s3.ap-northeast-2.amazonaws.com/exampledata/fludata09.csv")))
        ),
        uiOutput("inicond"),
        menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
        menuItem("Other", icon = icon("th"), tabName = "other",
                 badgeLabel = "new", badgeColor = "green")
      )
    ),
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(tabName = "dashboard",
                fluidRow(
                  withMathJax(),
                  valueBoxOutput("estbeta"),
                  valueBoxOutput("R0")
                  
                ),
                fluidRow(
                  box(title = "Estimation", status = "primary", solidHeader = TRUE, collapsible = T,
                      withLoader(highchartOutput(paste0("chart", "Estimation")), type="html", loader="loader6")),
                  uiOutput("chart_parametters")
                )
                
        )
      )
    )
  )
  
  
  
  server <- function(input, output, session){
    

    
    userFile <- eventReactive(input$file, {
      # If no file is selected, don't do anything
      #validate(need(input$file, message = FALSE))
      input$file
    })
    
    incidata <- eventReactive(input$file, {
      validate(need((grepl("csv", userFile()$name) == T), message = "Please upload zip file"))
      out <- fread(userFile()$datapath)
      return(out)
    })
    
    
    
    output$inicond <- renderUI({
      tagList(
        h2("Choice model"),
        withMathJax(),
        numericInput("ini_S", '감수성 (S) $$S_{0}$$', min = 0,  value = 49656776),
        radioButtons("E_yn", "잠복기 (E)", c("Yes" = T, "No" = F), selected = T, inline = T),
        conditionalPanel(condition = "input.E_yn == 'TRUE'",
                         fluidRow(
                           column(4, numericInput("ini_E", '$$E_{0}$$', min = 0, value = 0)),
                           column(4, numericInput("param_kappa", "$$\\kappa$$", value = 1.9, min = 0)),
                           column(4, sliderInput("param_le", "$$l_{E}$$", min =0, max = 1, value = 0.1))
                         )),
        
        numericInput("ini_I", '증상 감염 (I) $$I_{0}$$', min = 1, value = 1),
        
        radioButtons("A_yn", "무증상감염 (A)", c("Yes" = T, "No" = F),  selected = T, inline = T),
        conditionalPanel(condition = "input.A_yn == 'TRUE'",
                         fluidRow(
                           column(4, numericInput("ini_A", '$$A_{0}$$', min = 0, value = 0)),
                           column(4, numericInput("param_p", "$$p$$", value = 2/3, min = 0)),
                           column(4, sliderInput("param_la", "$$l_{A}$$", min =0, max = 1, value = 0.5))
                         )),
        
        radioButtons("Q_yn", "격리환자 (Q)", c("Yes" = T, "No" = F), selected = F, inline = T),
        conditionalPanel(condition = "input.Q_yn == 'TRUE'",
                         fluidRow(
                           column(4, numericInput("ini_Q", '$$Q_{0}$$', min = 0, value = 1)),
                           column(4, numericInput("param_alpha_q", "$$\\alpha_{q}$$", value = 1/6, min = 0)),
                           column(4, sliderInput("param_lq", "$$l_{q}$$", min =0, max = 1, value = 0.1)))
        ),
        
        
        fluidRow(
          column(6, numericInput("ini_R", '회복 (R) $$R_{0}$$', min = 0, value = 0)),
          column(6, numericInput("param_gamma", '회복률 $$\\gamma$$', value = 1/6, min = 0))
        ),
        
        
        radioButtons("D_yn", "사망 (D)", c("Yes" = T, "No" = F), selected = T, inline = T),
        conditionalPanel(condition = "input.D_yn == 'TRUE'",
                         fluidRow(
                           column(6, numericInput("ini_D", '$$D_{0}$$', min = 0, value = 1)),
                           column(6, numericInput("param_f", "$$f$$", value = 0, min = 0)))),
        
        fluidRow(
          column(6, numericInput("ini_Report", '$$Report_{0}$$', min = 0,  value = 1)),
          column(6, h4("Additional parameters"),
                 
                 shinyWidgets::dropdownButton(
                   tags$div(
                     style = "color: black !important;", # for text
                     numericInput("param_beta0", "$$\\beta_{ini}$$", value = 0.1, min = 0)
                   ),
                   circle = FALSE,
                   icon = icon("cog")
                 )
          )
        ))
    })
    
    
    
    
    
    out.ode <- reactive({
      req(incidata())
      params  <- c(sigma1 = ifelse(input$E_yn == T, 1, 0),
                   sigma2 = ifelse(input$E_yn == T, 0, 1),
                   sigma3 = ifelse(input$Q_yn == T, 1, 0),
                   sigma4 = ifelse(input$Q_yn == T, 0, 1),
                   l_E = ifelse(is.null(input$param_le), 0, input$param_le),
                   l_A = ifelse(is.null(input$param_la), 0, input$param_la),
                   l_Q = ifelse(is.null(input$param_lq), 0, input$param_lq),
                   p = ifelse(is.null(input$param_p), 0, input$param_p),
                   k = ifelse(is.null(input$param_kappa), 0, input$param_kappa),
                   f = ifelse(is.null(input$param_f), 0, input$param_f),
                   alpha_Q = ifelse(is.null(input$param_alpha_q), 0, input$param_alpha_q),
                   gamma = ifelse(is.null(input$param_gamma), 0, input$param_gamma)
      )
      
      
      init <-  c(S = input$ini_S, 
                 E = ifelse(input$E_yn == F, 0, input$ini_E), 
                 I = input$ini_I, 
                 A = ifelse(input$A_yn == F, 0, input$ini_A), 
                 Q = ifelse(input$Q_yn == F, 0, input$ini_Q),
                 R = input$ini_R,
                 D = ifelse(input$D_yn == F, 0, input$ini_D),
                 Report = input$ini_Report
      )
      
      times   <- seq(0, length(incidata()$V1) - 1, by = 1)
      
      ftemp <- function(x) {
        modelparams <- c(params, beta = x)
        tempoutput <- as.data.frame(ode(y = init, times = times, func = epiODE, parms = modelparams))
        return(abs(tempoutput$Report - incidata()$V1))
      }
      
      beta0 <- ifelse(is.null(input$param_beta0), 0.1, input$param_beta0)
      fitresult <- pracma::lsqnonlin(ftemp , beta0)
      beta <- fitresult$x   
      
      ## solve the model
      params <- c(params, beta = beta)
      output <- ode(y = init, times = times, func = epiODE, parms = params)
      output <- as.data.frame(output)
      
      R0 <- (beta*params["gamma"]^2*params["k"]*params["sigma2"] + beta*params["gamma"]^2*params["k"]*params["sigma1"]^2 - beta*params["gamma"]^2*params["k"]*params["p"]*params["sigma1"]^2 + beta*params["f"]*params["gamma"]*params["k"]*params["sigma2"] + beta*params["f"]*params["gamma"]*params["k"]*params["sigma1"]^2 + beta*params["f"]*params["gamma"]^2*params["l_E"]*params["sigma1"] + beta*params["f"]^2*params["gamma"]*params["l_E"]*params["sigma1"] - beta*params["gamma"]^2*params["k"]*params["p"]*params["sigma2"] + beta*params["gamma"]^3*params["l_E"]*params["sigma1"]*params["sigma4"] - beta*params["f"]*params["gamma"]*params["k"]*params["p"]*params["sigma1"]^2 + params["alpha_Q"]*beta*params["gamma"]^2*params["l_E"]*params["sigma1"]*params["sigma3"] + beta*params["f"]*params["gamma"]^2*params["l_E"]*params["sigma1"]*params["sigma4"] + beta*params["f"]^2*params["k"]*params["l_A"]*params["p"]*params["sigma2"] + beta*params["f"]^2*params["k"]*params["l_A"]*params["p"]*params["sigma1"]^2 - beta*params["f"]*params["gamma"]*params["k"]*params["p"]*params["sigma2"] + beta*params["gamma"]^2*params["k"]*params["l_A"]*params["p"]*params["sigma1"]^2*params["sigma4"] + params["alpha_Q"]*beta*params["f"]*params["gamma"]*params["l_E"]*params["sigma1"]*params["sigma3"] + params["alpha_Q"]*beta*params["gamma"]*params["k"]*params["l_Q"]*params["sigma2"]*params["sigma3"] + beta*params["f"]*params["gamma"]*params["k"]*params["l_A"]*params["p"]*params["sigma2"] + params["alpha_Q"]*beta*params["gamma"]*params["k"]*params["l_Q"]*params["sigma1"]^2*params["sigma3"] + beta*params["f"]*params["gamma"]*params["k"]*params["l_A"]*params["p"]*params["sigma1"]^2 + beta*params["gamma"]^2*params["k"]*params["l_A"]*params["p"]*params["sigma2"]*params["sigma4"] + params["alpha_Q"]*beta*params["f"]*params["k"]*params["l_A"]*params["p"]*params["sigma1"]^2*params["sigma3"] + params["alpha_Q"]*beta*params["gamma"]*params["k"]*params["l_A"]*params["p"]*params["sigma1"]^2*params["sigma3"] - params["alpha_Q"]*beta*params["gamma"]*params["k"]*params["l_Q"]*params["p"]*params["sigma1"]^2*params["sigma3"] + beta*params["f"]*params["gamma"]*params["k"]*params["l_A"]*params["p"]*params["sigma1"]^2*params["sigma4"] + params["alpha_Q"]*beta*params["f"]*params["k"]*params["l_A"]*params["p"]*params["sigma2"]*params["sigma3"] + params["alpha_Q"]*beta*params["gamma"]*params["k"]*params["l_A"]*params["p"]*params["sigma2"]*params["sigma3"] - params["alpha_Q"]*beta*params["gamma"]*params["k"]*params["l_Q"]*params["p"]*params["sigma2"]*params["sigma3"] + beta*params["f"]*params["gamma"]*params["k"]*params["l_A"]*params["p"]*params["sigma2"]*params["sigma4"])/(params["gamma"]^3*params["k"]*params["sigma4"] + params["f"]*params["gamma"]^2*params["k"] + params["f"]^2*params["gamma"]*params["k"] + params["alpha_Q"]*params["gamma"]^2*params["k"]*params["sigma3"] + params["f"]*params["gamma"]^2*params["k"]*params["sigma4"] + params["alpha_Q"]*params["f"]*params["gamma"]*params["k"]*params["sigma3"])
      
      
      return(list(beta = beta, output = output, R0 = R0, params = params, init=init))
    })
    
    output$estbeta <- renderValueBox({
      req(input$param_beta0)
      valueBox(
        out.ode()$beta, "Estimated beta", icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow"
      )
    })
    
    output$R0 <- renderValueBox({
      req(out.ode())
      valueBox(
        out.ode()$R0, "Estimated R0", icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green"
      )
    })
    
    output$chartEstimation <- renderHighchart({
      out.ode <- out.ode()$output
      highchart() %>% 
        hc_xAxis(time = out.ode$time) %>% 
        hc_add_series(name = "Data", data = incidata()$V1, type = "scatter") %>% 
        hc_add_series(name = "Cumulative incidence", data = out.ode$Report) %>% 
        hc_exporting(enabled = T) %>% 
        hc_tooltip(valueDecimals = 0)
      
    })
    
    
    output$chart_parametters <- renderUI({
      
      paramList <- names(which(c("S" = T, "E" = input$E_yn, "I" = T, "A"= input$A_yn, "Q" = input$Q_yn, "R" = T, "D" = input$D_yn) == T))
      
      do.call(tagList, lapply(paramList, function(i){
        box(title = i, status = "primary", solidHeader = TRUE,
            withLoader(highchartOutput(paste0("chart", i)), type="html", loader="loader6"))
      }))
    })
    
    
    observeEvent(c(out.ode()), {
      paramList <- names(which(c("S" = T, "E" = input$E_yn, "I" = T, "A"= input$A_yn, "Q" = input$Q_yn, "R" = T, "D" = input$D_yn) == T))
      
      for (i in paramList) {
        local({
          my_i <- i
          output[[paste0("chart",my_i)]] <- renderHighchart({
            out.ode <- out.ode()$output
            highchart() %>% 
              hc_xAxis(time = out.ode$time) %>% 
              hc_add_series(name = my_i, data = out.ode[[my_i]]) %>% 
              hc_exporting(enabled = T) %>% 
              hc_tooltip(valueDecimals = 0)
          })
          
        })
      }
    })
    
    
    
  }
  
  
  viewer <- browserViewer(browser = getOption("browser"))
  runGadget(ui, server, viewer = viewer)
}

