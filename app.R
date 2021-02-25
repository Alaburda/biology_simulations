library(shiny)
library(markdown)
library(deSolve)
library(plotly)
library(tidyverse)

time <- seq(0,50,by = 0.01)

competitive_inhibition <- function(time, substrate, params) {
    with(as.list(c(substrate, params)), {
        dE = k_r*ES - k_f*S*E + k_cat*ES - k_if*I*E + k_ir*EI
        dS = k_r*ES - k_f*S*E
        dI =  k_ir*EI - k_if*I*E
        dES = k_f*S*E - k_r*ES - k_cat*ES
        dP = k_cat*ES
        dEI = k_if*I*E - k_ir*EI
        
        return(list(c(dE, dS, dI, dES, dP, dEI)))
    })
}

uncompetitive_inhibition <- function(time, substrate, params) {
    with(as.list(c(substrate, params)), {
        dE = k_r*ES - k_f*S*E + k_cat*ES
        dS = k_r*ES - k_f*S*E
        dI =  k_ir*EIS - k_if*I*ES
        dES = k_f*S*E - k_r*ES - k_cat*ES - k_if*I*ES + k_ir*EIS 
        dP = k_cat*ES
        dEIS = k_if*I*ES - k_ir*EIS
        
        return(list(c(dE, dS, dI, dES, dP, dEIS)))
    })
}

noncompetitive_inhibition <- function(time, substrate, params) {
    with(as.list(c(substrate, params)), {
        dE = k_r*ES - k_f*S*E + k_cat*ES - k_if*I*E + k_ir*EI
        dS = k_r*ES - k_f*S*E - k_f*S*EI + k_r*EIS
        dI =  k_ir*EI - k_if*I*E + k_ir*EIS - k_if*ES*I
        dES = k_f*S*E - k_r*ES - k_cat*ES + k_ir*EIS - k_if*ES*I
        dP = k_cat*ES
        dEI = k_if*I*E - k_ir*EI - k_f*S*EI + k_r*EIS
        dEIS = k_f*S*EI - k_r*EIS - k_ir*EIS + k_if*ES*I
        return(list(c(dE, dS, dI, dES, dP, dEI, dEIS)))
    })
}

no_inhibition <- function(time, substrate, params) {
    with(as.list(c(substrate, params)), {
        dE = k_r*ES - k_f*S*E + k_cat*ES
        dS = k_r*ES - k_f*S*E
        dES = k_f*S*E - k_r*ES - k_cat*ES
        dP = k_cat*ES
        
        return(list(c(dE, dS, dES, dP)))
    })
}


inhibitions <- list("Competitive inhibition" = competitive_inhibition,
                    no_inhib = no_inhibition,
                    "Non-competitive inhibition" = noncompetitive_inhibition,
                    "Uncompetitive inhibition" = uncompetitive_inhibition)


ui <- navbarPage("Biology simulations",
                 navbarMenu("Ecology",
                            tabPanel("Lotka–Volterra",
                            titlePanel("Grobuonio ir plėšrūno modelis"), 
                            sidebarLayout(
                                sidebarPanel(
                                    sliderInput("rabbits",
                                                "Number of rabbits:",
                                                min = 0,
                                                max = 50,
                                                value = 30),
                                    sliderInput("foxes",
                                                "Number of foxes:",
                                                min = 0,
                                                max = 50,
                                                value = 30)
                                ),
                                
                                # Show a plot of the generated distribution
                                mainPanel(
                                    plotlyOutput("distPlot")
                                )
                            )
                 )
        ),
        navbarMenu("Biochemistry",
                   tabPanel("Enzyme kinetics",# Application title
                            titlePanel("Michaelis-Menten enzyme kinetics: types of inhibition"),
                            
                            # Sidebar with a slider input for number of bins 
                            sidebarLayout(
                                sidebarPanel(
                                    selectInput("type", "Type of inhibition:",
                                                c("Competitive inhibition",
                                                  "Non-competitive inhibition",
                                                  "Uncompetitive inhibition")),
                                    sliderInput("E",
                                                "Enzyme concentration:",
                                                min = 1,
                                                max = 50,
                                                value = 30),
                                    sliderInput("S",
                                                "Substrate concentration:",
                                                min = 1,
                                                max = 1000,
                                                value = 30),
                                    sliderInput("I",
                                                "Inhibitor concentration:",
                                                min = 0,
                                                max = 100,
                                                value = 1),
                                    sliderInput("k_f",
                                                "Rate of forward reaction:",
                                                min = 0,
                                                max = 1,
                                                value = 0.1),
                                    sliderInput("k_r",
                                                "Rate of reverse reaction:",
                                                min = 0,
                                                max = 1,
                                                value = 0.1),
                                    sliderInput("k_cat",
                                                "Rate of catalysis:",
                                                min = 0,
                                                max = 1,
                                                value = 0.1),
                                    sliderInput("k_3",
                                                "Rate of inhibitor-enzyme formation:",
                                                min = 0,
                                                max = 1,
                                                value = 0.1),
                                    sliderInput("k_3r",
                                                "Rate of inhibitor-enzyme dissociation:",
                                                min = 0,
                                                max = 1,
                                                value = 0.1)
                                ),
                                
                                # Show a plot of the generated distribution
                                mainPanel(
                                    plotOutput("kinetics_plot")
                                )
                            )))
)


server <- function(input, output) {
    
    ppred_func <- function (Time, State, Pars) {
        with(as.list(c(State, Pars)), {
            dx = x*(alpha - beta*y)
            dy = -y*(gamma - delta*beta*x)
            return(list(c(dx, dy)))
        })
    }
    Pars <- c(alpha = 1, beta = 0.1, gamma = 1.5, delta = .75)
    
    Time <- seq(0, 100, by = 1)
    
    generate_data <- reactive({
        State <- c(x = input$rabbits, y = input$foxes)
        out <- as.data.frame(ode(func = ppred_func, y = State, parms = Pars, times = Time)) %>% gather(key = "species", value = "pop", x, y)
        return(out)
    })
    
    output$distPlot <- renderPlotly({
        # draw the histogram with the specified number of bins
        print(
            ggplotly(ggplot(generate_data(), aes(x = time, y = pop, color = species)) + geom_line())
        )
    })
    
    enzyme_data <- reactive({
        
        
        params_inhib <- c(k_f = input$k_f, 
                          k_r = input$k_r, 
                          k_cat = input$k_cat, 
                          k_if = input$k_3, 
                          k_ir = input$k_3r
        )
        params_no_inhib <- c(k_f = input$k_f, 
                             k_r = input$k_r, 
                             k_cat = input$k_cat
        )
        
        substrate_inhib = c(E = input$E, 
                            S = input$S, 
                            I = input$I, 
                            ES = 0, 
                            P = 0,
                            EI = 0
        )
        
        substrate_no_inhib = c(E = input$E, 
                               S = input$S, 
                               ES = 0, 
                               P = 0
        )
        
        substrate_uncomp_inhib = c(E = input$E, 
                                   S = input$S, 
                                   I = input$I,
                                   ES = 0, 
                                   P = 0,
                                   EIS = 0
        )
        
        substrate_noncomp_inhib = c(E = input$E, 
                                    S = input$S, 
                                    I = input$I,
                                    ES = 0, 
                                    P = 0,
                                    EI = 0,
                                    EIS = 0
        )
        
        substrate <- list("Competitive inhibition" = substrate_inhib,
                          no_inhib = substrate_no_inhib,
                          "Non-competitive inhibition" = substrate_noncomp_inhib,
                          "Uncompetitive inhibition" = substrate_uncomp_inhib)
        
        comp_inhib <- as.data.frame(ode(func = inhibitions[[input$type]], 
                                        y = substrate[[input$type]], 
                                        parms = params_inhib, 
                                        times = time)) %>% mutate(Type = input$type)
        
        no_inhib <- as.data.frame(ode(func = no_inhibition, 
                                      y = substrate[["no_inhib"]], 
                                      parms = params_no_inhib, 
                                      times = time)) %>% mutate(Type = "Without inhibition")
        
        
        out <- bind_rows(comp_inhib,no_inhib) %>% gather(key = "ingre", value = "UI", -time,-Type) 
        return(out)
    })
    
    output$kinetics_plot <- renderPlot({
        ggplot(enzyme_data(), aes(x = time, y = UI, color = Type)) + 
            geom_line() + 
            facet_wrap(~ingre, scales = "free") + 
            theme_minimal() + 
            labs(x = "Time", y = "Units")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
