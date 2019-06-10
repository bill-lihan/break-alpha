# Authors: Victoria Savalei, Lihan (Bill) Chen, and Wolf Vanpaemel
library(shiny)
library(plotly) 

#------------------------------------------------------------------------------------------------------------------#
#Below are three R functions needed for this app: alpha.omega.pop, alpha.omega.pop.sim, and alpha.omega.pop.varyl

#alpha.omega.pop function computes alpha and omega (1-factor reliability) from the pop cov matrix with unit variances
#it takes only loadings (l) as input
alpha.omega.pop <- function(l){

	Sigmat <- l%*%t(l); Sigma <- Sigmat; diag(Sigma) <- 1
	totalvar <- sum(Sigma) #total variance
	p <- length(l) #number of items
	avecov <- (totalvar-p)/(p*(p-1))#average covariance
	alpha <- p^2*avecov/totalvar
	omega <- sum(Sigmat)/totalvar
	diff <- omega-alpha
	rels <- list(alpha,omega,diff)

	names(rels) <- c("alpha","omega","omina")

	return(rels)
}

#alpha.omega.pop.sim function computes alpha and omega (1-factor reliability) for a simulation based on the input U(a,b) distribution
#it takes only p (number of items), L (average loading) and R (loadings range) as input.
alpha.omega.pop.sim <- function(reps=1000,p=6,L=.7,R=.3){

	loadings <- runif(p*reps,min=L-.5*R,max=L+.5*R)

	rels <- matrix(0,nrow=reps,ncol=(4+p)) #

	for (i in (1:reps)) {

		l <- loadings[((i-1)*p+1):(i*p)]
		Sigmat <- l%*%t(l); Sigma<-Sigmat; diag(Sigma)<-1
		totalvar <- sum(Sigma) #total variance
		p <- length(l) #number of variables
		avecov <- (totalvar-p)/(p*(p-1))#average covariance

		alpha <- p^2*avecov/totalvar
		omega <- sum(Sigmat)/totalvar

		rels[i,1:2] <- c(alpha,omega)
		rels[i,3:(2+p)] <- l
		rels[i,(3+p)] <- max(l)-min(l)
		rels[i,(4+p)] <- omega-alpha

		colnames(rels) <- c("alpha","omega",paste("l",1:p,sep=""),
							"l.range","omega.minus.alpha")

	}

	return(rels)

}

#alpha.omega.pop.varyl function computes alpha and omega (1-factor reliability) for a simulation based on the user specified loadings
#it takes the p x 1 vector of user-specified loadings, l, as input
#it then selects one loading at a time and and varies it continuously from 0 to 1, keeping others the same.
alpha.omega.pop.varyl <- function(l){
	
	p <- length(l) #number of items
	reps <- seq(0,1,.05) #first loading will change #wolf: decreased granularity because the linetype became induistinguishable; was .01
	rels <- matrix(0,nrow=length(reps)*length(l),ncol=(6+p)) #stored reliabilities
	t <- 0 #counting index

	for (j in (1:p)) { #loops through all the loadings

		for (i in (1:length(reps))) { #loops through all the values of each loading
			t <- t+1
			ltemp <- l
			ltemp[j] <- reps[i] #change the jth loading
			Sigmat <- ltemp%*%t(ltemp); Sigma <- Sigmat; diag(Sigma) <- 1
			totalvar <- sum(Sigma) #total variance
			
			avecov <- (totalvar-p)/(p*(p-1))#average covariance
			alpha <- p^2*avecov/totalvar
			omega <- sum(Sigmat)/totalvar
			
			rels[t,1:2] <- c(alpha,omega)
			rels[t,3:(2+p)] <- ltemp
			rels[t,(3+p)] <- max(ltemp)-min(ltemp)
			rels[t,(4+p)] <- omega-alpha
			rels[t,(5+p)] <- j #the loading that is changing
			rels[t,(6+p)] <- ltemp[j] #the value of the changing loading
		} #end of i loop

	} #end of j loop

	colnames(rels) <- c("alpha","omega",paste("l",1:p,sep=""),"l.range","omega.minus.alpha","which.l","changing.loading")

	return(rels)
}

#--------------------------------------------------------------------------------------------------------------------#
# helper functions for repeated identical text formatting
renderLoadingRecap <- function(loadtxt, loadings, p, digits=3){
	return(renderText({
				tnoend <- paste(sapply(loadings[1:(p-1)], function(x){
								sprintf(paste0('%.', digits, 'f'), x)}), collapse=", ")
				tend <- paste(sapply(loadings[p], function(x){
								sprintf(paste0('%.', digits, 'f'), x)}), collapse=", ")
				paste(loadtxt, tnoend," and ", tend, ".", sep="")
			}))
}

renderRelOutput <- function(rel){
	return(renderText({ 
			 paste("Coefficient omega is ", round(rel$omega, 3), ". ",
			 			"Coefficient alpha is ", round(rel$alpha, 3), ". ",
			 			"Alpha underestimates omega by ", round(rel$omega-rel$alpha, 3), ", or ",
			 			ifelse(rel$omega-rel$alpha==0, 0, round((rel$omega-rel$alpha)*100/rel$omega, 2)),
			 			"%.", sep="")
	 		}))
}

#--------------------------------------------------------------------------------------------------------------------#
# UI code
M <- 1000 # constant; number of simultions for random generation; used in alpha.omega.pop.sim as reps input

ui <- fluidPage(

	h1("Break Coefficient Alpha!"),
	h4("How much can coefficient alpha differ from coefficient omega? This app will compute alpha and omega assuming a 1-factor 
		(congeneric, unidimensional) measurement model with user-specified or randomly generated factor loadings. In this scenario, alpha is a lower bound to internal consistency
		reliability (omega). But how big is the difference? Try to \"break\" alpha and find out!"),
	h5("Disclaimer: This app makes calculations using population covariance matrices 
      implied by the generated/specified loadings. Entering estimated loadings from your own dataset will not produce alpha 
      based on the observed covariance matrix. Do not use this app to obtain sample coefficient alpha."),
	h6("This app was developed by", HTML(paste0(a("Victoria Savalei", href="http://ubcsemlab.com/"), ",")), "Lihan (Bill) Chen and", a("Wolf Vanpaemel", href="https://ppw.kuleuven.be/okp/people/wolf_vanpaemel/byYearType/")),
	
	sidebarLayout(
		sidebarPanel(
			sliderInput("p", "Number of items on your scale:", 3, 30, 1, value=3), 
			
			radioButtons("custom", "How do you want to define factor loadings?", 
				c("Generate randomly"="rand","Specify manually"="spec")),
			
			conditionalPanel(
				condition = "input.custom == 'rand'",
				"Factor loadings will be drawn from a uniform distribution U(a,b) with mean L=(a+b)/2 and range R=(b-a). Setting R=0 implies 
				tau-equivalence.",br(),br(),
				sliderInput("aveloading", "Average Loading (L)", min=0, max=1,value=.5),
				uiOutput("sliderange"), actionButton("updateButton", "Break it once!"), br(), 
				actionButton("simButton",paste("Break it", M, "times and plot!",sep=" "))
				),
			
			conditionalPanel(
				condition = "input.custom == 'spec'",
				strong("Specify factor loadings for each variable:"),br(),
				uiOutput("slideload"),actionButton("updateButtonSpec", "Break it!")
				)
			),
		
		mainPanel(  
			conditionalPanel(
				condition = "input.custom == 'rand' && output.updateButton == 'GO'",
				h4(textOutput("printload")),
				h4(textOutput("alpha"))
				),
			
			conditionalPanel(
				condition = "input.custom == 'rand' && output.simButton == 'GO'",
				h4(textOutput("plottext")),
				plotlyOutput("plot1"),
				plotlyOutput("plot2")
				),
			
			conditionalPanel(
				condition = "input.custom == 'spec'",
				h4(textOutput("printloadspec")),
				h4(textOutput("alphaspec")), br(),
				h4(textOutput("plottextspec")),
				plotlyOutput("plotspec"),
				h4(textOutput("plottextspecadd"))
				)
			
			)
		)
	)

#--------------------------------------------------------------------------------------------------------------------#
# server code
server <- function(input, output) {

	nsdl = 3 #mumber of sign digits
	pm = 10 #cutoff for p to stop displaying some output to prevent clutter  
	
	#define input sliders   
	#for randomly generated factor loadings:
	output$sliderange <- renderUI({
		sliderInput("range", "Loadings Range (R)", min = 0, max = min(2*input$aveloading, 2*(1-input$aveloading)), value = min(input$aveloading, (1-input$aveloading))) 
	})
	
	#for user defined factor loadings: 
	numP <- reactive({input$p})  # get number of loadings from ui slider 
	
	output$slideload <- renderUI({   # create custom loading sliders
		lapply(1:(numP()), function(i) {
			sliderInput(inputId = paste0("load", i), # this creates input$load1, input$load2, etc.
				label = paste0("Loading ", i),
				min = 0, max = 1, value = .5)
		})
	})
	
	# button pressed to update randomly generated loadings
	observeEvent(input$updateButton,{      

		output$updateButton <- renderText({"GO"}) #for conditionalPanel to work
		output$simButton <- renderText({"NOGO"}) #for conditionalPanel to work
		outputOptions(output, "updateButton", suspendWhenHidden = FALSE)
		outputOptions(output, "simButton", suspendWhenHidden = FALSE)

		loadtxt <- "The generated loadings are " 
		genLoadingss <- runif(input$p, min=input$aveloading-.5*input$range, max=input$aveloading+.5*input$range)
		relss <- as.data.frame(alpha.omega.pop(genLoadingss))
		
		#define output text (recap)
		output$printload <- renderLoadingRecap(loadtxt, genLoadingss, isolate(input$p))
		output$alpha <- renderRelOutput(relss)

	}) # observeEvent(input$updateButton

	# button pressed to update specified loadings
	observeEvent(input$updateButtonSpec,{

		loadtxtspec <- "The specified loadings are " 
				
 		# specLoadings <- # defaults: we need this before numP() is reset
		all <- rep(.5, numP())
		# check for custom values
		for(i in 1:(numP())){
			# eval parse is ugly, but see comments on slideload
			eval(parse(text=paste0("all[", i, "] <- input$load", i)))
		}
		
		specLoadingss <- all # defaults: we need this before numP() is reset
		
		relss <- alpha.omega.pop(specLoadingss)
		
		output$printloadspec <- renderLoadingRecap(loadtxtspec, specLoadingss, isolate(input$p), 2)
		output$alphaspec <- renderRelOutput(relss)
				
		#define text for plot
		basictxt <- "The plot below varies each loading from 0 to 1 (value shown on the x-axis), while keeping the rest at the specified 
		values. The blue dots correspond to the specified value of each loading. You can hover over the curve to see specific values."

		output$plottextspec <- renderText({basictxt})

		addtxt <- "There are too many variables to include a legend; you can hover over each plot to find out which loading is changing."
		
		if(input$p<=pm){
			addtxt = paste0("Single-click on a linetype within the legend to remove/re-add the corresponding plot. ",
										"Double-click on a linetype to remove/re-add all other plots.")
		}  

		if(sum(duplicated(specLoadingss))){
			duptxt = paste0("The list of specified loadings contains at least one set of duplicate values. ",
										"Curves for duplicate loading values overlap. ")
		} else {
			duptxt = ""
		}

		output$plottextspecadd <- renderText({ 
			paste(duptxt, addtxt, sep="") 
		})

		#make plot
		relspecx <- as.data.frame(alpha.omega.pop.varyl(specLoadingss)) 

		#set up stuff for plotting
		xsta <- relspecx$changing.loading
		ysta <- relspecx$omega.minus.alpha
		xend <- c(xsta[2:length(xsta)],NA)
		yend <- c(ysta[2:length(ysta)],NA)

		gs <- nrow(relspecx)/input$p #gridsize, as set in alpha.omega.pop.varyl
		
		xend[seq(gs,nrow(relspecx),by=gs)]=xend[seq(gs,nrow(relspecx),by=gs)-1] 
		yend[seq(gs,nrow(relspecx),by=gs)]=yend[seq(gs,nrow(relspecx),by=gs)-1] 

		dimrel <- dim(relspecx)[1]
		
		for(i in 1:(numP())){
			addme=c(alpha.omega.pop(specLoadingss)$alpha,
							alpha.omega.pop(specLoadingss)$omega,
							t(specLoadingss),99,
							alpha.omega.pop(specLoadingss)$omina,i,
							t(specLoadingss[i]))
			relspecx=rbind(relspecx,addme)    
		}
		dimrel2 <- dim(relspecx)[1]

		#plot difference vs varing loading
		output$plotspec <- renderPlotly({  

			plt <- ggplot(data=relspecx[1:dimrel,],
								aes(text = paste0("omega = ",sprintf('%.3f', omega),
												"<br>alpha = ",sprintf('%.3f', alpha),
												"<br>alpha underestimates omega by ", sprintf('%.3f', omega-alpha),", or ",
												ifelse(omega-alpha==0, 0, sprintf('%.2f', (omega-alpha)*100/omega)), "%",
												"<br>changing loading value (loading ", which.l, ") = ",
												sprintf('%.2f', changing.loading)
							))) +
					geom_segment(data=relspecx[1:dimrel,],
							aes(x=-1, xend=-1, y=0, yend=0, linetype=as.factor(which.l)), color='black') + xlim(0,1) + #hack for solving ggplotly legend color issue
					geom_segment(data=relspecx[1:dimrel,],
							aes(x=changing.loading, xend=xend, y=omega.minus.alpha, yend=yend,
								linetype=as.factor(which.l), colour=omega)) +    
					geom_point(data=relspecx[(dimrel+1):dimrel2,],
							aes(x=changing.loading, y=omega.minus.alpha), colour= "blue") +
					xlab("Changing loading value") + ylab("Difference between omega and alpha") +
					scale_color_gradient2(name="omega",midpoint=.5,low="lightpink", mid="salmon",high="blue") +
					theme(legend.title = element_blank())
			
			if (isolate(input$p)<=pm){ # show linetype legend
				plt <- ggplotly(plt,tooltip = c("text")) %>%
						add_annotations(text="changing loading",  xref="paper", yref="paper",
							x=0.15, xanchor="right", 
							y=-.3, yanchor="bottom",    
							legendtitle=TRUE, showarrow=FALSE) %>% 
						layout(legend = list(orientation = "h", yanchor="bottom", xanchor="left",y = -.4, x =0.15))
			} else { # hide linetype legend
				plt <- ggplotly(plt,tooltip = c("text")) %>% 
						layout(showlegend = FALSE) 
			}
		})

	}) # closes observeEvent(input$updateButtonSpec

	observeEvent(input$simButton,{

		output$simButton <- renderText({"GO"}) #for conditionalPanel to work
		output$updateButton <- renderText({"NOGO"}) #for conditionalPanel to work
		outputOptions(output, "simButton", suspendWhenHidden = FALSE)
		outputOptions(output, "updateButton", suspendWhenHidden = FALSE)
		
		rellistx<-as.data.frame(alpha.omega.pop.sim(reps=M,p=input$p,L=input$aveloading,R=input$range))  
		
		#define output to appear on tooltip
		cols <- grep("^l[0-9]", colnames(rellistx)) #columns with loadings value
		loadlist <- colnames(rellistx)[cols]

		rellistx$all.loadings <- apply(round(rellistx[,loadlist],nsdl), 1, paste, collapse=" ")
		rellistx$min <- apply(round(rellistx[,loadlist],nsdl), 1, min)
		rellistx$max <- apply(round(rellistx[,loadlist],nsdl), 1, max)
		
		indloadtxt <- ""
		pSelect <- isolate(input$p)
		for(i in 1:(min(pSelect, pm))) {
			indloadtxt <- paste(indloadtxt, "<br>loading ", i, " = ",
						sprintf(paste0('%.', nsdl, 'f'), rellistx[,paste0("l",i)]))
		}

		if(pSelect > pm){
			indloadtxt <- paste0(indloadtxt, '\n...')
		}

		#define text
		output$plottext <- renderText({ 
				paste("In the plots below, each datapoint is based on alpha and omega for a 1-factor model with p = ",
					isolate(input$p), " loadings, randomly drawn from a uniform distribution with mean ",
					isolate(input$aveloading), " and range ", isolate(input$range),
					". You can hover over the curve to get specific values (only the first 10 loadings are shown).",
					sep="") 
			})
		
		#make alpha vs omega plot
		output$plot1 <- renderPlotly({
			pt <- ggplot(rellistx, 
						 aes(text = paste0("omega = ",sprintf('%.3f', omega),
							 	"<br>alpha = ", sprintf('%.3f', alpha),
							 	"<br>alpha underestimates omega by ",
						 		sprintf('%.3f', omega-alpha),
						 		", or ",
						 		ifelse(omega-alpha==0, 0, sprintf('%.2f', (omega-alpha)*100/omega)), "%",
								"<br>empirical loadings range = ",sprintf('%.3f', l.range),
								"<br>min loading = ",sprintf('%.3f', min),
								"<br>max loading = ",sprintf('%.3f', max),
								indloadtxt
							))) + 
				geom_point(mapping = aes(x=alpha, y=omega,color=l.range),size=1) +
				scale_color_gradient2("Empirical \n loadings range \n (largest minus \n smallest loading)",
					midpoint=.5,high="orange",mid="darkgreen",low="brown") +
				geom_abline(color="grey",slope=1, intercept=0) 
			
			pt <- ggplotly(pt,tooltip = c("text"))
		})

		#make difference vs loadings plot
		output$plot2 <- renderPlotly({
			pt <- ggplot(rellistx, 
					 aes(text = paste0("omega = ", sprintf('%.3f', omega),
					 				"<br>alpha = ", sprintf('%.3f', alpha),
								 	"<br>alpha underestimates omega by ",
							 		sprintf('%.3f', omega-alpha),
							 		", or ",
							 		ifelse(omega-alpha==0, 0, sprintf('%.2f', (omega-alpha)*100/omega)),"%",
									"<br>empirical loadings range = ", sprintf('%.3f', l.range),
									indloadtxt
						))) + 
			geom_point(mapping = aes(x=l.range, y=omega.minus.alpha,color=omega),size=1) +
			scale_color_gradient2(midpoint=.5,low="lightpink", mid="salmon",high="blue") +
			xlab("Empirical loadings range (largest minus smallest loading)")+
			ylab("Difference between omega and alpha")
			
			pt1 <- ggplotly(pt,tooltip = c("text") , dynamicTicks = TRUE)
		})
	}) #closes observeEvent(input$simButton

} #end of server

# start the app
shinyApp(ui = ui, server = server)