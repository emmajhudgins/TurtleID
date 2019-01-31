## Robust design analyses for Maldivian turtles - written by Emma J. Hudgins 2018. Adapted from Santostasi et al. 2016, 'A Robust Design Capture-Recapture Analysis of Abundance, Survival and Temporary Emigration of Three Odontocete Species in the Gulf of Corinth, Greece', available at: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166650

#Original script available at: https://dx.plos.org/10.1371/journal.pone.0166650.s005

library(RMark)

pseudo.input=read.csv("github_rcapture_sample.csv", header=TRUE, stringsAsFactors = F) 

# remove the lines of 0s (not relevant here)
mask <- apply(pseudo.input,1,sum)
pseudo.input <- pseudo.input[as.logical(mask),]

# format for RMark analysis, add character line of full counts
pseudo.input$ch <- apply(pseudo.input, 1 , paste , collapse = "" ) 

#group columns of data into open and closed population intervals, here we do 6 month blocks with individual months within each (must add up to number of columns in data -1), 0's are closed and 1's are open.
time.intervals<-c(rep(0, 5),1, rep(0, 5),1,rep(0, 5), 1, rep(0, 5),1,rep(0, 5),1, rep(0, 5),1,rep(0, 5), 1, rep(0, 5)) 

#ROBUST DESIGN ANALYSES#

# you can run this entire script at once and then look at the summary of the output to see which model fits best. Refer to Kendall et al. (2017) to interpret parameters (http://www.phidot.org/software/mark/docs/book/pdf/chap15.pdf)

### P(.)S(.)

# Markovian emigration
S.dot=list(formula=~1)
p.dot=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.01=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.dot),threads=2)


# Random emigration
S.dot=list(formula=~1)
p.dot=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.02=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.dot),threads=2)

# No emigration
S.dot=list(formula=~1)
p.dot=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.03=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.dot),threads=2)


### P(T)S(.)

# Markovian emigration
S.dot=list(formula=~1)
p.time=list(formula=~time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.04=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time),threads=2)


# Random emigration
S.dot=list(formula=~1)
p.time=list(formula=~time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.05=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time),threads=2)

# No emigration
S.dot=list(formula=~1)
p.time=list(formula=~time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.06=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time),threads=2)

### P(HET)S(.)

# Markovian emigration with heterogeneity (mixture) on detection
S.dot=list(formula=~1)
pi.dot=list(formula=~1)
p.het=list(formula=~mixture,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.07=mark(data = pseudo.input, model = "RDHet",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.het),threads=2)

# Random emigration with heterogeneity (mixture) on detection
S.dot=list(formula=~1)
pi.dot=list(formula=~1)
p.het=list(formula=~mixture,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE)
model.08=mark(data =pseudo.input, model = "RDHet",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.het),threads=2)

# No emigration with heterogeneity (mixture) on detection
S.dot=list(formula=~1)
pi.dot=list(formula=~1)
p.het=list(formula=~mixture,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0)
model.09=mark(data = pseudo.input, model = "RDHet",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.het),threads=2)


### P (SESSION) S(.)

# Markovian emigration
S.dot=list(formula=~1)
p.session=list(formula=~session,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.10=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.session),threads=2)


# Random emigration
S.dot=list(formula=~1)
p.session=list(formula=~session,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.11=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.session),threads=2)

# No emigration
S.dot=list(formula=~1)
p.session=list(formula=~session,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.12=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.session),threads=2)




### P (TIME.SESSION) S(.)

# Markovian emigration
S.dot=list(formula=~1)
p.time.session=list(formula=~-1+session:time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.13=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time.session),threads=2)


# Random emigration
S.dot=list(formula=~1)
p.time.session=list(formula=~-1+session:time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.14=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time.session),threads=2)

# No emigration
S.dot=list(formula=~1)
p.time.session=list(formula=~-1+session:time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.15=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time.session),threads=2)


### P(.)S(T)

# Markovian emigration
S.time=list(formula=~time)
p.dot=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.16=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.dot),threads=2)


# Random emigration
S.time=list(formula=~time)
p.dot=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.17=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.dot),threads=2)

# No emigration
S.time=list(formula=~time)
p.dot=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.18=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.dot),threads=2)

###P(T)S(T)

# Markovian emigration
S.time=list(formula=~time)
p.time=list(formula=~time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.19=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time),threads=2)


# Random emigration
S.time=list(formula=~time)
p.time=list(formula=~time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.20=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time),threads=2)

# No emigration
S.time=list(formula=~time)
p.time=list(formula=~1,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.21=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time),threads=2)

###P(HET)S(T)

# Markovian emigration with heterogeneity (mixture) on detection
S.time=list(formula=~time)
pi.dot=list(formula=~1)
p.het=list(formula=~mixture,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.22=mark(data = pseudo.input, model = "RDHet",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.het),threads=2)

# Random emigration with heterogeneity (mixture) on detection
S.time=list(formula=~time)
pi.dot=list(formula=~1)
p.het=list(formula=~mixture,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE)
model.23=mark(data = pseudo.input, model = "RDHet",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.het),threads=2)

# No emigration with heterogeneity (mixture) on detection
S.time=list(formula=~time)
pi.dot=list(formula=~1)
p.het=list(formula=~mixture,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0)
model.24=mark(data = pseudo.input, model = "RDHet",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.het),threads=2)



## P (SESSION) S(T)

# Markovian emigration
S.time=list(formula=~time)
p.session=list(formula=~session,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.25=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.session),threads=2)


# Random emigration
S.time=list(formula=~time)
p.session=list(formula=~session,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.26=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.session),threads=2)

# No emigration
S.time=list(formula=~time)
p.session=list(formula=~session,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.27=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.session),threads=2)

### P (TIME.SESSION) S(T)

# Markovian emigration
S.time=list(formula=~time)
p.time.session=list(formula=~-1+session:time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1)
GammaPrime.dot=list(formula=~1)
model.28=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaPrime=GammaPrime.dot,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time.session),threads=2)


# Random emigration
S.time=list(formula=~time)
p.time.session=list(formula=~-1+session:time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE) # gamma' = gamma''
model.29=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time.session),threads=2)

# No emigration
S.time=list(formula=~time)
p.time.session=list(formula=~-1+session:time,share=TRUE)
GammaDoublePrime.dot=list(formula=~1,share=TRUE,fixed=0) # gamma'=gamma''=0
model.30=mark(data = pseudo.input, model = "Robust",
              time.intervals=time.intervals,
              model.parameters=list(S=S.time,
                                    GammaDoublePrime=GammaDoublePrime.dot,
                                    p=p.time.session),threads=2)

results=collect.models()
results
saveRDS(results, file="github_models.rds")

#once you look at the results, see which model number has the lowest AIC, 
#then look at it directly with (summary(model.#))


