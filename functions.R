#===== FUNCTIONS TO SIMULATE OUTBREAKS, RECONSTRUCT THEM, ANALYZE AND PLOT RESULTS =====#

## Subset an outbreak object to keep only the first n.cases
subset.outbreak <- function(outbreak,to.keep){

  outbreak$n <- length(to.keep)
  outbreak$dna <- outbreak$dna[to.keep,]
  outbreak$onset <- outbreak$onset[to.keep]
  outbreak$id <- outbreak$id[to.keep]
  outbreak$ances <- outbreak$ances[to.keep]
  outbreak$group <- outbreak$group[to.keep]
  outbreak$nmut <- outbreak$nmut[to.keep]
  outbreak$ngen <- outbreak$ngen[to.keep]

  return(outbreak)
}

## Provide a mean R0, the overdispersion parameter k, and the upper limit
## beyond which you want to average the r value. The function returns
## the r value and the proportions they constitute
calc.offsp <- function(R0, k, cut) {
    
    ## Calculate the proportion of cases with that r
    calc.prop <- function(i, k, R0) {
        if(i == 0) return(pnbinom(0, size = k, mu = R0))
        return(pnbinom(i, size = k, mu = R0) -
               pnbinom(i - 1, size = k, mu = R0))
    }

    ## r and its proportion in the population
    df <- data.frame(r = 0:200,
                     prop = sapply(0:200, calc.prop, k, R0))

    ## split df into above and below cut
    ind <- which(df$r == cut)
    upper <- df[ind:nrow(df),]
    out <- df[1:ind,]
    out$prop[nrow(out)] <- sum(upper$prop)
    
    ## adjust prop to 1 so mean of subpop is accurate
    upper$prop <- upper$prop/sum(upper$prop)
    out$r[nrow(out)] <- sum(apply(upper, 1, prod))

    return(out)
    
}

## Provide a mean R0, the overdispersion parameter k, and the upper limit
## beyond which you want to average the r value. The function returns
## the r value and the proportions they constitute
calc.ss <- function(R0, k, cut) {
    
    ## Calculate the proportion of cases with that r
    calc.prop <- function(i, k, R0) {
        if(i == 0) return(pnbinom(0, size = k, mu = R0))
        return(pnbinom(i, size = k, mu = R0) -
               pnbinom(i - 1, size = k, mu = R0))
    }

    ## r and its proportion in the population
    df <- data.frame(r = 0:200,
                     prop = sapply(0:200, calc.prop, k, R0))

    ## split df into above and below cut
    ind <- which(df$r == cut)
    df1 <- df[1:(ind - 1),]
    df2 <- df[ind:nrow(df),]

    out <- data.frame(r = NA,
                      prop = c(sum(df1$prop), sum(df2$prop)))

    df1$prop <- df1$prop/out$prop[1]
    df2$prop <- df2$prop/out$prop[2]
    
    ## adjust prop to 1 so mean of subpop is accurate
    out$r[1] <- sum(apply(df1, 1, prod))
    out$r[2] <- sum(apply(df2, 1, prod))

    return(out)
    
}

## Returns a discretized gamma distribution
discr.gamma <- function(mean, sd) {
    w <- sapply(0:50, EpiEstim::DiscrSI, mean, sd)
    return(w)
}

## simulated contact tracing data
sim_ctd <- function(tTree, eps, lambda) {

    if(any(c(eps, lambda) < 0) | any(c(eps, lambda) > 1)) {
        stop('eps and lambda must be probabilities')
    }

    id <- unique(c(tTree[,1], tTree[,2]))
    id <- id[!is.na(id)]
    
    ## Sort tTree by value or alphabetically, This ensures A:B and B:A are both
    ## recognised when querying the contacts dataframe for transmission pairs
    tTree <- tTree %>%
        stats::na.omit() %>%
        apply(1, sort, FALSE) %>%
        t() %>%
        as.data.frame(stringsAsFactors = FALSE)

    tTree <- tTree[order(tTree[1]),]
    names(tTree) <- c('V1', 'V2')
    if(nrow(tTree) == 0) stop("No transmission observed")

    ## Create a dataframe of all potential contacts
    contacts <- as.data.frame(t(utils::combn(id, 2)))

    ## Create a column of logicals indicating whether a pair represents a
    ## transmission pair or not
    tTree$tPair <- TRUE

    ## Mark non-transmission pairs in the contacts dataframe. The merge function
    ## will mark pairs found in contacts but not in tTree as 'NA'. These are
    ## then converted to FALSE
    contacts <- merge(contacts, tTree, by = c('V1', 'V2'), all.x = TRUE)
    contacts$tPair[is.na(contacts$tPair)] <- FALSE

    ## Sample a number of rows given by a binomial distribution
    sampler <- function(x, prob) {
        x[sample(1:nrow(x), stats::rbinom(1, nrow(x), prob)), 1:3]
    }

    ## Sample transmission pairs with probability eps
    ## Sample non-transmission pairs with probability eps*lambda
    ctd <- rbind(sampler(contacts[contacts$tPair,], eps),
                 sampler(contacts[!contacts$tPair,], eps*lambda))

    rownames(ctd) <- NULL

    return(ctd)
}



#===== Functions to analyse outbreaks and store results =====#

## A function to compare accuracy and precision of CTD.outbreaker with DNA and with no DNA
compare.outbreakers <- function(runs, min.cases, max.cases, n.hosts, eps, xi,
                                w.dens, f.dens, R0, mu.transi, rate.import.case,
                                seq.length, group.freq, prior.eps,
                                print = FALSE, seed){
    
    analysis <- data.frame(variable = rep(c("accuracy","entropy"), 4*runs),
                           model=rep(c("t","tc","tg","tcg"), runs, each = 2),
                           value = NA)
    
    store.outbreak <- store.CTD <- store.t.result <- store.tc.result <-
    store.tg.result <- store.tcg.result <- list()

    ind  <- 0

    while((ind + 1) <= runs){

        browser()
        
        #tmp.seed <- (seed - 1)*runs + (ind + 1)
        
        #set.seed(tmp.seed)
        
        # simulate outbreak
        outbreak <- outbreaker::simOutbreak(R0, w.dens, n.hosts = n.hosts,
                                            max.cases, mu.transi = mu.transi,
                                            group.freq = group.freq,
                                            rate.import.case = rate.import.case,
                                            seq.length = seq.length)
        if(outbreak$n < min.cases) next
        if(outbreak$n > min.cases) {
            outbreak <- subset.outbreak(outbreak, 1:min.cases)
        }
       
        ## add delay between infection and sampling from f
        sample_delay  <- sample(0:(length(f.dens) - 1), outbreak$n,
                                prob = f.dens, replace = TRUE)
        outbreak$onset <- outbreak$onset + sample_delay
        
        print(ind + 1)

        ## simulated ctd
        tTree <- data.frame(i = outbreak$ances, j = outbreak$id)
        CTD <- sim_ctd(tTree, eps = eps, lambda = xi)[,1:2]
        if(nrow(CTD) == 0) CTD <- NULL
 
        ## time only
        data <- list(dates = outbreak$onset, w_dens = w.dens, f_dens = f.dens)
        config <- list(n_iter = 10000, sample_every = 50, find_import = FALSE,
                       move_kappa = FALSE, move_pi = FALSE, init_pi = 1,
                       init_tree = 'star')
        
        t.result <- outbreaker2::outbreaker(data = data, config = config)
        analysis$value[8*ind + 1]  <- get.acc(t.result, outbreak)
        analysis$value[8*ind + 2]  <- get.ent(t.result)

        ## time and contact
        data$ctd <- CTD
        
        tc.result <- outbreaker2::outbreaker(data = data, config = config)
        analysis$value[8*ind + 3]  <- get.acc(tc.result, outbreak)
        analysis$value[8*ind + 4]  <- get.ent(tc.result)
        
        ## time and genetic
        data$ctd <- NULL
        data$dna <- outbreak$dna
        config$init_tree  <- 'seqTrack'

        tg.result <- outbreaker2::outbreaker(data = data, config = config)
        analysis$value[8*ind + 5]  <- get.acc(tg.result, outbreak)
        analysis$value[8*ind + 6]  <- get.ent(tg.result)
        
        ## time, genetic and contact
        data$ctd <- CTD
       
        tcg.result <- outbreaker2::outbreaker(data = data, config = config)
        analysis$value[8*ind + 7]  <- get.acc(tcg.result, outbreak)
        analysis$value[8*ind + 8]  <- get.ent(tcg.result)
        
        ## remove genetic data to free up space
        outbreak$dna <- NULL

        ## store results        
        store.outbreak[[ind + 1]] <- outbreak
        store.CTD[[ind + 1]] <- CTD
        store.t.result[[ind + 1]] <- t.result
        store.tc.result[[ind + 1]] <- tc.result
        store.tg.result[[ind + 1]] <- tg.result
        store.tcg.result[[ind + 1]] <- tcg.result
        
        ind <- ind + 1

    }
    
    param <- list(n.hosts = n.hosts, min.cases = min.cases, R0 = R0,
                  w.dens = w.dens, n.hosts = n.hosts, max.cases = max.cases, 
                  mu.transi = mu.transi, rate.import.case = rate.import.case,
                  eps = eps, xi = xi, runs = runs, prior.eps = prior.eps)

    out <- list()
    out$param <- param
    out$outbreak <- store.outbreak
    out$CTD <- store.CTD
    out$t.result <- store.t.result
    out$tc.result <- store.tc.result
    out$tg.result <- store.tg.result
    out$tcg.result <- store.tcg.result
    out$analysis <- analysis

    return(out)
}

## Run outbreak comparisons across two SARS- and Ebola-like outbreaks
run.comparison <- function(param) {

    eps <- param$eps
    xi <- param$xi

    ## Defining global properties of run

    runs <- 20
    n.hosts <- 200
    min.cases <- 60
    max.cases <- min.cases + 1
    rate.import.case <- 0
    prior.eps <- c(0,1)

    ## Specifying ebola parameters
    ebola.w.dens <- discr.gamma(mean = 15.3, sd = 9.3)
    ebola.f.dens <- discr.gamma(mean = 5, sd = 4.7)

    ebola.R0 <- 1.5
    ebola.k <- 0.18
    #ebola.prop.ss <- 0.05
    #ebola.R0.ratio <- 5

    #ebola.small.R0 <- ebola.R0/(ebola.prop.ss*ebola.R0.ratio+(1-ebola.prop.ss))
    #ebola.R0.vec <- c(ebola.small.R0, ebola.R0.ratio*ebola.small.R0)
    #ebola.ss.vec <- c(1-ebola.prop.ss,ebola.prop.ss)
    ebola.R0.vec <- calc.offsp(ebola.R0, ebola.k, 100)$r
    ebola.ss.vec <- calc.offsp(ebola.R0, ebola.k, 100)$prop
    
    # Scale ebola.mu.transi by two thirds, as mu.transv is added to this
    ebola.mu.transi <- 1.24e-3/365*(2/3)
    ebola.seq.length <- 18958

    ## Specifying sars parameters
    sars.w.dens <- discr.gamma(mean = 8.4, sd = 3.8)
    sars.f.dens <- discr.gamma(mean = 4.85, sd = 3.89)    
    #sars.weib.shape <- (sars.w.sd/sars.w.mean)^-1.086
    #sars.weib.scale <- sars.w.mean/gamma(1+1/sars.weib.shape)
    #sars.w.dens <- dweibull(1:30,shape=sars.weib.shape,scale=sars.weib.scale)
    #shape <- sars.f.mean^2/sars.f.variance
    #rate <- sars.f.mean/sars.f.variance
    #sars.f.dens <- dgamma(1:30, shape = shape, rate = rate)

    sars.R0 <- 3.5
    sars.k <- 0.16
    #sars.prop.ss <- 0.05
    #sars.R0.ratio <- 10

    #sars.small.R0 <- sars.R0/(sars.prop.ss*sars.R0.ratio+(1-sars.prop.ss))
    #sars.R0.vec <- c(sars.small.R0,sars.R0.ratio*sars.small.R0)
    #sars.ss.vec <- c(1-sars.prop.ss,sars.prop.ss)
    sars.R0.vec <- calc.offsp(sars.R0, sars.k, 100)$r
    sars.ss.vec <- calc.offsp(sars.R0, sars.k, 100)$prop

    sars.mu.transi <- 1.14e-5*(2/3)
    sars.seq.length <- 29750

    ebola <- compare.outbreakers(runs = runs, min.cases = min.cases,
                                 n.hosts = n.hosts,eps = eps, xi = xi,
                                 w.dens = ebola.w.dens, f.dens = ebola.f.dens,
                                 R0 = ebola.R0.vec, mu.transi = ebola.mu.transi,
                                 seq.length = ebola.seq.length,
                                 max.cases = max.cases,
                                 group.freq = ebola.ss.vec,
                                 rate.import.case = rate.import.case,
                                 prior.eps = prior.eps, seed = param$seed)

    sars <- compare.outbreakers(runs = runs, min.cases = min.cases,
                                n.hosts = n.hosts,eps = eps, xi = xi,
                                w.dens = sars.w.dens, f.dens = sars.f.dens,
                                R0 = sars.R0.vec, mu.transi = sars.mu.transi,
                                seq.length = sars.seq.length,
                                max.cases = max.cases, group.freq = sars.ss.vec,
                                rate.import.case = rate.import.case,
                                prior.eps = prior.eps, seed = param$seed)

    out <- list(ebola = ebola, sars = sars)

    return(out)

}


#===== Connect to cluster and run =====#

## Connects to the cluster and returns the cluster object
set.up <- function(clust = 'hn') {

    setwd("/home/fc1915/mnt/fc1915/contacts")

    #didewin::didewin_config_global(credentials="fc1915",cluster="fi--dideclusthn")
    #didewin::web_login()

    options(didehpc.cluster = "fi--dideclusthn",
            didehpc.credentials = "~/.smbcredentials")

    if(clust == "mrc") options(didehpc.cluster = "fi--didemrchnb")
    
    didehpc::didehpc_config_global(temp = didehpc::path_mapping('tmp',
                                                                '/home/fc1915/mnt/tmp',
                                                                '//fi--didef3.dide.ic.ac.uk/tmp',
                                                                'T:'))
    
    didehpc::web_login()

    our.pkgs <- c("ggplot2","reshape2","outbreaker","dplyr","viridis","outbreaker2", "EpiEstim")
    pkg <- provisionr::package_sources(local="~/mnt/fc1915/contacts/outbreaker2_1.0-0.tar.gz")

    our.sources <- "~/mnt/fc1915/contacts/functions.R"

    ctx <- context::context_save("contexts", packages = our.pkgs,
                                 sources = our.sources,
                                 package_sources = pkg)

    obj <- didehpc::queue_didehpc(ctx)

    return(obj)
}

## Cluster connection for windows
set.up.wn <- function() {

    setwd("Q:\\outbreaker2")

    didewin::didewin_config_global(credentials="fc1915",cluster="fi--dideclusthn")
    didewin::web_login()

    our.pkgs <- c("ggplot2","reshape2","outbreaker2","outbreaker","dplyr","viridis")
    pkg <- context::package_sources(local="Q:\\outbreaker2")

    our.sources <- "Q:\\outbreaker2\\functions.R"
    
    ctx <- context::context_save("contexts_win",packages=our.pkgs,sources=our.sources,package_sources=pkg)

    obj <- didewin::queue_didewin(ctx)

    return(obj)
}

## Send the run to the cluster, with parameter values eps and xi
send.off <- function(obj, eps = NULL, xi = NULL, param = NULL) {

    if(is.null(eps)) eps <- seq(0,1,length=6)
    if(is.null(xi)) xi <- seq(0,1,length=6)
    if(is.null(param)) param <- expand.grid(xi = xi, eps = eps)

    param$seed <- seq_len(nrow(param))

    #queuer::enqueue_bulk(obj, param, run.comparison, do.call=FALSE, timeout=0)

    obj$enqueue_bulk(param, run.comparison, do_call = FALSE, timeout = 0)
    
}


#===== Functions to convert large stored data objects into smaller, useful R objects =====#

## Creates a new store object by loading previously saved results, and downloading cluster results
create.store <- function(obj, store = list(), bundle.name, dir, dl = TRUE,
                         load = FALSE, full = FALSE, saves = FALSE,
                         start = 1, addon = NULL){

    create.est <- function(disease,temp.runs){

        calc.est <- function(model, var){
            sapply(seq_len(temp.runs), function(i) mean(r[[disease]][[paste0(model,".result")]][[i]][[var]]))
        }

        calc.mu <- function(model, var) {

            genome.sizes <- data.frame(ebola = 18958, sars = 29750)

            sapply(seq_len(temp.runs), function(i) {
                outbreak <- r[[disease]]$outbreak[[i]]
                result <- r[[disease]][[paste0(model,".result")]][[i]]
                genome.size <- genome.sizes[[disease]]

                return(mean(get.mu(result, dna = outbreak$dna,
                                   burnin = 0, genome.size = genome.size)))
            })
        }

        data.frame(Model=rep(c("tcg","tc","tg","t"),3,each=temp.runs),
                   Variable=rep(c("eps","xi","mu"),each=4*temp.runs),
                   Value=c(calc.est("tcg","eps"),
                           calc.est("tc","eps"),
                           calc.est("tg","eps"),
                           calc.est("t","eps"),
                           calc.est("tcg","lambda"),
                           calc.est("tc","lambda"),
                           calc.est("tg","lambda"),
                           calc.est("t","lambda"),
                           calc.mu("tcg","mu"),
                           calc.mu("tc","mu"),
                           calc.mu("tg","mu"),
                           calc.mu("t","mu")))
    }

    create.CI <- function(disease,temp.runs){

        true.eps <- r[[disease]]$param$eps
        true.xi <- r[[disease]]$param$xi

        calc.CI <- function(model,var){
            true.val <- r[[disease]]$param[[var]]
            sapply(seq_len(temp.runs), function(i){
                distribution <- r[[disease]][[paste0(model,".result")]][[i]][[var]]
                quantile(distribution,0.025) <= true.val & true.val <= quantile(distribution,0.975)
            })
        }

        data.frame(Model=rep(c("tcg","tc"),2,each=temp.runs),
                   Variable=rep(c("eps","xi"),each=2*temp.runs),
                   Value=c(calc.CI("tcg","eps"),
                           calc.CI("tc","eps"),
                           calc.CI("tcg","xi"),
                           calc.CI("tc","xi")))
    }

    create.NTP <- function(i, r, disease) {

        outbreak <- r[[disease]]$outbreak[[i]]
        eps <- r[[disease]]$param$eps
        xi <- r[[disease]]$param$lambda
        n <- outbreak$n
        I <- sum(is.na(outbreak$ances))
        NTP <- analytical.NTP(eps, xi, n, I)
        return(NTP)

    }

    create.offsp <- function(r, disease) {
        as.vector(sapply(seq_len(r[[disease]]$param$runs),
                         function(i) get.offsp(r[[disease]]$outbreak[[i]])))
    }

    add.r <- function(store,r) {

        if(full) {
            store[[run]]$ebola <- r$ebola
            store[[run]]$sars <- r$sars
            return(store)

        } else {

            ## Adding Ebola
            temp.runs <- r$ebola$param$runs
            if(is.null(store[[run]]$ebola$param)) {
                store[[run]]$ebola$param  <- r$ebola$param
            } else {
                store[[run]]$ebola$param$runs  <- store[[run]]$ebola$param$runs + temp.runs
            }
            store[[run]]$ebola$analysis <- rbind(store[[run]]$ebola$analysis,r$ebola$analysis)

            ## Adding CI
            store[[run]]$ebola$CI <- rbind(store[[run]]$ebola$CI,create.CI("ebola",temp.runs=temp.runs))

            ## Adding est; the mean parameter estimate of the posterior distribution
            store[[run]]$ebola$est <- rbind(store[[run]]$ebola$est,create.est("ebola",temp.runs=temp.runs))

            ## Adding range; range is defined as n-2, and describes the range of the prior
            range <- mean(sapply(seq_len(temp.runs), function(i) r$ebola$outbreak[[i]]$n-2))
            store[[run]]$ebola$range <- mean(c(store[[run]]$ebola$range,range))

            ## Adding NTP; the number of potential non-transmission pairs
            NTP <- c(store[[run]]$ebola$NTP, sapply(1:temp.runs, create.NTP, r, "ebola"))
            store[[run]]$ebola$NTP <- NTP

            ## Adding the offspring distribution
            store[[run]]$ebola$offsp <- c(store[[run]]$ebola$offsp, create.offsp(r, 'ebola'))
            
            

            ## Adding SARS
            temp.runs <- r$sars$param$runs
            if(is.null(store[[run]]$sars$param)) store[[run]]$sars$param  <- r$sars$param
            else store[[run]]$sars$param$runs  <- store[[run]]$sars$param$runs + temp.runs
            store[[run]]$sars$analysis <- rbind(store[[run]]$sars$analysis,r$sars$analysis)

            ## Adding CI
            store[[run]]$sars$CI <- rbind(store[[run]]$sars$CI,create.CI("sars",temp.runs=temp.runs))

            ## Adding est
            store[[run]]$sars$est <- rbind(store[[run]]$sars$est,create.est("sars",temp.runs=temp.runs))

            ## Adding range
            range <- mean(sapply(seq_len(temp.runs), function(i) r$sars$outbreak[[i]]$n-2))
            store[[run]]$sars$range <- mean(c(store[[run]]$sars$range,range))

            ## Adding NTP; the number of potential non-transmission pairs
            NTP <- c(store[[run]]$sars$NTP, sapply(1:temp.runs, create.NTP, r, "sars"))
            store[[run]]$sars$NTP <- NTP

            ## Adding the offspring distribution
            store[[run]]$sars$offsp <- c(store[[run]]$sars$offsp, create.offsp(r, 'sars'))
            
            return(store)
        }
    }

    ## Load previously downloaded results into store
    if(load){
        files <- list.files(dir)
        files <- files[files!="store.RData"]
        pb <- txtProgressBar(min=1,max=length(files),style=3)
        for(file in files){
            setTxtProgressBar(pb,which(files==file))
            load(paste0(dir,file))
            run <- paste0("eps",r$ebola$param$eps,"xi",r$ebola$param$lambda)
            store <- add.r(store,r)
        }
    }

    ## Download results from the cluster
    if(dl) {
        task_bundle <- obj$task_bundle_get(bundle.name)
        ids <- task_bundle$ids
        pb <- txtProgressBar(min=1,max=length(ids),style=3)
        for(i in start:length(ids)){
            setTxtProgressBar(pb,i)
            task <- obj$task_get(ids[i])
            r <- task$result()
            run <- paste0("eps",r$ebola$param$eps,"xi",r$ebola$param$lambda)
            nm <- paste0(dir, run, addon, ".RData")
            if(nm %in% list.files(dir)) nm <- paste0(dir, run, addon,
                                                     runif(1, 0, 10000), ".RData")
            save(file=,r)
            store <- add.r(store,r)
            rm(r)
            gc()
        }
    }

    if(saves) save(store,file=paste0(dir,"store.RData"))

    return(store)
}

## Create store but without calculating mu estimate (takes ages)
create.nomu.store <- function(obj, store = list(), bundle.name, dir, dl = TRUE,
                              load = FALSE, full = FALSE,
                              start = 1, addon = NULL){

    create.est <- function(disease,temp.runs){

        calc.est <- function(model, var){
            sapply(seq_len(temp.runs),
                   function(i) mean(r[[disease]][[paste0(model,".result")]][[i]][[var]]))
        }

        data.frame(Model=rep(c("tcg","tc","tg","t"),2,each=temp.runs),
                   Variable=rep(c("eps","xi"),each=4*temp.runs),
                   Value=c(calc.est("tcg","eps"),
                           calc.est("tc","eps"),
                           calc.est("tg","eps"),
                           calc.est("t","eps"),
                           calc.est("tcg","lambda"),
                           calc.est("tc","lambda"),
                           calc.est("tg","lambda"),
                           calc.est("t","lambda")))
    }
    
    create.CI <- function(disease,temp.runs){

        true.eps <- r[[disease]]$param$eps
        true.xi <- r[[disease]]$param$lambda

        calc.CI <- function(model,var){
            true.val <- r[[disease]]$param[[var]]
            sapply(seq_len(temp.runs), function(i){
                distribution <- r[[disease]][[paste0(model,".result")]][[i]][[var]]
                quantile(distribution,0.025) <= true.val & true.val <= quantile(distribution,0.975)
            })
        }

        data.frame(Model=rep(c("tcg","tc"),2,each=temp.runs),
                   Variable=rep(c("eps","xi"),each=2*temp.runs),
                   Value=c(calc.CI("tcg","eps"),
                           calc.CI("tc","eps"),
                           calc.CI("tcg","xi"),
                           calc.CI("tc","xi")))
    }
    
    create.NTP <- function(i, r, disease) {

        outbreak <- r[[disease]]$outbreak[[i]]
        eps <- r[[disease]]$param$eps
        xi <- r[[disease]]$param$lambda
        n <- outbreak$n
        I <- sum(is.na(outbreak$ances))
        NTP <- analytical.NTP(eps, xi, n, I)
        return(NTP)

    }

    create.offsp <- function(r, disease) {
        as.vector(sapply(seq_len(r[[disease]]$param$runs),
                         function(i) get.offsp(r[[disease]]$outbreak[[i]])))
    }
    
    add.r <- function(store,r) {

        if(full) {
            store[[run]]$ebola <- r$ebola
            store[[run]]$sars <- r$sars
            return(store)

        } else {

            ## Adding Ebola
            temp.runs <- r$ebola$param$runs
            if(is.null(store[[run]]$ebola$param)) store[[run]]$ebola$param  <- r$ebola$param
            else store[[run]]$ebola$param$runs  <- store[[run]]$ebola$param$runs + temp.runs
            store[[run]]$ebola$analysis <- rbind(store[[run]]$ebola$analysis,r$ebola$analysis)

            ## Adding CI
            store[[run]]$ebola$CI <- rbind(store[[run]]$ebola$CI,create.CI("ebola",temp.runs=temp.runs))

            ## Adding est; the mean parameter estimate of the posterior distribution
            store[[run]]$ebola$est <- rbind(store[[run]]$ebola$est,create.est("ebola",temp.runs=temp.runs))

            ## Adding range; range is defined as n-2, and describes the range of the prior
            range <- mean(sapply(seq_len(temp.runs), function(i) r$ebola$outbreak[[i]]$n-2))
            store[[run]]$ebola$range <- mean(c(store[[run]]$ebola$range,range))

            ## Adding NTP; the number of potential non-transmission pairs
            NTP <- c(store[[run]]$ebola$NTP, sapply(1:temp.runs, create.NTP, r, "ebola"))
            store[[run]]$ebola$NTP <- NTP
            
            ## Adding the offspring distribution
            store[[run]]$ebola$offsp <- c(store[[run]]$ebola$offsp, create.offsp(r, 'ebola'))

            
            ## Adding SARS
            temp.runs <- r$sars$param$runs
            if(is.null(store[[run]]$sars$param)) store[[run]]$sars$param  <- r$sars$param
            else store[[run]]$sars$param$runs  <- store[[run]]$sars$param$runs + temp.runs
            store[[run]]$sars$analysis <- rbind(store[[run]]$sars$analysis,r$sars$analysis)

            ## Adding CI
            store[[run]]$sars$CI <- rbind(store[[run]]$sars$CI,create.CI("sars",temp.runs=temp.runs))

            ## Adding est
            store[[run]]$sars$est <- rbind(store[[run]]$sars$est,create.est("sars",temp.runs=temp.runs))

            ## Adding range
            range <- mean(sapply(seq_len(temp.runs), function(i) r$sars$outbreak[[i]]$n-2))
            store[[run]]$sars$range <- mean(c(store[[run]]$sars$range,range))

            ## Adding NTP; the number of potential non-transmission pairs
            NTP <- c(store[[run]]$sars$NTP, sapply(1:temp.runs, create.NTP, r, "sars"))
            store[[run]]$sars$NTP <- NTP

            ## Adding the offspring distribution
            store[[run]]$sars$offsp <- c(store[[run]]$sars$offsp, create.offsp(r, 'sars'))
            
            return(store)
        }
    }

    ## Load previously downloaded results into store
    if(load){
        files <- list.files(dir)
        #files <- files[-grep('store', files)]
        pb <- txtProgressBar(min=1,max=length(files),style=3)
        for(file in files){
            setTxtProgressBar(pb,which(files==file))
            load(paste0(dir,file))
            run <- paste0("eps",r$ebola$param$eps,"xi",r$ebola$param$xi)
            store <- add.r(store,r)
        }
    }
    
    ## Download results from the cluster
    if(dl) {
        task_bundle <- obj$task_bundle_get(bundle.name)
        ids <- task_bundle$ids
        pb <- txtProgressBar(min=1,max=length(ids),style=3)
        for(i in start:length(ids)){
            setTxtProgressBar(pb,i)
            task <- obj$task_get(ids[i])
            r <- task$result()
            run <- paste0("eps",r$ebola$param$eps,"xi",r$ebola$param$xi)
            nm <- paste0(dir, run,  runif(1, 0, 10000), addon, ".RData")
            save(file = nm, r)
            store <- add.r(store, r)
            rm(r)
            gc()
        }
    }
    
    return(store)
}

## Calculate the mutation rate mu in units per base per day
get.mu <- function(x, dna, burnin = 20000, genome.size) {

    if (!any(x$step > burnin))
        stop("burnin too high - no chain left")
    dat <- x[x$step > burnin, , drop = FALSE]

    ances <- dat[, grep("alpha", names(dat)), drop = FALSE]
    Tinf <- dat[, grep("t.inf", names(dat)), drop = FALSE]

    D <- as.matrix(ape::dist.dna(dna, model="N"))
    n <- ncol(ances)

    f1 <- function(vecAnces, vecTinf) {
        vecAnces <- as.integer(vecAnces)
        vecAnces[vecAnces < 1] <- NA
        nmut <- sapply(1:n, function(i) D[i, vecAnces[i]])
        vecTinf <- as.integer(vecTinf)
        deltaT <- vecTinf - vecTinf[vecAnces]
        out <- mean(nmut/deltaT, na.rm = TRUE)
        return(out)
    }

    out <- unlist(lapply(1:nrow(ances), function(i) f1(ances[i,
                                                             ], Tinf[i, ])))
    if (!is.null(genome.size))
        out <- out/genome.size

    return(out)

}

## Returns the tranmission tree of modal ancestries
get.tTree <- function(x) {

    ances <- grep("alpha", names(x))
    modal.ances <- sapply(ances, function(i) get.mode(x[[i]]))
    tTree <- data.frame(id = seq_along(ances), ances = modal.ances)

    return(tTree)

}

## Calculates the number of non-transmission pairs per person, from a tranmission tree and CTD
get.NTP <- function(tTree, CTD) {

    if(is.null(tTree) | is.null(CTD)) return(NULL)

    names(tTree) <- names(CTD) <- c("i", "j")

    # Merge returns the common rows between tTree and CTD, i.e. transmission pairs in the CTD
    TP.1 <- nrow(merge(tTree, CTD))

    # Swap columns, as contacts are symmetrical
    tTree <- data.frame(i = tTree$j, j = tTree$i)
    TP.2 <- nrow(merge(tTree, CTD))

    NTP <- 2*(nrow(CTD) - TP.1 - TP.2)/nrow(tTree)

    return(NTP)

}

## get the accuracy of outbreak reconstruction
get.acc <- function(result, outbreak) {

    id <- seq_len(outbreak$n)
    adder <- which(names(result)=="alpha_1") - 1
    samples <- length(result$step)

    #Determine the modal transmission network
    network <- data.frame(from=do.call(rbind, lapply(id,  function(i) ({
        modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
        if(length(modal.ances)==0) return(NA) else return(modal.ances)
    }))),  to=id)

    import <- which(is.na(network$from))

    transmission.id <- id[!sapply(result[id+adder], function(i) any(is.na(i)))]

    #Determine the proportion of correctly inferred ancestries
    num.correct <-  sum(outbreak$ances==network$from, na.rm=TRUE)
    num.correct <- num.correct + sum(is.na(outbreak$ances[is.na(network$from)]))
    acc <- round(num.correct/nrow(network), 2)

    return(acc)

}

## calculate the mean entropy of an outbreak run
get.ent <- function(res) {

    i <- grep('alpha', names(res))

    calc.ent <- function(x) {
        x <- as.character(x)
        fk <- table(x)/sum(table(x))
        return(-sum(log(fk)*fk))
    }
    
    ent <- sapply(res[i], calc.ent)

    return(round(mean(ent), 2))
    
}

## print errors of a cluster bundle
get.error <- function(bundle) {

    ids <- bundle$ids[bundle$status() == 'ERROR']
    for(id in ids) print(obj$task_get(id)$result()$message)
    
}

## return the offspring distribution of a simOutbreak objects
get.offsp <- function(sim) {

    offsp <- table(sim$ances)

    ## add zeroes which aren't listed by table
    offsp <- as.vector(c(rep(0, sim$n - length(offsp)), offsp))

    return(offsp)

}

## combine two r objects
combine.r <- function(tmp, r) {

    name <- names(r$ebola)[-c(1, 8)]
    for(dis in c('ebola', 'sars')){
        
        for(fac in name) {
            tmp[[dis]][[fac]] <- c(tmp[[dis]][[fac]], r[[dis]][[fac]])
        }

        tmp[[dis]]$analysis <- rbind(tmp[[dis]]$analysis, r[[dis]]$analysis)
        tmp[[dis]]$param$runs <- tmp[[dis]]$param$runs + r[[dis]]$param$runs
    
    }

    return(tmp)
}

## Calculates the theoretical NTP from CTD parameters
analytical.NTP <- function(eps, xi, n, I) {
    (n-3 + 2*I/n)*eps*xi
}

## Calculate xi from NTP
analytical.xi <- function(eps, ntp, n, I) {
    ntp/((n-3 + 2*I/n)*eps)
}

## Simple function to calculate the modal value of a vector
get.mode <- function(u) {
    u[which.max(tabulate(u))]
}

## Converts a comparison output into a wide dataframe from easy plotting, comparing accuracy and entropy
long2wide <- function(store,dat1="tcg",dat2="tg",subtract=FALSE){

  iterations <- length(store)
  par.names <- names(store)
  runs <- store[[1]]$ebola$param$runs

  wide <- as.data.frame(matrix(nrow=runs*iterations*4,ncol=7))
  names(wide) <- c("Disease","Eps","Xi","Variable","Value","Value2","Improvement")

  par.length <- runs*4

  counter <- 0

  for(run in par.names){

    eps.indices <- (counter*par.length+1):(par.length*(counter+1))

    wide$Eps[eps.indices] <- store[[run]]$ebola$param$eps
    wide$Xi[eps.indices] <- store[[run]]$ebola$param$xi

    eb.indices <- eps.indices[1:(length(eps.indices)/2)]

    wide$Disease[eb.indices] <- "Ebola"
    wide$Variable[eb.indices] <- c(rep("Accuracy",runs),rep("Entropy",runs))
    wide$Value[eb.indices] <- c(subset(store[[run]]$ebola$analysis,variable=="accuracy" & model==dat1)$value,
                                subset(store[[run]]$ebola$analysis,variable=="entropy" & model==dat1)$value)
    wide$Value2[eb.indices] <- c(subset(store[[run]]$ebola$analysis,variable=="accuracy" & model==dat2)$value,
                                 subset(store[[run]]$ebola$analysis,variable=="entropy" & model==dat2)$value)

    #N.store <- c()
    #for(i in 1:runs) N.store <- c(N.store,store[[run]]$ebola$outbreak[[i]]$n)
    #wide$N[eb.indices] <- rep(N.store,2)

    sars.indices <- tail(eps.indices,(length(eps.indices)/2))

    wide$Disease[sars.indices] <- "SARS"
    wide$Variable[sars.indices] <- c(rep("Accuracy",runs),rep("Entropy",runs))
    wide$Value[sars.indices] <- c(subset(store[[run]]$sars$analysis,variable=="accuracy" & model==dat1)$value,
                                  subset(store[[run]]$sars$analysis,variable=="entropy" & model==dat1)$value)
    wide$Value2[sars.indices] <- c(subset(store[[run]]$sars$analysis,variable=="accuracy" & model==dat2)$value,
                                   subset(store[[run]]$sars$analysis,variable=="entropy" & model==dat2)$value)

    #N.store <- c()
    #for(i in 1:runs) N.store <- c(N.store,store[[run]]$sars$outbreak[[i]]$n)
    #wide$N[sars.indices] <- rep(N.store,2)

    counter <- counter + 1
  }

  if(subtract){
    wide$Improvement <- wide$Value-wide$Value2
  } else{
    wide$Improvement <- 100*(wide$Value/wide$Value2 - 1)
  }

  wide$Improvement[wide$Variable=="Entropy"] <- -wide$Improvement[wide$Variable=="Entropy"]

  wide$Xi <- 1-wide$Xi

  return(wide)

}

## Work in progress to make the wide output even 'wider'
long2wwide <- function(store){

  out <- c()

  for(run in names(store)){
    for(disease in c("ebola","sars")){
      tmp <- store[[run]][[disease]]$analysis
      tmp <- cbind(tmp$variable,data.frame(split(tmp$value, tmp$model)))
      tmp$eps <- store[[run]][[disease]]$param$eps
      tmp$xi <- store[[run]][[disease]]$param$xi
      tmp$disease <- disease

      out <- rbind(out,tmp)
    }
  }

  names(out) <- c("Variable","T","TC","TG","TCG","Eps","Xi","Disease")

  out <- subset(out,Variable %in% c("accuracy","entropy"))

  out <- melt(out,id=c("Disease","Eps","Xi","Variable"),measure.vars=c("T","TC","TG","TCG"))

  names(out)[5] <- "Data"

  return(out)

}

## Simple function to capitalise a string
capitalise <- function(string) {
    substr(string, 1, 1) <- toupper(substr(string, 1, 1))
    return(string)
}

## Extracts a ggplot legend
extract.legend <- function(gplot){

    tmp <- ggplot_gtable(ggplot_build(gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)

}

## Arranges two plots horizontally with a join legend at the bottom (from p1)
joint.legend <- function(p1, p2, pos = 'bottom') {

    leg <- extract.legend(p1)

    if(pos == 'bottom') {
    p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = 'none'),
                                   p2 + theme(legend.position = 'none'),
                                   ncol = 2),
                      leg, nrow = 2, ncol = 1, heights = c(10, 0.5))
    }

    if(pos == 'right') {
        p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = 'none'),
                                      p2 + theme(legend.position = 'none'),
                                      nrow = 2),
                          leg, nrow = 1, ncol = 2, widths = c(10, 0.5))
    }

    return(p)

}

## Returns a named vector for colour/fill_manual; for consistent colour scheme
match.col <- function(var) {

    if(var == "data") {
        cols = rep(brewer.pal(4,"Set1"), 2)
        dats <- c("t", "tc", "tg", "tcg")
        names(cols) <- c(dats, toupper(dats))
        return(cols)
    }

    if(var == "disease") {
        pal <- wesanderson::wes_palette("Royal1",4)[c(3,2)]
        cols <- rep(pal, 2)
        names(cols) <- c("Ebola", "SARS", "ebola", "sars")
        return(cols)
    }
}

## Provide the disease name and its mutation rate will be returned
match.mu <- function(disease) {
    if(disease %in% c("ebola", "Ebola")) return(3.39726e-06)
    if(disease %in% c("sars", "SARS")) return(1.14e-5)
}

## Make a dataframe with mean accuracy values per grid point
mk.dwide <- function(store, mod = c("t", "tc", "tg", "tcg"), var = 'accuracy'){

    NTP <- analytical.NTP(1, 1, store[[1]]$ebola$param$min.cases, 1)

    dwide <- data.frame(matrix(nrow=length(store)*2*4,ncol=5,
                               dimnames=list(c(),c("Disease","Data","Eps","Xi","Accuracy"))))
    counter <- 1
    
    for(run in names(store)){
        for(disease in c("ebola","sars")){
            for(dat in c("t","tc","tg","tcg")){
                dwide$Disease[counter] <- disease
                dwide$Data[counter] <- toupper(dat)
                dwide$Eps[counter] <- store[[run]][[disease]]$param$eps
                dwide$Xi[counter] <- store[[run]][[disease]]$param$xi*NTP
                dwide$Accuracy[counter] <- mean(subset(store[[run]][[disease]]$analysis,
                                                       variable==var & model==dat)$value)
                counter <- counter + 1
            }
        }
    }

    df <- subset(dwide, Data %in% toupper(mod))
    df$Data <- factor(df$Data, levels = toupper(mod))

    return(df)

}

## Convert 'ebola' and 'sars' to 'EBOV' and 'SARS-CoV' respectively
mk.name <- function(x) {

    x[x == 'ebola'] <- 'EBOV'
    x[x == 'sars'] <- 'SARS-CoV'

    x
    
}


#===== Plotting functions ======#

## Plots a violin plot comparing t,tc,tg and tcg for Ebola and SARS, for a given parameter combination
vis.violin <- function(run){

    df <- rbind(subset(run$ebola$analysis,variable %in% c("accuracy","entropy")),
                subset(run$sars$analysis,variable %in% c("accuracy","entropy")))
    df$Disease <- rep(c("EBOV","SARS-CoV"),each=run$ebola$param$runs*8)
    names(df) <- c("Variable","Data","Value","Disease")

    df$Variable <- capitalise(as.character(df$Variable))
    df$Data <- factor(toupper(as.character(df$Data)),
                      levels=c("T","TC","TG","TCG"))
    df <- filter(df, Variable == 'Accuracy')

    p <- ggplot(df,aes(x=Data,y=Value,colour=Data)) + geom_sina() +
        facet_grid(. ~ Disease,scales="free") +
        ylab("Accuracy of outbreak reconstruction") +
        xlab("Outbreak data") +
        theme_set(theme_minimal(base_size = 18)) +
        scale_colour_brewer(name="Data",palette="Set1") +
        theme(legend.position  = 'none')

    return(p)
    
}

## Plots a line plot of accuracy and entropy for given xi, across all eps
vis.line <- function(store,xi=0){

    capitalise <- function(string) {
        substr(string, 1, 1) <- toupper(substr(string, 1, 1))
        return(string)
    }

    wwide <- long2wwide(store)

    col <- data.frame(dat=c("T","TC","TG","TCG"),col=brewer.pal(4,"Set1"))
    col$col <- as.character(col$col)

    wwide$Shift <- wwide$Eps + 0.0225
    wwide$Shift[wwide$Data=="TG"] <- wwide$Shift[wwide$Data=="TG"] - 0.015
    wwide$Shift[wwide$Data=="TC"] <- wwide$Shift[wwide$Data=="TC"] - 0.03
    wwide$Shift[wwide$Data=="T"] <- wwide$Shift[wwide$Data=="T"] - 0.045

    wwide$Disease <- capitalise(wwide$Disease)
    wwide$Disease[wwide$Disease=="Sars"] <- "SARS"
    wwide$Variable <- capitalise(as.character(wwide$Variable))

    wwide <- subset(wwide,Xi==xi)

    theme_set(theme_gray(base_size=15))

    p <- ggplot(wwide,aes(x=Eps,y=value,colour=Data)) + geom_point(aes(x=Shift),alpha=0.1,size=0.8) +
        stat_smooth(aes(fill=Data),level=0.99,method='loess') + facet_grid(Variable ~ Disease,scales="free") +
        scale_color_manual(values=col$col,labels=col$dat,guide=FALSE) +
        scale_fill_manual(values=col$col,guide=FALSE,labels = col$dat) +
        guides(fill   = guide_legend(override.aes=list(alpha=0.4,fill=col$col,
                                                       color=col$col))) +
        theme(legend.text.align=0.5) + xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +  ylab("Value")

    return(p)

}

## Plots the accuracy of transmission tree inference, with NTP as the y label
vis.acc <- function(store, mod = c("t", "tc", "tg", "tcg")){

    NTP <- analytical.NTP(1, 1, store[[1]]$ebola$param$min.cases, 1)
    
    df <- mk.dwide(store, mod)
    
    y_axs <- 1:length(store) %>%
        sapply(function(i) store[[i]]$ebola$param$xi) %>%
        unique %>%
        sort %>%
        "*"(NTP) %>%
        round(1)

    x_axs <- sort(unique(sapply(1:length(store), function(i) store[[i]]$ebola$param$eps)))

    mk.lab <- function(variable) {

        tmp <- c(sars = 'SARS-CoV', ebola = 'EBOV')
        tmp[variable]

    }


    p <- ggplot(df, aes(x = Eps, y = Xi, fill = Accuracy)) +
        geom_tile(color = "white") +
        facet_grid(Disease ~ Data, labeller = labeller(Disease = mk.lab)) +
        scale_fill_gradientn(colours=c("white","darkgreen"),limits=c(min(df$Accuracy),1)) +
        theme_set(theme_gray(base_size = 20)) +
        xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
        ylab(expression(paste("Number of non-transmission contacts (",psi,")"))) +
        scale_x_continuous(expand=c(0,0), breaks = x_axs) +
        scale_y_continuous(expand=c(0,0), breaks = y_axs) +
        theme(axis.ticks=element_blank(),panel.border = element_rect(fill=NA,colour="black",size=1)) +
        guides(fill = guide_colorbar(title.vjust = 0,title.hjust=0.5,reverse=TRUE,title="Accuracy"))

    p
}

## Plots the accuracy of transmission tree inference, with NTP as the y label
vis.ent <- function(store, mod = c("t", "tc", "tg", "tcg")){

    NTP <- analytical.NTP(1, 1, store[[1]]$ebola$param$min.cases, 1)
    
    df <- mk.dwide(store, mod, var = 'entropy')
    
    y_axs <- 1:length(store) %>%
        sapply(function(i) store[[i]]$ebola$param$xi) %>%
        unique %>%
        sort %>%
        "*"(NTP) %>%
        round(1)

    x_axs <- sort(unique(sapply(1:length(store), function(i) store[[i]]$ebola$param$eps)))

    mk.lab <- function(variable) {

        tmp <- c(sars = 'SARS-CoV', ebola = 'EBOV')
        tmp[variable]

    }

    p <- ggplot(df, aes(x = Eps, y = Xi, fill = Accuracy)) +
        geom_tile(color = "white") +
        facet_grid(Disease ~ Data, labeller = labeller(Disease = mk.lab)) +
        scale_fill_gradientn(colours=c("darkgreen","white"),
                             limits=c(min(df$Accuracy),max(df$Accuracy))) +
        theme_set(theme_gray(base_size = 20)) +
        xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
        ylab(expression(paste("Number of non-transmission contacts (",lambda,")"))) +
        scale_x_continuous(expand=c(0,0), breaks = x_axs) +
        scale_y_continuous(expand=c(0,0), breaks = y_axs) +
        theme(axis.ticks=element_blank(),
              panel.border = element_rect(fill=NA,colour="black",size=1)) +
        guides(fill = guide_colorbar(title.vjust = 0,title.hjust=0.5,reverse=TRUE,title="Entropy"))

    p
}

## Plot change in accuracy (absolute) between dat1 and dat2
## Set 'scaled' to TRUE if using the .scaled dataset
vis.rel <- function(store, dat1 = 'tc', dat2 = 'tg', var = 'accuracy'){

    NTP <- analytical.NTP(1, 1, store[[1]]$ebola$param$min.cases, 1)

    dwide <- mk.dwide(store, var = var)
    dwide$Disease <- mk.name(dwide$Disease)

    df1 <- filter(dwide, Data == toupper(dat1))
    df2 <- filter(dwide, Data == toupper(dat2))

    df1$Accuracy <- df1$Accuracy - df2$Accuracy

    if(var == 'entropy') {
        grad <- scale_fill_gradient2(low="darkgreen",mid="white",high="darkred")
    } else {
        grad <- scale_fill_gradient2(low="darkred",mid="white",high="darkgreen")
    }

    p <- ggplot(df1, aes(x = Eps, y = Xi, fill = Accuracy)) +
        geom_tile(color = "white") + facet_grid(Disease ~ ., labeller = label_parsed) +
        grad +
       #scale_fill_gradientn(colours=c("white","darkgreen"),limits=c(min(dwide$Accuracy),1)) +
        theme_set(theme_gray(base_size = 16)) +
        xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
        ylab(expression(paste("Number of non-infectious contacts per person (",psi,")"))) +
        #ylab(expression(paste("Number of non-transmission contacts (",lambda,")"))) +
        scale_x_continuous(expand=c(0,0),
                           breaks = round(sort(unique(df1$Eps)), 2)) +
        scale_y_continuous(expand=c(0,0),
                           breaks = round(sort(unique(df1$Xi)), 1)) +
        theme(axis.ticks   = element_blank(),
              panel.border = element_rect(fill=NA,colour="black",size=1)) +
        guides(fill = guide_colorbar(title.vjust = 0,
                                     title.hjust=0.5,
                                     reverse=TRUE,
                                     title=paste0("Change in\n", var)))

    p

}

## Plots the parameter estimate for mu
vis.est.mu <- function(store) {

    df1 <- NULL

    for(i in seq_along(store)) {
        for(disease in c("ebola","sars")) {
            for(mod in c("tg","tcg")) {
                temp <- subset(store[[i]][[disease]]$est, Model == mod & Variable == "mu")
                if(disease == "ebola") temp$Disease <- capitalise(disease)
                if(disease == "sars") temp$Disease <- toupper(disease)
                df1 <- rbind(df1, temp)
            }
        }
    }

    df1$Model <- toupper(df1$Model)

    df2 <- data.frame(Disease = c("Ebola","SARS"), Mu = c(3.39726e-06, 1.14e-5))

    p <- ggplot(df1, aes(Value)) + geom_density(aes(fill = Model), alpha=0.5) +
        facet_grid(. ~ Disease, scales = 'free') +
        ylab("Density") + xlab(expression(paste("Mutation rate ", mu, " (base",NULL^-1," day",NULL^-1,")"))) +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(name = "Data", values = match.col("data"))


    dat <- ggplot_build(p)
    yend <- 1.05 * max(dat$data[[1]]$density)
    p <- p + geom_segment(data = df2, aes(x = Mu, y = 0, xend = Mu, yend = yend),
                          linetype = "dotted", size = 1)


    p
}

## Plot the parameter estimates for eps and xi
## We average across dis and factor across tcg and tc
vis.grid.mu <- function(store, rm = c(0,1)) {

    df <- NULL

    for(i in seq_along(store)) {
        for(disease in c("ebola","sars")) {
            for(mod in c("tg","tcg")) {
                for(var in c("eps", "xi", "mu")){
                    temp <- subset(store[[i]][[disease]]$est, Model == mod & Variable == var)
                    temp$Eps <- store[[i]][[disease]]$param$eps
                    temp$Xi <- store[[i]][[disease]]$param$xi
                    if(disease == "ebola") temp$Disease <- "Ebola"
                    if(disease == "sars") temp$Disease <- "SARS"
                    df <- rbind(df, temp)
                }
            }
        }
    }

    df$Model <- toupper(df$Model)

    df1 <- subset(df, Variable == "mu" & Disease == "Ebola" & !Eps %in% rm & !Xi %in% rm)
    df2 <- subset(df, Variable == "mu" & Disease == "SARS"  & !Eps %in% rm & !Xi %in% rm)

    df1$lab.Eps <- factor(paste0("epsilon==",df1$Eps))
    df1$lab.Xi <- factor(paste0("lambda==",df1$Xi))

    df2$lab.Eps <- factor(paste0("epsilon==",df2$Eps))
    df2$lab.Xi <- factor(paste0("lambda==",df2$Xi))

    levels(df1$lab.Xi) <- rev(levels(df1$lab.Xi))
    levels(df2$lab.Xi) <- rev(levels(df2$lab.Xi))

    df1$Mu <- match.mu("ebola")
    df2$Mu <- match.mu("sars")

    create.histgrid <- function(df, main) {

        p <- ggplot(df1, aes(Value)) + geom_density(aes(fill = Model), alpha = 0.5) +
            scale_x_continuous(breaks = c(0.5)) +
            facet_grid(lab.Xi ~ lab.Eps, labeller = label_parsed, scales = 'free') +
            scale_fill_manual(values = match.col("data")) +
            ggtitle(main) + ylab("Density") +
            xlab(expression(paste("Mutation rate ", mu, " (base", NULL^-1, " day", NULL^-1, ")"))) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text.y      = element_blank(),
                  axis.ticks.y     = element_blank(),
                  legend.position  = "bottom")

        dat <- ggplot_build(p)
        yend <- max(dat$data[[1]]$density)
        p <- p + geom_segment(aes(x = Mu, y = 0, xend = Mu, yend = yend),
                              linetype = "dotted", size = 1)

        return(p)
    }


    p1 <- create.histgrid(df1,
                          expression(paste("Estimating the mutation rate of Ebola")))
    p2 <- create.histgrid(df2,
                          expression(paste("Estimating the mutation rate of SARS")))

    p <- joint.legend(p1, p2, 'bottom')

    return(p)
}

## Plot the paramters estimates for eps and xi
## We average across dis and factor across tcg and tc
vis.est <- function(store, dis = c("Ebola", "SARS"), rm = NULL) {

    df <- NULL

    for(i in seq_along(store)) {
        for(disease in c("ebola","sars")) {
            for(mod in c("tc","tcg")) {
                for(var in c("eps","xi")){
                    temp <- subset(store[[i]][[disease]]$est, Model == mod & Variable == var)
                    temp$Eps <- store[[i]][[disease]]$param$eps
                    temp$Xi <- store[[i]][[disease]]$param$xi
                    if(disease == "ebola") temp$Disease <- "EBOV"
                    if(disease == "sars") temp$Disease <- "SARS"
                    df <- rbind(df, temp)
                }
            }
        }
    }

    df$Model <- toupper(df$Model)
    names(df)[which(names(df) == 'Model')] <- 'Data'

    df1 <- subset(df, Variable == "eps" & Disease %in% dis & !Eps %in% rm & !Xi %in% rm)
    df2 <- subset(df, Variable == "xi" & Disease %in% dis & !Eps %in% rm & !Xi %in% rm)

    df1$lab.Eps <- factor(paste0("epsilon==",df1$Eps))
    df1$lab.Xi <- factor(paste0("lambda==",round(df1$Xi, 2)))

    df2$lab.Eps <- factor(paste0("epsilon==",df2$Eps))
    df2$lab.Xi <- factor(paste0("lambda==", round(df2$Xi, 2)))

    levels(df1$lab.Xi) <- rev(levels(df1$lab.Xi))
    levels(df2$lab.Xi) <- rev(levels(df2$lab.Xi))

    create.histgrid <- function(df, var, main, lab.x) {

        p <- ggplot(df, aes(Value)) + geom_density(aes(y = ..scaled.., fill = Data), alpha = 0.5) +
            scale_x_continuous(breaks = c(0.5)) +
            facet_grid(lab.Xi ~ lab.Eps, labeller = label_parsed, scales = 'free') +
            scale_fill_manual(values = match.col("data")) +
            ggtitle(main) + xlab(lab.x) + ylab("Density") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text.y      = element_blank(),
                  axis.ticks.y     = element_blank(),
                  legend.position  = "bottom")

        dat <- ggplot_build(p)
        yend <- max(dat$data[[1]]$density)
        p <- p + geom_segment(aes_string(x = var, y = 0, xend = var, yend = 1),
                              linetype = "dotted", size = 1)

        return(p)
    }

    p1 <- create.histgrid(df1, "Eps",
                          expression(paste("Estimating the reporting sensitivity (", epsilon, ")")),
                          expression(paste("Reporting sensitivity (", epsilon, ")")))
    p2 <- create.histgrid(df2, "Xi",
                          expression(paste("Estimating the contact specificity (", lambda, ")")),
                          expression(paste("Contact specificity (", lambda, ")")))

    p <- joint.legend(p1, p2, 'bottom')

    return(p)
}

## Plots the likelihood of the normal likelihood and the binomial likelihood
vis.mu.ll <- function() {

    our.f <- function(x, size, p) (p^x)*(1 - p)^(size - x)

    x <- 7
    seq <- 1:10

    mu.dist <- mean(sample(seq, size = 10000, replace = TRUE, prob = our.f(seq, 10000, x/10000)))
    count <- sum(sample(seq, size = 10000, replace = TRUE, prob = our.f(seq, 10000, x/10000)))

    df <- data.frame(x = seq,
                     Outbreaker2 = c(log(our.f(x, 10000, (seq)/10000)),
                                     log(our.f(seq, 10000, x/10000))),
                     Binomial = c(dbinom(x, 10000, (seq)/10000),
                                  dbinom(seq, 10000, x/10000)))

    df$func <- rep(c("L", "p"), 1, each = length(df$Outbreaker2)/2)
    mlt <- melt(df, id = c('x', 'func'))
    mlt$func <- factor(mlt$func, levels = c("p", "L"))

    levels(mlt$func) <- c("p", "L")
    levels(mlt$func) <- c(expression(paste("p(D | ", mu,")")),
                          expression(paste("p(", mu, " | D)")))

    p <- ggplot(mlt, aes(x, value, colour = variable)) +
        geom_line() + facet_grid(variable ~ func, scales = 'free', labeller = label_parsed) +
        ylab("Probability") + xlab(expression(paste("D or ", mu))) +
        geom_vline(xintercept = x, linetype = 'dotted') +
        theme(legend.position = 'none')

    p

}

## Plots the offspring distribution of a simOutbreak object
vis.offsp <- function(sim) {
    
    df <- data.frame(offsp = get.offsp(sim))

    ggplot(df, aes(offsp)) + geom_histogram(binwidth = 1)

}

## Plot the offspring distribution across all simulations
vis.all.offsp <- function(store) {

    df <- data.frame(disease = character(), offsp = numeric())
    
    for(i in seq_along(store)) {
        for(dis in c('ebola', 'sars')) {
            tmp <- data.frame(disease = dis,
                              offsp = store[[i]][[dis]]$offsp)
            df <- rbind(df, tmp)
        }
    }

    df2 <- data.frame(x = 0:max(df$offsp),
                      y = dnbinom(0:max(df$offsp), mu = 3.5, size = 0.16))
    
    ggplot(df, aes(x = offsp, fill = disease)) +
        geom_histogram(aes(y = ..count../sum(..count..)),
                       position = 'dodge',
                       binwidth = 1) +
        labs(x = 'Offspring', y = 'Proportion') +
        geom_point(data = df2, aes(x, y), fill = 'black', colour = 'black')

}


#===== Plot saving functions =====#

## Save your plot to the write up directory
g.save <- function(p, name, ext = 'svg', ...){

    ggsave(p,file=paste0('~/Dropbox/phd/contacts/figs/',
                         name,".", ext),...)

    return(NULL)
}

## Update all plots in the writeup directory
update.plots <- function(store){

    p <- vis.acc(store, mod = c("t", "tc", "tg", "tcg"))
    g.save(p, "acc", width = 25, height = 13, ext = 'png')

    p <- vis.ent(store, mod = c("t", "tc", "tg", "tcg"))
    g.save(p, "ent", width = 25, height = 13, ext = 'png')

    pp <- vis.est(store)
    g.save(p, "est", width = 20, height = 9, ext = 'png')

    p <- vis.rel(store)
    g.save(p, "rel", width = 10, height = 5, ext = 'png')

    p <- vis.grid.mu(store)
    g.save(p, "grid.mu", width = 20, height = 9)

}


#===== Deprecated functions =====#

## Plot the accuracy of transmission tree inference, for large and small
vis.acc.size <- function(store.s,store.l,mod=c("t", "tc", "tg", "tcg")){

    p1 <- plot.acc(store.s, mod = mod)
    p2 <- plot.acc(store.l, mod = mod)

    p <- joint.legend(p1, p2, 'right')

}

## Plots a heatmap of the accuracy in estimating eps and xi (plot.abs.grid without transmission accuracy)
vis.est.grid <- function(store,model="tcg"){

  capitalise <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }

  dwide <- data.frame(matrix(nrow=length(store)*2*2,ncol=5,dimnames=list(c(),c("Disease","Variable","Eps","Xi","Accuracy"))))
  counter <- 1

  for(run in names(store)){
    for(disease in c("ebola","sars")){
      for(variable in c("eps","xi")){
        if(disease=="sars") dwide$Disease[counter] <- toupper(disease)
        else dwide$Disease[counter] <- capitalise(disease)
        dwide$Variable[counter] <- capitalise(variable)
        dwide$Eps[counter] <- store[[run]][[disease]]$param$eps
        dwide$Xi[counter] <- store[[run]][[disease]]$param$xi
        true <- store[[run]][[disease]]$param[[variable]]
        est <- mean(subset(store[[run]][[disease]]$est,Model==model & Variable==variable)$Value)
        dwide$Accuracy[counter] <- 1-abs(true-est)
        counter <- counter + 1
      }
    }
  }

  dwide$Xi <- 1-dwide$Xi
  dwide$Variable <- factor(dwide$Variable)
  levels(dwide$Variable) <- c(expression(epsilon),expression(lambda))

  p <- ggplot(dwide,aes(x=Eps,y=Xi,fill=Accuracy)) + xlab("Reporting sensitivity") + ylab("False positive contact rate") +
    geom_tile(color="white") + facet_grid(Variable ~ Disease,labeller=label_parsed) +
    scale_fill_gradientn(name="Difference",colours=c("white","darkgreen"),limits=c(0,1)) +
    theme_set(theme_gray(base_size = 16)) +
    xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) + ylab(expression(paste("Contact specificity (",lambda,")"))) +
    scale_x_continuous(breaks=c(0.1,0.4,0.7,1),expand=c(0,0)) +
    scale_y_continuous(breaks=seq(0,1,0.2),expand=c(0,0)) +
    theme(axis.ticks           = element_blank(),
          panel.border         = element_rect(fill=NA,colour="black",size=1),
          legend.text.align    = 0.5,
          legend.title.align   = 0.5) +
    guides(fill=guide_colorbar(reverse=TRUE))

  p

}

## Plots a violin plot comparing t,tc,tg and tcg for Ebola and SARS, for a given parameter combination
vis.violin.2 <- function(par.run){

    total.plot <- rbind(subset(par.run$ebola$analysis,variable %in% c("accuracy","entropy")),
                        subset(par.run$sars$analysis,variable %in% c("accuracy","entropy")))
    total.plot$Disease <- rep(c("EBOV","SARS"),each=par.run$ebola$param$runs*8)
    names(total.plot) <- c("Variable","Model","Value","Disease")

    total.plot$Variable <- capitalise(as.character(total.plot$Variable))
    total.plot$Model <- factor(toupper(as.character(total.plot$Model)),levels=c("T","TC","TG","TCG"))

    p <- ggplot(total.plot,aes(x=Model,y=Value,colour=Model)) + geom_sina() +
        facet_grid(Variable ~ Disease,scales="free") + ylab("Value") + theme_set(theme_gray(base_size = 16)) +
        scale_colour_brewer(name="Data",palette="Set1") +
        theme(axis.title.x       = element_blank(),
              axis.text.x        = element_blank(),
              axis.ticks.x       = element_blank(),
              panel.grid.major     = element_blank(),
              legend.text.align  = 0.5,
              legend.title.align = 0)
    p

  return(p)
}

## Plots a histogram of the accuracy in estimating eps and xi
vis.est.hist <- function(store){

    capitalise <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }

    runs <- store[[1]][["ebola"]]$param$runs
    dwide <- data.frame(matrix(nrow=length(store)*2*2*2*runs,ncol=6,
                               dimnames=list(c(),c("Disease","Data","Variable","Eps","Xi","Accuracy"))))
    counter <- 1

    for(run in names(store)){
        for(disease in c("ebola","sars")){
            for(variable in c("eps","xi")){
                for(Data in c("tcg","tc")){
                    for(i in seq_len(runs)){
                        if(disease=="sars") dwide$Disease[counter] <- toupper(disease)
                        else dwide$Disease[counter] <- capitalise(disease)
                        dwide$Data[counter] <- toupper(Data)
                        dwide$Variable[counter] <- capitalise(variable)
                        dwide$Eps[counter] <- store[[run]][[disease]]$param$eps
                        dwide$Xi[counter] <- store[[run]][[disease]]$param$xi
                        true <- store[[run]][[disease]]$param[[variable]]
                        est <- subset(store[[run]][[disease]]$est,Model==Data & Variable==variable)$Value[i]
                        if(variable %in% c("eps","xi")) dwide$Accuracy[counter] <- 1-abs(true-est)
                        if(variable=="xi1") dwide$Accuracy[counter] <- 1-abs(true-est)/store[[run]][[disease]]$range
                        counter <- counter + 1
                    }
                }
            }
        }
    }

    dwide$Variable <- factor(dwide$Variable)
    levels(dwide$Variable) <- c(expression(epsilon),expression(lambda))

    cols <- c("firebrick","forestgreen")

    q <- ggplot(dwide,aes(x=Accuracy,colour=Data,fill=Data)) +
        facet_grid(Variable ~ Disease,labeller=label_parsed,scales='fixed') + geom_density(size=1,alpha=0.3) +
        scale_colour_manual(values=cols) + scale_fill_manual(values=cols) +
        theme_set(theme_gray(base_size = 16)) +
        scale_x_continuous(name="Accuracy",limits=c(0,1)) +
        scale_y_continuous(expand=c(0.05,0.0),name="Density") +
        theme(axis.ticks=element_blank(),legend.text.align=0.5,legend.title.align=0.5)

    q

}

## Plots a heat map of relative performance of two models
vis.rel.grid <- function(store,dat1="tcg",dat2="tg",subtract=FALSE,prop=FALSE,size=16){

    if(prop & !subtract) stop("subtract must be TRUE")
    wide <- long2wide(store,dat1,dat2,subtract=subtract)

    if(prop){
        wide$Improvement[wide$Variable=="Accuracy"] =
            wide$Improvement[wide$Variable=="Accuracy"]/
            (1-wide$Value2[wide$Variable=="Accuracy"])

        wide$Improvement[wide$Variable=="Entropy"] =
            wide$Improvement[wide$Variable=="Entropy"]/
            wide$Value2[wide$Variable=="Entropy"]
    }

    dwide <- summarise(group_by(wide, Disease, Variable, Eps,Xi),mean=mean(Improvement))
    if(subtract) leg.lab <- "Difference" else leg.lab <- "%\nChange"
    if(prop) leg.lab <- "Quality"

    p <- ggplot(dwide,aes(x=Eps,y=Xi,fill=mean)) + geom_tile(color="white") + facet_grid(Variable ~ Disease) +
        scale_fill_gradient2(name=leg.lab,low="darkred",mid="white",high="darkgreen") +
        theme_set(theme_gray(base_size = size)) +
        xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
        ylab(expression(paste("Contact specificity (",lambda,")"))) +
        scale_x_continuous(breaks=c(0.1,0.4,0.7,1),expand=c(0,0)) +
        scale_y_continuous(breaks=seq(0,1,0.2),expand=c(0,0)) +
        theme(axis.ticks=element_blank(),panel.border = element_rect(fill=NA,colour="black",size=1)) +
        guides(fill = guide_colorbar(title.vjust = 0,title.hjust=0.5,reverse=TRUE))
    p

    return(p)
}

## Plot.rel.grid but with seperate legends (took ages to make)
vis.rel.grid.2lg <- function(store, dat1 = "tcg", dat2 = "tg", subtract = FALSE, size=12){

    plot.margin = unit(c(1,0,-1,0.2), "line")

    create.grid.acc <- function(dwide,axis.title.x,size){
        p <- ggplot(dwide,aes(x=Eps,y=Xi,fill=mean)) +
            geom_tile(color="white") + facet_grid(. ~ Disease) +
            scale_fill_gradient2(name="Difference",low="darkred",mid="white",high="darkgreen") +
            xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
            theme_set(theme_gray(base_size = size)) +
            scale_x_continuous(breaks=c(0.1,0.4,0.7,1),expand=c(0,0)) +
            scale_y_continuous(breaks=seq(0,1,0.2),expand=c(0,0)) +
            theme(axis.ticks           = element_blank(),
                  legend.text.align    = 0.5,
                  legend.title.align   = 0.5,
                  axis.title.x         = axis.title.x,
                  axis.title.y         = element_blank(),
                  axis.text.x          = element_text(colour="transparent"),
                  plot.margin          = unit(c(1,0,-0.7,0.2), "line"),
                  panel.border         = element_rect(fill=NA,colour="black",size=1)) +
            guides(fill=guide_colorbar(reverse=TRUE))
        p
    }

    create.grid.ent <- function(store,axis.title.x,size){

        wide <-

        p <- ggplot(dwide,aes(x=Eps,y=Xi,fill=mean)) +
            geom_tile(color="white") + facet_grid(. ~ Disease) +
            scale_fill_gradient2(name="Difference",low="darkred",mid="white",high="darkgreen") +
            xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
            theme_set(theme_gray(base_size = size)) +
            scale_x_continuous(breaks=c(0.1,0.4,0.7,1),expand=c(0,0)) +
            scale_y_continuous(breaks=seq(0,1,0.2),expand=c(0,0)) +
            theme(axis.ticks           = element_blank(),
                  legend.text.align    = 0.5,
                  legend.title.align   = 0.5,
                  axis.title.x         = element_text(vjust=1),
                  axis.title.y         = element_blank(),
                  strip.text.x         = element_blank(),
                  plot.margin          = unit(c(0.3,0,0,0.2), "line"),
                  panel.border         = element_rect(fill=NA,colour="black",size=1)) +
            guides(fill=guide_colorbar(reverse=TRUE))
        p
    }

  wide <- long2wide(store,dat1,dat2,subtract=subtract)
  dwide <- summarise(group_by(wide, Disease, Variable, Eps,Xi),mean=mean(Improvement))

  p <- plot.rel.grid(store,size=size)

  p1 <- create.grid.acc(subset(dwide,Variable=="Accuracy"),element_blank(),size=size)
  p1$facet$rows <- p$facet$rows

  p2 <- create.grid.ent(subset(dwide,Variable=="Entropy"),NULL,size=size)
  p2$facet$rows <- p$facet$rows

  label1 = textGrob(expression(paste("Contact specificity (",lambda,")")), rot = 90, vjust = 0.5,
                    gp=gpar(fontsize=size))

  grid.arrange(label1,arrangeGrob(p1, p2, nrow = 2), nrow=1,
               widths=unit.c(unit(1.5, "lines"), unit(1, "npc") - unit(1.5, "lines")))
}

## Plots the absolute values of transmission tree, eps and xi inference accuracy
vis.abs.grid <- function(store,mod="tcg"){

    capitalise <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }

    dwide <- data.frame(matrix(nrow=length(store)*2*3,ncol=5,
                               dimnames=list(c(),c("Disease","Variable","Eps","Xi","Accuracy"))))
    counter <- 1

    for(run in names(store)){
        for(disease in c("ebola","sars")){
            for(variable in c("transmission","eps","xi")){
                if(disease=="sars") dwide$Disease[counter] <- toupper(disease)
                else dwide$Disease[counter] <- capitalise(disease)
                dwide$Variable[counter] <- capitalise(variable)
                dwide$Eps[counter] <- store[[run]][[disease]]$param$eps
                dwide$Xi[counter] <- store[[run]][[disease]]$param$xi
                if(variable %in% c("eps","xi")){
                    true <- store[[run]][[disease]]$param[[variable]]
                    est <- mean(subset(store[[run]][[disease]]$est,Model==mod & Variable==variable)$Value)
                    if(variable=="eps") dwide$Accuracy[counter] <- 1-abs(true-est)
                    if(variable=="xi") dwide$Accuracy[counter] <- 1-abs(true-est)/store[[run]][[disease]]$range
                }else if(variable=="transmission"){
                    dwide$Accuracy[counter] <- mean(subset(store[[run]][[disease]]$analysis,
                                                           variable=="accuracy" & model==mod)$value)
                }
                counter <- counter + 1
            }
        }
    }

    dwide$Variable <- factor(dwide$Variable,levels=c("Transmission","Eps","Xi"))
    levels(dwide$Variable) <- c("Transmission",expression(epsilon),expression(lambda))

    p <- ggplot(dwide,aes(x=Eps,y=Xi,fill=Accuracy)) +
        geom_tile(color="white") + facet_grid(Variable ~ Disease,labeller=label_parsed) +
        scale_fill_gradientn(colours=c("white","darkgreen"),limits=c(min(dwide$Accuracy),1)) +
        theme_set(theme_gray(base_size = 16)) +
        xlab(expression(paste("Reporting sensitivity (",epsilon,")"))) +
        ylab(expression(paste("Number of non-transmission contacts (",lambda,")"))) +
        scale_x_continuous(breaks=c(0.1,0.4,0.7,1),expand=c(0,0)) +
        scale_y_continuous(breaks=seq(0,10,2),expand=c(0,0)) +
        theme(axis.ticks=element_blank(),panel.border = element_rect(fill=NA,colour="black",size=1)) +
        guides(fill = guide_colorbar(title.vjust = 0,title.hjust=0.5,reverse=TRUE,title="Accuracy"))

    p
}

## Plot the paramters estimates for eps and xi
## We use use only only dat (i.e. tcg or tg) and factor across diseases
vis.est.dis <- function(store, dat = "tcg", rm = c(0,1)) {

    df <- NULL

    for(i in seq_along(store)) {
        for(disease in c("ebola","sars")) {
            for(mod in c("tc","tcg")) {
                for(var in c("eps","xi")){
                    temp <- subset(store[[i]][[disease]]$est, Model == mod & Variable == var)
                    temp$Eps <- store[[i]][[disease]]$param$eps
                    temp$Xi <- store[[i]][[disease]]$param$xi
                    if(disease == "ebola") temp$Disease <- "Ebola"
                    if(disease == "sars") temp$Disease <- "SARS"
                    df <- rbind(df, temp)
                }
            }
        }
    }

    df1 <- subset(df, Variable == "eps" & Model == dat & !Eps %in% rm & !Xi %in% rm)
    df2 <- subset(df, Variable == "xi" & Model == dat & !Eps %in% rm & !Xi %in% rm)

    df1$lab.Eps <- paste0("epsilon==",df1$Eps)
    df1$lab.Xi <- paste0("xi==",df1$Xi)

    df2$lab.Eps <- paste0("epsilon==",df2$Eps)
    df2$lab.Xi <- paste0("xi==",df2$Xi)

    create.histgrid <- function(df, var, main) {

        p <- ggplot(df, aes(Value)) + geom_density(aes(fill = Disease), alpha = 0.5) +
            scale_x_continuous(breaks = c(0.5)) +
            facet_grid(lab.Xi ~ lab.Eps, labeller = label_parsed, switch = 'both') +
            scale_fill_manual(values = match.col("disease"), labels = c("Ebola ", "SARS")) +
            ggtitle(main) + xlab("Value") + ylab("Density") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text.y      = element_blank(),
                  axis.ticks.y     = element_blank(),
                  legend.position  = "bottom")

        dat <- ggplot_build(p)
        yend <- max(dat$data[[1]]$density)
        p <- p + geom_segment(aes_string(x = var, y = 0, xend = var, yend = yend),
                              linetype = "dotted", size = 1)

        return(p)
    }

    cols <- wesanderson::wes_palette("Royal1", 4)

    p1 <- create.histgrid(df1, "Eps",
                          expression(paste("Estimating the reporting sensitivity (", epsilon, ")")))
    p2 <- create.histgrid(df2, "Xi",
                          expression(paste("Estimating the contact specificity (", xi, ")")))

    p <- joint.legend(p1, p2, 'bottom')

    return(p)
}

## A function to simulate contact tracing data (CTD) from a simOutbreak object
## NTP.simCTD takes the number of NTP as the argument for xi
NTP.simCTD <- function(outbreak,eps=1,xi=0,full=FALSE){

  if(outbreak$n==1) stop("No transmission observed")
  if(xi > (outbreak$n-2)) stop(paste0("Can't have this many contacts (",xi+1,
                                      ") per person, outbreak is too small (",outbreak$n," individuals)"))

  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
    
    accept.reject <- function(pair,xi,eps){
    if(pair[3]) return(stats::runif(1,0,1) < eps)
    else return(stats::runif(1,0,1) < eps*xi/(outbreak$n-2))
  }

  import <- which(is.na(outbreak$ances))

  infec.contact <- cbind(outbreak$ances[-import],outbreak$id[-import])

  potent.CTD <- as.data.frame(t(utils::combn(outbreak$id,2)))
  colnames(potent.CTD) = c("i","j")

  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject,xi=xi,eps=eps))

  CTD <- potent.CTD[potent.CTD$accept,1:2]
  rownames(CTD) <- NULL

  if(full) return(potent.CTD)

  return(CTD)
}

## A function to simulate contact tracing data (CTD) from a simOutbreak object
## clust.simCTD takes p(c(NTP)) as the argumennt for xi
clust.simCTD <- function(outbreak, eps = 1,xi = 0 ){
    
    if(outbreak$n==1) {
        warning("No transmission observed, returning NULL")
        return(NULL)
    }

    is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))

    accept.reject <- function(pair){
        if(pair[3]) return(stats::runif(1,0,1) < eps)
        else return(stats::runif(1,0,1) < xi*eps)
    }

    import <- which(is.na(outbreak$ances))

    infec.contact <- cbind(outbreak$ances[-import],outbreak$id[-import])

    potent.CTD <- as.data.frame(t(utils::combn(outbreak$id,2)))
    colnames(potent.CTD) = c("i","j")

    potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
    potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))

    CTD <- potent.CTD[potent.CTD$accept,1:2]
    rownames(CTD) <- NULL

    return(CTD)
}

## A function to analyse the outbreaker output for accuracy, entropy and confidence
result.analysis <- function(result, true.outbreak, print = FALSE, kappa = FALSE){
  id <- seq_len(true.outbreak$n)
  adder <- which(names(result)=="alpha_1")-1
  samples <- length(result$step)
  
  #Determine the modal transmission network
  network <- data.frame(from=do.call(rbind,lapply(id, function(i) ({
    modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
    if(length(modal.ances)==0) return(NA) else return(modal.ances)
  }))), to=id)

  if(kappa){

      modal.kappa <- function(i){
      index <- i + kappa.adder
      mode <- as.integer(names(which.max(table(result[[index]]))))
      if(length(mode)==0) return(NA)
      else return(mode)
  }

      kappa.adder <- which(names(result)=="kappa.1")-1
      kappa.result <- sapply(id,modal.kappa)
      kappa.acc <- sum(kappa.result==true.outbreak$kappa,na.rm=TRUE)/sum(!is.na(true.outbreak$kappa))
  } else kappa.acc <- NULL

  import <- which(is.na(network$from))

  #Define the indices to call the times of infection
  t.inf.index <- which(names(result)=="t.inf.1")-1+id

  #Determine the median posterior time of onset for plotting
  onset <- unlist(lapply(result[t.inf.index],median))

  #Scale onset to begin at 0
  onset <- onset - min(onset)

  #Determine confidence in our results

  transmission.id <- id[!sapply(result[id+adder],function(i) any(is.na(i)))]

  entropy <- round(mean(sapply(result[transmission.id+adder],
                               function(i) {fk <- table(i)/sum(table(i))
  -sum(log(fk)*fk)})),2)

  #Determine the proportion of correctly inferred ancestries
  num.correct <-  sum(true.outbreak$ances==network$from,na.rm=TRUE)
  num.correct <- num.correct + sum(is.na(true.outbreak$ances[is.na(network$from)]))
  accuracy <- round(num.correct/nrow(network),2)

  out <- list(analysis=data.frame(accuracy=accuracy,entropy=entropy))
  out$analysis$kappa.acc <- kappa.acc
  return(out)
}

## Calculate NTP in four ways, across available data
est.NTP <- function(dir, plot = FALSE) {

    # Calculate the analytical NTP from true parameter values
    f1 <- function(r) {
        act.eps <- r$ebola$param$eps
        act.xi <- r$ebola$param$xi
        n <- r$ebola$param$min.cases
        I <- 1
        analytical.NTP(act.eps, act.xi, n, I)
    }

    # Calculate the analytical NTP from estimated parameter values
    f2 <- function(i, r) {
        est.eps <- mean(r$ebola$tcg.result[[i]]$eps)
        est.xi <- mean(r$ebola$tcg.result[[i]]$xi)
        n <- r$ebola$param$min.cases
        I <- 1
        analytical.NTP(est.eps, est.xi, n, I)
    }

    # Count the NTP from the actual tTree and CTD
    f3 <- function(i, r) {
        result <- r$ebola$outbreak[[i]]
        tTree <- data.frame(i = result$id, j = result$ances)
        CTD <- r$ebola$CTD[[i]]
        get.NTP(tTree, CTD)
    }

    # Count the NTP from the inferred tTree and CTD
    f4 <- function(i, r) {
        result <- r$ebola$tcg.result[[i]]
        tTree <- get.tTree(result)
        CTD <- r$ebola$CTD[[i]]
        get.NTP(tTree, CTD)
    }

    files <- paste0(dir, list.files(dir)[-grep("store", list.files(dir))])

    out <- data.frame(matrix(nrow = length(files), ncol = 4))
    names(out) <- c("an.true", "an.est", "count.true", "count.est")

    pb <- txtProgressBar(min = 1, max = length(files), style = 3)

    for(i in seq_along(files)) {

        setTxtProgressBar(pb, i)

        load(files[i])
        runs <- 1:r$ebola$param$runs

        out$an.true[i] <- f1(r)
        out$an.est[i] <- mean(sapply(runs, f2, r))

        if(length(r$ebola$CTD) == 0) {
            out$count.true[i] <- 0
            out$count.est[i] <- 0
        } else {
            out$count.true[i] <- mean(unlist(sapply(runs, f3, r)))
            out$count.est[i] <- mean(unlist(sapply(runs, f4, r)))
        }
    }

    out$run <- factor(1:nrow(out))

    mlt <- reshape2::melt(out, id = 'run')

    p <- ggplot(mlt, aes(x = run, y = value, fill = variable)) +
        geom_bar(stat = 'identity', position = 'dodge')

    if(plot) return(p)
    else return(out)

}
