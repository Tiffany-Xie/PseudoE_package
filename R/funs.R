##Time
timeSeq <- function(ts, T, mid=TRUE) {
        if (mid) return(seq(ts/2, T, by=ts))
        return(seq(0, T, by=ts))
}

midpoints <- function(ts){
        return((ts[-1]+ts[-length(ts)])/2)
}

######################################################################

erlang <- function(x, n, γ) {
        (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

InR_geom <- function(t, states, params) {
        with(as.list(c(params)), {
                I <- states[1:n]
                R <- states[[n+1]]

                Iprev <- c(0, I[1:(n-1)])
                dI <- a*r^(c(0,0:(n-2)))*Iprev - (a*r^(0:(n-1))+μ)*I
                dR <- a*r^(n-1)*I[[n]] - μ*R
                return(list(c(dI, dR)))
        })
}


SIR <- function(time, states, params) {
        with(params, {
                stopifnot(length(outrate)==n)
                S <- states[[1]]
                I <- states[2:(n+1)]
                R <- states[[n+2]]

                inc <- β*S*sum(I)/N

                outflow <- outrate*I
                inrate <- outrate-μ
                inflow <- c(inc, (I*inrate)[0:(n-1)])

                dS <- μ*N - inc - μ*S
                dI <- inflow - outflow
                dR <- inrate[n]*I[[n]] - μ*R
                dC <- inc

                return(list(c(dS, dI, dR, dC)))
        })
}

######################################################################

SInRFlow <- function(β, D, n, μ, S0, I0, ts, T) {
        gamma <- 1/D
        outrate <- rep(n*gamma + μ, times=n)

        params <- list(β=β, n=n, μ=μ, N = S0+I0, outrate=outrate)
        states <- c(S0, I0, numeric(n-1), 0, 0)
        names(states) <- c("S", paste0("I", 1:n), "R", "inc")
        return(Integration2(params, states, ts, T, SIR))
}

sinnerFlow <- function(β, D, kappa, n, μ, S0, I0, ts, T) {
        r <- kappa2r(kappa, n)
        a <- (1-1/r^n)/(D*(1-1/r))
        #print(c(r, a))
        outrate <- a*r^(0:(n-1)) + μ

        params <- list(β=β, n=n, μ=μ, N=S0+I0, outrate=outrate)
        states <- c(S0, I0, numeric(n-1), 0, 0)
        names(states) <- c("S", paste0("I", 1:n), "R", "inc")
        return(Integration2(params, states, ts, T, SIR))
}

######################################################################

r2kappa <- function(r, n, offset=0){
        delta <- (1:n) - (n+1)/2
        res <- exp(delta*log(r))
        kappa <- sum(res^2)/(sum(res)^2)
        return(kappa - offset)
}

kappa2r <- function(kappa, n){
        rmax <- 2*(1+kappa)/(1-kappa)
        if(kappa>=1) return(NA)
        if(kappa<1/n) return(NA)
        u <- uniroot(r2kappa, interval=c(1, rmax), n=n, offset=kappa)
        return(u$root)
}

######################################################################

parGenerator <- function(n, mu, kappa, μ=0) {
        r <- kappa2r(kappa, n)
        a <- (1-1/r^n)/(mu*(1-1/r))
        return(c(n = n, r = r, a = a, μ = μ))
}

parCheck <- function(flow, ts, T) {
        flow <- diff(flow)
        time <- timeSeq(ts, T)
        tot <- sum(flow)
        mu <- sum(flow*time)
        S <- sum(flow*time^2)
        kappa <- S/mu^2 - 1
        return(c(tot = tot, mu = mu, kappa = kappa))
}

######################################################################

Integration <- function(n, mu, kappa, ts, T, model=InR_geom) {
        time <- timeSeq(ts, T, mid=FALSE)
        params <- parGenerator(n, mu, kappa)
        states <- c(1, numeric(n))
        names(states) <- c(paste0("I", 1:n), "R")
        soln <- ode(y = states,
                    times = time,
                    func = model,
                    parms = params)
        return(soln[,"R"])
}

Integration2 <- function(params, states, ts, T, model) {
        time <- timeSeq(ts, T, FALSE)
        soln <- ode(y = states,
                    times = time,
                    func = model,
                    parms = params)
        return(soln)
}
######################################################################

compPlot <- function(gamm, ode, ts, T) {
        df <- data.frame(Time = timeSeq(ts, T, FALSE), Gamma=gamm, ODE = ode)
        ggplot(df, aes(x=Time)) + geom_line(aes(y=Gamma, color = "Geometric")) +
                geom_line(aes(y=ODE, color = "ODE")) +
                labs(title = "ODE & Gamma (CDF)", y = "Cumulative Density")
}

compPlotDens <- function(gamm, ode, ts, T) {
        df <- data.frame(Time = timeSeq(ts, T), Gamma=diff(gamm), ODE = diff(ode))
        ggplot(df, aes(x=Time)) + geom_line(aes(y=Gamma, color = "Geometric")) +
                geom_line(aes(y=ODE, color = "ODE")) +
                labs(title = "ODE & Gamma (PDF*ts)", y = "Probability Density (*ts)")
}


######################################################################

gammaFlowDens <- function (mu, kappa, ts, T){
        b <- boundaries(ts, T)
        cum <- dgamma(b, 1/kappa, 1/(mu*kappa))
        return(ts*cum)
}

######################################################################
