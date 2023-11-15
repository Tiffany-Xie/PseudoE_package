# R Packages Try

erlang <- function(x, n, γ) {
        (n*γ)^n*x^(n-1)*exp(-n*γ*x)/factorial(n-1)
}

SInR_geom <- function(t, states, params) {
        with(as.list(c(params)), {
                I <- states[1:n]
                R <- states[[n+1]]

                Iprev <- c(0, I[1:(n-1)])
                dI <- a*r^(c(0,0:(n-2)))*Iprev - (a*r^(0:(n-1))+μ)*I
                dR <- a*r^(n-1)*I[[n]] - μ*R
                return(list(c(dI, dR)))
        })
}

Edens <- function(time, γ, nE) {
        df <- data.frame(Time = time)
        df$PE <- erlang(time, nE, γ)
        return(df)
}

r2kappa <- function(r, n, offset=0){
        delta <- (1:n) - (n+1)/2
        res <- exp(delta*log(r))
        kappa <- sum(res^2)/(sum(res)^2)
        return(kappa - offset)
}

kappa2r <- function(kappa, n){
        rmax <- 2*(1+kappa)/(1-kappa) ## BAD CODE, XNR please fix
        if(kappa>=1) return(NA)
        if(kappa<1/n) return(NA)
        u <- uniroot(r2kappa, interval=c(1, rmax), n=n, offset=kappa)
        return(u$root)
}

r2a <- function(r, n, mean) {
        return((1-1/r^n)/(mean*(1-1/r)))
}

PEdens <- function(time, r, a, nPE, model) {
        df <- expand.grid(Time = time, a = a, r = r, nPE = nPE)
        states <- c(1, numeric(nPE))
        names(states) <- c(paste0("I", 1:nPE), "R")
        params <- c(n = nPE, a = a, r = r)
        soln <- ode(y = states,
                    times = time,
                    func = model,
                    parms = params)
        df$PPE <- c(diff(soln[,"R"])/diff(time), NA)
        return(na.omit(df))
}

parComp <- function(df, mean, kappa, a, r, nPE) {
        numM <- sum(df$Time * df$PPE * ts)
        numV <- sum(df$Time^2 * df$PPE * ts) - numM^2
        numK <- numV/numM^2
        fM <- 1/a *(1-1/r^nPE)/(1-1/r)
        fV <- 1/a^2 *(1-1/r^(2*nPE))/(1-1/r^2)
        fK <- fV/fM^2
        df <- data.frame(TheoMean = mean,
                         FormulaMean = fM,
                         ActualMean = numM,
                         TheoKappa = kappa,
                         FormulaKappa = fK,
                         ActualKappa = numK)
        return(df)
}

Cfplot <- function(dfE, dfPE, mean) {
        ggplot(dfE, aes(x=Time, y=PE, color = "Erlang")) + geom_line(linewidth=1.5) +
                geom_vline(xintercept = mean, color = "red", linetype = "dashed", linewidth = 1) +
                geom_line(data=dfPE, aes(x=Time, y=PPE, color="ODE")) +
                labs(title = paste0("Erlang vs. Pseudo-Erlang: Same Mean and Kappa (nE=", dfE$nE[1], ", nPE=", dfPE$nPE[1], ", Mean=", mean, ")"))
}


