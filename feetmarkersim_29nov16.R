

feet_marker_sim <- function(
    tmax=1000,
    d=0.1, #coordination benefit
    g=0.1, #extra coordination benefit among mutualists
    h=0.1, #mis-coordination cost for mutualists
    a=0.9, #probability of assorting on marker
    m0=0.01, #proportion of population considering migration
    mu=0.5, #proportion of migrants that are payoff-biased
    s=c(0.5, 0.5), #proportion of total population in each group
    init_p=c(0.55,0.45), #intial proportion of behavior 0 in each group
    init_q=c(0.51,0.49), #intial proportion of marker 0 in each group 
    draw=TRUE ) {

    # define and check group number
    ngroups=length(s)
    if ( length(init_p)!=ngroups | length(init_q)!=ngroups ) {
      stop("!! s, p and q vectors must be equal length !!")
    }
    if (sum(s)!=1){
      stop("!! s must sum to 1 !!")
    }

    # init matrix of states
    # x[[k]][i,j] gives frequency of xijk
    xblank <- matrix(0,nrow=2,ncol=2)
    x <- list(xblank)
    for ( i in 1:ngroups ) {
        x[[i]] <- xblank
    }

    # init vector of group-mean fitnesses
    fit <- rep(0, times=ngroups)
    fitp1 <- rep(0, times=ngroups)
    fitq1 <- rep(0, times=ngroups)
    fitq0 <- rep(0, times=ngroups)
    
    # fitness function
    # w0jk = 1-h+(d+g+h)((1-a)*p0k+a*x0jk/qjk)
    # w1jk = 1+d*((1-a)*p1k+a*x1jk/qjk)
    wijk <- function( i , j , k ) {

        i <- i+1 #turn behav and marker indices into matrix indices
        j <- j+1
        
        if (i==1) {
            w <- 1 - h + (d + g + h)*( (1 - a)*(x[[k]][1,1] + x[[k]][1,2]) +
                                        a*x[[k]][1,j]/(x[[k]][1,j] + x[[k]][2,j]) )
        } else {
            w <- 1 + d*( (1 - a)*(x[[k]][2,1] + x[[k]][2,2]) +
                            a*x[[k]][2,j]/(x[[k]][1,j] + x[[k]][2,j]) )
        } #else

        w
    }

    
    # payoff-biased migration function 1 outgroup
    mij <- function( i , j ) {
        #mij <- m0( 1 + mu*(fit[j] - fit[i]) ) #original R&B payoff-biased migration
        mij <- m0*(1-mu)/2 + m0*mu*(fit[j]/(fit[i]+fit[j]))   
        mij 
    }
    
    # payoff-biased migration function when compare to 2 neighboring groups
        m_kto <- function( k , t , o ) { #migrants from ingroup k, to target t, given the opposite side group o
        
        mkt <- m0*(1-mu)/3 + m0*mu*(fit[t]/(fit[k]+fit[t]+fit[o]))   #those that move to target group
        mkt
    }


    # group index bounding function
    groupindex <- function( k ) {
        k <- ifelse( k > ngroups , 1 , k )
        k <- ifelse( k < 1 , ngroups , k )
        k
    }
    
    #initial frequencies with no linkage
    for ( k in 1:ngroups ) {
        if ( missing(init_p) || missing(init_q) ) {
            p <- runif(1)  
            q <- runif(1)
            x[[k]][1,1] <- p*q
            x[[k]][1,2] <- p*(1-q)
            x[[k]][2,1] <- (1-p)*q
            x[[k]][2,2] <- (1-p)*(1-q)
        } else {
            p <- init_p[k]  
            q <- init_q[k]
            x[[k]][1,1] <- p*q
            x[[k]][1,2] <- p*(1-q)
            x[[k]][2,1] <- (1-p)*q
            x[[k]][2,2] <- (1-p)*(1-q) 
        } #else
    }


    #initialize matrices to hold p, q, D, and s for all groups
    p_mat <- NULL     #behavior freq
    q_mat <- NULL     #marker freq
    D_mat <- NULL     #linkage
    s_mat <- s
    fit_mat <- fit
    for ( k in 1:ngroups ) {
        p_mat <- cbind(p_mat, x[[k]][1,1] + x[[k]][1,2])
        q_mat <- cbind(q_mat, x[[k]][1,1] + x[[k]][2,1])
        D_mat <- cbind(D_mat, x[[k]][1,1]*x[[k]][2,2] - x[[k]][1,2]*x[[k]][2,1])
    }
    
    print("START")
    print(x)
    
    for ( t in 1:tmax ) {

        xnew <- x 

        # learning
        for( k in 1:ngroups ) {
            w11k <- wijk( 1 , 1 , k )
            w10k <- wijk( 1 , 0 , k )
            w01k <- wijk( 0 , 1 , k )
            w00k <- wijk( 0 , 0 , k )
            wbar <- x[[k]][1,1]*w00k + x[[k]][1,2]*w01k + x[[k]][2,1]*w10k + x[[k]][2,2]*w11k
            xnew[[k]][1,1] <- x[[k]][1,1] * w00k / wbar
            xnew[[k]][1,2] <- x[[k]][1,2] * w01k / wbar
            xnew[[k]][2,1] <- x[[k]][2,1] * w10k / wbar
            xnew[[k]][2,2] <- x[[k]][2,2] * w11k / wbar
            
            fitp1[k] <- (w10k+w11k)/ (w10k+w11k+w01k+w00k)
            fitq1[k] <- (w01k+w11k)/ (w10k+w11k+w01k+w00k)
            fitq0[k] <- (w00k+w10k)/ (w10k+w11k+w01k+w00k)
            fit[k] <- wbar #save group-mean fitness for later
        } #k

        xnewnew <- xnew
        snew <- s

        # migrate
        for( k in 1:ngroups ) {
            n1 <- groupindex(k+1)   #index for right side group
            n2 <- groupindex(k-1)   #index for left side group
            m_k_n1 <- m_kto( k , n1, n2) #from k (ingroup) to n1 (target)
            m_n1_k <- m_kto( n1, k, n2 )
            m_k_n2 <- m_kto( k , n2, n1 ) 
            m_n2_k <- m_kto( n2, k , n1)

            for ( i in 1:2 ) {
                for ( j in 1:2 ) {

                    #for now assume only two groups
                    xnewnew[[k]][i,j] <- ( xnew[[k]][i,j]*s[k]*( 1 - m_k_n1 - m_k_n2) +   # those staying in k
                                            xnew[[n1]][i,j]*s[n1]*m_n1_k +                # those coming from n1
                                            xnew[[n2]][i,j]*s[n2]*m_n2_k)/                # those coming from n2
                                            ( s[k]*( 1 - m_k_n1 ) + s[n1]*m_n1_k  + s[n2]*m_n2_k )  #total new population in group k
                                            
                } #j
            } #i

            #recursion for group size    
            snew[k] <- s[k]*(1 - m_k_n1 - m_k_n2) + s[n1]*m_n1_k + s[n2]*m_n2_k

        } #k

        x <- xnewnew
        s <- snew
        
        
        # draw
        if ( draw==TRUE ) {
            p <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][1,2] )
            q <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][2,1] )
            linkage <- sapply( 1:ngroups , function(z) x[[z]][1,1]*x[[z]][2,2] - x[[z]][1,2]*x[[z]][2,1] )
            linkage <- linkage + 0.5
            
            #update results matrices
            # p_mat <- rbind(p_mat, p)
            # q_mat <- rbind(q_mat, q)
            # D_mat <- rbind(D_mat, linkage-0.5)
            s_mat <- rbind(s_mat, s)
            # fit_mat <-rbind(fit_mat, fit/sum(fit))
            
            if ( t==1 ) {
                plot( 1:ngroups , p , ylim=c(0,1) , lwd=2 )
            } else {
                rect( 0,-.02,ngroups+.1,1.1 , col="white" , border="white" )
                points( 1:ngroups , p , lwd=2 )
            }
            lines( 1:ngroups , q , col="red" )          #frequency of marker
            lines( 1:ngroups , linkage , col="blue" )   #linkage btwn bhvr & marker per grp
            lines( 1:ngroups , s , col="green")         #proportion of population in each group
            lines( 1:ngroups , fit/sum(fit) , col="cyan")         #prop of total fitness of each group
            points( 1:ngroups , fit/sum(fitp1) , col="cyan")         #prop of group fitness of behavior 1
            lines( 1:ngroups , fit/sum(fitq1) , col="lightblue", lty=2)         #prop of group fitness of marker 1
            lines( 1:ngroups , fit/sum(fitq0) , col="lightblue", lty=3)         #prop of group fitness of marker 1
            legend (ngroups*.6, 1, legend=c("prop. bhvr", "prop. marker", "linkage", "prop in grp", "fitness", "payoff p1", "payoff q1"), 
                    col=c("black", "red", "blue", "green", "cyan", "cyan", "lightblue"), pch=c(1,19,19,19,19,1,1))
        }
        
    } # t
    print("END")
    invisible(x)

    print(s[1:ngroups])
    print(sum(s[1:ngroups])) #for testing = 1
    print(fit)
    
    # plot(0:tmax, p_mat[,1], ylim=c(0,1) , lwd=2, col="black" )
    # points(0:tmax, p_mat[,2], ylim=c(0,1) , lwd=2, col="blue")
    # lines(0:tmax, q_mat[,1], ylim=c(0,1) , lwd=2, col="black")
    # lines(0:tmax, q_mat[,2], ylim=c(0,1) , lwd=2, col="blue")
    # lines(0:tmax, D_mat[,1], ylim=c(0,1) , lwd=2, col="black", lty=3)
    # lines(0:tmax, D_mat[,2], ylim=c(0,1) , lwd=2, col="blue", lty=3)
    # lines(0:tmax, s_mat[,1], ylim=c(0,1) , lwd=3, col="black", lty=4)
    # lines(0:tmax, s_mat[,2], ylim=c(0,1) , lwd=3, col="blue", lty=2)
    # lines(0:tmax, fit_mat[,1], ylim=c(0,1) , lwd=3, col="darkgrey", lty=3)
    # lines(0:tmax, fit_mat[,2], ylim=c(0,1) , lwd=3, col="lightblue", lty=3)
    # 
    print(x)

} #function

feet_marker_sim(
    tmax=100,
    d=0.5, #coordination benefit
    g=2, #extra coordination benefit among mutualists
    h=0, #mis-coordination cost for mutualists
    a=0, #probability of assorting on marker
    m0=0.03, #proportion of each group that migrates
    mu=0, #of all migrants, proportion that engage in payoff biased migration
    s=rep(.1,10), #proportion of total population in each group
    init_p=c(rep(.9,4),rep(.1,6)), #intial proportion of behavior 0 in each group
    init_q=c(rep(.6,4),rep(.4,6)), #intial proportion of marker 0 in each group 
    draw=TRUE )


