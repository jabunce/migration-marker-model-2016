

feet_marker_sim <- function( tmax=1000 , ngroups=2 ,
    d=0.1, #coordination benefit
    g=0.1, #extra coordination benefit among mutualists
    h=0.1, #mis-coordination cost for mutualists
    a=0.9, #probability of assorting on marker
    m0=0.01, #base proportion of each group that migrates
    mu=0.5, #mean-payoff bias in migration decisions
    s=c(0.5, 0.5), #proportion of total population in each group 
    draw=TRUE ) {


    # init matrix of states
    # x[[k]][i,j] gives frequency of xijk
    xblank <- matrix(0,nrow=2,ncol=2)
    x <- list(xblank)
    for ( i in 1:ngroups ) {
        x[[i]] <- xblank
    }

    # init vector of group-mean fitnesses
    fit <- rep(0, times=ngroups)

    
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

    
    # payoff-biased migration function
    mij <- function( i , j ) {
        mij <- m0*( 1 + mu*(fit[j] - fit[i]) )
        mij
    }


    # group index bounding function
    groupindex <- function( k ) {
        k <- ifelse( k > ngroups , 1 , k )
        k <- ifelse( k < 1 , ngroups , k )
        k
    }
    
    #initial frequencies with no linkage
    for ( k in 1:ngroups ) {
        p <- k*0.3 #runif(1)  
        q <- k*0.3 #runif(1)
        x[[k]][1,1] <- p*q
        x[[k]][1,2] <- p*(1-q)
        x[[k]][2,1] <- (1-p)*q
        x[[k]][2,2] <- (1-p)*(1-q)
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

            fit[k] <- wbar #save group-mean fitness for later
        } #k

        xnewnew <- xnew
        snew <- s

        # migrate
        for( k in 1:ngroups ) {
            n1 <- groupindex(k+1)
            n2 <- groupindex(k-1)
            m_k_n1 <- mij( k , n1 )
            m_n1_k <- mij( n1 , k )
            #m_k_n2 <- mij( k , n2 ) #for now assume only two groups
            #m_n2_k <- mij( n2 , k )

            for ( i in 1:2 ) {
                for ( j in 1:2 ) {

                    #for now assume only two groups
                    xnewnew[[k]][i,j] <- ( xnew[[k]][i,j]*s[k]*( 1 - m_k_n1 ) +
                                            xnew[[n1]][i,j]*s[n1]*m_n1_k )/
                                            ( s[k]*( 1 - m_k_n1 ) + s[n1]*m_n1_k )
                } #j
            } #i

            #recursion for group size    
            snew[k] <- s[k]*(1 - m_k_n1) + s[n1]*m_n1_k

        } #k

        x <- xnewnew
        s <- snew
        
        
        # draw
        if ( draw==TRUE ) {
            p <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][1,2] )
            q <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][2,1] )
            linkage <- sapply( 1:ngroups , function(z) x[[z]][1,1]*x[[z]][2,2] - x[[z]][1,2]*x[[z]][2,1] )
            linkage <- linkage + 0.5
            
            if ( t==1 ) {
                plot( 1:ngroups , p , ylim=c(0,1) , lwd=2 )
            } else {
                rect( 0,-.02,ngroups+.1,1.1 , col="white" , border="white" )
                points( 1:ngroups , p , lwd=2 )
            }
            lines( 1:ngroups , q , col="red" )
            lines( 1:ngroups , linkage , col="blue" )
            lines( 1:ngroups , s , col="green") #proportion of population in each group
        }
        
    } # t
    print("END")
    invisible(x)

    print(s[1:ngroups])
    print(fit)

    
} #function

feet_marker_sim( ngroups=2, tmax=500,
    d=0.5, #coordination benefit
    g=0, #extra coordination benefit among mutualists
    h=0, #mis-coordination cost for mutualists
    a=0.9, #probability of assorting on marker
    m0=0.01, #base proportion of each group that migrates
    mu=0, #mean-payoff bias in migration decisions
    s=c(0.5, 0.5), #proportion of total population in each group 
    draw=TRUE )


