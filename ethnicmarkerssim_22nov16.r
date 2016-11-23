

ethnicsim <- function( tmax=1000 ,
                        ngroups=2 ,
                        delta=0.1,
                        e=0.1,
                        m=0.01,
                        init_p=c(0.55,0.45), #intial proportion of behavior 0 in each group
                        init_q=c(0.51,0.49), #intial proportion of marker 0 in each group
                        draw=TRUE ) {
    
    # init matrix of states
    # x[[k]][i,j] gives frequency of xijk
    xblank <- matrix(0,nrow=2,ncol=2)
    x <- list(xblank)
    for ( i in 1:ngroups ) {
        x[[i]] <- xblank
    }
    
    # fitness function
    # wijk = 1 + d*( e*pik + (1-e)*xijk/qjk )
    wijk <- function( i , j , k ) {
        w <- 1 + delta*( e*( x[[k]][i,1] + x[[k]][i,2] ) + (1-e)*x[[k]][i,j]/( x[[k]][1,j] + x[[k]][2,j] ) )
        w
    }
    
    # group index bounding function
    groupindex <- function( k ) {
        k <- ifelse( k > ngroups , 1 , k )
        k <- ifelse( k < 1 , ngroups , k )
        k
    }
    
    # random initial frequencies with no linkage
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
    p_mat <- NULL
    q_mat <- NULL
    D_mat <- NULL
    for ( k in 1:ngroups ) {
        p_mat <- cbind(p_mat, x[[k]][1,1] + x[[k]][1,2])
        q_mat <- cbind(q_mat, x[[k]][1,1] + x[[k]][2,1])
        D_mat <- cbind(D_mat, x[[k]][1,1]*x[[k]][2,2] - x[[k]][1,2]*x[[k]][2,1])
    }
    
    print("START")
    print(x)
    
    for ( t in 1:tmax ) {
        # migrate
        xnew <- x
        for( k in 1:ngroups ) {
            n1 <- groupindex(k+1)
            n2 <- groupindex(k-1)
            for ( i in 1:2 ) {
                for ( j in 1:2 ) {
                    xnew[[k]][i,j] <- (1-m)*x[[k]][i,j] + m/2 * x[[n1]][i,j] + m/2 * x[[n2]][i,j]
                } #j
            } #i
        } #k
        x <- xnew
        
        # learning
        for( k in 1:ngroups ) {
            w11k <- wijk( 1 , 1 , k )
            w10k <- wijk( 1 , 2 , k )
            w01k <- wijk( 2 , 1 , k )
            w00k <- wijk( 2 , 2 , k )
            wbar <- x[[k]][1,1]*w11k + x[[k]][1,2]*w10k + x[[k]][2,1]*w01k + x[[k]][2,2]*w00k
            xnew[[k]][1,1] <- x[[k]][1,1] * w11k / wbar
            xnew[[k]][1,2] <- x[[k]][1,2] * w10k / wbar
            xnew[[k]][2,1] <- x[[k]][2,1] * w01k / wbar
            xnew[[k]][2,2] <- x[[k]][2,2] * w00k / wbar
        } #k
        x <- xnew
        
        # draw
        if ( draw==TRUE ) {
            p <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][1,2] )
            q <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][2,1] )
            linkage <- sapply( 1:ngroups , function(z) x[[z]][1,1]*x[[z]][2,2] - x[[z]][1,2]*x[[z]][2,1] )
            linkage <- linkage + 0.5


            #update results matrices
            p_mat <- rbind(p_mat, p)
            q_mat <- rbind(q_mat, q)
            D_mat <- rbind(D_mat, linkage-0.5)
            

            if ( t==1 ) {
                plot( 1:ngroups , p , ylim=c(0,1) , lwd=2 )
            } else {
                rect( 0,-.02,ngroups+.1,1.1 , col="white" , border="white" )
                points( 1:ngroups , p , lwd=2 )
            }
            lines( 1:ngroups , q , col="red" )
            lines( 1:ngroups , linkage , col="blue" )
        }
        
    } # t
    print("END")
    invisible(x)

        print(s[1:ngroups])
    print(fit)

    plot(0:tmax, p_mat[,1], ylim=c(0,1) , lwd=2, col="black" )
    points(0:tmax, p_mat[,2], ylim=c(0,1) , lwd=2, col="blue")
    lines(0:tmax, q_mat[,1], ylim=c(0,1) , lwd=2, col="black")
    lines(0:tmax, q_mat[,2], ylim=c(0,1) , lwd=2, col="blue")
    lines(0:tmax, D_mat[,1], ylim=c(0,1) , lwd=2, col="black", lty=3)
    lines(0:tmax, D_mat[,2], ylim=c(0,1) , lwd=2, col="blue", lty=3)
    
}

ethnicsim( ngroups=2,
            tmax=500,
            delta=0.5,
            e=0.1,
            m=0.01,
            init_p=c(0.55,0.45), #intial proportion of behavior 0 in each group
            init_q=c(0.51,0.49), #intial proportion of marker 0 in each group
            draw=TRUE )










