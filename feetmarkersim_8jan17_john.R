rm(list = ls())

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
    draw=TRUE,
    plots=TRUE ) {

 
    # define and check group number
    init_s=s
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
    s_mat <- s        #proportion of total pop in each group

    for ( k in 1:ngroups ) {
        p_mat <- cbind(p_mat, x[[k]][1,1] + x[[k]][1,2])
        q_mat <- cbind(q_mat, x[[k]][1,1] + x[[k]][2,1])
        D_mat <- cbind(D_mat, x[[k]][1,1]*x[[k]][2,2] - x[[k]][1,2]*x[[k]][2,1])
    }

    # init vector of group-mean fitnesses
    fit <- rep(0, times=ngroups)
    fitp0 <- rep(0, times=ngroups)
    fitq0 <- rep(0, times=ngroups)
    fit_mat <- fit    #mean fitness in each group



    # fitness function
    # w0jk = 1-h+(d+g+h)((1-a)*p0k+a*x0jk/qjk)
    # w1jk = 1+d*((1-a)*p1k+a*x1jk/qjk)
    wijk <- function( i , j , k ) {

        i <- i+1 #turn behav and marker indices into matrix indices
        j <- j+1
        
        if (i==1) {
            w <- 1 - h + (d + g + h)*( (1 - a)*(x[[k]][1,1] + x[[k]][1,2]) +            #stag payoff (i=0, index=1)
                                        a*x[[k]][1,j]/(x[[k]][1,j] + x[[k]][2,j]) )
        } else {
            w <- 1 + d*( (1 - a)*(x[[k]][2,1] + x[[k]][2,2]) +                          #hare payoff (i=1, index=2)
                            a*x[[k]][2,j]/(x[[k]][1,j] + x[[k]][2,j]) )
        } #else

        w
    }

    
    # payoff-biased migration function with one outgroup
    mij <- function( i , j ) {
        
        #mij <- m0( 1 + mu*(fit[j] - fit[i]) ) #original R&B payoff-biased migration
        
        #new payoff-biased migration:

        mij <- m0*(1-mu)/2 + m0*mu*(fit[j]/(fit[i]+fit[j]))

        #Expanded, this is mij <- m0*( 1/2 - mu/2 + mu*(fit[j]/(fit[i]+fit[j]) )
        #If fit[i] = fit[j], then mij equals m0/2
        #There are only two groups, so this migration happens twice, and the total
        #migration from i to j is 2*m0/2 = m0.

        #If fit[j] >> fit[i], the most migration from i to j we'll get is 2*m0*(1/2 + mu/2) = m0*(1 + mu)
        #and if fit[j] << fit[i], the least migration we'll get is m0*(1 - mu)
        #which means the effect of mu will be less than in the original R&B model.
        #As in the original model, still have to keep m0 and mu small enough such that
        #migration rates are between 0 and 1.

        mij 
    }
    
    # payoff-biased migration function with two neighboring groups, one on each side
    m_kto <- function( k , t , o ) { #migrants from ingroup k, to target t, given the opposite side group o
        
        #mkt <- m0*(1-mu)/3 + m0*mu*(fit[t]/(fit[k]+fit[t]+fit[o]))   #those that move to target group
        #Expanded, this is mkt <- m0*( 1/3 - mu/3 + mu*(fit[t]/(fit[k]+fit[t]+fit[o])) )
        #If fit[t] = fit[k] = fit[o], then mij equals m0/3
        #Problem: This out-migration happens twice (one for each neighbor). If all three
        #groups have equal fitness, then migration out of group k is only 2/3*m0, but it should be m0.

        #Maybe instead use:

        mkt <- m0*( 1/2 - mu/3 + mu*(fit[t]/(fit[k]+fit[t]+fit[o])) )

        #Now, if all three groups have equal fitnessm mkt = m0/2, which means the
        #total out-migration to the two neighbors is m0.

        #If fit[t] >> fit[k] = fit[o], the most migration from k to t we'll get is m0*(1/2 + 2/3*mu), and
        #from k to o is m0*(1/2 - mu/3). Summing these, total out-migration from k is
        #m0*(1/2 - mu/3) + m0*(1/2 + 2/3*mu) = m0*(1/2 - mu/3 + 1/2 + 2/3*mu) = m0*(1 + mu/3)

        #If fit[t] = fit[o] >> fit[k], the most migration from k to t we'll get is m0*(1/2 + (1/2-1/3)*mu), and
        #from k to o is m0*(1/2 + (1/2-1/3)*mu). Summing these, total out-migration from k is
        #2*m0*(1/2 + (1/2-1/3)*mu) = m0*(1 + (1-2/3)*mu) = m0*(1 + mu/3)

        #If fit[k] >> fit[o] = fit[t], the most migration from k to t we'll get is m0*(1/2 - mu/3), and
        #from k to o is m0*(1/2 - mu/3). Summing these, total out-migration from k is
        #2*m0*(1/2 - mu/3) = m0*(1 - 2/3*mu)

        mkt
    }


    # group index bounding function
    groupindex <- function( k ) {
        k <- ifelse( k > ngroups , 1 , k )
        k <- ifelse( k < 1 , ngroups , k )
        k
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
            
            fitp0[k] <- (w00k+w01k)/ (w10k+w11k+w01k+w00k)
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

                    xnewnew[[k]][i,j] <- ( xnew[[k]][i,j]*s[k]*( 1 - m_k_n1 - m_k_n2 ) +   # those staying in k
                                            xnew[[n1]][i,j]*s[n1]*m_n1_k +                # those coming from n1
                                            xnew[[n2]][i,j]*s[n2]*m_n2_k )/                # those coming from n2
                                            ( s[k]*(1 - m_k_n1 - m_k_n2) + s[n1]*m_n1_k  + s[n2]*m_n2_k )  #new prop total population in group k
                                            ######Note: an error in this last line is now fixed
                                            
                } #j
            } #i

            #recursion for group size    
            snew[k] <- s[k]*(1 - m_k_n1 - m_k_n2) + s[n1]*m_n1_k + s[n2]*m_n2_k

        } #k

        x <- xnewnew
        s <- snew
        
        
        p <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][1,2] )
        q <- sapply( 1:ngroups , function(z) x[[z]][1,1] + x[[z]][2,1] )
        linkage <- sapply( 1:ngroups , function(z) x[[z]][1,1]*x[[z]][2,2] - x[[z]][1,2]*x[[z]][2,1] )
        linkage <- linkage + 0.5
            
        #update results matrices
        p_mat <- rbind(p_mat, p)
        q_mat <- rbind(q_mat, q)
        D_mat <- rbind(D_mat, linkage-0.5)
        s_mat <- rbind(s_mat, s)
        fit_mat <-rbind(fit_mat, fit/sum(fit))


        # draw
        if ( draw==TRUE ) {
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
            points( 1:ngroups , fitp0 , col="cyan")         #prop of group fitness of behavior 0
            lines( 1:ngroups , fitq0 , col="lightblue", lty=3)         #prop of group fitness of marker 0
            legend (ngroups*.6, 1, legend=c("prop. bhvr 0", "prop. marker 0", "linkage", "prop in grp", "fitness", "payoff p0", "payoff q0"), 
                    col=c("black", "red", "blue", "green", "cyan", "cyan", "lightblue"), pch=c(1,19,19,19,19,1,1))
        } # /if draw
        
    } # /t

   
    print("END")
    invisible(x)


    # Trajectories Plot ------------------------------------------------------------
    if ( plots==TRUE ) {
        plot(0:tmax, p_mat[,1], ylim=c(0,1) , lwd=2, col="black" )
        points(0:tmax, p_mat[,2], ylim=c(0,1) , lwd=2, col="blue")
        lines(0:tmax, q_mat[,1], ylim=c(0,1) , lwd=2, col="black")
        lines(0:tmax, q_mat[,2], ylim=c(0,1) , lwd=2, col="blue")
        lines(0:tmax, D_mat[,1], ylim=c(0,1) , lwd=2, col="black", lty=3)
        lines(0:tmax, D_mat[,2], ylim=c(0,1) , lwd=2, col="blue", lty=3)
        lines(0:tmax, s_mat[,1], ylim=c(0,1) , lwd=3, col="black", lty=4)
        lines(0:tmax, s_mat[,2], ylim=c(0,1) , lwd=3, col="blue", lty=2)
        lines(0:tmax, fit_mat[,1], ylim=c(0,1) , lwd=3, col="darkgrey", lty=3)
        lines(0:tmax, fit_mat[,2], ylim=c(0,1) , lwd=3, col="lightblue", lty=3)
    } # if plots

    #print(x)

    #output matrices
    bound_out <- list(p_mat=p_mat,
                      q_mat=q_mat,
                      D_mat=D_mat,
                      s_mat=s_mat,
                      fit_mat=fit_mat)
    return(bound_out)

    
} #function


res <- feet_marker_sim(
  tmax=100,
  d=0.5, #coordination benefit (hare hunting together is better than alone)
  g=0.5, #extra coordination benefit among mutualists (stag is better than hare)
  h=0.2, #mis-coordination cost for mutualists (lost hare opp)
  a=0.1, #probability of assorting on marker
  m0=0.1, #proportion of each group that migrates
  mu=0.1, #of all migrants, proportion that engage in payoff biased migration
  s=rep(0.5,2), #proportion of total population in each group
  init_p=c(rep(.55,1),rep(.45,1)), #intial proportion of behavior 0 (stag) in each group
  init_q=c(rep(.55,1),rep(.45,1)), #intial proportion of marker 0 in each group 
  draw=F,    #make ending graphs?
  plots=T )  #make time-series plots?

#res$p_mat




##################################################################
#Explore equilibria

a <- seq(0, 1, 0.1)
m0 <- c(0.01, 0.1, 0.5)

tmax <- 300

#initialize vectors
grp1a_p_vec <- NULL
grp2a_p_vec <- NULL

grp1a_q_vec <- NULL
grp2a_q_vec <- NULL

grp1a_D_vec <- NULL
grp2a_D_vec <- NULL

grp1a_s_vec <- NULL
grp2a_s_vec <- NULL

grp1m0_list <- vector("list", length(m0))
grp2m0_list <- vector("list", length(m0))

counter <- 0

#Initial conditions: Group 1 is a minority (small s) with high payoff equilibirum (high p = stag)
#   Group 2 is majority with low payoff equilibrium. Markers start out nearly even in both groups

#R&B predictions:
#1) when markers not important (low a), minority should be converted to all-hare (p=0 in both groups).
#2) if you can assort on marker (high a), minorities with high stag should evolve marked boundaries to protect
#   themselves from effects of payoff-biased migration. In other words, minority with high stag
#   stays a minority with high stag, majority with high hare stays majority with high hare, and
#   population-level norm-marker covariance evolves.


for ( m0_inc in m0 ) {

    counter <- counter + 1

    for ( a_inc in a ) {

        res <- feet_marker_sim(
                    tmax=tmax,
                    d=0.2, #coordination benefit (hare hunting together is better than alone)
                    g=0.2, #extra coordination benefit among mutualists (stag is better than hare)
                    h=0.2, #mis-coordination cost for mutualists (lost hare opp)
                    a=a_inc, #probability of assorting on marker
                    m0=m0_inc, #proportion of each group that migrates
                    mu=0.2, #of all migrants, proportion that engage in payoff biased migration
                    s=c(0.2, 0.8), #proportion of total population in each group
                    init_p=c(0.9, 0.1), #intial proportion of behavior 0 (stag) in each group
                    init_q=c(0.5, 0.4), #intial proportion of marker 0 in each group 
                    draw=F,    #make ending graphs?
                    plots=F )  #make time-series plots?

        grp1a_p_vec <- rbind(grp1a_p_vec, res$p_mat[tmax,1])
        grp2a_p_vec <- rbind(grp2a_p_vec, res$p_mat[tmax,2])

        grp1a_q_vec <- rbind(grp1a_q_vec, res$q_mat[tmax,1])
        grp2a_q_vec <- rbind(grp2a_q_vec, res$q_mat[tmax,2])

        grp1a_D_vec <- rbind(grp1a_D_vec, res$D_mat[tmax,1])
        grp2a_D_vec <- rbind(grp2a_D_vec, res$D_mat[tmax,2])

        grp1a_s_vec <- rbind(grp1a_s_vec, res$s_mat[tmax,1])
        grp2a_s_vec <- rbind(grp2a_s_vec, res$s_mat[tmax,2])

    } #a

    #make lists of matrices containing p, q, D, and s time series for each group
    grp1m0_list[[counter]] <- cbind(grp1a_p_vec, grp1a_q_vec, grp1a_D_vec, grp1a_s_vec)
    grp2m0_list[[counter]] <- cbind(grp2a_p_vec, grp2a_q_vec, grp2a_D_vec, grp2a_s_vec)

    #reset vectors
    grp1a_p_vec <- NULL
    grp2a_p_vec <- NULL
    grp1a_q_vec <- NULL
    grp2a_q_vec <- NULL
    grp1a_D_vec <- NULL
    grp2a_D_vec <- NULL
    grp1a_s_vec <- NULL
    grp2a_s_vec <- NULL

} #m0



pdf(file="./mig_equilibs.pdf",
  #file="C:/Users/jabunce/Desktop/mig_equilibs.pdf", 
  height=10, width=12)
par(mfrow=c(3,2), oma=c(5,7,5,5), mar=c(3,2,2,2))


#group 1 plot, m0[1]

plot( a , grp1m0_list[[1]][,1] , ylim=c(0,1) , lwd=8 , type="n", col="black", cex.axis=1.5, #plot set-up
    #main="gnk=0.3, mk=mnk=0.1, dnk=2, beta=0.2",
    xlab="Prob of assortment on norm (a)", ylab=c("Long-run freqs") )

lines( a , grp1m0_list[[1]][,1] , ylim=c(0,1) , lwd=9 , type="l", col="black" )            #p m0[1]
lines( a , grp1m0_list[[1]][,2] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,0,1,0.5) )     #q
lines( a , grp1m0_list[[1]][,3] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,1,0,0.5) )     #D
lines( a , grp1m0_list[[1]][,4] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp1m0_list[[2]][,1] , ylim=c(0,1) , lwd=5 , type="l", col="black" )            #p m0[2]
# lines( a , grp1m0_list[[2]][,2] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp1m0_list[[2]][,3] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp1m0_list[[2]][,4] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp1m0_list[[3]][,1] , ylim=c(0,1) , lwd=1 , type="l", col="black" )            #p m0[3]
# lines( a , grp1m0_list[[3]][,2] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp1m0_list[[3]][,3] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp1m0_list[[3]][,4] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(1,0,0,0.5) )     #s

text(x=0.2, y=0.85, label="Grp 1, m0=0.01", cex=2)

legend(x="right", c("p","q","D", "s"),lty=1,
    lwd=c(10), col=c("black", rgb(0,0,1,0.5), rgb(0,1,0,0.5), rgb(1,0,0,0.5)), cex=1.5 )


#group 2 plot, m0[1]

plot( a , grp2m0_list[[1]][,1] , ylim=c(0,1) , lwd=8 , type="n", col="black", cex.axis=1.5, #plot set-up
    #main="gnk=0.3, mk=mnk=0.1, dnk=2, beta=0.2",
    xlab="Prob of assortment on norm (a)", ylab=c("Long-run freqs") )

lines( a , grp2m0_list[[1]][,1] , ylim=c(0,1) , lwd=9 , type="l", col="black" )            #p m0[1]
lines( a , grp2m0_list[[1]][,2] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,0,1,0.5) )     #q
lines( a , grp2m0_list[[1]][,3] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,1,0,0.5) )     #D
lines( a , grp2m0_list[[1]][,4] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp2m0_list[[2]][,1] , ylim=c(0,1) , lwd=5 , type="l", col="black" )            #p m0[2]
# lines( a , grp2m0_list[[2]][,2] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp2m0_list[[2]][,3] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp2m0_list[[2]][,4] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp2m0_list[[3]][,1] , ylim=c(0,1) , lwd=1 , type="l", col="black" )            #p m0[3]
# lines( a , grp2m0_list[[3]][,2] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp2m0_list[[3]][,3] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp2m0_list[[3]][,4] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(1,0,0,0.5) )     #s

text(x=0.2, y=0.85, label="Grp 2, m0=0.01", cex=2)



#group 1 plot, m0[2]

plot( a , grp1m0_list[[1]][,1] , ylim=c(0,1) , lwd=8 , type="n", col="black", cex.axis=1.5, #plot set-up
    #main="gnk=0.3, mk=mnk=0.1, dnk=2, beta=0.2",
    xlab="Prob of assortment on norm (a)", ylab=c("Long-run freqs") )

# lines( a , grp1m0_list[[1]][,1] , ylim=c(0,1) , lwd=9 , type="l", col="black" )            #p m0[1]
# lines( a , grp1m0_list[[1]][,2] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp1m0_list[[1]][,3] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp1m0_list[[1]][,4] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(1,0,0,0.5) )     #s

lines( a , grp1m0_list[[2]][,1] , ylim=c(0,1) , lwd=5 , type="l", col="black" )            #p m0[2]
lines( a , grp1m0_list[[2]][,2] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,0,1,0.5) )     #q
lines( a , grp1m0_list[[2]][,3] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,1,0,0.5) )     #D
lines( a , grp1m0_list[[2]][,4] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp1m0_list[[3]][,1] , ylim=c(0,1) , lwd=1 , type="l", col="black" )            #p m0[3]
# lines( a , grp1m0_list[[3]][,2] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp1m0_list[[3]][,3] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp1m0_list[[3]][,4] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(1,0,0,0.5) )     #s

text(x=0.2, y=0.85, label="Grp 1, m0=0.1", cex=2)

#legend(x="left", c("p","q","D", "s"),lty=1,
#    lwd=c(10), col=c("black", rgb(0,0,1,0.5), rgb(0,1,0,0.5), rgb(1,0,0,0.5)), cex=1.5 )


#group 2 plot, m0[1]

plot( a , grp2m0_list[[1]][,1] , ylim=c(0,1) , lwd=8 , type="n", col="black", cex.axis=1.5, #plot set-up
    #main="gnk=0.3, mk=mnk=0.1, dnk=2, beta=0.2",
    xlab="Prob of assortment on norm (a)", ylab=c("Long-run freqs") )

# lines( a , grp2m0_list[[1]][,1] , ylim=c(0,1) , lwd=9 , type="l", col="black" )            #p m0[1]
# lines( a , grp2m0_list[[1]][,2] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp2m0_list[[1]][,3] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp2m0_list[[1]][,4] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(1,0,0,0.5) )     #s

lines( a , grp2m0_list[[2]][,1] , ylim=c(0,1) , lwd=5 , type="l", col="black" )            #p m0[2]
lines( a , grp2m0_list[[2]][,2] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,0,1,0.5) )     #q
lines( a , grp2m0_list[[2]][,3] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,1,0,0.5) )     #D
lines( a , grp2m0_list[[2]][,4] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp2m0_list[[3]][,1] , ylim=c(0,1) , lwd=1 , type="l", col="black" )            #p m0[3]
# lines( a , grp2m0_list[[3]][,2] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp2m0_list[[3]][,3] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp2m0_list[[3]][,4] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(1,0,0,0.5) )     #s

text(x=0.2, y=0.85, label="Grp 2, m0=0.1", cex=2)



#group 1 plot, m0[3]

plot( a , grp1m0_list[[1]][,1] , ylim=c(0,1) , lwd=8 , type="n", col="black", cex.axis=1.5, #plot set-up
    #main="gnk=0.3, mk=mnk=0.1, dnk=2, beta=0.2",
    xlab="Prob of assortment on norm (a)", ylab=c("Long-run freqs") )

# lines( a , grp1m0_list[[1]][,1] , ylim=c(0,1) , lwd=9 , type="l", col="black" )            #p m0[1]
# lines( a , grp1m0_list[[1]][,2] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp1m0_list[[1]][,3] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp1m0_list[[1]][,4] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp1m0_list[[2]][,1] , ylim=c(0,1) , lwd=5 , type="l", col="black" )            #p m0[2]
# lines( a , grp1m0_list[[2]][,2] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp1m0_list[[2]][,3] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp1m0_list[[2]][,4] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(1,0,0,0.5) )     #s

lines( a , grp1m0_list[[3]][,1] , ylim=c(0,1) , lwd=1 , type="l", col="black" )            #p m0[3]
lines( a , grp1m0_list[[3]][,2] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,0,1,0.5) )     #q
lines( a , grp1m0_list[[3]][,3] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,1,0,0.5) )     #D
lines( a , grp1m0_list[[3]][,4] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(1,0,0,0.5) )     #s

text(x=0.2, y=0.85, label="Grp 1, m0=0.5", cex=2)

#legend(x="left", c("p","q","D", "s"),lty=1,
#    lwd=c(10), col=c("black", rgb(0,0,1,0.5), rgb(0,1,0,0.5), rgb(1,0,0,0.5)), cex=1.5 )


#group 2 plot, m0[3]

plot( a , grp2m0_list[[1]][,1] , ylim=c(0,1) , lwd=8 , type="n", col="black", cex.axis=1.5, #plot set-up
    #main="gnk=0.3, mk=mnk=0.1, dnk=2, beta=0.2",
    xlab="Prob of assortment on norm (a)", ylab=c("Long-run freqs") )

# lines( a , grp2m0_list[[1]][,1] , ylim=c(0,1) , lwd=9 , type="l", col="black" )            #p m0[1]
# lines( a , grp2m0_list[[1]][,2] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp2m0_list[[1]][,3] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp2m0_list[[1]][,4] , ylim=c(0,1) , lwd=10 , type="l", col=rgb(1,0,0,0.5) )     #s

# lines( a , grp2m0_list[[2]][,1] , ylim=c(0,1) , lwd=5 , type="l", col="black" )            #p m0[2]
# lines( a , grp2m0_list[[2]][,2] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,0,1,0.5) )     #q
# lines( a , grp2m0_list[[2]][,3] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(0,1,0,0.5) )     #D
# lines( a , grp2m0_list[[2]][,4] , ylim=c(0,1) , lwd=6 , type="l", col=rgb(1,0,0,0.5) )     #s

lines( a , grp2m0_list[[3]][,1] , ylim=c(0,1) , lwd=1 , type="l", col="black" )            #p m0[3]
lines( a , grp2m0_list[[3]][,2] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,0,1,0.5) )     #q
lines( a , grp2m0_list[[3]][,3] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(0,1,0,0.5) )     #D
lines( a , grp2m0_list[[3]][,4] , ylim=c(0,1) , lwd=2 , type="l", col=rgb(1,0,0,0.5) )     #s

text(x=0.2, y=0.85, label="Grp 2, m0=0.5", cex=2)


mtext(text="Long-run values", side=2, outer=TRUE, line=2, cex=2, las=3, adj=0.5)
mtext(text="Prob of assortment on norm (a)", side=1, outer=TRUE, line=1, cex=2, las=1, adj=0.5)
mtext(text="s=c(0.2, 0.8), init_p=c(0.9, 0.1), init_q=c(0.5, 0.4), d=g=h=mu=0.2",
    side=3, outer=TRUE, line=1, cex=2, las=1, adj=0.5)

graphics.off()

#conclusions:
#If migration is even moderately high, the minority group gets converted to all-hare, regardless of
#assortment on markers.
#If migration and assortment are both low, you kind-of get a marked ethnic boundary.
#If migration is low and assortment is high, the majority adopts the minority norms and markers.
#The minority never stays a minority. Population prop always reaches ~0.5 in each group.