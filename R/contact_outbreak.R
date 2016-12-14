#' A contact-based outbreak simulator
#'
#' Under development, please do not use without contacting the author.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' 
#' @export
#' 
#' @param lambda The average number of contacts per infectious case, taken to be
#'     the rate of a Poisson distribution.
#' 
#' @param R The average reproduction number, i.e. number of secondary cases per
#'     infected individual. Note that this number cannot exceed the average
#'     number of contacts \code{lambda}.
#' 
#' @param r_SI A function generating random delays from the serial interval
#'     distribution, i.e. the delay between primary and secondary symptoms. See
#'     details for more information on how to specify this distribution.
#'
#' @param r_contact A function generating random delays from onset of symptom
#'     to non-infectious contacts.
#'
#' @param max_time The maximum number of time units the simulation can be run
#'     for.
#'
#' @param max_cases The maximum number of new cases the simulation can be run
#'     for.
#'
#' @details
#' We recommend using the package \code{distcrete} for generating the serial
#'     interval distribution. See: \url{https://github.com/reconhub/distcrete}.
#'
#'
#' @examples
#' 
contact_outbreak <- function(lambda, R, r_SI, r_contact,
                         max_time = 100, max_cases = 1e3){
    if (lambda < 0) {
        stop("lambda cannot be negative")
    }
    if (R < 0) {
        stop("R cannot be negative")
    }
    if (max_time < 0) {
        stop("max_time cannot be negative")
    }
    if (max_cases < 1) {
        stop("max_cases cannot be less than 1")
    }
    if (R > lambda) {
        stop("R cannot exceed lambda")
    }
    if (!inherits(r_si, "function")) {
        stop("r_si is not a function")
    }
    if (!inherits(r_contact, "function")) {
        stop("r_contact is not a function")
    }
    
    
    ## Rc is the contact-reprodution number, i.e. the probability that a new
    ## contact will lead to a secondary infection.
    
    Rc <- R / lambda

    
    ## The simulator works as a loop iterating trough new infectious cases; the
    ## steps are:

    ## i) determine the number of contacts ~ Poisson(lambda)

    ## ii) determine which contacts create a new case ~ Rc

    ## iii) determine when the contacts took place, depending on whether they
    ## led to a secondary case (r_SI) or not (r_contact)

    ## As the number of new cases will typically grow exponentially, we define
    ## two exit conditions: by exceeding a given simulation time, and by
    ## exceeding a maximum number of cases.


    ## Initialize the simulations
    n_cases <- 1L
    t <- 0L
    id_new_cases <- 1L
    continue <- (n_cases <= max_cases) && (t <= max_time)

    ## Create output material
    ## for contacts
    out_from  <- NA_integer_ # contact from
    out_to  <- 1 # contact to
    out_infectious  <- TRUE # TRUE = contact led to transmission
    out_date  <- 0 # date of contact

    ## for linelist
    out_case_id  <- 1L # ID of case
    out_onset  <- 0 # onset of symptom
    names(out_onset) <- out_case_id
    
    
    while (continue) {
        continue <- (n_cases <= max_cases) && (t <= max_time)

        ## i is the ID of the current infectious case
        i <- new_cases[1L]
        
        for (i in id_new_cases) {
            i_onset <- out_onset[as.character(i)] 
            
            ## i) how many contacts?
            n_contacts <- stats::rpois(1L, lambda)
            new_ids <- seq(max(out_to)+1, length = n_contacts)

            ## ii) are contacts leading to new cases?
            are_cases <- sample(c(TRUE, FALSE), n_contacts, replace = TRUE,
                                prob = c(Rc, 1 - Rc))
            new_cases <- new_ids[are_cases]
            
            ## iii) when does this happen?
            contacts_dates <- t + integer(n_contacts)
            contacts_dates[!are_cases] <- r_contact(sum(!are_cases))
            new_onsets <- i_onset + r_si(sum(are_cases))
            contacts_dates[are_cases] <- as.integer(runif(n_contacts,
                                                          min = t,
                                                          max = new_onsets))
            
            ## store contact information
            out_from <- c(out_from, i)
            out_to <- c(out_to, new_ids)
            out_infectious <- c(out_infectious, are_cases)
            out_dates <- c(out_dates, contacts_dates)

            ## store cases in linelist
            out_case_id  <- c(out_case_id, new_cases)
            out_onset  <- c(out_onset, new_onsets)
   
        }
    }

    
}

