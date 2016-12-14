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
                             max_cases = 1e3) {
    if (lambda < 0) {
        stop("lambda cannot be negative")
    }
    if (R < 0) {
        stop("R cannot be negative")
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
    ## exit conditions: by exceeding a maximum number of cases.
    

    ## Initialize the simulations
    n_cases <- 1L
    t <- 0L
    new_cases <- 1L
    continue <- (n_cases <= max_cases)

    ## Create output material
    ## for contacts
    out_ct <- list()
    out_ct$from  <- NULL # contact from
    out_ct$to  <- NULL # contact to
    out_ct$infectious  <- NULL # TRUE = contact led to transmission
    out_ct$date  <- NULL # date of contact

    ## for linelist
    out_ll <- list()
    out_ll$case_id  <- 1L # ID of case
    out_ll$onset  <- 0L # onset of symptom
    names(out_ll$onset) <- out_ll$case_id
    
    
    while (continue) {

        ## i is the ID of the current infectious case from which contacts
        ## (infectious or not) are being created; it is taken from a vector of
        ## case indices which grows every time new cases are created; an exit
        ## condition is embedded in 'continue' if this vector has a length of 0,
        ## so that as new iterations start, there is always a new 'i' to use
        
        i <- new_cases[1L]
        new_cases <- new_cases[-1L]
        
        ## i) how many contacts?
        n_contacts <- stats::rpois(1L, lambda)
        new_ids <- seq(max(1L, out_ct$to)+1, length.out = n_contacts)

        ## ii) are contacts leading to new cases?
        are_cases <- sample(c(TRUE, FALSE), n_contacts, replace = TRUE,
                            prob = c(Rc, 1 - Rc))
        new_cases <- c(new_cases, new_ids[are_cases])
        
        ## iii) when does this happen?
        ## note: the date of onset of 'i' serves as a basis 
        i_onset <- out_ll$onset[as.character(i)] 
        contacts_dates <- rep(i_onset, n_contacts)

        ## contact doesn't lead to a case
        contacts_dates[!are_cases] <- r_contact(sum(!are_cases))

        ## contact leads to a new case
        new_onsets <- i_onset + r_si(sum(are_cases))
        names(new_onsets) <- new_ids[are_cases]
        contacts_dates[are_cases] <- as.integer(round(runif(sum(are_cases),
                                                      min = i_onset,
                                                      max = new_onsets)))
        names(contacts_dates) <- new_ids
        
        ## store contact information
        out_ct$from <- c(out_ct$from, rep(i, n_contacts))
        out_ct$to <- c(out_ct$to, new_ids)
        out_ct$infectious <- c(out_ct$infectious, are_cases)
        out_ct$dates <- c(out_ct$dates, contacts_dates)

        ## store cases in linelist
        out_ll$case_id  <- c(out_ll$case_id, new_ids[are_cases])
        out_ll$onset  <- c(out_ll$onset, new_onsets)


        ## exit conditions: exceed maximum number of cases, or no new
        ## infectious cases left
        n_cases <- length(out_ll$case_id)
        continue <- (n_cases <= max_cases) && (length(new_cases) > 0)
        
    }

    ## shape the output: the returned object will be an epicontacts object,
    ## implemented in the epicontacts package; the loop above will classically
    ## overshoot the maximum number of cases, so we need to prune the resulting
    ## object to make sure there is not more cases than requested.
    
    out_ct <- data.frame(out_ct)
    out_ll <- data.frame(out_ll)

    out <- epicontacts::make_epicontacts(linelist = out_ll,
                                         contacts = out_ct,
                                         directed = TRUE)

    ## remove excess cases in linelist
    out <- out[i = seq_len(max_cases)]

    ## remove corresponding contacts
    id_linelist <- epicontacts::get_id(out, "linelist")
    
    ct_to_keep <- (!out$contacts$infectious) |
        (out$contacts$from %in% id_linelist) |
            (out$contacts$to %in% id_linelist)
 
    out <- out[j = ct_to_keep]
   
    return(out)
    
}

