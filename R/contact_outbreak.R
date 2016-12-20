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
#' if (require(distcrete)) {
#' ## get distributions
#'  SI <- distcrete("gamma", 1L, w = 0, 10, 0.65)
#'  d_contacts <- distcrete("exp", 1L, w=0, 0.05)
#'
#' ## simulate outbreak and contacts
#'  set.seed(1)
#'  x <- contact_outbreak(3, 1.8, SI$r, d_contacts$r)
#'  x
#'  plot(x, group = "case_def")
#'
#'  if (require(incidence)) {
#'    plot(incidence(x$linelist$onset, 7))
#'  }
#' }
#'
#'
contact_outbreak <- function(lambda, R, r_SI, r_contact,
                             max_cases = 100) {
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
    if (!inherits(r_SI, "function")) {
        stop("r_SI is not a function")
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
    out_ll$id  <- 1L # ID of case
    out_ll$case_def <- "confirmed"
    out_ll$onset  <- 0L # onset of symptom
    names(out_ll$onset) <- out_ll$id


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
        new_ids <- seq(max(out_ll$id) + 1, length.out = n_contacts)

        ## ii) are contacts leading to new cases?
        are_cases <- sample(c(TRUE, FALSE), n_contacts, replace = TRUE,
                            prob = c(Rc, 1 - Rc))
        new_cases <- c(new_cases, new_ids[are_cases])
        new_def <- rep("non-case", n_contacts)
        new_def[are_cases] <- "confirmed"

        ## iii) when does this happen?
        ## note: the date of onset of 'i' serves as a basis
        i_onset <- out_ll$onset[as.character(i)]
        contacts_dates <- rep(i_onset, n_contacts)

        ## contact doesn't lead to a case
        contacts_dates[!are_cases] <- r_contact(sum(!are_cases))

        ## contact leads to a new case
        new_onsets <- rep(NA_integer_, n_contacts)
        new_onsets[are_cases] <- i_onset + r_SI(sum(are_cases))
        names(new_onsets) <- new_ids
        contacts_dates[are_cases] <- as.integer(
            round(stats::runif(sum(are_cases),
                               min = i_onset,
                               max = new_onsets[are_cases])))
        names(contacts_dates) <- new_ids

        ## store contact information
        out_ct$from <- c(out_ct$from, rep(i, n_contacts))
        out_ct$to <- c(out_ct$to, new_ids)
        out_ct$infectious <- c(out_ct$infectious, are_cases)
        out_ct$date_of_contact <- c(out_ct$date_of_contact, contacts_dates)

        ## store cases in linelist
        out_ll$id  <- c(out_ll$id, new_ids)
        out_ll$case_def <- c(out_ll$case_def, new_def)
        out_ll$onset  <- c(out_ll$onset, new_onsets)
        out_ll$date_of_isolation <- NA

        ## exit conditions: exceed maximum number of cases, or no new
        ## infectious cases left
        n_cases <- sum(out_ll$case_def == "confirmed")
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

    ## ## remove excess cases in linelist
    ## out <- out[i = seq_len(max_cases)]

    ## ## remove corresponding contacts
    ## id_linelist <- epicontacts::get_id(out, "linelist")

    ## ct_to_keep <- (!out$contacts$infectious) |
    ##     (out$contacts$from %in% id_linelist) |
    ##         (out$contacts$to %in% id_linelist)

    ## out <- out[j = ct_to_keep]

    return(out)

}

#' Replays an outbreak up until a specific point in time
#'
#' @param full_outbreak the full outbreak datasets
#'
#' @return a list with named components having functions (replay/isolate_at) operating on
#'         the outbreak
#'
#' @examples
#'
#' set.seed(1235)
#' SI <- distcrete("gamma", 1L, w = 0, 10, 0.65)
#' d_contacts <- distcrete("exp", 1L, w=0, 0.05)
#' x <- contact_outbreak(3, 1.8, SI$r, d_contacts$r, max_cases = 10)
#' replay <- replay_outbreak(x)
#' plot(replay$replay(99))
#' outbreak_2 <- replay$isolate_at(2, 24)
#' plot(outbreak_2$replay(99))
#' outbreak_3 <- outbreak_2$isolate_at(3, 25)
#' plot(outbreak_3$replay(99))
#'
#' @export
replay_outbreak <- function(full_outbreak) {
  base_outbreak <- full_outbreak
  stopifnot("date_of_contact" %in% colnames(base_outbreak$contacts))
  stopifnot("onset" %in% colnames(base_outbreak$linelist))
  stopifnot("date_of_isolation" %in% colnames(base_outbreak$linelist))
  stopifnot("case_def" %in% colnames(base_outbreak$linelist))
  stopifnot(sum(base_outbreak$linelist$onset == 0, na.rm = TRUE) == 1) # only one index case
  replay <- function(until = pmax(max(base_outbreak$linelist$onset, na.rm = TRUE),
                        max(base_outbreak$contacts$date_of_contact, na.rm = TRUE))) {
    new_contacts <- dplyr::filter_(base_outbreak$contacts, ~date_of_contact <= until)
    new_linelist <- dplyr::filter_(base_outbreak$linelist,
                                   ~id %in% unique(c(new_contacts$from, new_contacts$to)) | onset == 0)
    new_linelist <- dplyr::mutate_(new_linelist,
                                   case_def = ~dplyr::if_else(onset > until, "non-case", "confirmed", "non-case"),
                                   onset = ~dplyr::if_else(onset > until, NA_real_, onset))

    # remove all links from people that have been isolated
    new_contacts <- dplyr::select_(dplyr::filter_(dplyr::left_join(dplyr::mutate_(new_contacts, from = ~as.character(from)),
                                                    dplyr::select_(new_linelist, .dots = c("id", "date_of_isolation")),
                                                    by = c("from" = "id")),
                                   ~is.na(date_of_isolation) | date_of_contact < date_of_isolation), ~-date_of_isolation)

    # filter out contacts that go from a non-case to a non-case
    non_cases <- dplyr::filter_(new_linelist, ~case_def == "non-case")
    new_contacts <- dplyr::filter_(new_contacts, ~!(from %in% non_cases$id & to %in% non_cases$id))

    # remove non-cases that have no links
    new_linelist <- dplyr::filter_(new_linelist,
                                   ~ case_def != "non-case" |
                                     id %in% unique(c(new_contacts$from, new_contacts$to)) |
                                     onset == 0)

    epicontacts::make_epicontacts(linelist = new_linelist,
                                  contacts = new_contacts,
                                  directed = TRUE)
  }
  isolate_at <- function(case_id, date) {
    new_linelist <- dplyr::mutate_(base_outbreak$linelist, date_of_isolation = ~ifelse(id == case_id, date, date_of_isolation))

    # now remove prune the outbreak
    new_contacts <- base_outbreak$contacts

    # remove invalid contacts
    new_contacts <- dplyr::left_join(dplyr::mutate_(new_contacts, from = ~as.character(from)),
                                       dplyr::select_(new_linelist, .dots = c("id", "date_of_isolation")),
                                       by = c("from" = "id"))
    new_contacts <- dplyr::select_(
      dplyr::filter_(new_contacts, ~is.na(date_of_isolation) | date_of_contact < date_of_isolation),
      ~-date_of_isolation)
    new_outbreak <- epicontacts::make_epicontacts(linelist = new_linelist,
                                                  contacts = new_contacts,
                                                  directed = TRUE)

    # now remove anything that is not connected
    g <- epicontacts::as.igraph.epicontacts(new_outbreak)
    root_id <- dplyr::filter_(new_outbreak$linelist, ~onset == 0)$id
    connected_component <- igraph::bfs(g, root = root_id, neimode = "out", unreachable = FALSE)$order
    connected_component <- connected_component[!is.na(connected_component)]
    replay_outbreak(new_outbreak[as.character(connected_component), as.character(connected_component)])
  }
  list(replay = replay, isolate_at = isolate_at, full_outbreak = base_outbreak)
}
