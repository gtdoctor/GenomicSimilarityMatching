# load or install function

load_or_install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

load_or_install(c("data.table", 
                "parallel"))


# PARAMETERS #

GSM_B <- function(wdt=ccdistances,
                  datadir=DATADIR, 
                  outdirp=OUTDIRP, 
                  SAMPLEFILE=SAMPLEFILE, 
                  SAMPLEFILE_IDCOL=2 ,
                  SAMPLEFILE_PHENOCOL=6,
                  SAMPLEFILE_BATCHCOL=NULL,
                  MATCHSELF=FALSE,
                  COMMONBATCH=NULL,
                  FORBIDDEN=NULL,      
                  fileoutstem = FILEOUTSTEM,
                  maxthresh=1,
                  minthresh=0,
                  nc = NC,
                  ncmin = NCMIN,
                  useMILP=MILP,
                  STRATA_REASSIGN=FALSE,
                  REFINE_TOPK=20,
                  MILPverbose=FALSE
                  )
{
  setwd(datadir)
  cohort=fread(SAMPLEFILE)

if (is.null(SAMPLEFILE_BATCHCOL)) {
  print("No batch information given: any case can be matched to any controls")
} else {
  print(sprintf("Batch column is %s", SAMPLEFILE_BATCHCOL))
  if (MATCHSELF == TRUE) {
    if (!is.null(COMMONBATCH)) {
      print(sprintf("Matching self batch and common batch %s", COMMONBATCH))
    } else {
      print("Matching self batch only (no common batch given)")
    }
  } else {
    if (!is.null(COMMONBATCH)) {
      print(sprintf("Matching only to common batch %s", COMMONBATCH))
    } else {
      stop("Error: not self-matching and no common batch given")
    }
  }
}



  ## guard against missing data, duplicates:
  if (!is.null(COMMONBATCH)) {
  req_cols <- c(SAMPLEFILE_IDCOL, SAMPLEFILE_PHENOCOL, SAMPLEFILE_BATCHCOL)
  } else { 
  req_cols <- c(SAMPLEFILE_IDCOL, SAMPLEFILE_PHENOCOL)
  }
bad <- vapply(req_cols, function(i) {
  any(is.na(cohort[[i]]) | cohort[[i]] == "")
}, logical(1))

if (any(bad)) {
  stop(
    sprintf(
      "Missing or NA values detected in cohort columns: %s",
      paste(req_cols[bad], collapse = ", ")
    ),
    call. = FALSE
  )
}

id <- cohort[[SAMPLEFILE_IDCOL]]

if (any(duplicated(id))) {
  dup_ids <- unique(id[duplicated(id)])
  stop(
    sprintf(
      "Duplicated sample IDs detected (%d): %s",
      length(dup_ids),
      paste(head(dup_ids, 10), collapse = ", ")
    ),
    call. = FALSE
  )
}


  outdirp=file.path(outdirp) # output directory parent. No final slash 
  

  fileoutkept=paste0("kept.", fileoutstem)
  fileoutremoved=paste0("removed.",fileoutstem)
  fileoutallmatched=paste0("allmatches.",fileoutstem)
  fileoutfiltmatched=paste0("filteredmatches.",fileoutstem)
  

  ac = function(x){as.character(x)}
  an = function(x){as.numeric(ac(x))}
  
    ##define first_match function "localbestworst"
  localbestworst <- function(wdt, nc, ncmin){
    
   
   
    
    # An empty data.table
    matched_casescontrols <- data.table(
      i_index = numeric(),
      j_index = numeric(),
      d_value = numeric()
    )
    
    # Lists for keeping and removing
    cases_kept <- character()
    cases_removed <- character()
    controls_kept <- character()  
    
    # Create "shortest_distances" matrix with i rows and three cols (d, j, count - distance to current control, count of matches)
    shortest_distances <- matrix(NA, nrow = nrow(wdt), ncol = 3)
    colnames(shortest_distances) <- c("d", "j", "count")
    
    # Assign NA to cases where there aren't enough controls left.
    eligible_counts <- rowSums(!is.na(wdt))
    feasible_case_indices<-which(eligible_counts >= ncmin)
    insufficient_possible <- which(eligible_counts < ncmin)
    
    wdt[insufficient_possible , ] <- NA
    
    print(sprintf(" %s cases removed for insufficient possible controls. %s feasible cases", length(insufficient_possible), length(feasible_case_indices)))
    print("Running the local worst-best algorithm")
    
    # For each case i, find smallest D and its j (control)_index
    min_d = numeric()
    j_index = numeric()
    index_largest_d = numeric()
    replace_i = numeric()
    
    for (i in 1:nrow(wdt)) {
      # Find the minimum distance and its index if there are non-NA values
      j_index <- which.min(wdt[i, ])
      if (length(j_index)>0) {
        min_d <- wdt[i, j_index]
        shortest_distances[i, ] <- c(min_d, j_index, 0) # update the shortest_distances matrix
      } else {
        # Assign NA to d and j but keep the row consistent - ie no matching controls at the start. 
        shortest_distances[i, ] <- c(NA, NA, 0) # Assign 0 or NA to 'count' 
      }
    }
    
    # priority match are those with <ncmin and still matchable (relevant where nc != ncmin )
    priority_match <- rep_len(1,length.out = nrow(shortest_distances))
    priority_match[which(is.na(shortest_distances[,1]))] <-0
    
    
    # Loop until conditions are met
    repeat {
      
      
      # Check if the process should continue
      # if there are no matches left, index_largest_d will be NA. 
      # If in a previous round a case has matched nc times, it will have been left as NA
      if (all(is.na(shortest_distances[, 1]))) break
      
      if(any(priority_match==1)){
        index_largest_d <- which(priority_match==1)[which.max(shortest_distances[,1][priority_match==1])]        
      }else{
        index_largest_d <- which.max(shortest_distances[,1])
      }
      
      j_index <- an(shortest_distances[index_largest_d, 2])
      d_value <- an(shortest_distances[index_largest_d, 1])
      
    # add the largest_d to the matched_casecontrol list and increment the count for that case
      
      matched_casescontrols <- rbind(
        matched_casescontrols,
        data.table(
          j_index = j_index,
          i_index = index_largest_d,
          d_value = d_value
        )
      )
      shortest_distances[index_largest_d, 3] <- shortest_distances[index_largest_d, 3] + 1
      
      # block the chosen control
      wdt[, j_index] <- NA
      
      
      # identify the control used in prev step, and all the cases matched to it. 
      replace_i <- which(shortest_distances[, 2 ] == j_index )
      
      # add NA to the cases and matched controls that are to be replaced. These will remain NA if no further match.
      shortest_distances[replace_i, 1:2] <- NA
      #set priority_match value to 0 temporarilty for cases that are being reassesed.
      priority_match[replace_i]<-0
      
      # if case has been fully matched , do not include it in search for new controls
      if (shortest_distances[index_largest_d,3] == nc) {
        # Remove index_largest_d from the replace_i if nc has been met
        replace_i <- setdiff(replace_i, index_largest_d)
      }  
      
      # nb if there are no other cases sharing controls and nc has been reached, replace_i will be empty 
      # and the loop restarts
      # find replacement controls j and d, if needed:
      if (length(replace_i) > 0) {
        has_non_na <- numeric()
        # Identify rows in replace_i with any non-NA value; has_non_na is a logical vector of replace_i length
        has_non_na <- apply(!is.na(wdt[replace_i, , drop = FALSE]), 1, any)
        
        # Check if there are any rows with non-NA values
        if (any(has_non_na)) {
          # Step 2: For rows with any non-NA values, find the minimum and its index
          for (i in replace_i[has_non_na]) {
            r.j_index <- which.min(wdt[i, ])
            r.min_d <- wdt[i,r.j_index]
            #     j.rmNA = c(j.rmNA, r.j_index)
            shortest_distances[i, 1] <- r.min_d
            shortest_distances[i, 2] <- r.j_index
            
            #for cases that do find  a new match, only give priority if they are undermatched. 
            if(shortest_distances[i,3]<ncmin){priority_match[i]<- 1}
          }
          # Block out these controls (cols) for further use
          #wdt[,unique(j.rmNA)] <- NA
        }
      }
    }    # end repeat block 
    
    
    return(list(
      matched_casescontrols = matched_casescontrols,
      insufficient_possible = insufficient_possible,
      feasible_case_indices = feasible_case_indices
    ))
  } # end localbestworst

if( useMILP==TRUE){
  load_or_install(c("ompr", 
                    "ompr.roi", 
                    "ROI", 
                    "igraph", 
                    "ROI.plugin.highs"))

### define function to convert matrix to edges_dt
mat_to_edgedt <- function(mat, case_ids = rownames(mat), ctrl_ids = colnames(mat)) 
{
  idx <- which(!is.na(mat), arr.ind = TRUE)
  data.table(
    i_index = case_ids[idx[, 1]],
    j_index = ctrl_ids[idx[, 2]],
    d_value = mat[idx]
  )
}

### define milp match function
MILP_solver <- function(casectrledges, ncmin, solver = "highs" )
{
  sel_lists <- list()
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  
  g   <- igraph::graph_from_data_frame(casectrledges[, .(i_index, j_index)], directed = FALSE)
  cmp <- igraph::components(g)$membership
  case_comp   <- cmp[match(as.character(casectrledges$i_index), names(cmp))]
  edges_split <- split(casectrledges, case_comp, drop = TRUE)
  
  print("MILP part1_start")
  sprintf("[%s] Subgraphs=%d\n", format(Sys.time(), "%H:%M:%S"), length(edges_split))
  
  part1 <- function(DT, ncmin, solver, idx) {
    
    DT <- data.table::copy(DT); if (nrow(DT) == 0L) return(list(cases = character(0), matches = DT[0]))
    cases <- sort(unique(DT$i_index))
    ctrls <- sort(unique(DT$j_index))
    DT[, i := match(i_index, cases)]
    DT[, j := match(j_index, ctrls)]
    DT[, eid := .I]
    n_cases <- length(cases)
    n_ctrls <- length(ctrls)
    case_edges <- split(DT$eid, DT$i)
    ctrl_edges <- split(DT$eid, DT$j); Eids <- DT$eid
    
    model1 <- ompr::MIPModel() %>%
      ompr::add_variable(x[e], e = Eids, type = "binary") %>%
      ompr::add_variable(y[i], i = 1:n_cases, type = "binary") %>%
      ompr::add_constraint(ompr::sum_expr(x[e], e = ctrl_edges[[j]]) <= 1, j = 1:n_ctrls) %>%
      # enforce exactly ncmin edges if y[i]==1 and 0 if y[i]==0
      ompr::add_constraint(ompr::sum_expr(x[e], e = case_edges[[i]]) >= ncmin * y[i], i = 1:n_cases) %>%
      ompr::add_constraint(ompr::sum_expr(x[e], e = case_edges[[i]]) <= ncmin * y[i], i = 1:n_cases) %>%
      ompr::set_objective(ompr::sum_expr(y[i], i = 1:n_cases), "max")
    
    sol1 <- ompr::solve_model(
      model1,
      ompr.roi::with_ROI(
        solver   = solver,
        verbose  = MILPverbose,
        control  = list(
          threads    = 1,
          time_limit = 10800,
          presolve   = "on"
        )
      )
    )
    
    y_sel <- ompr::get_solution(sol1, y[i])
    if (!nrow(y_sel)) return(list(cases = character(0), matches = DT[0]))
    
    chosen_cases <- cases[y_sel$i[y_sel$value == 1]]
    x_sel <- ompr::get_solution(sol1, x[e])
    chosen_eids <- x_sel$e[x_sel$value == 1]
    matches <- DT[eid %in% chosen_eids, .(i_index, j_index, d_value)][order(i_index, j_index)]
    
    list(cases = chosen_cases, matches = matches)
    
  } # end part 1 define 

# run part 1
  res_list<- lapply(seq_along(edges_split), function(ii) {
    part1(edges_split[[ii]], ncmin = ncmin, solver = solver, idx = ii)
  })
  
  res_list <- list(
    cases = unlist(lapply(res_list, `[[`, "cases")),
    matches = rbindlist(lapply(res_list, `[[`, "matches"))
  )
  
# output of part 1 is a filtered df of matcheable cases, and the edges of these matches. If running part 2, all needed is the list of matchable cases

  # combine
  selected_cases <- unique(res_list$cases)
  if (length(selected_cases) == 0) {  print("No deficit_cases matchable by MILP");   return(list()) 
  } else { 
    print(paste0(length(selected_cases), " deficit cases matchable. Optimising..."))}
  #all_matches <- rbindlist(lapply(res_list, `[[`, "matches"), use.names = TRUE, fill = TRUE)
  #all_matches <- unique(all_matches)
  
    

# define part 2 (min cost max flow)
  part2 <- function(DT, sc) {
    if (is.null(sc) || !length(sc)) return(list(cases = character(0), matches = DT[0][, .(i_index, j_index, d_value)]))
    DT <- data.table::copy(DT)[i_index %in% sc]
    if (!nrow(DT)) return(list(cases = character(0), matches = DT[0][, .(i_index, j_index, d_value)]))
    
    cases <- sort(unique(DT$i_index))
    ctrls <- sort(unique(DT$j_index))
    DT[, i := match(i_index, cases)]
    DT[, j := match(j_index, ctrls)]
    DT[, eid := .I]
    
    n_cases <- length(cases)
    n_ctrls <- length(ctrls)
    case_edges <- split(DT$eid, DT$i)
    ctrl_edges <- split(DT$eid, DT$j)
    Eids <- DT$eid
    costs <- DT$d_value[match(Eids, DT$eid)]
    
    model2 <- ompr::MIPModel() %>%
      ompr::add_variable(x[e], e = Eids, type = "continuous", lb = 0, ub = 1) %>%  # TU ⇒ integral for b-matching
      ompr::add_constraint(ompr::sum_expr(x[e], e = ctrl_edges[[j]]) <= 1, j = 1:n_ctrls) %>%
      ompr::add_constraint(ompr::sum_expr(x[e], e = case_edges[[i]]) == ncmin, i = 1:n_cases) %>%
      ompr::set_objective(ompr::sum_expr(costs[match(e, Eids)] * x[e], e = Eids), "min")
    
    sol2 <- ompr::solve_model(
      model2,
      ompr.roi::with_ROI(solver = solver, verbose = MILPverbose,
                         control = list(threads = 1, presolve = "on"))
    )
    
    x_sel <- ompr::get_solution(sol2, x[e])
    chosen_eids <- x_sel$e[x_sel$value > 1e-9]
    matches <- DT[eid %in% chosen_eids, .(i_index, j_index, d_value)][order(i_index, j_index)]
    
    # Derive which cases actually got ncmin edges
    got_counts <- matches[, .N, by = i_index]
    chosen_cases <- got_counts[N == ncmin, i_index]
    
    list(cases = chosen_cases, matches = matches)
  
  

   }  # end part 2 define 


#  subset deficit_cases to solvable cases
     sc_by_comp <- split(as.character(selected_cases),
                        cmp[match(as.character(selected_cases), names(cmp))])

    keys <- intersect(names(edges_split), names(sc_by_comp))

    res_list <- lapply(keys, function(k) part2(edges_split[[k]], sc_by_comp[[k]]))

  # combine
  all_cases <- unique(unlist(lapply(res_list, `[[`, "cases"), use.names = FALSE))
  all_matches <- rbindlist(lapply(res_list, `[[`, "matches"), use.names = TRUE, fill = TRUE)
  all_matches <- unique(all_matches)


  list(cases = all_cases, matches = all_matches)
} # end MILP_solver
} # end if(useMILP==TRUE)

### BATCH BLOCKING ~~~
preblock_wdt <- function(wdt,
                         case_batch,
                         control_batch,
                         matchself ,
                         common_batches) {

  if (is.null(common_batches)) {
    common_batches <- character(0)
  }
  non_common_ctrl <- which(!control_batch %in% common_batches)
  # Rule 1: If not self-matching, only allow controls in common_batches
  if (!matchself) {
    if (length(non_common_ctrl)) {
      wdt[, non_common_ctrl] <- NA_real_
    }
    return(wdt)
  }

  # Rule 2: If self-matching
  if (!length(non_common_ctrl)) return(wdt)
# for each non-common control batch , block out cases that are from different batch
  for (b in unique(control_batch[non_common_ctrl])) {
    ctrl_idx <- which(control_batch == b)
    case_idx <- which(case_batch != b)
    if (length(ctrl_idx) && length(case_idx)) {
      wdt[case_idx, ctrl_idx] <- NA_real_
    }
  }
  wdt
}
### Part 3 Refinement ###

## helper: build sparse admissible subgraph on kept cases x kept controls.
## edges kept:
##   1) nearest top_k kept controls per kept case — a control may be in several cases' top_k,
##      giving it multiple candidate edges and so freedom to be reassigned
##   2) every incumbent (currently assigned) edge, re-read from the full wdt
## Consequences:
##   - no currently selected control can be discarded by pruning
##   - the incumbent assignment is itself feasible, so refinement can only lower total distance
##   - a control in no case's top_k has its incumbent edge only, so stays with its current case
## Indices: i_index / j_index are ORIGINAL wdt indices; ok/ord/keep_j are newwdt column positions
build_refine_edges_sparse <- function(wdt,
                                      kept_case_idx,
                                      kept_ctrl_idx = NULL,
                                      current_matches = NULL,
                                      top_k) {
  if (is.null(current_matches) || !nrow(current_matches)) {
    stop("current_matches must be provided and non-empty")
  }

  if (is.null(kept_ctrl_idx)) {
    kept_ctrl_idx <- sort(unique(as.integer(current_matches$j_index)))
  } else {
    kept_ctrl_idx <- sort(unique(as.integer(kept_ctrl_idx)))
  }

  kept_case_idx <- sort(unique(as.integer(kept_case_idx)))

  newwdt <- wdt[kept_case_idx, kept_ctrl_idx, drop = FALSE] ## a new smaller matrix of all included strata
  rownames(newwdt) <- as.character(kept_case_idx)   # preserve ORIGINAL wdt row indices
  colnames(newwdt) <- as.character(kept_ctrl_idx)   # preserve ORIGINAL wdt col indices

  edges_list <- vector("list", length = length(kept_case_idx) + 1L)

  ## top_k controls per kept case
  for (ii in seq_along(kept_case_idx)) {
    d <- newwdt[ii, ]
    ok <- which(!is.na(d))
    if (!length(ok)) next

    ord <- ok[order(d[ok], kept_ctrl_idx[ok])]
    keep_j <- head(ord, top_k)

    edges_list[[ii]] <- data.table(
      i_index = kept_case_idx[ii],
      j_index = kept_ctrl_idx[keep_j],
      d_value = as.numeric(d[keep_j])
    )
  }

  ## always keep incumbent matched edges
  incumbent <- data.table(
    i_index = as.integer(current_matches$i_index),
    j_index = as.integer(current_matches$j_index)
  )
  incumbent <- unique(incumbent)

  incumbent[, d_value := wdt[cbind(i_index, j_index)]]
  incumbent <- incumbent[!is.na(d_value)]

  edges_list[[length(edges_list)]] <- incumbent

  edgesfinal <- unique(rbindlist(edges_list, use.names = TRUE, fill = TRUE))
  edgesfinal <- edgesfinal[order(i_index, d_value, j_index)]

  list(
    newwdt = newwdt,
    edgesfinal = edgesfinal,
    kept_case_idx = kept_case_idx,
    kept_ctrl_idx = kept_ctrl_idx
  )
}


## refinement MILP:
## rules:
## - every selected control must be used exactly once
## - every case must have at least ncmin and at most nc controls
## - minimise total distance
##
## IMPORTANT:
## i_index / j_index in returned matches remain ORIGINAL indices from wdt
refine_fixed_casectrl_assignment <- function(casectrledges,
                                             ncmin,
                                             nc,
                                             solver = "highs",
                                             MILPverbose = FALSE) {
  if (!nrow(casectrledges)) stop("casectrledges is empty")

  DT <- data.table::copy(casectrledges)

  cases <- sort(unique(as.integer(DT$i_index)))   # original case indices
  ctrls <- sort(unique(as.integer(DT$j_index)))   # original ctrl indices

  n_cases <- length(cases)
  n_ctrls <- length(ctrls)

  ## necessary count check
  if (n_ctrls < n_cases * ncmin) {
    stop(sprintf("Infeasible by counts: %s controls < %s cases * ncmin=%s",
                 n_ctrls, n_cases, ncmin))
  }
  if (n_ctrls > n_cases * nc) {
    stop(sprintf("Infeasible by counts: %s controls > %s cases * nc=%s",
                 n_ctrls, n_cases, nc))
  }

  ## local indices for solver
  DT[, i := match(as.integer(i_index), cases)]
  DT[, j := match(as.integer(j_index), ctrls)]
  DT[, eid := .I]

  case_edges <- split(DT$eid, DT$i)
  ctrl_edges <- split(DT$eid, DT$j)
  Eids <- DT$eid
  costs <- DT$d_value

  model <- ompr::MIPModel() %>%
    ompr::add_variable(x[e], e = Eids, type = "continuous", lb = 0, ub = 1) %>%

    ## every selected control used exactly once
    ompr::add_constraint(
      ompr::sum_expr(x[e], e = ctrl_edges[[j]]) == 1,
      j = 1:n_ctrls
    ) %>%

    ## each selected case between ncmin and nc
    ompr::add_constraint(
      ompr::sum_expr(x[e], e = case_edges[[i]]) >= ncmin,
      i = 1:n_cases
    ) %>%
    ompr::add_constraint(
      ompr::sum_expr(x[e], e = case_edges[[i]]) <= nc,
      i = 1:n_cases
    ) %>%

    ## minimise total distance
    ompr::set_objective(
      ompr::sum_expr(costs[e] * x[e], e = Eids),
      sense = "min"
    )

  sol <- ompr::solve_model(
    model,
    ompr.roi::with_ROI(
      solver = solver,
      verbose = MILPverbose,
      control = list(
        threads = 1,
        presolve = "on",
        time_limit = 10800
      )
    )
  )

  x_sel <- ompr::get_solution(sol, x[e])
  if (!nrow(x_sel)) stop("No solution returned by solver")

  ## because this is a flow/b-matching LP, solution should be integral
    nonint <- x_sel[abs(x_sel$value - round(x_sel$value)) > 1e-8, ]
  if (nrow(nonint)) {
    stop(sprintf("LP returned %s non-integer variables; expected integral solution", nrow(nonint)))
  }

  chosen_eids <- x_sel$e[x_sel$value > 0.5]

  matches <- DT[eid %in% chosen_eids, .(i_index, j_index, d_value)][order(i_index, d_value, j_index)]

  ## safety checks
  ctrl_counts <- matches[, .N, by = j_index]
  if (any(ctrl_counts$N != 1L)) stop("Refinement error: some controls were not used exactly once")

  case_counts <- matches[, .N, by = i_index]
  if (any(case_counts$N < ncmin | case_counts$N > nc)) {
    stop("Refinement error: some cases outside ncmin..nc")
  }

  list(
    cases = sort(unique(as.integer(matches$i_index))),
    controls = sort(unique(as.integer(matches$j_index))),
    matches = matches
  )
}

### RUN THE FUNCTIONS
#####------------- Operations on distance matrix ------- #########

text_choice=sprintf("Max controls = %s, Min controls = %s, minthresh = %s centile, maxthresh = %s of median",
nc, ncmin,
minthresh,   maxthresh)
print(text_choice)

# calculate cutoff based on median of all case-control differences
print(paste0("Sense check: nrow of sample IDs loaded: ", nrow(cohort)))
print(paste0("Sense check: nrow of distancematrix loaded: ", nrow(wdt), " (should equal number cases)"))

# update matrix to remove impermissible matches based on batch criteria and maxdisting.
ncommonbatch=0
case_ids <- cohort[[SAMPLEFILE_IDCOL]][cohort[[SAMPLEFILE_PHENOCOL]] == 2]
control_ids <- cohort[[SAMPLEFILE_IDCOL]][cohort[[SAMPLEFILE_PHENOCOL]] == 1]

if (nrow(wdt) != length(case_ids)) {
  stop(sprintf("nrow(wdt)=%s but number of cases=%s", nrow(wdt), length(case_ids)))
}
if (ncol(wdt) != length(control_ids)) {
  stop(sprintf("ncol(wdt)=%s but number of controls=%s", ncol(wdt), length(control_ids)))
}

if (!is.null(SAMPLEFILE_BATCHCOL)) {
  batch <- cohort[[SAMPLEFILE_BATCHCOL]]
  if (!is.null(COMMONBATCH)) {
    ncommonbatch <- sum(batch %in% COMMONBATCH)
  }

  case_batch <- batch[match(case_ids, cohort[[SAMPLEFILE_IDCOL]])]
  control_batch <- batch[match(control_ids, cohort[[SAMPLEFILE_IDCOL]])]

  if (anyNA(case_batch)) stop("Some case_ids not found in cohort$IID (batch missing).")
  if (anyNA(control_batch)) stop("Some control_ids not found in cohort$IID (batch missing).")

  wdt <- preblock_wdt(
    wdt,
    case_batch = case_batch,
    control_batch = control_batch,
    common_batches = COMMONBATCH,
    matchself = MATCHSELF
  )
}

if (!is.null(FORBIDDEN)) {
  cat(sprintf("Forbidden matches file %s given. Head of file:\n", FORBIDDEN))
  forbidden <- fread(FORBIDDEN, header = TRUE, colClasses = "character") # header row "ctrl case" is expected.
  setnames(forbidden, c("ctrl", "case")) # in case typos
  print(head(forbidden))
  cat(sprintf("Number of forbidden provided: %s\n", nrow(forbidden)))
  case_index    <- setNames(seq_along(case_ids), case_ids)
  control_index <- setNames(seq_along(control_ids), control_ids)
  fi <- case_index[forbidden$case]
  fj <- control_index[forbidden$ctrl]
  ok <- !is.na(fi) & !is.na(fj)
  cat(sprintf("Pruning %s forbidden matches\n", sum(ok)))
  wdt[cbind(fi[ok], fj[ok])] <- NA_real_
}

# calculate median distance between allowed case-ctrl pairs
threshs <- quantile(wdt, probs = c(minthresh, 0.5), na.rm = TRUE, names = FALSE)
mindist <- threshs[1]
median_distance <- threshs[2]
maxdist <- maxthresh * median_distance

if (is.null(SAMPLEFILE_BATCHCOL)) { ncommonbatch <- length(control_ids) }
print(sprintf("Seeking to match %s cases to pool of %s controls, of which %s are matchable for all cases. Calculated maxdist based on cutoff of %s of median distance is %s",
              length(case_ids), length(control_ids), ncommonbatch, maxthresh, round(maxdist, 2)))
print(sprintf("Calculated mindist based on minthresh of %s is %s", minthresh, round(mindist, 2)))

# Create a maxdisted version of the ccdistances matrix:
wdt[wdt > maxdist] <- NA
if (minthresh > 0) {
  wdt[wdt < mindist] <- NA
}

###### Run the first match
firstmatch <- localbestworst(wdt = wdt, nc = nc, ncmin = ncmin)

###### Count cases and then apply nc filter.
firstmatch_case_counts <- firstmatch$matched_casescontrols[, .N, by = i_index]

nmatch1 <- sum(firstmatch_case_counts$N >= ncmin)
print(sprintf("Local worst-best algorithm complete. %s cases matched", nmatch1))

##### identify deficit matches
firstmatch_cases_kept_indices <- firstmatch_case_counts[N >= ncmin, i_index]

deficit_cases_indices <- setdiff(firstmatch$feasible_case_indices, firstmatch_cases_kept_indices)
deficit_matches <- firstmatch$matched_casescontrols[firstmatch$matched_casescontrols$i_index %in% deficit_cases_indices]

if (length(deficit_cases_indices) <= 1 || length(unique(deficit_matches$j_index)) < ncmin) {
  print("One or fewer deficit cases or no excess controls")
  final_matched_casecontrols <- firstmatch$matched_casescontrols[
    firstmatch$matched_casescontrols$i_index %in% firstmatch_cases_kept_indices, ]
} else {
  if (useMILP == TRUE) {
    print("Performing MILP matching on deficit cases")

    wdt_def <- wdt[deficit_cases_indices, deficit_matches$j_index, drop = FALSE]
    rownames(wdt_def) <- as.character(deficit_cases_indices)
    colnames(wdt_def) <- as.character(deficit_matches$j_index)

    edges_def <- mat_to_edgedt(mat = wdt_def)

    print(sprintf("Running MILP on %s admissible edges", nrow(edges_def)))
    secondmatch <- MILP_solver(casectrledges = edges_def, ncmin = ncmin, solver = "highs")

    final_matched_casecontrols <- firstmatch$matched_casescontrols[
      firstmatch$matched_casescontrols$i_index %in% firstmatch_cases_kept_indices, ]

    if (!is.null(secondmatch$matches) && nrow(secondmatch$matches) > 0) {
      print(sprintf("MILP returns %s additional fully matched cases",
                    length(unique(secondmatch$matches$i_index))))
      final_matched_casecontrols <- rbind(final_matched_casecontrols, secondmatch$matches)
    } else {
      print("No additional MILP matches")
    }
  } else {
    print("MILP not performed for deficit cases")
    final_matched_casecontrols <- firstmatch$matched_casescontrols[
      firstmatch$matched_casescontrols$i_index %in% firstmatch_cases_kept_indices, ]
  }
}

cases_kept_indices <- an(unique(final_matched_casecontrols$i_index))
print(paste0("Number of cases matched = ", length(cases_kept_indices)))

################# final MILP: optional reassignment ###############
# Reallocates the already-selected controls among the kept cases to minimise
# total distance. Introduces no new controls and drops no cases.
# STRATA_REASSIGN=FALSE (default) skips this stage, reproducing
# GSM_B_local_MILP_batch_mindist.forbid.R.
if (STRATA_REASSIGN == TRUE && useMILP == TRUE && nrow(final_matched_casecontrols) > 0) {
  print(sprintf("Refining final matched set by sparse bounded assignment (top %s controls per kept case + incumbent edges)",
                REFINE_TOPK))

  refine_input <- build_refine_edges_sparse(
    wdt = wdt,
    kept_case_idx = cases_kept_indices,
    current_matches = final_matched_casecontrols,
    top_k = REFINE_TOPK
  )

  print(sprintf(
    "Refinement graph: %s kept cases, %s kept controls, %s admissible edges",
    length(refine_input$kept_case_idx),
    length(refine_input$kept_ctrl_idx),
    nrow(refine_input$edgesfinal)
  ))

  refined <- refine_fixed_casectrl_assignment(
    casectrledges = refine_input$edgesfinal,
    ncmin = ncmin,
    nc = nc,
    solver = "highs",
    MILPverbose = MILPverbose
  )

  final_matched_casecontrols <- refined$matches
  cases_kept_indices <- an(unique(final_matched_casecontrols$i_index))
} else {
  print("Reassignment of kept controls not performed")
}

# names of controls kept and removed
controls_kept <- control_ids[an(final_matched_casecontrols$j_index)]
controls_removed <- setdiff(control_ids, controls_kept)
cases_kept <- case_ids[cases_kept_indices]
cases_removed <- setdiff(case_ids, cases_kept)

# convert indices back to names
filtered_matched_names <- data.frame(
  control_id = controls_kept,
  case_id = case_ids[an(final_matched_casecontrols$i_index)],
  distance = final_matched_casecontrols$d_value,
  stringsAsFactors = FALSE
)
filtered_matched_names <- filtered_matched_names[order(filtered_matched_names$case_id,
                                                       filtered_matched_names$distance), ]

# create merged list of IIDs to keep.
outputtable1 <- data.frame(IID = ac(c(cases_kept, controls_kept)))
outputtable2 <- data.frame(IID = ac(c(cases_removed, controls_removed)))

#### RESULTS #####
current_date <- Sys.Date()
filename_append = sprintf("maxt%smint%s_minm%smaxm%s_%s",
                          maxthresh, minthresh,
                          ncmin, nc,
                          current_date)

for (name in c("fileoutkept", "fileoutremoved", "fileoutallmatched", "fileoutfiltmatched")) {
  value <- get(name)
  new_var_name <- paste0("o", name)
  assign(new_var_name, paste0(value, filename_append, ".tsv"))
}

nctrlskept = length(controls_kept)
ncaseskept = length(cases_kept)
ncasesremoved = length(cases_removed)
print(text_choice)
text_results = paste0("# Results: ", length(control_ids), " controls and ", length(case_ids),
                      " cases processed. ", nctrlskept , " controls matched to ",
                      ncaseskept, " cases, and ", ncasesremoved, " cases removed.")
print(text_results)

#### --- saving out
outdir = paste0(outdirp, "/", fileoutstem, "_", filename_append)
dir.create(outdir, recursive = TRUE)
setwd(outdir)

write(paste0("# Ancestry Matching: ", text_choice, "\n", text_results), file = ofileoutkept)
write.table(outputtable1, file = ofileoutkept, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
write.table(outputtable2, file = ofileoutremoved, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = FALSE)
write.table(filtered_matched_names, file = ofileoutfiltmatched, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = FALSE)

print(sprintf("Files saved to %s", outdir))

} # end GSM_B