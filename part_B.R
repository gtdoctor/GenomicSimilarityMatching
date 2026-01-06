
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
#  threshcutoff  # number 0-1 # smaller == tighter
#  nc  #  max number of controls that permitted per case within distance t for case to be included
#  ncmin # <=nc case excluded if fewer controls found. 
# useMILP#  increase set of possible matches using mixed integer linear programming implementation of min cost max flow algorithm 
GSM_B <- function(wdt=ccdistances,
                  datadir=DATADIR, 
                  outdirp=OUTDIRP, 
                  psam=PSAM, 
                  fileoutstem = FILEOUTSTEM,
                  threshcutoff=THRESHCUTOFF,
                  nc = NC,
                  ncmin = NCMIN,
                  useMILP=MILP, 
                  MILPverbose=FALSE)
{
  setwd(datadir)
  cohort=fread(psam, col.names = c("FID", "IID", "F", "M", "Sex", "PHENO1"))
  outdirp=file.path(outdirp) # output directory parent. No final slash 
  
  setwd(outdirp)
  
  text_choice=paste0("Max controls = ", nc, ", Min controls = ", ncmin, ", Thresholding = ", threshcutoff, " of median")
  print(text_choice)
  fileoutkept=paste0("kept.", fileoutstem)
  fileoutremoved=paste0("removed.",fileoutstem)
  fileoutallmatched=paste0("allmatches.",fileoutstem)
  fileoutfiltmatched=paste0("filteredmatches.",fileoutstem)
  
  ac = function(x){as.character(x)}
  an = function(x){as.numeric(ac(x))}
  
    ##define first_match function "localbestworst"
  localbestworst <- function(wdt, nc, ncmin){
    
    
    
    
    # Initialize an empty data.table
    matched_casescontrols <- data.table(
      i_index = numeric(),
      j_index = numeric(),
      d_value = numeric()
    )
    
    # Initialize lists for keeping and removing
    cases_kept <- character()
    cases_removed <- character()
    controls_kept <- character()  # If using lists, initialize as vector() and then convert to unique values at the end
    controls_removed <- character()  # Similarly, initialize appropriately
    
    # Create "shortest_distances" matrix with i rows and three cols (d, j, count)
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
    
    # priority match are those with <ncmin and still matchable
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
      
      #add the largest_d to the matched_casecontrol list and increment the count for that case
      
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


### RUN THE FUNCTIONS
#####------------- Operations on distance matrix ------- #########
  
  
  # calculate cutoff based on median of all case-control differences
  median_distance <- median(wdt, na.rm = TRUE)
  
  threshold=threshcutoff * median_distance
  
  case_ids <- cohort$IID[cohort$PHENO1==2]  # Assuming the case IDs are stored in 'IID' column of 'cases_use'
  control_ids <- cohort$IID[cohort$PHENO1==1]  # Assuming the case IDs are stored in 'IID' column of 'cases_use'
  print(text_choice)
  print(sprintf("Seeking to match %s cases to pool of %s controls, Calculated threshold based on cutoff of %s of median distance is %s", length(case_ids), length(control_ids), threshcutoff, round(threshold,2)))
  
  
  # Create a thresholded version of the ccdistances matrix: 
  # Replace NA for any D > threshold in ccdistances:
  
  wdt[wdt > threshold] <- NA
  
## Run the first match
firstmatch<-localbestworst(wdt=wdt, nc=nc, ncmin=ncmin)
# Count cases and then apply nc filter.
firstmatch_case_counts <- firstmatch$matched_casescontrols[, .N, by = i_index]

nmatch1<-sum(firstmatch_case_counts$N>=ncmin)
print(sprintf("Local worst-best algorithm complete. %s cases matched", nmatch1))
                    
## identify deficit matches - # Deficit cases are feasible cases 
##make deficit matrix - deficit cases and previously reserved but unused ctrls.

all_case_indices <- seq_len(nrow(wdt))
firstmatch_cases_kept_indices <- firstmatch_case_counts[N >= ncmin, i_index] # firstmatch_cases_kept_indices is a DT

  # Exclude cases with insufficient possible matches from the universe first
  # deficit cases include feasilble cases which were not matched at all in first_match, 
  # but if the set exists at all, it must mean that some other cases part-matched

deficit_cases_indices<-setdiff(firstmatch$feasible_case_indices,firstmatch_cases_kept_indices)

if(length(deficit_cases_indices)<=1){
  print("One or fewer deficit cases")
final_matched_casecontrols<-firstmatch$matched_casescontrols[firstmatch$matched_casescontrols$i_index %in% firstmatch_cases_kept_indices,]
}else{
  if(useMILP==TRUE){
    print("Performing MILP mathcing on deficit cases")
    # Filtered matches for deficit edges and all deficit cases 
    deficit_matches <- firstmatch$matched_casescontrols[firstmatch$matched_casescontrols$i_index %in% deficit_cases_indices]
    wdt_def <- wdt[deficit_cases_indices, deficit_matches$j_index ]
    rownames(wdt_def)<-deficit_cases_indices
    colnames(wdt_def)<- deficit_matches$j_index


    ## Run remainder match
    # convert to edges_datatable and run milp 
    edges_def <- (mat_to_edgedt(mat=wdt_def))

    sprintf("Running MILP on %s cases and %s controls", nrow(edges_def), ncol(edges_def) )
    secondmatch <- MILP_solver(casectrledges = edges_def, ncmin = ncmin, solver = "highs")

    ## combine matches
    final_matched_casecontrols<-firstmatch$matched_casescontrols[firstmatch$matched_casescontrols$i_index %in% firstmatch_cases_kept_indices,]
    if(!is.null(secondmatch$matches)) {
      print(sprintf("MILP returns %s additional fully matched cases", length(unique(secondmatch$matches$i_index))))
      final_matched_casecontrols <- rbind(final_matched_casecontrols,secondmatch$matches)
    } else {print("No additional MILP matches")}
  } 
 else {
  print("MILP not performed for deficit cases")
  final_matched_casecontrols<-firstmatch$matched_casescontrols[firstmatch$matched_casescontrols$i_index %in% firstmatch_cases_kept_indices,]
  }

}

cases_kept_indices <- an(unique(final_matched_casecontrols$i_index))
print(paste0("Number of cases matched = ", length(cases_kept_indices)))



# names of controls kept and removed 

controls_kept <- control_ids[an(final_matched_casecontrols$j_index)]
controls_removed <- setdiff(control_ids, controls_kept)
cases_kept <- case_ids[cases_kept_indices]
cases_removed <- setdiff(case_ids, cases_kept)

# convert indices back to names
filtered_matched_names <- data.frame(control_id = controls_kept, 
                                     case_id = case_ids[an(final_matched_casecontrols$i_index)], 
                                     distance = final_matched_casecontrols$d_value,
                                     stringsAsFactors = FALSE)
filtered_matched_names <- filtered_matched_names[order(filtered_matched_names$case_id, filtered_matched_names$distance), ]


# create merged list of IIDs to keep. 
outputtable1 <- data.frame(IID = ac(c(cases_kept, controls_kept)))
outputtable2 <- data.frame(IID = ac(c(cases_removed, controls_removed)))

#### RESULTS #####
#naming outfiles

current_date <- Sys.Date()
filename_append = paste0("threshold", threshcutoff, "_minmatches", ncmin,"_maxmatches", nc, "_", current_date)

for (name in c("fileoutkept", "fileoutremoved", "fileoutallmatched", "fileoutfiltmatched")) {
  # Get the value of the object identified by the name
  value <- get(name)
  # Create the new variable name
  new_var_name <- paste0("o", name)
  # Assign the modified value to the new variable
  assign(new_var_name, paste0(value, filename_append, ".tsv"))
}


nctrlskept = length(controls_kept)
ncaseskept = length(cases_kept)
ncasesremoved= length(cases_removed)
print(text_choice)
text_results=(paste0("# Results: ", length(control_ids), " controls and ", length(case_ids), " cases processed. ", nctrlskept , " controls matched to ", ncaseskept, " cases, and ",ncasesremoved," cases removed."))
print(text_results)
#### --- saving out
#

# writing out files
outdir=paste0(outdirp,"/", fileoutstem, "_",filename_append)
dir.create(outdir, recursive = F)
setwd(outdir)

write(paste0("# Ancestry Matching: ", text_choice, "\n", text_results), file = ofileoutkept )
write.table(outputtable1,file=ofileoutkept, quote=F,row.names=F,col.names=F,sep="\t", append = TRUE)
write.table(outputtable2,file=ofileoutremoved,quote=F,row.names=F,col.names=F,sep="\t", append = FALSE)
write.table(filtered_matched_names, file=ofileoutfiltmatched, quote=F, row.names=F,col.names=F,sep="\t", append = FALSE)




} # end GSM_B



