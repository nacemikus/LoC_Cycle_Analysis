library(rmatio)
library(tidyverse)
# library(cowplot) # for plotting grids
# library(lme4)
# library(xtable)
library(lmerTest)
source("theme_functions.r")
library(readr)
library(stringr)
library(jsonlite)
library(GGally)
rm(list = ls())
getwd()

per_trial_data_baseline <- read_csv("per_trial_data_baseline.csv", show_col_types = FALSE)
per_trial_data_weekly <- read_csv("per_trial_data_weekly.csv", show_col_types = FALSE)
per_trial_data_final <- read_csv("per_trial_data_final.csv", show_col_types = FALSE)

participant_data = per_trial_data_baseline %>% select(prolific_id, phq, shaps, ius, gad, )
# sort dates fu:

parse_any_date <- function(date_vec) {
  cleaned <- date_vec %>%
    # Replace backslashes, Unicode dashes, dots, commas, etc.
    str_replace_all("\\\\+", "/") %>%     # double backslashes to slash
    str_replace_all("–", "-") %>%         # unicode dash
    str_replace_all("—", "-") %>%         # em dash
    str_replace_all("[.]", "-") %>%
    str_replace_all("[,]", "") %>%
    str_replace_all("\\s+", "-")          # spaces to dash
  
  # Define common date formats (like the ones you mentioned)
  orders <- c(
    "ymd",        # 2025-07-08
    "dmy",        # 08-07-2025
    "mdy",        # 07/08/2025
    "Ymd",        # 20250708
    "d-b-Y",      # 08-jul-2025
    "Y-B-d",      # 2025-July-08
    "Y/m/d",      # 2025/07/08
    "d/m/Y",      # 08/07/2025
    "m/d/y",      # 07/08/25
    "Y.m.d"       # 2025.07.08
  )
  
  parsed <- parse_date_time(cleaned, orders = orders, exact = FALSE)
  return(as_date(parsed))  # ensure Date class
}
per_trial_data_weekly$date_of_bleeding_dt   <- parse_any_date(per_trial_data_weekly %>% pull(date))
per_trial_data_baseline$date_of_bleeding_dt   <- parse_any_date(per_trial_data_baseline %>% pull(date1))
per_trial_data_final$date_of_bleeding_dt <- NA

per_trial_data_baseline$record_date_dt = as.Date(per_trial_data_baseline$record_date, format = f)
per_trial_data_weekly$record_date_dt = as.Date(per_trial_data_weekly$record_date, format = f)
per_trial_data_final$record_date_dt <- NA

per_trial_data_weekly$day_of_cycle =   as.integer(per_trial_data_weekly$record_date_dt -per_trial_data_weekly$date_of_bleeding_dt)
per_trial_data_baseline$day_of_cycle =   as.integer(per_trial_data_baseline$record_date_dt -per_trial_data_baseline$date_of_bleeding_dt)
per_trial_data_final$day_of_cycle <- NA

per_trial_data_baseline = per_trial_data_baseline %>% filter(trial_number > 0)

# get questionnaire data #####


# Define the ordered options
phq_diff_levels <- c(
  "Not difficult at all",
  "Somewhat difficult",
  "Very difficult",
  "Extremely difficult"
)

# Function to parse the JSON string and map to 0–3
parse_phq_diff <- function(json_str) {
  parsed <- fromJSON(json_str)
  response <- parsed$Colors[[1]]
  match(response, phq_diff_levels) - 1  # match() is 1-based, so subtract 1
}

# Apply to your column
per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(phq_diff_num = sapply(phq_diff, parse_phq_diff),
         gad_diff_num = sapply(gad_diff, parse_phq_diff))
# compute sum scors of other scales ####
json_str = per_trial_data_baseline$shaps[1]
compute_scale_score <- function(json_str,
                                reverse_vector = NULL,
                                attention_check_position = NULL,
                                attention_check_target = NULL,
                                scale_prefix = "shaps",
                                max_score = 3,
                                min_score = 0) {
  # Parse answers
  if (is.na(json_str) || json_str == "") {
    return(setNames(as.list(rep(NA_real_, length(reverse_vector))), 
                    paste0(scale_prefix, "_item", seq_along(reverse_vector))))
  }
  
  answers <- tryCatch(fromJSON(json_str), error = function(e) return(rep(NA_real_, length(reverse_vector))))
  scores <- as.numeric(answers)
  n_items <- length(scores)
  # attention_check_flags_vect = rep(0, n_items)
  # attention_check_flags_vect[attention_check_position] = 1
  # Reverse scoring
  if (!is.null(reverse_vector)) {
    scores <- ifelse(reverse_vector[seq_len(n_items)],
                     max_score - scores + min_score,
                     scores)
  }
  if (!is.null(attention_check_position))  {
    attention_check_passed = scores[attention_check_position] == attention_check_target
    scores = scores[-attention_check_position]
    
  }
  
  # Build named list of items
  item_scores <- setNames(as.list(scores), paste0(scale_prefix, "_item", seq_along(scores)))
  
  if (!is.null(attention_check_position))  {
    item_scores[[paste0("attention_check_",scale_prefix)]] <- attention_check_passed
    
  }
  
  
  item_scores
  # Add attention check flags
  # if (!is.null(attention_check_flags))    item_scores[attention_check_position] = paste0(item_scores[attention_check_position], "_is_attention_check")
  
}


# Function to compute SHAPS score for one JSON string
shaps_reverse <- c(
  FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
  FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE
)
max_score_shaps =3
min_score_shaps =0

phq8_attention_check_question = 7
phq8_attention_check_target = 2

ius_attention_check_question = 8
ius_attention_check_target = 0
per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(shaps_item_data = lapply(shaps, compute_scale_score,
                                  reverse_vector = shaps_reverse,
                                  max_score = 3,
                                  scale_prefix = "shaps")) %>%
  unnest_wider(shaps_item_data)
per_trial_data_baseline <- per_trial_data_baseline %>%
  rowwise() %>%
  mutate(shaps_total = sum(c_across(starts_with("shaps_item")), na.rm = TRUE)) %>%
  ungroup()

per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(item_data = lapply(phq, compute_scale_score,
                            attention_check_position = phq8_attention_check_question,
                            attention_check_target = phq8_attention_check_target,
                            scale_prefix = "phq")) %>%
  unnest_wider(item_data)

per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(item_data = lapply(gad, compute_scale_score,
                            scale_prefix = "gad")) %>%
  unnest_wider(item_data)
per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(item_data = lapply(ius, compute_scale_score,
                            attention_check_position = ius_attention_check_question,
                            attention_check_target = ius_attention_check_target,
                            scale_prefix = "ius")) %>%
  unnest_wider(item_data)


per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(item_data = lapply(stress, compute_scale_score,
                            scale_prefix = "stress")) %>%
  unnest_wider(item_data)

per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(item_data = lapply(pss, compute_scale_score,
                            scale_prefix = "pss")) %>%
  unnest_wider(item_data)

per_trial_data_baseline <- per_trial_data_baseline %>%
  mutate(item_data = lapply(drsp, compute_scale_score,
                            scale_prefix = "drsp")) %>%
  unnest_wider(item_data)

prospective_items <- paste0("ius_item", c(1, 2, 4, 8, 9, 11))
inhibitory_items  <- paste0("ius_item", c(3, 5, 6, 7, 10, 12))

per_trial_data_baseline <- per_trial_data_baseline %>%
  rowwise() %>%
  mutate(shaps_total = sum(c_across(starts_with("shaps_item")), na.rm = TRUE),
         phq_total = sum(c_across(starts_with("phq_item")), na.rm = TRUE),
         gad_total = sum(c_across(starts_with("gad_item")), na.rm = TRUE),
         pss_total = sum(c_across(starts_with("pss_item")), na.rm = TRUE),
         drsp_total = sum(c_across(starts_with("drsp_item")), na.rm = TRUE),
         stress_total = sum(c_across(starts_with("stress_item")), na.rm = TRUE),
         ius_total = sum(c_across(starts_with("ius_item")), na.rm = TRUE),
         ius_prospective = sum(c_across(all_of(prospective_items)), na.rm = TRUE),
         ius_inhibitory  = sum(c_across(all_of(inhibitory_items)), na.rm = TRUE)) %>%
  ungroup()



per_trial_data_baseline = per_trial_data_baseline %>%
  mutate(
    assigned_sex_clean = sapply(assigned_sex, function(x) fromJSON(x)$assigned_sex),
    menstrual_cycle_clean = sapply(menstrual_cycle, function(x) fromJSON(x)$menstrual_cycle)
  ) 
per_trial_data_baseline %>%
  select(prolific_id,  assigned_sex_clean, menstrual_cycle_clean, shaps_total, phq_total:ius_inhibitory) %>%
  group_by(prolific_id) %>%
  summarize(
    assigned_sex_clean = first(assigned_sex_clean), 
    menstrual_cycle_clean = first(menstrual_cycle_clean),
    shaps_total   = first(shaps_total),
    phq_total       = first(phq_total),
    gad_total       = first(gad_total),
    pss_total       = first(pss_total),
    drsp_total      = first(drsp_total),
    stress_total    = first(stress_total),
    ius_total       = first(ius_total),
    ius_prospective = first(ius_prospective),
    ius_inhibitory  = first(ius_inhibitory),
    .groups = "drop"
  ) %>%
  mutate(
    shaps_total_s       = ave(shaps_total, FUN = scale),
    phq_total_s       = ave(phq_total, FUN = scale),
    gad_total_s       = ave(gad_total, FUN = scale),
    pss_total_s       = ave(pss_total, FUN = scale),
    drsp_total_s      = ave(drsp_total, FUN = scale),
    stress_total_s    = ave(stress_total, FUN = scale),
    ius_total_s       = ave(ius_total, FUN = scale),
    ius_prospective_s = ave(ius_prospective, FUN = scale),
    ius_inhibitory_s  = ave(ius_inhibitory, FUN = scale)
  ) -> per_participant_scaled


per_participant_scaled %>% glimpse

per_participant_scaled$group = case_when(per_participant_scaled$assigned_sex_clean == "male" ~ "male",
                                         per_participant_scaled$assigned_sex_clean == "female"  & per_participant_scaled$menstrual_cycle_clean=="Yes" ~ "female (cycling)",
                                         TRUE ~ "female (non-cycling)"
                                         )



per_trial_data_baseline_reduced = per_trial_data_baseline %>% select(prolific_id:MotivationQ, record_date : day_of_cycle)
per_trial_data_weekly_reduced = per_trial_data_weekly %>% select(prolific_id:MotivationQ, record_date : day_of_cycle)
per_trial_data_final_reduced = per_trial_data_final %>% select(prolific_id:MotivationQ, record_date : day_of_cycle)
length(per_trial_data_baseline_reduced)
length(per_trial_data_weekly_reduced)
length(per_trial_data_final_reduced)


per_trial_data_baseline_reduced$session = "baseline"
per_trial_data_weekly_reduced$session = "weekly"
per_trial_data_final_reduced$session = "finaly"

per_trial_data_all = rbind(per_trial_data_baseline_reduced,per_trial_data_weekly_reduced)
per_trial_data_all = rbind(per_trial_data_all,per_trial_data_final_reduced)
dim(per_trial_data_all)



per_trial_data_all = per_trial_data_all %>%  mutate(across(AttribQ : MotivationQ, ~ as.numeric(ave(., FUN = scale)), .names = "{.col}_s"))
per_trial_data_all = per_trial_data_all %>%  mutate(across(AttribQ : MotivationQ, ~ as.numeric(ave(., group =per_trial_data_all$SESSION_ID, FUN = scale)), .names = "{.col}_sc"))

per_trial_data_all = per_trial_data_all %>% merge(per_participant_scaled, by ="prolific_id")


per_trial_data_all = per_trial_data_all %>% 
  merge(per_trial_data_all %>% 
          group_by(prolific_id) %>% 
          summarize(baseline_record_date_dt = record_date_dt[session == "baseline"] %>% unique) %>% ungroup , 
        by ="prolific_id")
per_trial_data_all$day_from_baseline = per_trial_data_all$record_date_dt - per_trial_data_all$baseline_record_date_dt

per_trial_data_all$day_from_baseline %>% as.integer() %>% hist

per_trial_data_all$day_of_cycle %>% unique

per_trial_data_all$cycle_phase = case_when( (per_trial_data_all$day_of_cycle >0) & (per_trial_data_all$day_of_cycle <14) ~ "Follicular",
                                            (per_trial_data_all$day_of_cycle >13) & (per_trial_data_all$day_of_cycle <28) ~ "Luteal",
                                            TRUE ~ NA_character_)

# get move by move data;


move_data_temp = c()


sid = per_trial_data_all$SESSION_ID[1]
# collect move by move data 
for (sid in (per_trial_data_all$SESSION_ID %>% unique)) {
  per_trial_data_all_id = per_trial_data_all %>% filter(SESSION_ID == sid)
  # per_trial_data_weekly_id = per_trial_data_baseline %>% filter(prolific_id == pid)
  for (t in per_trial_data_all_id$trial_number) { 
    json_data<- 
      per_trial_data_all_id$key_presses[t] %>% 
      fromJSON( )  # there should be 8 variables in key_presses 
    
    
    if(length(json_data) >0) {
      # json_data = remove_two_in_triplets(json_data)
      # if (length(json_data)%% 8 != 0) {
      #   
      #   warning(paste("Uneven number of elements — maybe a corrupt or incomplete trial", t,"in pid", pid, "index i =", i))
      #   break
      # }
      # json_data = json_data[json_data != "stayed_still"]
      json_data = json_data %>%
        matrix(ncol = 9, byrow = TRUE)
      # json_data = json_data[json_data[,6] != "error",] # remove errors
      move_data_temp_tr <- as_tibble(json_data, .name_repair = "minimal") %>%
        set_names(c(
          "time_from_trial_start", "trial_duration", "yellow_duration", "turn_to_yellow_time", "RT",
          "outcome", "chosen_direction", "actual_direction", "outcome_type"
        )) %>%
        mutate(
          time_from_trial_start = as.numeric(time_from_trial_start),
          trial_duration = as.numeric(trial_duration),
          yellow_duration = as.numeric(yellow_duration),
          turn_to_yellow_time = as.numeric(turn_to_yellow_time),
          RT = as.numeric(RT),
          outcome = as.character(outcome),
          chosen_direction = as.character(chosen_direction),
          actual_direction = as.character(actual_direction),
          outcome_type = as.character(outcome_type)
        )
      move_data_temp_tr$SESSION_ID = sid
      move_data_temp_tr$trial = t
      move_data_temp_tr$move_no = dim(move_data_temp_tr)[1]
      move_data_temp_tr$move =1: dim(move_data_temp_tr)[1]
      
    }else {
      move_data_temp_tr <- tibble(
        time_from_trial_start = NA_real_,
        trial_duration = NA_real_,
        yellow_duration = NA_real_,
        turn_to_yellow_time =  NA_real_,
        RT = NA_real_,
        outcome = NA_character_,
        chosen_direction = NA_character_,
        actual_direction = NA_character_,
        outcome_type = NA_character_,
        SESSION_ID = sid,
        trial = t,
        move_no = NA_integer_,
        move = NA_integer_
      )
    }
    move_data_temp = rbind(move_data_temp, move_data_temp_tr)
  }
}

saveRDS(move_data_temp, "per_move_data_baseline.rds")



saveRDS(per_trial_data_baseline, "per_trial_data_baseline_4r.rds")
##### rest ###################################
########################################







# Extract parts with regex
files <- list.files(data_folder, pattern = "\\.csv$", full.names = TRUE)
info <- tibble(
  filepath = files,
  filename = basename(files)
) %>%
  separate(filename, into = c("prolific_pid", "session_id", "study_id_ext"), sep = "_") %>%
  mutate(study_id = str_remove(study_id_ext, "\\.csv$")) %>%
  select(-study_id_ext)


excl_ID = c();
trial_data = data.frame()
move_data = data.frame()
pid_prev = -999
for (i in 1:length(info$csv_files)) {
  if (!(i %in% excl_ID)) {
    
    print(i)
    pid <- info$prolific_pid[i]
    
    
    data_id <- read_csv(files[i], show_col_types = FALSE)
    no_trials = sum(!is.na(data_id$trial_start_time))
    # 
    # Initialize an empty data frame with the correct number of rows (8) and relevant columns
    data_id_temp <- data.frame(
      index = numeric(no_trials)#,
      # ID = numeric(8),
      # trial = numeric(8),
      # stressQ = numeric(8),
      # controlQ = numeric(8),
      # attribQ = numeric(8),
      # arousalQ = numeric(8),
      # valenceQ = numeric(8),
      # control_level = numeric(8),
      # outcome_total = numeric(8),
      # outcome = numeric(8)
    )
    
    # Filter rows with 'html-slider-response' and 'image-slider-response'
    html_slider_data<- subset(data_id, trial_type %in% c("html-slider-response", "image-slider-response"))
    if (dim(html_slider_data)[1] >  0 ) {
      
      
      # Filter rows for 'ctl-task-cpt', and limit the selection to specific rows if applicable
      trial_data_temp <- subset(data_id, trial_type == "ctl-task-cpt")
      
      if (dim(trial_data_temp)[1] == 8) {
        
        data_id_temp$index =  i
        
        data_id_temp$ID = pid
        
        data_id_temp$datetime = info$datetime[i]
        data_id_temp$trial  =1:8
        data_id_temp$session = session
        data_id_temp$session_type = NA
        
        
        
        
        
        
        # Extract and match responses for stress, control, merit, arousal, and valence
        data_id_temp$stressQ <- as.numeric(html_slider_data$response[grepl("stressed", html_slider_data$stimulus)])
        data_id_temp$controlQ <- as.numeric(html_slider_data$response[grepl("follow", html_slider_data$stimulus)])
        data_id_temp$attribQ <- as.numeric(html_slider_data$response[grepl("result", html_slider_data$stimulus)])
        data_id_temp$arousalQ <- as.numeric(html_slider_data$response[grepl("arousal", html_slider_data$stimulus)])
        data_id_temp$valenceQ <- as.numeric(html_slider_data$response[grepl("valence", html_slider_data$stimulus)])
        data_id_temp$control_level = trial_data_temp$control_level%>% as.numeric()
        data_id_temp$outcome = trial_data_temp$outcome%>% as.numeric()
        data_id_temp$outcome_total = trial_data_temp$outcome_total%>% as.numeric()
        data_id_temp$goal_achieved =  data_id_temp$outcome>0
        
        data_id_temp$goal_selected = trial_data_temp$goal_selected %>% as.numeric()
        data_id_temp$distance_large = trial_data_temp$goal_1_distance == 4
        data_id_temp$goal_1_value = trial_data_temp$goal_1_value
        data_id_temp$goal_2_value = trial_data_temp$goal_2_value
        # data_id_temp$ = trial_data_tempcounts
        
        move_data_temp = c()
        for (t in 1:8) { 
          json_data<- 
            trial_data_temp$key_presses[t] %>% 
            fromJSON( )  # there should be 8 variables in key_presses 
          
          
          if(length(json_data) >0) {
            json_data = remove_two_in_triplets(json_data)
            if (length(json_data)%% 8 != 0) {
              
              warning(paste("Uneven number of elements — maybe a corrupt or incomplete trial", t,"in pid", pid, "index i =", i))
              break
            }
            # json_data = json_data[json_data != "stayed_still"]
            json_data = json_data %>%
              matrix(ncol = 8, byrow = TRUE)
            json_data = json_data[json_data[,5] != "error",] # remove errors
            move_data_temp_tr <- as_tibble(json_data, .name_repair = "minimal") %>%
              set_names(c(
                "time_from_trial_start", "trial_duration", "yellow_duration", "RT",
                "outcome", "chosen_direction", "actual_direction", "outcome_type"
              )) %>%
              mutate(
                time_from_trial_start = as.numeric(time_from_trial_start),
                trial_duration = as.numeric(trial_duration),
                yellow_duration = as.numeric(yellow_duration),
                RT = as.numeric(RT),
                outcome = as.character(outcome),
                chosen_direction = as.character(chosen_direction),
                actual_direction = as.character(actual_direction),
                outcome_type = as.character(outcome_type)
              )
            move_data_temp_tr$trial = t
            move_data_temp_tr$session = session
            move_data_temp_tr$move_no = dim(move_data_temp_tr)[1]
            move_data_temp_tr$move =1: dim(move_data_temp_tr)[1]
            
          }else {
            move_data_temp_tr <- tibble(
              time_from_trial_start = NA_real_,
              trial_duration = NA_real_,
              yellow_duration = NA_real_,
              RT = NA_real_,
              outcome = NA_character_,
              chosen_direction = NA_character_,
              actual_direction = NA_character_,
              outcome_type = NA_character_,
              trial = t,
              session = session,
              move_no = NA_integer_,
              move = NA_integer_
            )
          }
          move_data_temp = rbind(move_data_temp, move_data_temp_tr)
          
          
        }
        
        move_data_2_trial_data = move_data_temp %>% group_by(trial) %>% summarize(mean_logRT= mean(log(RT), na.rm =T),
                                                                                  move_no = move_no %>% first,
                                                                                  prop_good_moves = sum(outcome_type == "good_move")/move_no,
                                                                                  prop_shook = sum(outcome == "oops")/move_no,
                                                                                  no_lucky_good_moves = sum(outcome == "oops" & outcome_type == "good_move")
        ) %>% ungroup()
        
        
        
        
        extracted_data <- t(sapply(trial_data_temp$counts, function(x) {
          json_data <- fromJSON(x)
          return(c(
            num_key_presses = json_data$num_key_presses,
            num_perfect = json_data$num_perfect,
            num_errors = json_data$num_errors,
            num_randdir = json_data$num_randdir,
            num_whilered = json_data$num_whilered
          ))
        }))
        # Convert the matrix to a data frame and ensure columns are numeric
        extracted_data <- as.data.frame(extracted_data)
        extracted_data[] <- lapply(extracted_data, as.numeric) 
        row.names(extracted_data) <- NULL
        extracted_data$number_valid_key_presses = extracted_data$num_key_presses - extracted_data$num_errors
        extracted_data$num_proportion_random =  extracted_data$num_randdir/extracted_data$number_valid_key_presses
        extracted_data$num_proportion_good_by_chance =  extracted_data$num_randdir/extracted_data$number_valid_key_presses
        
        
        data_id_temp = cbind(data_id_temp,extracted_data)
        data_id_temp = cbind(data_id_temp,move_data_2_trial_data %>% select(-trial))
        
        move_data= rbind(move_data,move_data_temp)
        trial_data = rbind(trial_data, data_id_temp)
      }
    }
  }
}

# write to csv;
write.csv(file = "trial_data.csv", x =trial_data)

trial_data %>% ggplot(aes(x = num_proportion_random, y = stressQ)) + geom_point(alpha = 0.2) + geom_smooth(method = "lm") + theme_mininal()
trial_data
# prepare data ####
trial_data$control_level = trial_data$control_level %>% round(4)
trial_data = trial_data %>% filter(control_level != 0.8286, goal_selected !=-9999 )

# do some data checks #####
# successes and fails across control levels
trial_data %>% group_by(control_level, goal_selected, distance_large) %>% summarize(mean_goal_achieved= mean(goal_achieved)) %>% View

trial_data %>% group_by(index) %>% summarize(mean_goal_achieved= mean(goal_achieved)) %>% pull(mean_goal_achieved) %>% hist

trial_data %>% group_by(index) %>% summarize(sum_outcome= sum(outcome)) %>% pull(sum_outcome) %>% hist()

trial_data %>% group_by(index) %>% summarize(sum_outcome= sum(outcome)) %>% pull(sum_outcome) %>% max()
trial_data %>% group_by(index) %>% summarize(sum_outcome= sum(outcome)) %>% pull(sum_outcome) %>% mean()
trial_data %>% group_by(index) %>% summarize(sum_outcome= sum(outcome)) %>% pull(sum_outcome) %>% median()

trial_data %>% group_by(index) %>% summarize(sum_outcome= sum(outcome)) %>% pull(sum_outcome) %>% median()

ggplot(data = trial_data, aes(x=control_level, y = num_proportion_random)) +  geom_point() + geom_smooth(method= "lm")
ggplot(data = trial_data, aes(x=controlQ, y = num_proportion_random)) +  geom_point() + geom_smooth(method= "lm")


ggpairs(trial_data %>% filter(goal_achieved == TRUE) %>% select(controlQ,stressQ,attribQ,arousalQ,valenceQ), title = "Pairs Plot of Outcome Variables")

ggpairs(trial_data %>% filter(goal_achieved == FALSE) %>% select(controlQ,stressQ,attribQ,arousalQ,valenceQ), title = "Pairs Plot of Outcome Variables")


ggplot(data = trial_data %>% filter( num_proportion_random != 0), aes(x=num_proportion_random, y = stressQ)) +  geom_point() + geom_smooth(method= "lm")
trial_data %>% filter( num_proportion_random == 0) %>% pull(stressQ)  %>% hist() 
ggplot(data = trial_data %>% filter( num_proportion_random != 0), aes(x=num_proportion_random, y = stressQ)) +  geom_point() + geom_smooth(method= "lm")

# rescale variables;

trial_data$stressQ_sc = trial_data$stressQ %>% ave(group = trial_data$ID, FUN = scale)
trial_data$stressQ_s = trial_data$stressQ %>% ave( FUN = scale)
trial_data$attribQ_s = trial_data$attribQ %>% ave(FUN = scale)
trial_data$attribQ_sc = trial_data$attribQ %>% ave(group = trial_data$ID, FUN = scale)
trial_data$num_proportion_random_s = trial_data$num_proportion_random %>% ave(FUN = scale)
trial_data$controlQ_s = trial_data$controlQ %>% ave(FUN = scale)

trial_data %>% ggplot(aes(x=attribQ_sc, y = stressQ_sc)) + facet_wrap(~goal_achieved)+ geom_point() + geom_smooth(method= "lm")

trial_data %>% ggplot(aes(x=attribQ_sc, y = stressQ_sc)) + geom_point() + geom_smooth(method= "lm")


trial_data %>% ggplot(aes(x=trial, y = stressQ_sc)) +  geom_point() + geom_smooth(method= "lm")

trial_data %>% ggplot(aes(x=trial, y = stressQ_sc)) +  geom_point() + geom_smooth(method= "lm")
trial_data %>% filter() %>% pull(stressQ)  %>% hist() 
trial_data %>% ggplot(aes(x=goal_achieved, y = stressQ_sc)) +  geom_point() + geom_boxplot(alpha = 0.2)
trial_data %>% ggplot(aes(x=goal_achieved, y = attribQ)) +  geom_point() + geom_boxplot(alpha = 0.2)
trial_data %>% ggplot(aes(x=goal_achieved, y = attribQ)) +  geom_point() + geom_boxplot(alpha = 0.2)
contrasts(trial_data$goal_achieved)
model_attrib = lmer(data = trial_data, attribQ_s~ goal_achieved + (goal_achieved | ID))
model_attrib %>% summary()
model2 = lmer(data = trial_data, stressQ_s~ num_proportion_random_s + (num_proportion_random_s | ID))
model2 %>% summary()

model2 = lmer(data = trial_data, stressQ_s~ num_proportion_random_s +goal_achieved*attribQ_s +(num_proportion_random_s + +goal_achieved*attribQ_s| ID))
model2 %>% summary()

model2 = lmer(data = trial_data, stressQ_s~ goal_achieved*attribQ_s +(goal_achieved*attribQ_s| ID))
model2 %>% summary()

trial_data = trial_data %>% mutate(failed  = goal_achieved == F)
model2 = lmer(data = trial_data, stressQ_s~ failed*controlQ_s +(failed*controlQ_s| ID))
model2 %>% summary()

model2 = lmer(data = trial_data, stressQ_s~ failed*attribQ_s +(failed*attribQ_s| ID))
model2 %>% summary()

model2 = lmer(data = trial_data, controlQ_s~ failed*attribQ_s +(failed*attribQ_s| ID))
model2 %>% summary()



current_score = 2
total_possible_score = 130
0.03
x #trial_data.participant_final_score
y #trial_data.participant_final_score
x0 = 1
x = 5
# max(min(
y = 
  x/ - 0.03 * (x - (x0+ 100*(current_score/total_possible_score)))

# , 100), 0);


# REST #################
# Inspect the first few rows to check the structure of your data
head(data)
trialdata = as.data.frame(trialdata$the_matrix);
names(trialdata) = c('subjectID', 'session', 'trial', 'control', 'gain1', 'loss1', 'gain2', 'loss2', 
                     'dist1', 'dist2', 'canvasRotation', 'goalChosen', 'outcome', 'choiceRT', 'lastPressTime',
                     'lastPressDistRem', 'numOfTaps', 'numNotIntended', 'numGood', 'numBad', 'numWrong',
                     'stressQ', 'controlQ', 'attribQ', 'regretQ');
trialdata <- trialdata %>% select(subjectID:regretQ) # remove rows with no names
trialdata$trialtype = factor(trialdata$gain1 > 0, labels = c("LoseOnly","WinOnly"));
# correct numBad and numGood
# first trial likely misclassified
trialdata$numBad = trialdata$numBad -1 
trialdata$numGood = trialdata$numGood +1

# filter trials that are useless:
trialdata[trialdata == -99] = NA;

# trialdata <- trialdata %>% filter(numOfTaps > 0, numNotIntended >= numBad )



trialdata = trialdata %>% mutate(goal_achieved = case_when(loss1 >0 & outcome == 0 ~ 1,
                                                           loss1==0 & outcome > 0 ~1,
                                                           TRUE ~ 0),
                                 goal_achieved_f = factor(goal_achieved, 
                                                          levels = c(0,1), labels= c('unsuccessful', 'successful')),
                                 session = factor(session),
                                 subjectID = as.integer(subjectID),
                                 trial = as.integer(trial))


trialdata <- trialdata %>% mutate(stressQ_s = scale(stressQ),
                                  attribQ_s = scale(attribQ),
                                  controlQ_s = scale(controlQ),
                                  control_s = scale(control))

# correct numBad and define numGood accordinly (to match better the press by press data )


trialdata$dist_diff =  trialdata$dist2 - trialdata$dist1 
trialdata$dist_diff = trialdata$dist_diff/max(trialdata$dist_diff, na.rm = T)
trialdata$gain_diff = trialdata$gain2 - trialdata$gain1 
trialdata$gain_diff = trialdata$gain_diff/max(trialdata$gain_diff, na.rm = T)
trialdata$loss_diff = trialdata$loss2 - trialdata$loss1 
trialdata$loss_diff = trialdata$loss_diff/min(trialdata$loss_diff, na.rm = T)



trialdata <- trialdata %>% mutate(goal_achieved_lag =lag(goal_achieved),
                                  numGoodNotIntended = numNotIntended - numBad,
                                  numIntended = numOfTaps - numNotIntended)



trialdata <- trialdata %>% mutate(prop_followed_command = numIntended/numOfTaps,
                                  prop_good_moves = numGood/numOfTaps,
                                  prop_good_move_by_chance =numGoodNotIntended/numOfTaps)
# trialdata$prop_good_move_by_chance %>% unique()

trialdata <- trialdata %>% mutate(prop_followed_command_s = prop_followed_command %>% ave(FUN = scale),
                                  prop_good_moves_s = prop_good_moves %>% ave(FUN = scale),
                                  prop_good_move_by_chance_s = prop_good_move_by_chance %>% ave(FUN = scale),
                                  prop_followed_command_sc = prop_followed_command %>% ave(group = subjectID, FUN = scale),
                                  prop_good_move_by_chance_sc = prop_good_move_by_chance %>% ave(group = subjectID, FUN = scale),
                                  prop_good_moves_sc = prop_good_moves %>% ave(group = subjectID, FUN = scale))

# pdi data #####

pdi_scores <- read.csv ('Data/questionnaire data trt/pdi_scores.csv', header = T, sep="," )
pdi_scores <- pdi_scores %>% mutate(subjectID = participant)
pdi_scores <- pdi_scores %>% mutate(pdi_t = as.vector(ave(logit((pdi_total-0.5)*0.95+ 0.5), FUN = scale )),
                                    ID = factor(participant))
trialdata <- merge(trialdata, pdi_scores, by = "subjectID") 





# loc data
load('Data/questionnaire data trt/LOC/LOC.Rda')
loc_scores <- df4
trialdata$Internality = 0;
trialdata$Chance = 0;
trialdata$PowerfulOthers = 0;
for (sub in unique(trialdata$subjectID)) {
  # print(sub)
  trialdata[trialdata$subjectID == sub,]$Internality = df4[df4$LC02_01 == sub,]$Internal_Locus_of_Control;
  trialdata[trialdata$subjectID == sub,]$Chance = df4[df4$LC02_01 == sub,]$Chance;
  trialdata[trialdata$subjectID == sub,]$PowerfulOthers = df4[df4$LC02_01 == sub,]$Powerful_Others;
}

trialdata$intCluster = 1*(trialdata$Internality > 35) + 1*(trialdata$Internality > 29);

ID_scores <- loc_scores[,"LC02_01"]

load('Data/questionnaire data trt/PSC/PSC.Rda')
psc_scores <- cbind(ID_scores,df3)
# View(psc_scores)
# load('../../Data/questionnaire data trt/BIS_11/BIS_11.Rda')
# bis11_scores <-cbind(ID_scores, df2)

# load('../../Data/questionnaire data trt/BIS_BAS/BIS_BAS.Rda')
# bis_bas_scores <- cbind(ID_scores,df1)

load('Data/questionnaire data trt/CTI/CTI.Rda')
cti_scores <- cbind(ID_scores,df)

load('Data/questionnaire data trt/AS/AS.Rda')
as_scores <- cbind(ID_scores,df6)
# pdi_scores$participant = as.factor(pdi_scores$participant)
# View(data_all)
# psc_scores %>% glimpse()
# cti_scores %>% glimpse
# data_all <- left_join(data_all, arrange(pdi_scores,participant), by = c("ID"= "participant")) 
trialdata <- left_join(trialdata, arrange(psc_scores %>% select(ID_scores, Total),ID_scores), by = c("subjectID"= "ID_scores")) 
trialdata <- left_join(trialdata, arrange(cti_scores %>% select(ID_scores, total_score),ID_scores), by = c("subjectID"= "ID_scores")) %>% rename(cti_total = total_score)

trialdata <- left_join(trialdata, arrange(as_scores %>% select(ID_scores, Attachment_Anxiety, Attachment_Avoidance),ID_scores), by = c("subjectID"= "ID_scores")) 

# trialdata %>% glimpse

# psc cluster
vec_temp <- (trialdata %>% group_by(subjectID) %>% summarize(psc = Total[1]) %>% select(psc)) # %>%quantile(probs = c(0.33, 0.66))
# quantile(vec_temp[["psc"]], probs = c(0.33,0.5, 0.66))



trialdata <- trialdata %>% rename(psc_total = Total)
trialdata$pscCluster = factor(1*(trialdata$psc_total > 32) + 1*(trialdata$psc_total > 25) ,levels = c(0,1,2), labels = c("low", "mid", "high"))
trialdata$pscCluster2 = factor(trialdata$psc_total > 30, levels = c(FALSE,TRUE), labels = c("low", "high"))
trialdata$goalChosen01 =  trialdata$goalChosen -1

shrink <- function(a, s =0.05) (a -0.5)*(1-s) + 0.5
# trialdata$control
trialdata <- trialdata %>% mutate(controlQ_prev = lag(shrink(controlQ/100)%>% logit) ,
                                  control_prev = lag(control), 
                                  goal_achieved_prev_f = lag(goal_achieved_f),
                                  attribQ_prev = lag(shrink(attribQ/100)%>% logit) ,
                                  goalChosen_next = as.factor(lead(goalChosen01)),
                                  dist_diff =dist2 -dist1,
                                  controlQ_prev_s = controlQ_prev %>% ave(FUN = scale),
                                  attribQ_prev_s = attribQ_prev %>% ave(FUN = scale),
                                  dist_diff_s = dist_diff%>% ave(FUN = scale),
                                  controlQ_s = controlQ%>% ave(FUN = scale),
                                  control_s = control%>% ave(FUN = scale),
                                  attribQ_s = attribQ%>% ave(FUN = scale),
                                  stressQ_s = stressQ%>% ave(FUN = scale),
                                  controlQ_prev_sc = controlQ_prev %>% ave(FUN = scale),
                                  attribQ_prev_sc = attribQ_prev %>% ave(group = subjectID, FUN = scale),
                                  dist_diff_sc = dist_diff%>% ave(group = subjectID,FUN = scale),
                                  controlQ_sc = controlQ%>% ave(group = subjectID,FUN = scale),
                                  control_sc = control%>% ave(group = subjectID,FUN = scale),
                                  stressQ_sc = stressQ%>% ave(group = subjectID,FUN = scale),
                                  attribQ_sc = attribQ%>% ave(group = subjectID,FUN = scale))

saveRDS(trialdata, file = "Data/trialdata.rds")
# prepare trialdata with effort data
pbyp_dataset = read.csv("Data/pbyp_dataset_v1.csv") # read.csv("/Users/Fede/Desktop/Research/__projects/ctl-task-2.0/Data/pbyp_dataset_v1.csv")
pbyp_dataset$rt[pbyp_dataset$rt < 200] <- NA
pbyp_dataset <- pbyp_dataset %>% filter(!is.na(move)) 
pbyp_dataset$log_rt = log10(pbyp_dataset$rt);
pbyp_dataset$log_rt_zscored = 0;
pbyp_dataset$rt_zscored = 0;
pbyp_dataset$sub_index = NA;
pbyp_dataset$session = as.factor(pbyp_dataset$session)

pbyp_dataset$followed_command = (pbyp_dataset$sel_dir == pbyp_dataset$true_dir)
pbyp_dataset$good_move = (pbyp_dataset$current_distance - pbyp_dataset$distance_rem) >0 
pbyp_dataset$bad_move =(pbyp_dataset$current_distance - pbyp_dataset$distance_rem) <= 0 
pbyp_dataset$last_move = c(pbyp_dataset$move %>% diff() < 0,TRUE)

pbyp_dataset %>% filter(move == 2, trial == 2) %>% group_by(session) %>% summarize(N=n())
trialdata %>% filter(trial == 2) %>% group_by(session) %>% summarize(N=n())

# pbyp_dataset$good_move[pbyp_dataset$last_move == 1] = pbyp_dataset$followed_command[pbyp_dataset$last_move == 1]
pbyp_dataset$good_move_by_chance = (pbyp_dataset$good_move == 1 & (!pbyp_dataset$followed_command))


pbyp_dataset_sum <- pbyp_dataset %>% group_by(subjectID, trial, session) %>% summarize(N = n(),
                                                                                       mean_log_rt = mean(log_rt, na.rm = T))

# numOfTaps_pbyp = N,
# numBad_pbyp = sum(bad_move, na.rm = T),
# numGood_pbyp = sum(good_move, na.rm = T),
# numGood_byChance_pbyp = sum(good_move_by_chance, na.rm = T),
# prop_good_moves_pbyp = sum(good_move, na.rm = T)/N,
# prop_good_move_by_chance_pbyp = sum(good_move_by_chance, na.rm = T)/N,
# prop_followed_command_pbyp = sum(followed_command, na.rm = T)/N)






trialdata_pbyp <- merge(pbyp_dataset_sum, trialdata)
trialdata_pbyp %>% glimpse

(trialdata_pbyp$numOfTaps == trialdata_pbyp$numOfTaps_pbyp ) %>% mean()
trialdata_pbyp <- trialdata_pbyp %>% mutate(mean_log_rt_s = mean_log_rt %>% ave(FUN = scale),
                                            mean_log_rt_sc = mean_log_rt %>% ave(group = subjectID, FUN = scale))
# prop_followed_command_s = prop_followed_command %>% ave(FUN = scale),
# prop_good_moves_s = prop_good_moves %>% ave(FUN = scale),
# prop_good_move_by_chance_s = prop_good_move_by_chance %>% ave(FUN = scale),
# prop_followed_command_sc = prop_followed_command %>% ave(group = subjectID, FUN = scale),
# prop_good_move_by_chance_sc = prop_good_move_by_chance %>% ave(group = subjectID, FUN = scale),
# prop_good_moves_sc = prop_good_moves %>% ave(group = subjectID, FUN = scale))


# trialdata_pbyp %>% glimpse

# trialdata_pbyp %>% glimpse
# ((trialdata_pbyp$numBad_pbyp - (trialdata_pbyp$numBad) )== 0) %>% mean
# ((trialdata_pbyp$numGood_pbyp - (trialdata_pbyp$numGood) )== 0) %>% mean

saveRDS(trialdata, file = "Data/trialdata.rds")
saveRDS(trialdata_pbyp, file = "Data/trialdata_pbyp.rds")
saveRDS(pbyp_dataset, file = "Data/pbyp_dataset.rds")

# add qs to stan data 


data_for_stan <- readRDS("~/mnt/p/userdata/mikusn22/data/LOC_Amisul_naltrexone_2019_2020/Control Task/control-task2/Scripts/ctl-task2-R/data4stan_winloss_sa02_nosa02_lowsa02.rds")
# data_for_stan <- readRDS("~/mnt/p/userdata/mikusn22/data/LOC_Amisul_naltrexone_2019_2020/Control Task/control-task2/Scripts/ctl-task2-R/data4stan_rho.rds")

# # loc data

load('Data/questionnaire data trt/LOC/LOC.Rda')
loc_scores <- df4

data_for_stan$Internality = 0;
data_for_stan$Chance = 0;
data_for_stan$PowerfulOthers = 0;
for (sub in unique(data_for_stan$subjectID)) {
  # print(sub)
  data_for_stan[data_for_stan$subjectID == sub,]$Internality = df4[df4$LC02_01 == sub,]$Internal_Locus_of_Control;
  data_for_stan[data_for_stan$subjectID == sub,]$Chance = df4[df4$LC02_01 == sub,]$Chance;
  data_for_stan[data_for_stan$subjectID == sub,]$PowerfulOthers = df4[df4$LC02_01 == sub,]$Powerful_Others;
}

data_for_stan$intCluster = 1*(data_for_stan$Internality > 35) + 1*(data_for_stan$Internality > 29);

ID_scores <- loc_scores[,"LC02_01"]

data_for_stan$Internality_s = data_for_stan$Internality %>% ave(FUN = scale)
data_for_stan$Chance_s = data_for_stan$Chance %>% ave(FUN = scale)
data_for_stan$PowerfulOthers_s = data_for_stan$PowerfulOthers %>% ave(FUN = scale)


# perceived stress score
load('Data/questionnaire data trt/PSC/PSC.Rda')
psc_scores <- cbind(ID_scores,df3)
psc_scores$psc_total_s = psc_scores$Total %>% ave(FUN = scale)
data_for_stan <- left_join(data_for_stan, arrange(psc_scores %>% select(ID_scores, psc_total_s),ID_scores), by = c("subjectID"= "ID_scores")) 


# bis 11

load('Data/questionnaire data trt/BIS_11/BIS_11.Rda')

df2 = df2 %>% filter(CASE %in% df3$CASE)
bis11_scores <-cbind(ID_scores, df2)

bis11_scores$Bis_total_s = bis11_scores$Bis_total %>% ave(FUN = scale)
bis11_scores$Bis_Self_Control_s = bis11_scores$Self_Control_1st_order  %>% ave(FUN = scale)

data_for_stan <- left_join(data_for_stan, arrange(bis11_scores %>% select(ID_scores, Bis_total_s, Bis_Self_Control_s),ID_scores), by = c("subjectID"= "ID_scores"))


#bis bas


load('Data/questionnaire data trt/BIS_BAS/BIS_BAS.Rda')

df1 = df1 %>% filter(CASE %in% df3$CASE)
bis_bas_scores <- cbind(ID_scores,df1)
bis_bas_scores$BAS_drive_s = bis_bas_scores$BAS_drive %>% ave(FUN = scale)
bis_bas_scores$BAS_funseeking_s = bis_bas_scores$BAS_funseeking %>% ave(FUN = scale)
bis_bas_scores$BAS_rewardresponsivness_s = bis_bas_scores$BAS_rewardresponsivness %>% ave(FUN = scale)
bis_bas_scores$BIS_s = bis_bas_scores$BIS %>% ave(FUN = scale)

data_for_stan <- left_join(data_for_stan, arrange(bis_bas_scores %>% select(ID_scores, BAS_drive_s, BAS_funseeking_s, BAS_rewardresponsivness_s, BIS_s),ID_scores), by = c("subjectID"= "ID_scores")) 


# cti traums
load('Data/questionnaire data trt/CTI/CTI.Rda')
cti_scores <- cbind(ID_scores,df)
cti_scores$cti_total_s = cti_scores$total_score %>% ave(FUN = scale)
data_for_stan <- left_join(data_for_stan, arrange(cti_scores %>% select(ID_scores, cti_total_s),ID_scores), by = c("subjectID"= "ID_scores")) 

# attachment scores

load('Data/questionnaire data trt/AS/AS.Rda')
as_scores <- cbind(ID_scores,df6)
as_scores$Attachment_Anxiety_s = as_scores$Attachment_Anxiety %>% ave(FUN = scale)
as_scores$Attachment_Avoidance_s = as_scores$Attachment_Avoidance %>% ave(FUN = scale)

data_for_stan <- left_join(data_for_stan, arrange(as_scores %>% select(ID_scores, Attachment_Anxiety_s, Attachment_Avoidance_s),ID_scores), by = c("subjectID"= "ID_scores")) 


# pdi_scores$participant = as.factor(pdi_scores$participant)
# View(data_all)
# psc_scores %>% glimpse()
# cti_scores %>% glimpse
# data_all <- left_join(data_all, arrange(pdi_scores,participant), by = c("ID"= "participant")) 

saveRDS(data_for_stan, file = "Data/data_for_stan_wQs.rds")
saveRDS(data_for_stan, file = "Data/data_for_stan_wQs_loss_bias.rds")

