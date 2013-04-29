popmenu job_title_pm                file jobs.ms
popdown job_title_ok                file jobs.ms
popmenu marc_input_style_pm         file job_common.ms
popmenu marc_input_style_run_adv_pm file job_common.ms
popmenu marc_version_run_pm         file job_common.ms
group job_solver_gr                 file job_common.ms
group user_domains_gr               file domain_decomposition.ms
group user_domains_generate_gr      file domain_decomposition.ms
group user_domains_tail_gr          file domain_decomposition.ms
group job_solver_opts_gr            file job_common.ms
popmenu ddm_options                 file job_common.ms

popmenu job_run_popmenu {

 group {

  label {
   position 0 0
   size 50 4
   text "RUN JOB"
  }

  button {
   position 1 +5
   size 24 4
   text "USER SUBROUTINE FILE"
   popmenu usersub_file_browser_pm
   set $usersub_file_browser_label "SELECT USER SUBROUTINE FILE"
   set $usersub_file_browser_command "*job_usersub_file"
   help job_usersub_file
  }

  toggle {
   position +26 =
   size 22 4
   text "SELECTED USER SUBS"
   toggle job_usersubs
   help job_run_seluser
   visible job_usersubs
   popmenu job_usersub_pm
  }

  display {
   position 1 +4
   size 48 4
   display job_usersub_file
  }

  button {
   position 1 +4
   size 12 4
   text "EDIT"
   command "*job_edit_usersub_file"
   visible job_usersub_file
  }

  button {
   position +12 =
   size 12 4
   text "CLEAR"
   command "*job_clear_usersub_file"
   visible job_usersub_file
  }

  roller {
   position +12 =
   size 24 4
   nvalues 3
   texts "COMPILE / NO SAVE"
         "COMPILE AND SAVE"
         "RUN SAVED EXECUTABLE"
   help job_run_compile
   roller "job_option user_source"
   visible job_usersub_file
   popmenu job_user_source_pm
  }

  label {
   position 1 +6
   size 48 10
   border_width 1
   border_color black
  }

  toggle {
   position +1 +1
   size 23 4
   text "PARALLELIZATION"
   help job_run_parallel
   toggle "or(*job_option parallel:on, \
              solver_multi_procs, \
              solver_multi_threads)"
   popmenu job_run_parallelization_pm
   set $popname2 "job_run_parallelization_pm"
  }

  label {
   position +23 =
   size 23 8
   border_width 1
   border_color black
  }

  display {
   position = =
   size 23 4
   display job_ddm_domains
  }

  display {
   position = +4
   size 23 4
   display job_solver_procs
   visible solver_allows_multi_procs
  }

  display {
   position = =
   size 23 4
   display job_solver_threads
   visible "and(solver_allows_multi_threads,\
                not(and(job_solver_mfront_sparse,\
                        *job_option parallel:on)))"
  }

  button {
   position 1 +7
   size 8 4
   text "TITLE"
   popmenu job_title_pm
   command "*job_title"
  }

# see also job_common.ms
# see also the ADVANCED JOB SUBMISSION popmenu in this file

  label {
   position +10 =
   size 7 4
   text "STYLE"
   border_width 1
   border_color black
  }

  roller {
   position +7 =
   size 16 4
   nvalues 3
   texts "TABLE-DRIVEN"
         "MULTI-PHYSICS"
         "OLD"
   rollers "job_input_style_table_driven"
           "job_input_style_multi_physics"
           "job_input_style_old"
   popmenu marc_input_style_pm
   help job_option_input_style
  }

  button {
   position +18 =
   size 13 4
   text "SAVE MODEL"
   command "*save_model"
  }

  button {
   position 1 +6
   size 16 6
   text "SUBMIT (1)"
   command "*submit_job 1 *monitor_job"
  }

  button {
   position +16 =
   size 16 6
   text "ADV. JOB SUBM."
   popmenu job_submit_adv_pm
  }

  button {
   position +16 =
   size 16 6
   text "DAMASK"
   popmenu damask
  }

  button {
   position 1 +6
   size 16 6
   text "UPDATE"
   command "*update_job"
  }

  button {
   position +16 =
   size 16 6
   text "MONITOR"
   command "*monitor_job"
  }

  button {
   position +16 =
   size 16 6
   text "KILL"
   command "*kill_job *monitor_job"
  }

  label {
   position 1 +7
   size 32 4
   text "STATUS"
   border_width 1
   border_color black
  }

  display {
   position +32 =
   size 16 4
   display "job_status"
  }

  label {
   position -32 +4
   size 32 4
   text "CURRENT INCREMENT (CYCLE)"
   border_width 1
   border_color black
  }

  display {
   position +32 =
   size 16 4
   display "job_increment"
  }

  label {
   position -32 +4
   size 32 4
   text "SINGULARITY RATIO"
   border_width 1
   border_color black
  }

  float {
   position +32 =
   size 16 4
   display "job_singularity_ratio"
  }

  label {
   position -32 +4
   size 32 4
   text "CONVERGENCE RATIO"
   border_width 1
   border_color black
  }

  float {
   position +32 =
   size 16 4
   display "job_convergence_ratio"
  }


  label {
   position 1 +4
   size 32 4
   text "ANALYSIS TIME"
   border_width 1
   border_color black
  }

  float {
   position +32 =
   size 16 4
   display "job_analysis_time"
  }

  label {
   position -32 +4
   size 32 4
   text "WALL TIME"
   border_width 1
   border_color black
  }

  float {
   position +32 =
   size 16 4
   display "job_time"
  }

  frame {
   position 1 +4
   size 48 8

   group {
    label {
     position 0 0
     size 48 8 
     text "TOTAL"
     border_width 1
     border_color black
    }
    
    label {
     position +6 =
     size 12 4 
     text "CYCLES"
     border_width 1
     border_color black
    }
    
    integer {
     position +12 =
     size 10 4
     display "acc_job_cycles"
     border_width 1
     border_color black
    }

    label {
     position -12 +4
     size 12 4 
     text "SEPARATIONS"
     border_width 1
     border_color black
    }
 
    integer {
     position +12 =
     size 10 4
     display "acc_job_separations"
     border_width 1
     border_color black
    }

    label {
     position +10 -4 
     size 10 4 
     text "CUT BACKS"
     border_width 1
     border_color black
    }
 
    integer {
     position +10 =
     size 10 4
     display "acc_job_cut_backs"
     border_width 1
     border_color black
    }

    label {
     position -10 +4
     size 10 4 
     text "REMESHES"
     border_width 1
     border_color black
    }

    integer {
     position +10 =
     size 10 4
     display "acc_job_remeshes"
     border_width 1
     border_color black
    }
   }
  }

  label {
   position 1 +8
   size 18 4
   text "EXIT NUMBER"
   border_width 1
   border_color black
  }

  integer {
   position +18 =
   size 10 4
   display "job_exit"
  }

  button {
   position +10 =
   size 20 4
   text "EXIT MESSAGE"
   popmenu job_exit_msg_pm
   help exit_message
  }

  label {
   position 1 +6
   size 5 4
   text "EDIT"
   border_width 1
   border_color black
  }

  button {
   position +5 = 
   size 12 4
   text "OUTPUT FILE"
   command "*job_edit_output"
  }

  button {
   position +12 =
   size 9 4
   text "LOG FILE"
   command "*job_edit_log_file"
  }

  button {
   position +9 =
   size 12 4
   text "STATUS FILE"
   command "*job_edit_status_file"
  }

  button {
   position +12 =
   size 10 4
   text "ANY FILE"
   set $edit_browser_label "EDIT FILE"
   set $edit_browser_command "*edit_file"
   popmenu edit_browser_popmenu
   help edit_file
  }

  popdown {
   position 1 +6
   size 32 4
   text "OPEN POST FILE (RESULTS MENU)"
   screen results
   command "*post_open_default"
  }

  button {
   position 1 +6
   size 12 8
   text "RESET"
   command "*job_submit_reset"
  }

  popdown {
   position +36 =
   size 12 8
   text "OK"
  }
 }

 window job_run_wi {
  parent mentat
  origin 35 4
  size 50 109
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

popmenu usersub_file_browser_pm {

 group {

  browser usersub_file_browser {
   position 0 0
   size 82 72
   text "$usersub_file_browser_label" ($usersub_file_browser_label)
   filter "*.f *.f90"
   nvisible 10
   select_files true
   cancel_action popdown
   ok_action popdown
   command "$usersub_file_browser_command" ($usersub_file_browser_command)
  }
 }

 window {
  title " Select MARC User Subroutine File"
  origin 30 30
  size 82 72
 }

 mode dialog
 disable_right_mouse_button
}

popmenu host_file_browser_pm {

 group {

  browser host_file_browser {
   position 0 0
   size 82 72
   text "$host_file_browser_label" ($host_file_browser_label)
   filter "*"
   nvisible 10
   select_files true
   cancel_action popdown
   ok_action popdown
   command "$host_file_browser_command" ($host_file_browser_command)
  }
 }

 window {
  title " Select MARC Host File"
  origin 30 30
  size 82 72
 }

 mode dialog
 disable_right_mouse_button
}

popmenu job_usersub_pm {

 group {

  label {
   position 0 0
   size 66 4
   text "CURRENTLY SELECTED USER SUBROUTINES"
  }

  display {
   position 1 +5
   size 64 86
   display "job_usersubs"
  }

  popdown {
   position 27 +88
   size 12 8
   text "OK"
  }
 }

 window {
  parent mentat
  origin 38 8
  size 66 102
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

popmenu job_network_pm {
 group {
  label {
   position 0 0
   size 50 4
   text "NETWORK SETTINGS"
  }

  button {
   position 1 +5
   size 24 4
   text "HOST FILE"
   popmenu host_file_browser_pm
   set $host_file_browser_label "SELECT HOST FILE"
   set $host_file_browser_command "*job_host_file"
   help job_host_file
  }

  button {
   position +24 =
   size 12 4
   text "EDIT"
   command "*job_edit_host_file"
   visible job_host_file
   help job_edit_host_file
  }

  button {
   position +12 =
   size 12 4
   text "CLEAR"
   command "*job_clear_host_file"
   visible job_host_file
   help job_clear_host_file
  }

  display {
   position 1 +4
   size 48 4
   display job_host_file
  }

  toggle {
   position 1 +5
   size 24 4
   text "COPY INPUT FILE"
   toggle "*job_option copy_input_file:on"
   true_command "*job_option copy_input_file:on"
   false_command "*job_option copy_input_file:off"
   help job_host_copy_inputfile
  }

  toggle {
   position +24 =
   size 24 4
   text "COPY POST FILE"
   toggle "*job_option copy_post_file:on"
   true_command "*job_option copy_post_file:on"
   false_command "*job_option copy_post_file:off"
   help job_host_copy_inputfile
  }

  label {
   position 1 +6
   size 48 4
   text "HOSTS"
   border_width 1
   border_color black
   visible job_usersub_file
  }

  oneonly {
   position +12 =
   size 18 4
   text "COMPATIBLE"
   oneonly "*job_option network_hosts:compatible"
   command "*job_option network_hosts:compatible"
   visible job_usersub_file
   help job_host_comp
  }

  oneonly {
   position +18 =
   size 18 4
   text "INCOMPATIBLE"
   oneonly "*job_option network_hosts:incompatible"
   command "*job_option network_hosts:incompatible"
   visible job_usersub_file
   help job_host_comp
  }

  popdown {
   position 19 +6
   size 12 8
   text "OK"
  }
 }

 window {
  parent mentat
  origin 38 8
  size 50 35
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

popmenu job_submit_adv_pm {
 group {
  label {
   position 0 0
   size 50 4
   text "ADVANCED JOB SUBMISSION"
  }

  label {
   position 1 +5
   size 19 4
   text "INITIAL ALLOCATION"
   border_width 1
   border_color black
  }
   
  button {
   position +19 =
   size 15 4
   text "GENERAL MEMORY"
   help job_param_general_init_allocation
   command "*job_param general_init_allocation"
  }

  text {
   position +15 =
   size 10 4
   display "job_param_value_general_init_allocation"
   command "*job_param general_init_allocation"
  }

  label {
   position +10 =
   size 4 4
   text "MB"
   border_width 1
   border_color black
  }

  toggle {
   position 1 +5
   size 32 4
   text "OUT-OF-CORE ELEMENT STORAGE"
   help job_param_elsto
   toggle "*job_option elsto:on"
   true_command "*job_option elsto:on"
   false_command "*job_option elsto:off"
  }

  toggle {
   position 1 +4
   size 32 4
   text "OUT-OF-CORE INCREMENTAL BACKUP"
   help job_param_inc_backup_storage
   toggle "*job_option inc_backup_storage:out_of_core"
   true_command "*job_option inc_backup_storage:out_of_core"
   false_command "*job_option inc_backup_storage:in_core"
  }

  toggle {
   position +34 =
   size 14 4
   text "CHECK SIZES"
   help job_run_check
   toggle "*job_option check:on"
   true_command "*job_option check:on"
   false_command "*job_option check:off"
  }

  label {
   position 1 +6
   size 16 4
   text "INTEGER SIZE"
   border_width 1
   border_color black
  }

  roller {
   position +16 =
   size 32 4
   nvalues 3
   texts "DEFAULT"
         "4-BYTE INTEGERS"
         "8-BYTE INTEGERS"
   popmenu marc_integer_size_pm
   help job_run_intsize
   roller "job_option marc_integer_size"
  }

  label {
   position 1 +6
   size 48 4
   text "MARC INPUT FILE"
   border_width 1
   border_color black
  }

  label {
   position 1 +4
   size 9 4
   text "VERSION"
   border_width 1
   border_color black
  }

  roller {
   position +9 =
   size 14 4
   nvalues 14
   texts
    "DEFAULT"
    "2010"
    "2008"
    "2007"
    "2005R3"
    "2005"
    "2003"
    "2001"
    "2000"
    "8"
    "K7"
    "K6.2"
    "K5"
    "K4"
   help job_param_version
   rollers
    "job_input_version_default"
    "job_input_version_2010"
    "job_input_version_2008"
    "job_input_version_2007"
    "job_input_version_2005r3"
    "job_input_version_2005"
    "job_input_version_2003"
    "job_input_version_2001"
    "job_input_version_2000"
    "job_input_version_8"
    "job_input_version_k7"
    "job_input_version_k6"
    "job_input_version_k5"
    "job_input_version_k4"
   popmenu marc_version_run_pm
  }

# see also job_common.ms
# see also the RUN JOB popmenu in this file

  label {
   position +14 =
   size 7 4
   text "STYLE"
   border_width 1
   border_color black
  }

  roller {
   position +7 =
   size 18 4
   nvalues 3
   texts "TABLE-DRIVEN"
         "MULTI-PHYSICS"
         "OLD"
   rollers "job_input_style_table_driven"
           "job_input_style_multi_physics"
           "job_input_style_old"
   popmenu marc_input_style_run_adv_pm
   help job_option_input_style
  }

  toggle {
   position 1 +4
   size 48 4
   text "EXTENDED PRECISION"
   help job_run_precision
   toggle "*job_option inp_file_prec:extended"
   true_command "*job_option inp_file_prec:extended"
   false_command "*job_option inp_file_prec:normal"
  }

  button {
   position 1 +6
   size 24 4
   text "SCRATCH DIRECTORY"
   set $directory_browser_label "JOB SCRATCH DIRECTORY"
   set $directory_browser_command "*job_scratch_directory"
   popmenu directory_browser_popmenu
   help job_scratch_directory
  }

  button {
   position +24 =
   size 24 4
   text "CLEAR"
   command "*job_clear_scratch_directory"
   visible job_scratch_directory
  }

  text {
   position 1 +4
   size 48 4
   display job_scratch_dir
   command "*job_scratch_directory"
  }

#ifdef DCOM
  toggle {
   position 1 +6
   size 8 4
   text "DCOM"
   toggle "*job_option dcom:on"
   help job_run_dcom
   true_command "*job_option dcom:on"
   false_command "*job_option dcom:off"
   visible win32_available
  }

  button {
   position +8 =
   size 12 4
   text "HOSTNAME"
   command "*job_dcom_hostname"
   visible "and(win32_available, *job_option dcom:on)"
  }
   
  text job_dcom_hostname {
   position +12 =
   size 28 4
   display "job_dcom_hostname"
   command "*job_dcom_hostname"
   visible "and(win32_available, *job_option dcom:on)"
  }
#endif

  button {
   position 1 +6
   size 24 4
   text "TITLE"
   popmenu job_title_pm
   command "*job_title"
  }

  button {
   position +24 =
   size 24 4
   text "SAVE MODEL"
   command "*save_model"
  }

  button {
   position +4 +6
   size 20 4
   text "WRITE INPUT FILE"
   command "*job_write_input"
  }

  button {
   position = +4
   size 20 4
   text "EDIT INPUT FILE"
   command "*job_edit_input"
  }

  popdown {
   position 1 +5
   size 20 6
   text "SUBMIT 1"
   command "*submit_job 1 *monitor_job"
  }

  popdown {
   position +28 =
   size 20 6
   text "EXECUTE 1"
   command "*execute_job 1 *monitor_job"
  }

  popdown {
   position -28 +6
   size 20 6
   text "SUBMIT 2"
   command "*submit_job 2 *monitor_job"
  }

  popdown {
   position +28 =
   size 20 6
   text "EXECUTE 2"
   command "*execute_job 2 *monitor_job"
  }

  popdown {
   position -28 +6
   size 20 6
   text "SUBMIT 3"
   command "*submit_job 3 *monitor_job"
  }

  popdown {
   position +28 =
   size 20 6
   text "EXECUTE 3"
   command "*execute_job 3 *monitor_job"
  }

  popdown {
   position 19 +8
   size 12 8
   text "OK"
  }
 }

 window {
  parent mentat
  origin 38 8
#ifdef DCOM
  size 50 100
#else
  size 50 94
#endif
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

popmenu damask {

#ifdef QT_MENTAT
 text "DAMASK.MPIE.DE"
#endif

 group {
#ifndef QT_MENTAT
  label {
   position 0 0
   size 50 4
   text "DAMASK.MPIE.DE"
  }
#endif

  label {
   position 1 6
   size 13   6
   text "Optimzation"
   border_width 1
   border_color black
  }

  label {
   position +13 =
   size 20 6
   text "write Input"
   border_width 1
   border_color black
  }

  label {
   position +18 =
   size 30 6
   text "do not write Inp."
   border_width 1
   border_color black
  }

  label {
   position -29 +6
   size 9 6
   text "O2"
   border_width 1
   border_color black
  }

  popdown {
   position +9 =
   size 20 6
   text "Submit"
   command "*submit_job 4 *monitor_job"
  }

  popdown {
   position +20 =
   size 20 6
   text "Execute"
   command "*execute_job 4 *monitor_job"
  }

  label {
   position -29 +6
   size 9 6
   text "O1"
   border_width 1
   border_color black
  }

  popdown {
   position +9 =
   size 20 6
   text "Submit"
   command "*submit_job 5 *monitor_job"
  }

  popdown {
   position +20 =
   size 20 6
   text "Execute"
   command "*execute_job 5 *monitor_job"
  }

  label {
   position -29 +6
   size 9 6
   text "O0"
   border_width 1
   border_color black
  }

  popdown {
   position +9 =
   size 20 6
   text "Submit"
   command "*submit_job 6 *monitor_job"
  }

  popdown {
   position +20 =
   size 20 6
   text "Execute"
   command "*execute_job 6 *monitor_job"
  }

  label {
   position -32 +6
   size 12 6
   text ""
   border_width 1
   border_color black
  }

  label {
   position = +6
   size 12 6
   text "O2 / OpenMP"
   border_width 1
   border_color black
  }

  popdown {
   position +12 =
   size 20 6
   text "Submit"
   command "*submit_job 7 *monitor_job"
  }

  popdown {
   position +20 =
   size 20 6
   text "Execute"
   command "*execute_job 7 *monitor_job"
  }

  label {
   position -32 +6
   size 12 6
   text "O1 / OpenMP"
   border_width 1
   border_color black
  }

  popdown {
   position +12 =
   size 20 6
   text "Submit"
   command "*submit_job 8 *monitor_job"
  }

  popdown {
   position +20 =
   size 20 6
   text "Execute"
   command "*execute_job 8 *monitor_job"
  }

  label {
   position -32 +6
   size 12 6
   text "O0 / OpenMP"
   border_width 1
   border_color black
  }

  popdown {
   position +12 =
   size 20 6
   text "Submit"
   command "*submit_job 9 *monitor_job"
  }

  popdown {
   position +20 =
   size 20 6
   text "Execute"
   command "*execute_job 9 *monitor_job"
  }

  label {
   position -32 +6
   size 12 6
   text ""
   border_width 1
   border_color black
  }

  popdown {
   position 19 +8
   size 12 8
   text "CANCEL"
  }
}


 window {
  parent mentat
  origin 38 8
#ifdef DCOM
  size 50 100
#else
  size 50 94
#endif
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

popmenu job_exit_msg_pm {

 group {

  label {
   position 0 0
   size 86 4
   text "EXIT MESSAGE"
  }

  display {
   position +1 +5
   size 84 74
   display "job_exit_msg"
  }

  popdown {
   position 37 +76
   size 12 8
   text "OK"
  }
 }

 window {
  parent mentat
  origin 38 8
  size 86 90
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

popmenu job_user_source_pm {
 window {
  parent job_run_wi
  origin 25 17
  size 24 17
  border_width 1
  border_color border
  background_color body
 }
 mode dialog

 group job_user_source_gr
}


group job_user_source_gr {
 popdown {
  position 0 0
  size 24 4
  text "COMPILE / NO SAVE"
  command "*job_option user_source:compile_nosave"
 }

 popdown {
  position = +4
  size 24 4
  text "COMPILE AND SAVE"
  command "*job_option user_source:compile_save"
 }

 popdown {
  position = +4
  size 24 4
  text "RUN SAVED EXECUTABLE"
  command "*job_option user_source:run_saved"
 }

 popdown {
  position = +5
  size 24 4
  text "CANCEL"
 }
}

popmenu job_run_parallelization_pm {

 group {

  label {
   position 0 0
   size 50 4
   text "PARALLELIZATION"
  }

  frame {
   position 1 +5
   size 48 28
   group job_ddm_gr
  }

  frame {
   position 1 +29
   size 48 26
   group job_parallel_matrix_solver_gr
  }

  frame {
   position 1 +27
   size 48 28
   group job_parallel_env_gr
   visible "or(*job_option parallel:on, \
               solver_multi_procs)"
  }

  popdown {
   position 19 +30
   size 12 8
   text "OK"
  }
 }

 window {
  parent mentat
  origin 38 8
  size 50 100
  background_color body
  border_width 1
  border_color border
  buffering single
 }
 mode dialog
}

group job_ddm_gr {
 label {
  position 0 0
  size 48 28
  text "DOMAIN DECOMPOSITION"
  border_width 1
  border_color black
 }

 toggle {
  position +1 +4
  size 30 4
  text "USE DDM"
  toggle "*job_option parallel:on"
  help job_run_ddm_use
  true_command "*job_option parallel:on"
  false_command "*job_option parallel:off"
 }

 frame {
  position +2 +5
  size 44 18
  group job_ddm_use_gr
  visible "*job_option parallel:on"
 }
}

group job_ddm_use_gr {
 roller {
  position 0 0
  size 30 4
  nvalues 2
  texts
   "DECOMPOSITION IN MARC"
   "DECOMPOSITION IN MENTAT"
  roller "job_option ddm_generator"
  popmenu job_ddm_generator_pm
  help job_run_ddm_generator
 }

 frame {
  position 0 +5
  size 44 8
  group job_ddm_fea_solver_gr
  visible "*job_option ddm_generator:fea_solver"
 }

 frame {
  position = =
  size 44 8
  group job_ddm_preprocessor_gr
  visible "*job_option ddm_generator:preprocessor"
 }

 text {
  position 0 +9
  size 22 4
  text "SINGLE INPUT FILE"
  readonly
  visible "*job_option ddm_generator:fea_solver"
 }

 roller {
  position = =
  size 22 4
  nvalues 2
  texts
   "MULTIPLE INPUT FILES"
   "SINGLE INPUT FILE"
  roller "job_option ddm_single_input"
  commands "*job_option ddm_single_input:off"
           "*job_option ddm_single_input:on"
  visible "*job_option ddm_generator:preprocessor"
  help job_run_ddm_single_input
 }

 roller {
  position +23 =
  size 21 4
  nvalues 2
  texts
   "MULTIPLE POST FILES"
   "SINGLE POST FILE"
  roller "job_option ddm_single_post"
  commands "*job_option ddm_single_post:off"
           "*job_option ddm_single_post:on"
  help job_run_ddm_single_post
 }
}

popmenu job_ddm_generator_pm {

 window {
  parent mentat
  origin 42 26
  size 30 13
  border_width 1
  border_color border
  background_color body
 }

 mode dialog

 group job_ddm_generator_gr
}

group job_ddm_generator_gr {
 popdown {
  position 0 0
  size 30 4
  text "DECOMPOSITION IN MARC"
  command "*job_option ddm_generator:fea_solver"
  help job_run_ddm_generator
 }

 popdown {
  position = +4
  size 30 4
  text "DECOMPOSITION IN MENTAT"
  command "*job_option ddm_generator:preprocessor"
  help job_run_ddm_generator
 }

 popdown {
  position = +5
  size 30 4
  text "CANCEL"
  help job_run_ddm_generator
 }
}

group job_ddm_fea_solver_gr {
 button {
  position 0 0
  size 10 4
  text "# DOMAINS"
  command "*job_param ndomains"
 }

 text {
  position +10 =
  size 4 4
  display "job_param_value_ndomains"
  command "*job_param ndomains"
 }

 label {
  position +4 =
  size 7 4
  text "METHOD"
  border_width 1
  border_color black
 }

 roller {
  position +7 =
  size 23 4
  nvalues 7
  texts
   "METIS BEST"
   "METIS ELEMENT"
   "METIS NODE"
   "REC. COORD. BISECTION"     
   "VECTOR"
   "RADIAL"
   "ANGULAR"
  help set_job_decomp_type
  rollers
   "*job_option ddm_method:metis_best"
   "*job_option ddm_method:metis_element_based"
   "*job_option ddm_method:metis_node_based"
   "*job_option ddm_method:recur_coord_bisect"     
   "*job_option ddm_method:vector"
   "*job_option ddm_method:radial"
   "*job_option ddm_method:angular"
  popmenu job_ddm_method_pm
 }

 button {
  position +4 +4
  size 19 4
  text "ADVANCED SETTINGS"
  screen job_ddm_fea_solver_sc
  command "@popdown($popname2) @popdown($popname1)" ($popname2, $popname1)
 }
}

popmenu job_ddm_method_pm {

 window {
  parent mentat
  origin 56 31
  size 30 33
  border_width 1
  border_color border
  background_color body
 }

 mode dialog

 group job_ddm_method_gr
}

group job_ddm_method_gr {

 popdown {
  position 0 0
  size 30 4
  text "METIS BEST"
  command "*job_option ddm_method:metis_best"
 }

 popdown {
  position = +4
  size 30 4
  text "METIS ELEMENT BASED"
  command "*job_option ddm_method:metis_element_based"
 }

 popdown {
  position = +4
  size 30 4
  text "METIS NODE BASED"
  command "*job_option ddm_method:metis_node_based"
 }

 popdown {
  position = +4
  size 30 4
  text "RECURSIVE COORDINATE BISECTION"
  command "*job_option ddm_method:recur_coord_bisect"
 }

 popdown {
  position = +4
  size 30 4
  text "VECTOR"
  command "*job_option ddm_method:vector"
 }

 popdown {
  position = +4
  size 30 4
  text "RADIAL"
  command "*job_option ddm_method:radial"
 }

 popdown {
  position = +4
  size 30 4
  text "ANGULAR"
  command "*job_option ddm_method:angular"
 }

 popdown {
  position = +5
  size 30 4
  text "CANCEL"
 }
}

screen job_ddm_fea_solver_sc {

 menu {

  group {
   units 32 120

   label {
    position 0 0
    size 32 4
    text "JOB PARALLELIZATION"
   }

   label {
    position 0 5
    size 6 4
    text "NAME"
   }

   display {
    position +6 =
    size 26 4
    display "job_name"
   }

   label {
    position 0 +5
    size 32 56
    text "ADVANCED DECOMPOSITION IN MARC"
    border_width 1
    border_color black
   }

   button {
    position 1 +4
    size 16 4
    text "# DOMAINS"
    command "*job_param ndomains"
   }

   text {
    position +16 =
    size 14 4
    display "job_param_value_ndomains"
    command "*job_param ndomains"
   }

   label {
    position 1 +4
    size 7 4
    text "METHOD"
    border_width 1
    border_color black
   }

   roller {
    position +7 =
    size 23 4
    nvalues 7
    texts
     "METIS BEST"
     "METIS ELEMENT"
     "METIS NODE"
     "REC. COORD. BISECTION"     
     "VECTOR"
     "RADIAL"
     "ANGULAR"
    help set_job_decomp_type
    rollers
     "*job_option ddm_method:metis_best"
     "*job_option ddm_method:metis_element_based"
     "*job_option ddm_method:metis_node_based"
     "*job_option ddm_method:recur_coord_bisect"     
     "*job_option ddm_method:vector"
     "*job_option ddm_method:radial"
     "*job_option ddm_method:angular"
    popmenu job_ddm_method2_pm
   }

   frame {
    position 1 +5
    size 30 20
    group job_ddm_direction_gr
    visible "or(*job_option ddm_method:vector, \
                *job_option ddm_method:radial, \
                *job_option ddm_method:angular)"
   }

   toggle {
    position 1 +21
    size 30 4
    text "DOMAIN ISLAND REMOVAL"
    toggle "*job_option ddm_island_removal:on"
    true_command "*job_option ddm_island_removal:on"
    false_command "*job_option ddm_island_removal:off"
   }

   toggle {
    position = +4
    size 30 4
    text "COARSE GRAPH"
    toggle "*job_option ddm_graph:coarse"
    true_command "*job_option ddm_graph:coarse"
    false_command "*job_option ddm_graph:fine"
    visible "or(*job_option ddm_method:metis_best, \
                *job_option ddm_method:metis_element_based, \
                *job_option ddm_method:metis_node_based)"
   }
   
   toggle {
    position = +4
    size 30 4
    text "DETECT CONTACT"
    toggle "*job_option ddm_detect_contact:on"
    true_command "*job_option ddm_detect_contact:on"
    false_command "*job_option ddm_detect_contact:off"
    visible "or(*job_option ddm_method:metis_best, \
                *job_option ddm_method:metis_element_based)"
   }

   button {
    position 1 +5
    size 30 4
    text "ELEMENT WEIGHT COEFFICIENT"
    command "*job_param ddm_elem_weight_coeff"
   }

   text {
    position = +4
    size 30 4
    display "job_param_value_ddm_elem_weight_coeff"
    command "*job_param ddm_elem_weight_coeff"
   }

   frame {
    position 0 +6
    size 32 8
    group ok2_gr
   }

   frame select_frame
   frame inactive_return
  }
  window main_window
  disable_right_mouse_button

 }

 menu transform

} # job_ddm_fea_solver_sc

popmenu job_ddm_method2_pm {

 window {
  parent mentat
  origin 1 22
  size 30 33
  border_width 1
  border_color border
  background_color body
 }

 mode dialog

 group job_ddm_method_gr
}

group job_ddm_direction_gr {
 button {
  position 0 0
  size 15 4
  text "DIRECTION"
  command "*job_vector ddm_sort_direction_x"
  visible "*job_option ddm_method:vector"
 }

 button {
  position = =
  size 15 4
  text "DIRECTION"
  command "*job_vector ddm_sort_direction_x"
  visible "not(*job_option ddm_method:vector)"
 }

 button {
  position +15 =
  size 15 4
  text "FROM / TO"
  command "*job_vector_from_to ddm_sort_direction_x"
 }

 text {
  position -15 +4
  size 10 4
  command "*job_param ddm_sort_direction_x"
  display "job_param_value_ddm_sort_direction_x"
  help ddm_job_decomp_user_direction_x
 }

 text {
  position +10 =
  size 10 4
  command "*job_param ddm_sort_direction_y"
  display "job_param_value_ddm_sort_direction_y"
 }

 text {
  position +10 =
  size 10 4
  command "*job_param ddm_sort_direction_z"
  display "job_param_value_ddm_sort_direction_z"
 }

 frame {
  position 0 +4
  size 30 8
  group job_ddm_sort_point_gr
  visible "not(*job_option ddm_method:vector)"
 }
}

group job_ddm_sort_point_gr {
 label {
  position 0 0
  size 14 4
  border_width 1
  border_color black
  text "POINT ON AXIS"
 }

 roller {
  position +14 =
  size 10 4
  nvalues 2
  texts
   "DEFAULT"
   "USER"
  roller "job_option ddm_sort_point"
  popmenu job_ddm_sort_point_meth_pm
 }

 button {
  position +10 =
  size 6 4
  text "SET"
  command "*job_position ddm_sort_point_x"
  visible "*job_option ddm_sort_point:user"
 }

 text {
  position 0 +4
  size 10 4
  command "*job_param ddm_sort_point_x"
  display "job_param_value_ddm_sort_point_x"
  visible "*job_option ddm_sort_point:user"
 }

 text {
  position +10 =
  size 10 4
  command "*job_param ddm_sort_point_y"
  display "job_param_value_ddm_sort_point_y"
  visible "*job_option ddm_sort_point:user"
 }

 text {
  position +10 =
  size 10 4
  command "*job_param ddm_sort_point_z"
  display "job_param_value_ddm_sort_point_z"
  visible "*job_option ddm_sort_point:user"
 }
}

popmenu job_ddm_sort_point_meth_pm {

 window {
  parent mentat
  origin 15 35
  size 14 13
  border_width 1
  border_color border
  background_color body
 }

 mode dialog

 group job_ddm_sort_point_meth_gr
}

group job_ddm_sort_point_meth_gr {
 popdown {
  position 0 0
  size 14 4
  text "DEFAULT"
  command "*job_option ddm_sort_point:default"
 }

 popdown {
  position = +4
  size 14 4
  text "USER DEFINED"
  command "*job_option ddm_sort_point:user"
 }

 popdown {
  position = +5
  size 14 4
  text "CANCEL"
 }
}

group job_ddm_preprocessor_gr {
 label {
  position 0 0
  size 10 4
  text "# DOMAINS"
  border_width 1
  border_color black
 }

 integer {
  position +10 =
  size 4 4
  display valid_domains
 }

 button {
  position +4 =
  size 30 4
  text "USER DOMAINS"
  screen job_ddm_user_domains_sc
  command "@popdown($popname2) @popdown($popname1)" ($popname2, $popname1)
  help job_run_ddm_user_domains
 }
}

screen job_ddm_user_domains_sc {

 menu {

  group {
   units 32 120

   label {
    position 0 0
    size 32 4
    text "USER DOMAINS"
   }

   frame {
    position 0 5
    size 32 28
    group user_domains_gr
   }

   frame {
    position 0 +29
    size 32 49
    group user_domains_generate_gr
   }

   frame {
    position 0 +50
    size 32 8
    group user_domains_tail_gr
   }

   frame {
    position 0 92
    size 32 8
    group ok2_gr
   }

   frame select_frame
   frame inactive_return
  }
  window main_window
  disable_right_mouse_button

 }

 menu transform

} # job_ddm_user_domains_sc

group job_parallel_matrix_solver_gr {
 label {
  position 0 0
  size 48 26
  text "MATRIX SOLVER"
  border_width 1
  border_color black
 }

 label {
  position +3 +4
  size 12 4
  text "SOLUTION"
  border_width 1
  border_color black
 }

 roller {
  position +12 =
  size 18 4
  nvalues 2
  texts
    "SYMMETRIC"
    "NONSYMMETRIC"
  rollers job_nonsym_off
          job_nonsym_on
  commands "*job_option solver_nonsym:off"
           "*job_option solver_nonsym:on"
  help job_param_solver_method
 }

 label {
  position -12 +4
  size 6 4
  text "TYPE"
  border_width 1
  border_color black
 }

 roller {
  position +6 =
  size 24 4
  nvalues 9
  help job_param_solver_method
  texts
    "MULTIFRONTAL SPARSE"
    "MIXED DIRECT-ITERATIVE"
    "CASI ITERATIVE"
    "PARDISO DIRECT SPARSE"
    "MUMPS PARALLEL DIRECT"
    "HARDWARE SPARSE"
    "ITERATIVE SPARSE"
    "DIRECT PROFILE"
    "DIRECT SPARSE"
  rollers job_solver_mfront_sparse
          job_solver_mixed_direct_iterative
          job_solver_it_ext
          job_solver_pardiso
          job_solver_mumps
          job_solver_sparse
          job_solver_it_sparse
          job_solver_dir_profile
          job_solver_dir_sparse
  popmenu job_set_parallel_matrix_solver_pm
 }

 button {
  position +24 =
  size 14 4
  text "OPTIONS"
  popmenu job_parallel_matrix_solver_opt_pm
 }

 frame {
  position = +4
  size 14 4
  group job_run_solver_ddm_opts_gr
  visible "*job_option parallel:on"
 }

 frame {
  position 1 +5
  size 46 8
  group job_solver_multi_procs_gr
  visible solver_allows_multi_procs
 }

 frame {
  position = =
  size 46 8
  group job_solver_multi_threads_gr
  visible solver_allows_multi_threads
 }
}

popmenu job_set_parallel_matrix_solver_pm {
 window {
  parent mentat
  origin 48 54
  size 24 42
  border_width 1
  border_color border
  background_color body
 }
 mode dialog

 group job_solver_gr
}

popmenu job_parallel_matrix_solver_opt_pm {
 group {
  label {
   position 0 0
   size 38 4
   text "MATRIX SOLVER OPTIONS"
  }

  frame {
   position 1 +5
   size 36 23
   group job_solver_opts_gr
  }

  popdown {
   position 13 +25
   size 12 8
   text "OK"
  }
 }

 window {
  parent mentat
  origin 41 12
  size 38 39
  border_width 1
  border_color border
  background_color body
 }
 mode dialog
}

group job_run_solver_ddm_opts_gr {
 button {
  position 0 0
  size 14 4
  text "DDM OPTIONS"
  popmenu ddm_options
# see also job_common.ms!
  visible "not(or(job_solver_it_sparse, \
                  job_solver_it_ext, \
                  job_solver_mixed_direct_iterative, \
                  job_solver_pardiso, \
                  job_solver_mumps))"
 }
}

group job_solver_multi_procs_gr {
 frame {
  position 0 0
  size 46 8
  group job_solver_multi_procs_parallel_off_gr
  visible "*job_option parallel:off"
 }

 frame {
  position = =
  size 46 8
  group job_solver_multi_procs_parallel_on_gr
  visible "*job_option parallel:on"
 }
}

group job_solver_multi_procs_parallel_off_gr {
 toggle {
  position 0 0
  size 30 4
  text "MULTIPLE SOLVER PROCESSES"
  true_command "*job_option nsolver_procs_serial:on"
  false_command "*job_option nsolver_procs_serial:off"
  toggle "*job_option nsolver_procs_serial:on"
  help job_run_multithreading
 }

 button {
  position +2 +4
  size 14 4
  text "# PROCESSES"
  command "*job_param nsolver_procs"
  visible "*job_option nsolver_procs_serial:on"
 }

 text {
  position +14 =
  size 14 4
  display "job_param_value_nsolver_procs"
  command "*job_param nsolver_procs"
  visible "*job_option nsolver_procs_serial:on"
 }
}

group job_solver_multi_procs_parallel_on_gr {
 toggle {
  position 0 0
  size 30 4
  text "MULTIPLE SOLVER PROCESSES"
  help job_run_multithreading
  toggle true
 }

 label {
  position +2 +4
  size 14 4
  text "# PROCESSES"
  border_width 1
  border_color black
 }

 roller {
  position +14 =
  size 14 4
  nvalues 2
  texts "AUTOMATIC"
        "USER"
  commands "*job_option nsolver_procs_ddm:automatic"
           "*job_option nsolver_procs_ddm:user"
  help job_run_multithreading
  rollers "*job_option nsolver_procs_ddm:automatic"
          "*job_option nsolver_procs_ddm:user"
 }

 frame {
  position +14 =
  size 16 4
  group job_nsolver_procs_ddm_automatic_gr
  visible "*job_option nsolver_procs_ddm:automatic"
 }

 frame {
  position = =
  size 16 4
  group job_nsolver_procs_ddm_user_gr
  visible "*job_option nsolver_procs_ddm:user"
 }
}

group job_nsolver_procs_ddm_automatic_gr {
 label {
  position 0 0
  size 8 4
  text "VALUE"
  border_width 1
  border_color black
 }

 integer {
  position +8 =
  size 8 4
  display valid_domains
  visible "*job_option ddm_generator:preprocessor"
 }

 integer {
  position = =
  size 8 4
  display "job_param_ndomains"
  visible "*job_option ddm_generator:fea_solver"
 }
}

group job_nsolver_procs_ddm_user_gr {
 button {
  position 0 0
  size 8 4
  text "VALUE"
  command "*job_param nsolver_procs"
 }

 text {
  position +8 =
  size 8 4
  display "job_param_value_nsolver_procs"
  command "*job_param nsolver_procs"
 }
}

group job_solver_multi_threads_gr {
 frame {
  position 0 0
  size 30 8
  group job_solver_multi_threads_mfront_sparse_gr
  visible "and(job_solver_mfront_sparse,*job_option parallel:off)"
 }

 frame {
  position = =
  size 46 8
  group job_solver_multi_threads_pardiso_parallel_off_gr
  visible "and(job_solver_pardiso,*job_option parallel:off)"
 }

 frame {
  position = =
  size 46 8
  group job_solver_multi_threads_pardiso_parallel_on_gr
  visible "and(job_solver_pardiso,*job_option parallel:on)"
 }
}

group job_solver_multi_threads_mfront_sparse_gr {
 toggle {
  position 0 0
  size 30 4
  text "MULTIPLE SOLVER PROCESSES"
  toggle "*job_option mfront_sparse_multi_threading:on"
  true_command "*job_option mfront_sparse_multi_threading:on"
  false_command "*job_option mfront_sparse_multi_threading:off"
  help job_run_multithreading
 }

 button {
  position +2 +4
  size 14 4
  text "# PROCESSES"
  command "*job_param nthreads"
  visible "*job_option mfront_sparse_multi_threading:on"
 }

 text {
  position +14 =
  size 14 4
  display "job_param_value_nthreads"
  command "*job_param nthreads"
  visible "*job_option mfront_sparse_multi_threading:on"
 }
}

group job_solver_multi_threads_pardiso_parallel_off_gr {
 toggle {
  position 0 0
  size 30 4
  text "MULTIPLE SOLVER PROCESSES"
  toggle "*job_option pardiso_multi_threading:on"
  true_command "*job_option pardiso_multi_threading:on"
  false_command "*job_option pardiso_multi_threading:off"
  help job_run_multithreading
 }

 button {
  position +2 +4
  size 14 4
  text "# PROCESSES"
  command "*job_param nthreads"
  visible "*job_option pardiso_multi_threading:on"
 }

 text {
  position +14 =
  size 14 4
  display "job_param_value_nthreads"
  command "*job_param nthreads"
  visible "*job_option pardiso_multi_threading:on"
 }
}

group job_solver_multi_threads_pardiso_parallel_on_gr {
 toggle {
  position 0 0
  size 30 4
  text "MULTIPLE SOLVER PROCESSES"
  help job_run_multithreading
  toggle true
 }

 label {
  position +2 +4
  size 14 4
  text "# PROCESSES"
  border_color black
  border_width 1
 }

 roller {
  position +14 =
  size 14 4
  nvalues 2
  texts "AUTOMATIC"
        "USER"
  commands "*job_option pardiso_multi_threading_ddm:automatic"
           "*job_option pardiso_multi_threading_ddm:user"
  help job_run_multithreading
  rollers "*job_option pardiso_multi_threading_ddm:automatic"
          "*job_option pardiso_multi_threading_ddm:user"
 }

 frame {
  position +14 =
  size 16 4
  group job_pardiso_multi_threads_ddm_automatic_gr
  visible "*job_option pardiso_multi_threading_ddm:automatic"
 }

 frame {
  position = =
  size 16 4
  group job_pardiso_multi_threads_ddm_user_gr
  visible "*job_option pardiso_multi_threading_ddm:user"
 }
}

group job_pardiso_multi_threads_ddm_automatic_gr {
 label {
  position 0 0
  size 8 4
  text "VALUE"
  border_width 1
  border_color black
 }

 integer {
  position +8 =
  size 8 4
  display valid_domains
  visible "*job_option ddm_generator:preprocessor"
 }

 integer {
  position = =
  size 8 4
  display "job_param_ndomains"
  visible "*job_option ddm_generator:fea_solver"
 }
}

group job_pardiso_multi_threads_ddm_user_gr {
 button {
  position 0 0
  size 8 4
  text "VALUE"
  command "*job_param nthreads"
 }

 text {
  position +8 =
  size 8 4
  display "job_param_value_nthreads"
  command "*job_param nthreads"
 }
}

group job_parallel_env_gr {
 label {
  position 0 0
  size 48 28
  text "PARALLELIZATION ENVIRONMENT"
  border_width 1
  border_color black
 }

 roller {
  position 1 +4
  size 30 4
  nvalues 2
  texts "SINGLE MACHINE"
        "NETWORK"
  help job_run_ddm_setup
  roller "job_option parallel_setup"
  commands "*job_option parallel_setup:single"
           "*job_option parallel_setup:network"
 }

 frame {
  position +2 +5
  size 44 18
  group job_parallel_env_network_gr
  visible "*job_option parallel_setup:network"
 }
}

group job_parallel_env_network_gr {
 button {
  position 0 0
  size 28 4
  text "HOST FILE"
  popmenu host_file_browser_pm
  set $host_file_browser_label "SELECT HOST FILE"
  set $host_file_browser_command "*job_host_file"
  help job_host_file
 }

 button {
  position +28 =
  size 8 4
  text "EDIT"
  command "*job_edit_host_file"
  help job_edit_host_file
  visible job_host_file
 }

 button {
  position +8 =
  size 8 4
  text "CLEAR"
  command "*job_clear_host_file"
  help job_clear_host_file
  visible job_host_file
 }

 display {
  position 0 +4
  size 44 4
  display job_host_file
 }

 frame {
  position 0 +5
  size 44 9
  group job_parallel_env_network_ddm_gr
  visible "*job_option parallel:on"
 }
}

group job_parallel_env_network_ddm_gr {
 toggle {
  position 0 0
  size 22 4
  text "COPY INPUT FILE"
  toggle "*job_option copy_input_file:on"
  true_command "*job_option copy_input_file:on"
  false_command "*job_option copy_input_file:off"
  help job_host_copy_inputfile
 }

 toggle {
  position +23 =
  size 21 4
  text "COPY POST FILE"
  toggle "*job_option copy_post_file:on"
  true_command "*job_option copy_post_file:on"
  false_command "*job_option copy_post_file:off"
  help job_host_copy_inputfile
 }

 label {
  position 0 +5
  size 10 4
  text "HOSTS"
  border_width 1
  border_color black
  visible job_usersub_file
 }

 roller {
  position +10 =
  size 18 4
  nvalues 2
  texts "COMPATIBLE"
        "INCOMPATIBLE"
  roller "job_option network_hosts"
  commands "*job_option network_hosts:compatible"
           "*job_option network_hosts:incompatible"
  help job_host_comp
  visible job_usersub_file
 }
}

popmenu marc_integer_size_pm {
 window {
  parent mentat
  origin 55 32
  size 32 17
  border_width 1
  border_color border
  background_color body
 }
 mode dialog

 group marc_integer_size_gr
}

group marc_integer_size_gr {
 popdown {
  position 0 0
  size 32 4
  text "DEFAULT"
  command "*job_option marc_integer_size:default"
  help job_run_intsize
 }

 popdown {
  position = +4
  size 32 4
  text "4-BYTE INTEGERS"
  command "*job_option marc_integer_size:i4"
  help job_run_intsize
 }

 popdown {
  position = +4
  size 32 4
  text "8-BYTE INTEGERS"
  command "*job_option marc_integer_size:i8"
  help job_run_intsize
 }

 popdown {
  position = +5
  size 32 4
  text "CANCEL"
  help job_run_intsize
 }
}
