## Add continental population level of 1KGP participants

label1KGPpop <- function(id_arr){
  super_arr <- rep(NA, length(id_arr))
  
  # 1KGP subpopulations super population categories
  eas_pop <- c('CHB','JPT','CHS','CDX','KHV')
  eur_pop <- c('CEU','TSI','FIN','GBR','IBS')
  afr_pop <- c('YRI','LWK','GWD','MSL','ESN') # no admixed Africans
  amr_pop <- c('MXL','PUR','CLM','PEL')
  sas_pop <- c('GIH','PJL','BEB','STU','ITU')
  
  super_arr[which(id_arr %in% eas_pop)] <- 'EAS'
  super_arr[which(id_arr %in% eur_pop)] <- 'EUR'
  super_arr[which(id_arr %in% afr_pop)] <- 'AFR'
  super_arr[which(id_arr %in% amr_pop)] <- 'AMR'
  super_arr[which(id_arr %in% sas_pop)] <- 'SAS'
    
  return(super_arr)
}