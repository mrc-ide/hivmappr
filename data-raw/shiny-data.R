devtools::load_all("~/Documents/Research/spatial-estimates/ART-coverage/manuscript/prevartcov/")

library(rgdal)

#' Load Malawi data
data(mwsh)

distord <- order(factor(mwsh$region, c("Northern", "Central", "Southern")), -coordinates(mwsh)[,2])

sh <- mwsh[distord, c("region", "zone", "district")]
dir.create("../inst/extdata/mwsh/")
writeOGR(sh, "../inst/extdata/mwsh/", "districts", driver="ESRI Shapefile")


rita <- structure(list(district = c("Chitipa", "Karonga", "Rumphi", "Mzimba", 
"Kasungu", "Nkhotakota", "Ntchisi", "Dowa", "Salima", "Lilongwe", 
"Mchinji", "Dedza", "Ntcheu", "Mangochi", "Machinga", "Zomba", 
"Chiradzulu", "Blantyre", "Mwanza", "Thyolo", "Mulanje", "Phalombe", 
"Chikwawa", "Nsanje", "Balaka", "Neno", "Nkhata Bay", "Likoma"
), nsamp = c(218L, 322L, 209L, 1045L, 764L, 332L, 229L, 674L, 
351L, 2462L, 596L, 712L, 531L, 811L, 529L, 769L, 317L, 1271L, 
114L, 601L, 541L, 376L, 478L, 276L, 335L, 131L, 247L, 14L), npos = c(18L, 
24L, 5L, 80L, 34L, 25L, 10L, 28L, 30L, 208L, 28L, 30L, 50L, 98L, 
57L, 86L, 68L, 220L, 16L, 106L, 85L, 48L, 50L, 35L, 48L, 18L, 
21L, 1L), nrecent = c(1L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 3L, 
0L, 0L, 0L, 2L, 0L, 1L, 1L, 4L, 1L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L)), .Names = c("district", "nsamp", "npos", "nrecent"), row.names = c(NA, 
-28L), class = "data.frame")


df <- dplyr::inner_join(as.data.frame(mwsh[distord,]), rita)
df$npos <- with(df, round(nsamp * prev_survey))

df <- df[c("district", "pop15pl", "pop15to49", "prev_survey", "prev_survey_se",
           "nsamp", "npos", "nrecent", "adultart", "anc_clients", "ancrt_n",
           "ancrt_pos", "ancrt_art")]

write.csv(df, "../inst/extdata/mwdf.csv", row.names=FALSE)
