load(file = here::here("./data/c_set_w_act_pentamer.RData"))

ebaa_extra <- ebaa_extra[ebaa_extra$antigen %in% c("DT","FHA","FIM","IPV1","IPV2","IPV3","PRN","PT","TT","WHOLE"), ]

save(ebaa_extra, file = here::here("./data/c_set.RData"))
