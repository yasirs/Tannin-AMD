tryCatch({
    library(msigdbr)
    m_df <- msigdbr(species = "Homo sapiens", category = "C2")
    if ("SAUL_SEN_MAYO" %in% m_df$gs_name) {
        print("FOUND: SAUL_SEN_MAYO")
        # Save it to a file directly if found
        senmayo <- m_df[m_df$gs_name == "SAUL_SEN_MAYO", ]
        write.csv(senmayo, "results/three-axes/gene-sets/senmayo_temp.csv", row.names=FALSE)
    } else {
        print("NOT FOUND in C2")
    }
}, error = function(e) {
    print(paste("Error:", e$message))
})
