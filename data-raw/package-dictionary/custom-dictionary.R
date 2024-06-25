custom_user_spelling_dictionary = function(custom_words, pkg) {

  # Adapted from tmsalab / edmdata
  # Create an R-level dictionary using an RDS file.
  saveRDS(custom_words, file = file.path(".aspell", paste0(pkg, ".rds")), version = 2)

  # For the {spelling} package, write out a simple WORDLIST.
  writeLines(custom_words, file.path("inst", "WORDLIST"))

  # Provide backward compatibility support for base R
  # Note, must be run from package root.
  # For more details, see: ?"aspell-utils"
  utils::aspell_write_personal_dictionary_file(
    custom_words,
    out = file.path(".aspell", "words.pws")
  )

  invisible(custom_words)
}

custom_words = c("recalibration", "Recalibration", "MLE", "LLO", "Frequentist")

custom_user_spelling_dictionary(custom_words, "BRcal")
