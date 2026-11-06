# Saving and loading pipeline objects.
# Objects are saved as .rds (one object, read back with whatever name you like). Files written by the old
# scripts with save() can still be read, and the old format can also be written for downstream scripts.

#Save an object as .rds. If legacy_name is given, also write the old save() format (file name without
#".rds") holding the list's elements as separate objects, e.g. COMB_mats and filter_params.
save_object = function(object, path, legacy_name = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  if (!is.null(legacy_name)) {
    legacy_path = sub("\\.rds$", "", path)
    env = new.env()
    if (identical(legacy_name, TRUE)) list2env(object, env) else assign(legacy_name, object, envir = env)
    save(list = ls(env), envir = env, file = legacy_path)
  }
  invisible(path)
}

#Read a mats_and_params object: either an .rds list(COMB_mats, filter_params) or an old save() file
read_mats_and_params = function(path) {
  if (grepl("\\.rds$", path)) {
    x = readRDS(path)
  } else {
    env = new.env()
    load(path, envir = env)
    x = list(COMB_mats = env$COMB_mats, filter_params = env$filter_params)
  }
  if (is.null(x$COMB_mats) || is.null(x$filter_params)) stop(path, " does not contain COMB_mats and filter_params")
  x
}
