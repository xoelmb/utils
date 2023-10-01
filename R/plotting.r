fig <- function (w = NA, h = NA) 
{
    options(repr.plot.width = w, repr.plot.height = h)
}

save_functions <- function (basename = "functions", functions = NULL, envir.get = NULL, 
    R_script = T, Rda = T) 
{
    require(glue)
    if (is.null(envir.get)) {
        envir.get <- globalenv()
    }
    if (is.null(functions)) {
        env.vars <- sort(ls(envir = envir.get))
        functions <- env.vars[sapply(sapply(env.vars, get, envir = envir.get), 
            class) == "function"]
    }
    if (Rda) {
        save(list = functions, envir = envir.get, file = glue("{basename}.rda"))
    }
    if (R_script) {
        script <- paste(sapply(1:length(functions), function(i, 
            e = envir.get) {
            paste(functions[i], paste(deparse(get(x = functions[i], 
                envir = e)), collapse = "\n"), sep = " <- ")
        }), collapse = "\n\n")
        cat(script, file = glue("{basename}.r"))
        return(script)
    }
}

set_parallel <- function (slurm_check = T, max_cores = Inf, future = T, future.strategy = "multicore", 
    future.mem = Inf) 
{
    n.cores <- NA
    if (slurm_check) {
        n.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    }
    n.cores <- if (!is.na(n.cores) & n.cores > 1) 
        n.cores
    else parallel::detectCores()
    n.cores <- min(max_cores, n.cores)
    if (future) {
        if (!is.null(future.mem)) {
            options(future.globals.maxSize = future.mem)
        }
        future::plan(strategy = future.strategy, workers = n.cores)
    }
    return(n.cores)
}