restore.factors <- function (x, reference, factor.vars = NULL) 
{
    if (is.null(factor.vars)) {
        factor.vars <- colnames(reference)[unlist(lapply(reference, 
            is.factor))]
    }
    for (fct.var in factor.vars) {
        if (!fct.var %in% colnames(x)) {
            (next)()
        }
        if (!all(unique(x[, fct.var]) %in% levels(reference[, 
            fct.var]))) {
            (next)()
        }
        x[, fct.var] <- factor(as.character(x[, fct.var]), levels = levels(reference[, 
            fct.var]))
    }
    return(x)
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

summarise.metadata <- function (meta.data, group.vars, annotation.vars, keep_factors = T, 
    num.fun = mean, char.fun = function(x) {
        paste(unique(as.character(x)), collapse = ",")
    }) 
{
    d <- meta.data
    require(dplyr)
    sum.d <- d %>% group_by_at(group.vars) %>% select_at(annotation.vars) %>% 
        summarise_all(function(x) {
            if (is.numeric(x)) {
                num.fun(x)
            }
            else {
                char.fun(x)
            }
        }) %>% as.data.frame()
    if (keep_factors) {
        sum.d <- restore.factors(x = sum.d, reference = meta.data)
    }
    return(sum.d)
}