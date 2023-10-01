dimred.auto <- function (object, lognorm = T, scale.factor = 10000, plot.vars = NULL, 
    base.size = 5, umap = T, tsne = T, pca = T) 
{
    require(Seurat)
    suppressMessages(suppressWarnings(expr = {
        if (lognorm) {
            object <- NormalizeData(object, scale.factor = scale.factor)
        }
        object <- ScaleData(object)
        object <- FindVariableFeatures(object)
        npcs <- min(ncol(object) - 1, 50)
        object <- RunPCA(object, npcs = npcs)
        options(repr.plot.height = 8, repr.plot.width = 16)
    }))
    if (!is.null(plot.vars) & pca) {
        pca.plots <- NULL
        for (pvar in plot.vars) {
            p <- LabelClusters(DimPlot(object, pt.size = 4, reduction = "pca", 
                group.by = pvar, ), id = pvar)
            pca.plots <- if (is.null(pca.plots)) {
                p
            }
            else {
                pca.plots + p
            }
        }
        plot(pca.plots)
    }
    if (umap) {
        object <- RunUMAP(object, dims = 1:npcs, verbose = F, 
            spread = 4, min.dist = 0)
        if (!is.null(plot.vars)) {
            umap.plots <- NULL
            for (pvar in plot.vars) {
                p <- LabelClusters(DimPlot(object, pt.size = 4, 
                  reduction = "umap", group.by = pvar, ), id = pvar)
                umap.plots <- if (is.null(umap.plots)) {
                  p
                }
                else {
                  umap.plots + p
                }
            }
            plot(umap.plots)
        }
    }
    if (tsne) {
        object <- RunTSNE(object, dims = 1:npcs, perplexity = 5)
        if (!is.null(plot.vars)) {
            tsne.plots <- NULL
            for (pvar in plot.vars) {
                p <- LabelClusters(DimPlot(object, pt.size = 4, 
                  reduction = "tsne", group.by = pvar, ), id = pvar)
                tsne.plots <- if (is.null(tsne.plots)) {
                  p
                }
                else {
                  tsne.plots + p
                }
            }
            plot(tsne.plots)
        }
    }
    return(object)
}

example.pseudobulk_generate <- function (seurat.dataname = "ifnb", group.vars = "stim", cluster.var = "seurat_annotations", 
    sel.method = "simspec", agg.method = "aggregate", simspec.ratio = 1/300, 
    samples.per.group = 3, cells.per.sample = NULL, lognorm = T, 
    scale.factor = 10000, verbose = T, return. = F, ...) 
{
    require(Seurat)
    require(SeuratData)
    require(glue)
    print(glue("First we load the object {seurat.dataname}"))
    object <- LoadData(seurat.dataname)
    print(head(object@meta.data))
    psbk <- generate.pseudobulk(object, group.vars = group.vars, 
        cluster.var = cluster.var, simspec.ratio = simspec.ratio, 
        sel.method = sel.method, agg.method = agg.method, samples.per.group = samples.per.group, 
        cells.per.sample = cells.per.sample, verbose = verbose, 
        ...)
    psbk <- dimred.auto(psbk, scale.factor = scale.factor, lognorm = lognorm, 
        plot.vars = c(cluster.var, group.vars))
    if (return.) {
        return(list(data = object, pseudobulk = psbk))
    }
    return()
}

generate.pseudobulk <- function (object, group.vars = NULL, cluster.var = NULL, annotation.vars = "all", 
    slot = "counts", features = NULL, samples.per.group = NULL, 
    cells.per.sample = NULL, simspec.ratio = 1/100, sel.method = "recursive", 
    agg.method = "aggregate", verbose = T, ...) 
{
    clgroup.vars <- if (is.null(cluster.var)) {
        group.vars
    }
    else {
        c(cluster.var, group.vars)
    }
    if (verbose) {
        message("Grouping cells")
    }
    grouped.cells <- get.grouped.cells(object, clgroup.vars)
    if (sel.method == "simspec") {
        if (verbose) {
            message("Sampling cells in neighbors space (simspec)")
        }
        sampled.cells <- lapply(grouped.cells, function(x) {
            suppressMessages(suppressWarnings(sample.simspec(object[, 
                x], ratio = simspec.ratio, ...)))
        })
    }
    else {
        if (verbose) {
            message("Sampling cells recursively")
        }
        sampled.cells <- lapply(grouped.cells, sample.recursively, 
            n_samples = samples.per.group, sample_size = cells.per.sample, 
            value = T, verbose = F)
    }
    if (verbose) {
        message("Naming samples")
    }
    sampled.cells <- lapply(names(sampled.cells), function(x) {
        y <- sampled.cells[[x]]
        names(y) <- paste(x, as.character(1:length(y)), sep = "_rep")
        y
    })
    sampled.cells <- do.call(sampled.cells, what = "c")
    if (verbose) {
        message("Aggregating expression")
    }
    psbk.objects <- lapply(names(sampled.cells), function(psbk.lab) {
        if (verbose) {
            i <- which(names(sampled.cells) == psbk.lab)
            p.i <- round(length(sampled.cells)/10 * (1:10), 0)
            if (i %in% p.i) {
                message(glue::glue("+ {i}/{length(sampled.cells)} ({round(i/length(sampled.cells)*100,2)}%)"))
            }
        }
        suppressMessages(suppressWarnings(expr = {
            psbk <- get.average.expression(get.seurat.from.cells(object, 
                cells = sampled.cells[[psbk.lab]]), group.vars = group.vars, 
                cluster.var = cluster.var, slot = slot, features = features, 
                verbose = verbose, method = agg.method, annotation.vars = annotation.vars, 
                ...)
            psbk <- rename.seurat.cells(psbk, psbk.lab)
            lab.split <- strsplit(psbk.lab, split = "_rep", fixed = T)[[1]]
            psbk$pseudobulk.id <- psbk.lab
            psbk$pseudobulk.group <- lab.split[1]
            psbk$pseudobulk.rep <- lab.split[2]
            psbk
        }))
        return(psbk)
    })
    if (verbose) {
        message(glue::glue("Merging {length(psbk.objects)} objects"))
    }
    psbk <- merge.seurat(psbk.objects)
    if (verbose) {
        message(glue::glue("Restoring factors"))
    }
    psbk@meta.data <- restore.factors(x = psbk@meta.data, reference = object@meta.data)
    return(psbk)
}

get.average.expression <- function (object, group.vars = NULL, cluster.var = NULL, annotation.vars = "all", 
    slot = "counts", features = NULL, method = "aggregate", verbose = T, 
    ...) 
{
    require(Seurat)
    warn <- if (verbose) {
        function(x) {
            x
        }
    }
    else {
        suppressMessages
    }
    group.vars <- c(cluster.var, group.vars)
    groups.df <- unique(FetchData(object = object, vars = group.vars))
    rownames(groups.df) <- apply(groups.df, 1, paste, collapse = "_")
    if (!is.null(annotation.vars)) {
        if (annotation.vars[1] == "all") {
            annotation.vars <- colnames(object@meta.data)[!colnames(object@meta.data) %in% 
                group.vars]
        }
        groups.df <- cbind(groups.df[, -(1:length(group.vars))], 
            summarise.metadata(object@meta.data, group.vars = group.vars, 
                annotation.vars = annotation.vars))
    }
    groups.df$pseudobulk.id <- rownames(groups.df)
    avg.vars <- group.vars[sapply(apply(groups.df[, group.vars, 
        drop = F], 2, unique), length) != 1]
    if (length(avg.vars) == 0) {
        object$fake_group <- "t"
        avg.vars <- c("fake_group")
    }
    if (method == "average") {
        agg.exp <- warn(Seurat::AverageExpression(object, group.by = avg.vars, 
            slot = slot, features = features, return.seurat = T, 
            verbose = verbose, ...))
    }
    else if (method == "aggregate") {
        agg.exp <- warn(Seurat::AggregateExpression(object, group.by = avg.vars, 
            slot = slot, features = features, return.seurat = T, 
            verbose = verbose, ...))
    }
    else {
        stop("`method` must be \"aggregate\" or \"average\"")
    }
    if (ncol(agg.exp) == 1) {
        agg.exp@meta.data <- groups.df
    }
    else {
        agg.exp@meta.data[, colnames(agg.exp@meta.data) %in% 
            colnames(groups.df)] <- NULL
        agg.exp@meta.data <- cbind(agg.exp@meta.data, groups.df[rownames(agg.exp@meta.data), 
            ])
    }
    return(agg.exp)
}

get.grouped.cells <- function (object, group.vars) 
{
    grouped.cells <- split(colnames(object), apply(object@meta.data[, 
        group.vars], 1, paste, collapse = "_"))
    return(grouped.cells)
}

get.seurat.from.cells <- function (object, cells) 
{
    require(dplyr)
    if (!any(duplicated(cells))) {
        return(object[, cells])
    }
    else {
        cells <- data.frame(cells) %>% group_by(cells) %>% mutate(rep = 1:length(cells), 
            new.cell = paste(cells, rep, sep = "_"))
        new.obj <- lapply(split(cells, cells$rep), function(celldata) {
            new.obj <- object[, celldata$cells]
            new.obj <- rename.seurat.cells(object = new.obj, 
                new.names = celldata$new.cell)
            new.obj
        }) %>% merge.seurat(merge.data = T)
        return(new.obj)
    }
}

merge.seurat <- function (seurat_list, merge.data = T) 
{
    return(purrr::reduce(seurat_list, function(x, y) {
        merge(x = x, y = y, merge.data = merge.data)
    }))
}

rename.seurat.cells <- function (object, new.names, assays = c("RNA"), slots = c("counts", 
    "data", "scale.data")) 
{
    require(Seurat)
    meta.data <- object@meta.data
    rownames(meta.data) <- new.names
    counts <- object[["RNA"]]@counts
    colnames(counts) <- new.names
    newobj <- CreateSeuratObject(counts = counts, row.names = rownames(object), 
        min.cells = 0, min.features = 0, assay = "RNA", meta.data = meta.data)
    for (assay in names(object)) {
        if (!assay %in% assays) {
            (next)()
        }
        for (sl in slotNames(object[[assay]])) {
            if (!sl %in% slots) {
                (next)()
            }
            slotdata <- slot(object[[assay]], sl)
            if (0 %in% dim(slotdata)) {
                (next)()
            }
            colnames(slotdata) <- new.names
            slot(newobj[[assay]], sl) <- slotdata
        }
    }
    return(newobj)
}

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

sample.recursively <- function (choices, n_samples = NULL, sample_size = NULL, value = T, 
    seed = NULL, verbose = F, duplicates = T, ...) 
{
    samples <- list(...)$samples
    if (all(is.null(n_samples), is.null(sample_size))) {
        stop("Must provide one of `n_samples` or `sample_size`")
    }
    if (is.null(samples)) {
        if (verbose) {
            message(glue::glue("Initializing for n_samples:{n_samples}, sample_size:{sample_size}"))
        }
        samples <- list(remaining = choices)
    }
    if (is.null(seed)) {
        seed <- sample(1:99999, 1)
    }
    if (is.null(sample_size)) {
        sample_size <- round(length(choices)/n_samples, 0)
        if (verbose) {
            message(glue::glue("`sample_size` set to {sample_size}"))
        }
    }
    else if (sample_size < 1) {
        sample_size <- round(length(choices) * sample_size, 0)
        if (verbose) {
            message(glue::glue("`sample_size` set to {sample_size}"))
        }
    }
    if (is.null(n_samples)) {
        n_samples <- round(length(choices)/sample_size, 0)
        if (verbose) {
            message(glue::glue("`n_samples` set to {n_samples}"))
        }
    }
    set.seed(seed)
    if (length(samples) - 1 == n_samples) {
        if (verbose) {
            message("Finished")
        }
        if (value) {
            return(samples[-1])
        }
        else {
            masks <- lapply(samples[-1], function(elements) {
                choices %in% elements
            })
            return(masks)
        }
    }
    else {
        n.rem <- length(samples$remaining)
        if (verbose) {
            message(glue::glue("Available for sample: {n.rem}"))
        }
        sid <- length(samples) + 1
        if (sample_size > n.rem) {
            if (verbose) {
                message(glue::glue("It is lower than {sample_size}. Appending remaining and sampling {sample_size-n.rem}"))
            }
            extra.s <- samples$remaining
            new.choices <- choices[!choices %in% extra.s]
            while (length(new.choices) < sample_size - n.rem) {
                new.choices <- c(new.choices, choices)
            }
            s <- sample(x = 1:length(new.choices), size = sample_size - 
                n.rem, replace = F)
            samples$remaining <- c(new.choices[-s], extra.s)
            s <- c(extra.s, new.choices[s])
            samples[[sid]] <- sort(s)
        }
        else {
            new.choices <- samples$remaining
            s <- sample(x = 1:length(new.choices), size = sample_size, 
                replace = F)
            samples[[sid]] <- sort(new.choices[s])
            samples$remaining <- new.choices[-s]
        }
        return(sample.recursively(choices = choices, n_samples = n_samples, 
            sample_size = sample_size, seed = seed + 1, value = value, 
            verbose = verbose, duplicates = duplicates, samples = samples))
    }
}

sample.simspec <- function (object, ratio = 1/10, ...) 
{
    require(dplyr)
    require(Seurat)
    npcs <- min(50, ncol(object) - 1)
    snn <- (ScaleData(object, verbose = F) %>% FindVariableFeatures(verbose = F) %>% 
        RunPCA(verbose = F, npcs = npcs) %>% FindNeighbors(verbose = F, 
        dims = 1:npcs))@graphs$RNA_snn
    sel_df <- simspec::construct_pseudocells(snn, ratio = ratio, 
        ...)
    sampled.cells <- lapply(sort(unique(sel_df[sel_df[, 2] != 
        -1, 2])), function(ix) {
        colnames(object)[sel_df[, 2] == ix]
    })
    return(sampled.cells)
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