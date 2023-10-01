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