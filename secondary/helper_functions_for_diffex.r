# Helper functions
library(edgeR)

# Mark cells that are positive ONLY for a given construct/vector/plasmid,
# and negative for a target.
mark_vector_pos_target_neg <- function(
    seurat_obj,
    perturbed_cells_by_guide,
    df_guide,
    guides_on_vector,
    guides4target,
    print_counts = T,
    pos_label = "vector_positive",
    neg_label = "target_negative") {
    all_cells <- Cells(seurat_obj)
    target_positive_cells <- c()
    dummy_perturbed <- c()

    Idents(seurat_obj) <- "vector_B"
    # Select all cells in which the target of the guide has been perturbed
    # by any guide in guides4target list.
    for (guide in guides4target) {
        dummy_perturbed <- unlist(perturbed_cells_by_guide[[guide]])
        target_positive_cells <- union(target_positive_cells, dummy_perturbed)
    }

    target_negative_cells <- unlist(setdiff(all_cells, target_positive_cells))
    if (length(target_negative_cells) > 0) {
        seurat_obj <-
            SetIdent(
                seurat_obj,
                cells = target_negative_cells, value = neg_label
            )
    }

    # Now find the cells with the guides on the vector, then mark them
    #
    vector_positive_cells <- c()
    dummy_perturbed <- c()

    for (guide in guides_on_vector) {
        dummy_perturbed <- unlist(perturbed_cells_by_guide[[guide]])
        vector_positive_cells <- union(vector_positive_cells, dummy_perturbed)
        if (print_counts) {
            cat(guide, " ", length(dummy_perturbed), "\n")
        }
    }

    if (length(vector_positive_cells) > 0) {
        seurat_obj <-
            SetIdent(
                seurat_obj,
                cells = vector_positive_cells, value = pos_label
            )
    }



    print(guides_on_vector)
    print(guides4target)
    df_tally <- as.data.frame(table(Idents(seurat_obj)))

    n_vec_positive <- length(vector_positive_cells)
    n_target_neg <- length(target_negative_cells)
    n_all <- length(all_cells)

    if (print_counts) {
        cat(
            blue(
                "Vector+ =",
                n_vec_positive,
                "; Target- =",
                n_target_neg,
                "; All =",
                n_all,
                "\n"
            )
        )
    }
    seurat_obj
}

# Mark cells that are positive ONLY for a given guide
# (not construct), and negative for a target
# The purpose is to run diffex analysis for each guide seperately.
mark_guide_pos_target_neg <- function(
    seurat_obj, perturbed_cells_by_guide, the_guide, guides, print_counts = T) {
    all_cells <- Cells(seurat_obj)
    perturbed_cells <- c()
    dummy_perturbed <- c()

    # Select all cells in which the target of the
    # guide has been successfully perturbed
    # by either the specified guide or others that target the same gene/DE.
    for (guide in guides) {
        dummy_perturbed <- unlist(perturbed_cells_by_guide[[guide]])
        perturbed_cells <- union(perturbed_cells, dummy_perturbed)
    }

    unperturbed_cells <- unlist(setdiff(all_cells, perturbed_cells))
    Idents(seurat_obj) <- "target_negative"

    # Select the cells perturbed by that guide
    perturbed_by_guide <- unlist(perturbed_cells_by_guide[[the_guide]])
    # Select the cells pertubed by other guides targeting the same gene/DE
    perturbed_by_other_guides <- c()
    dummy_perturbed <- c()

    for (guide in guides) {
        if (guide != the_guide) {
            # This is one of the other guides
            dummy_perturbed <- unlist(perturbed_cells_by_guide[[guide]])
            perturbed_by_other_guides <-
                union(perturbed_by_other_guides, dummy_perturbed)
            cat(guide, " ", length(perturbed_by_other_guides), "\n")
        }
    }

    # Now select the cell only perturbed by the specified guide
    # but not the others targeting the same gene/DE
    perturbed_only_by_specified_guide <-
        unlist(setdiff(perturbed_by_guide, perturbed_by_other_guides))
    if (length(perturbed_only_by_specified_guide) > 0) {
        seurat_obj <- SetIdent(
            seurat_obj,
            cells = perturbed_only_by_specified_guide,
            value = "guide_positive"
        )
    }

    if (print_counts) {
        n_gplus <- length(perturbed_only_by_specified_guide)
        n_gother <- length(perturbed_by_other_guides)
        n_tminus <- length(unperturbed_cells)
        n_all <- length(all_cells)
        cat(
            blue(
                "Guide+ =",
                n_gplus,
                "SisterG+ =",
                n_gother,
                "; Target- =",
                n_tminus, ";
                All =",
                n_all,
                "\n"
            )
        )
        if (n_gplus + n_gother + n_tminus != n_all) {
            cat(red("Sum not equal to cell count!!!"))
            break
        }
    }
    seurat_obj
}


# Mark cells that are positive for given guides as target_positive, all others as target_negative
mark_target_pos_neg_TEST <- function(
    seurat_obj,
    perturbed_cells_by_guide,
    guides,
    print_counts = T,
    pos_label = "target_positive",
    neg_label = "target_negative",
    guide_labels = NULL) {
    all_cells <- Cells(seurat_obj)
    perturbed_cells <- c()
    dummy_perturbed <- c()

    Idents(seurat_obj) <- neg_label
    for (i in seq_along(guides)) {
        guide <- guides[i]

        dummy_perturbed <- unlist(perturbed_cells_by_guide[[guide]])
        cat(guide, " ", length(perturbed_cells), "\n")

        if (!is.null(guide_labels)) {
            label <- guide_labels[i]
        } else {
            label <- pos_label
        }

        if (length(perturbed_cells) > 0) {
            seurat_obj <- SetIdent(
                seurat_obj,
                cells = perturbed_cells,
                value = label
            )
        }
    }
    # unperturbed_cells = unlist(setdiff(all_cells, perturbed_cells))

    if (print_counts) {
        n_gplus <- length(perturbed_cells)
        n_gminus <- length(unperturbed_cells)
        cat(blue("Guide+ =", n_gplus, "; Guide- =", n_gminus, "\n"))
    }

    seurat_obj
}


# Mark cells that are positive for given guides as target_positive, all others as target_negative
mark_target_pos_neg <- function(
    seurat_obj,
    perturbed_cells_by_guide,
    guides,
    print_counts = T,
    pos_label = "target_positive",
    neg_label = "target_negative") {
    all_cells <- Cells(seurat_obj)
    perturbed_cells <- c()
    dummy_perturbed <- c()

    for (guide in guides) {
        dummy_perturbed <- unlist(perturbed_cells_by_guide[[guide]])
        perturbed_cells <- union(perturbed_cells, dummy_perturbed)
        if (print_counts) {
            cat(guide, " ", length(perturbed_cells), "\n")
        }
    }

    unperturbed_cells <- unlist(setdiff(all_cells, perturbed_cells))
    Idents(seurat_obj) <- neg_label

    if (length(perturbed_cells) > 0) {
        seurat_obj <- SetIdent(
            seurat_obj,
            cells = perturbed_cells,
            value = pos_label
        )
    }

    if (print_counts) {
        n_gplus <- length(perturbed_cells)
        n_gminus <- length(unperturbed_cells)
        if (print_counts) {
            cat(blue("Guide+ =", n_gplus, "; Guide- =", n_gminus, "\n"))
        }
    }

    seurat_obj
}

# Mark cells that are positive for given guides as
# 'target_positive', unperturbed as 'unperturbed',
# all others as 'target_negative'. The 'unperturbed'
# are the cells with either no guides at all or only
# positive for non-targeting guides. I wanted use them as
# a gold-standard set of controls, but their number is
# pretty low.
mark_pert_unpert <- function(
    seurat_obj, perturbed_cells_by_guide, df_guides, guides, print_counts = F) {
    guide_classes <-
        unique(filter(df_guide, guide1 %in% guides | guide2 %in% guides)$class)
    if ("nt_ctrl" %in% guide_classes) {
        stop("You can't pass non-targeting guides to mark_pert_unpert()")
    }
    select_targeting_guides <- df_guides$class == "targeting"
    df_targeting <- df_guides[select_targeting_guides, ]
    targeting_guides <- c(df_targeting$guide1, df_targeting$guide2)
    seurat_dummy <-
        mark_target_pos_neg(
            seurat_obj, perturbed_cells_by_guide, targeting_guides,
            print_counts = F
        )
    unperturbed_cells <- Cells(subset(seurat_dummy, idents = "target_negative"))

    seurat_dummy <- mark_target_pos_neg(
        seurat_obj, perturbed_cells_by_guide, guides,
        print_counts = F
    )
    seurat_dummy <- SetIdent(
        seurat_dummy,
        cells = unperturbed_cells, value = "unperturbed"
    )
    if (print_counts) {
        n_unpert <- length(Cells(subset(seurat_dummy, idents = "unperturbed")))
        n_pert <- length(
            Cells(subset(seurat_dummy, idents = "target_positive"))
        )
        cat(blue("Guide+ =", n_pert, "; Unperturbed =", n_unpert, "\n"))
    }
    seurat_dummy
}

# Mark two different set of guides, A and B
mark_a_b <- function(
    seurat_obj, arg_perturbed_cells_by_guide, guides_a, guides_b) {
    all_cells <- Cells(seurat_obj)
    a_cells <- c()
    b_cells <- c()
    dummy <- c()

    Idents(seurat_obj) <- ""
    for (guide in guides_a) {
        dummy <- unlist(arg_perturbed_cells_by_guide[[guide]])
        a_cells <- union(a_cells, dummy)
    }

    dummy <- c()
    for (guide in guides_b) {
        dummy <- unlist(arg_perturbed_cells_by_guide[[guide]])
        b_cells <- union(b_cells, dummy)
    }

    # Discard overlapping cells
    overlapping_cells <- intersect(a_cells, b_cells)
    a_cells <- setdiff(a_cells, overlapping_cells)
    b_cells <- setdiff(b_cells, overlapping_cells)

    n_a <- length(a_cells)
    n_b <- length(b_cells)
    cat(blue("A =", n_a, "; B =", n_b, "\n"))

    seurat_obj <- SetIdent(seurat_obj, cells = a_cells, value = "A")
    seurat_obj <- SetIdent(seurat_obj, cells = b_cells, value = "B")
    seurat_obj
}

# Find guides that belong to a class/subclass
get_guides_by_subclass <- function(df_guides, column, value) {
    # column variable can be either 'class' or 'subclass'
    select_subclass <- df_guides[, column] == value
    df_subclass <- df_guides[select_subclass, ]
    subclass_guides <- c(df_subclass$guide1, df_subclass$guide2)
    subclass_guides
}

# Produce Vlnplot for given targets
vlnplot_for_targets <- function(
    seurat_obj, df_guide, perturbed_cells_by_guide, targets) {
    plt_list <- list()
    for (i in seq_along(targets)) {
        target <- targets[i]
        print(target)
        guides <- get_guides_by_subclass(df_guide, "alias", target)
        seurat_dummy <- mark_target_pos_neg(
            seurat_obj,
            perturbed_cells_by_guide,
            guides,
            print_counts = T
        )

        options(repr.plot.width = 5, repr.plot.height = 4)
        plt <- VlnPlot(
            object = seurat_dummy,
            features = target,
            idents = NULL,
            pt.size = 0.,
            sort = F,
            ncol = 1,
        ) + geom_boxplot(
            width = .2, color = "black", alpha = 0.2
        ) + theme(legend.position = "none")
        plt_list[[i]] <- plt
    }
    plt_list
}



# Produce Vlnplot for plasmids
vlnplot_for_plasmids <- function(
    seurat_obj, df_guides, perturbed_cells_by_guide) {
    plt_list <- list()

    for (i in 1:nrow(df_guides)) {
        target <- df_guides[i, "alias"]
        guides_on_plasmid <- unlist(
            as.list(t(df_guides[i, c("guide1", "guide2")]))
        )

        cat(blue(target, ":"), paste(guides_on_plasmid, collapse = ","), "\n")
        seurat_dummy <- mark_target_pos_neg(
            seurat_rna,
            perturbed_cells_by_guide,
            guides_on_plasmid,
            print_counts = T,
            pos_label = "plasmid_positive",
            neg_label = "plasmid_negative"
        )

        options(repr.plot.width = 5, repr.plot.height = 4)
        plt <- VlnPlot(
            object = seurat_dummy,
            features = target,
            idents = NULL,
            pt.size = 0.,
            sort = F,
            ncol = 1,
        ) +
            geom_boxplot(width = 0.2, color = "black", alpha = 0.2) +
            theme(
                legend.position = "none",
                plot.subtitle = element_text(hjust = 0.5)
            ) +
            labs(subtitle = paste0("Guides: ", guides_on_plasmid[1], ",b"))

        plt_list[[i]] <- plt
    }
    plt_list
}



# Produce RidgePlot for given targets
ridgeplot_for_targets <- function(
    seurat_obj, df_guide, perturbed_cells_by_guide, targets) {
    plt_list <- list()
    for (i in seq_along(targets)) {
        target <- targets[i]
        print(target)
        guides <- get_guides_by_subclass(df_guide, "alias", target)
        seurat_rna <- mark_target_pos_neg(
            seurat_obj,
            perturbed_cells_by_guide,
            guides,
            print_counts = T
        )

        options(repr.plot.width = 20, repr.plot.height = 4)
        plt <- RidgePlot(
            object = seurat_rna,
            features = targets,
            idents = NULL,
            sort = F,
            ncol = length(targets),
        ) + +theme(legend.position = "none")
        plt_list[[i]] <- plt
    }
    plt_list
}

# Get perturbed cells
get_perturbed_cells <- function(
    seurat_obj = NULL, df_thresholds = NULL, donor_id = NULL) {
    if (!is.null(donor_id)) {
        seurat_obj <- subset(seurat_obj, subset = donor == donor_id)
    }
    libraries <- unique(seurat_obj$library)
    seurat_libs <- list()

    for (i in seq_along(libraries)) {
        lib <- libraries[i]
        seurat_libs[[i]] <- subset(seurat_obj, subset = library == lib)
    }
    names(seurat_libs) <- libraries
    perturbed_cells_by_guide <- list()

    for (i in 1:nrow(df_thresholds)) {
        perturbed_cells_in_all_libs <- list()
        guide <- df_thresholds$guide[i]
        # Loop over libraries
        for (lib in libraries) {
            seurat_lib <- seurat_libs[[lib]]
            threshold <- df_thresholds[i, lib]
            cells_in_lib <- Cells(seurat_lib)
            sgrna_counts <- seurat_lib[["sgRNA"]]@counts
            select_perturbed <- sgrna_counts[guide, cells_in_lib] >= threshold
            perturbed_cells_in_library <- cells_in_lib[select_perturbed]
            if (!is.na(threshold)) {
                perturbed_cells_in_all_libs <-
                    append(
                        perturbed_cells_in_all_libs,
                        perturbed_cells_in_library
                    )
            }
        }
        perturbed_cells_by_guide[[i]] <- perturbed_cells_in_all_libs
    }
    names(perturbed_cells_by_guide) <- df_thresholds$guide
    perturbed_cells_by_guide
}


# This function is somewhat redundant. It runs on a list of Seurat objects
get_all_perturbed_cells_by_guide <- function(
    seurat_obj_libs = NULL, df_thresholds = NULL) {
    perturbed_cells_by_guide <- list()

    for (i in 1:nrow(df_thresholds)) {
        perturbed_cells_in_all_libs <- list()
        guide <- df_thresholds$guide[i]
        libraries <- names(seurat_obj_libs)
        # Loop over libraries
        for (lib in libraries) {
            seurat_lib <- seurat_obj_libs[[lib]]
            threshold <- df_thresholds[i, lib]
            cells_in_lib <- Cells(seurat_lib)
            sgrna_counts <- seurat_lib[["sgRNA"]]@counts
            select_perturbed <- sgrna_counts[guide, cells_in_lib] >= threshold
            perturbed_cells_in_library <- cells_in_lib[select_perturbed]
            if (!is.na(threshold)) {
                perturbed_cells_in_all_libs <-
                    append(
                        perturbed_cells_in_all_libs,
                        perturbed_cells_in_library
                    )
            }
        }
        perturbed_cells_by_guide[[i]] <- perturbed_cells_in_all_libs
    }
    names(perturbed_cells_by_guide) <- df_thresholds$guide
    perturbed_cells_by_guide
}




get_neighboring_genes <- function(
    bm, df_target_coords,
    genes_in_assay,
    target_upstream_range = 1000000,
    target_downstream_range = 1000000) {
    attributes <- c(
        "ensembl_gene_id", "chromosome_name", "start_position", "end_position",
        "strand", "hgnc_symbol", "refseq_ncrna", "refseq_ncrna_predicted"
    )
    filters <- c("chromosome_name", "start", "end")
    target_neigbors_list <- list()

    for (i in 1:nrow(df_target_coords)) {
        target <- df_target_coords[i, ]
        chr <- target$chromosome
        start <- target$start_position - target_upstream_range
        end <- target$end_position + target_downstream_range
        cat(blue(target$hgnc_symbol, chr, start, end, "\n"))

        values <- list(chromosome = chr, start = start, end = end)
        neighbor_genes <- getBM(
            attributes = attributes,
            filters = filters,
            values = values,
            mart = bm
        )
        n_neighbors <- length(unique(neighbor_genes$ensembl_gene_id))
        select_non_na <- !is.na(
            neighbor_genes$hgnc_symbol
        ) & (neighbor_genes$hgnc_symbol != "")
        neighbor_genes <- neighbor_genes[select_non_na, ]
        n_neighbors_w_names <- length(unique(neighbor_genes$hgnc_symbol))

        if (n_neighbors_w_names > 0) {
            dummy <- as.vector(unique(neighbor_genes$hgnc_symbol))
            # some neighbors are not in the cellranger reference. Remove those.
            select_genes_in_assay <- dummy %in% genes_in_assay
            target_neigbors_list[[target$hgnc_symbol]] <- dummy[select_genes_in_assay]
        } else {
            print("")
        }
    }
    target_neigbors_list
}


# Split perturbed (or control) cells into random replicate
# groups for diffex testing
create_random_replicates <- function(
    seurat_obj, perturbed_id_list, control_id_list) {
    perturbed_cells <- WhichCells(seurat_obj, idents = "guide_positive")
    control_cells <- WhichCells(seurat_obj, idents = "guide_negative")
    n_pos <- length(perturbed_cells)
    n_neg <- length(control_cells)
    perturbed_labels <- sample(perturbed_id_list, n_pos, replace = TRUE)
    control_labels <- sample(control_id_list, n_neg, replace = TRUE)
    seurat_obj <- SetIdent(
        seurat_obj,
        cells = perturbed_cells,
        value = perturbed_labels
    )
    seurat_obj <- SetIdent(
        seurat_obj,
        cells = control_cells,
        value = control_labels
    )

    # double check replicate labels
    if (
        any(
            WhichCells(
                seurat_obj,
                idents = perturbed_id_list
            ) != perturbed_cells
        ) |
            any(
                WhichCells(
                    seurat_obj,
                    idents = control_id_list
                ) != control_cells
            )
    ) {
        cat(red("Failed assignment!!! "))
    } else {
        cat(green("Random replicate creation OK"))
    }
    seurat_obj
}


