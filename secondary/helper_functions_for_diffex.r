# Helper functions

# Mark cells that are positive ONLY for a given guide (not construct), and negative for a target
# The purpose is to run diffex analysis for each guide seperately.
mark_guide_pos_target_neg <- function(seurat_obj, perturbed_cells_by_guide, the_guide, guides, print_counts = T) {
    all_cells = Cells(seurat_obj)
    perturbed_cells = c()
    dummy_perturbed = c()

    # Select all cells in which the target of the guide has been successfully perturbed    
    # by either the specified guide or others that target the same gene/DE.
    for(guide in guides) {
        dummy_perturbed = unlist(perturbed_cells_by_guide[[guide]])
        perturbed_cells = union(perturbed_cells, dummy_perturbed)    
        #cat(guide," ", length(perturbed_cells), "\n")
    }

    unperturbed_cells = unlist(setdiff(all_cells, perturbed_cells))
    Idents(seurat_obj) <- "target_negative"

    # Select the cells perturbed by that guide
    perturbed_by_guide = unlist(perturbed_cells_by_guide[[the_guide]])
    # Select the cells pertubed by other guides targeting the same gene/DE
    perturbed_by_other_guides = c()
    dummy_perturbed = c()    
    
    for(guide in guides) {
        if(guide != the_guide) {
            # This is one of the other guides
            dummy_perturbed = unlist(perturbed_cells_by_guide[[guide]])
            perturbed_by_other_guides = union(perturbed_by_other_guides, dummy_perturbed)
            cat(guide," ", length(perturbed_by_other_guides), "\n")
        }
    }
    
    # Now select the cell only perturbed by the specified guide
    # but not the others targeting the same gene/DE
    perturbed_only_by_specified_guide = unlist(setdiff(perturbed_by_guide, perturbed_by_other_guides))
    if(length(perturbed_only_by_specified_guide) > 0) {
        seurat_obj <- SetIdent(
            seurat_obj, 
            cells = perturbed_only_by_specified_guide, 
            value = "guide_positive"
        )
    }

    if(print_counts) {
        n_gplus   = length(perturbed_only_by_specified_guide)
        n_gother  = length(perturbed_by_other_guides)    
        n_tminus  = length(unperturbed_cells)    
        n_all     = length(all_cells)
        cat(blue("Guide+ =",n_gplus, "SisterG+ =",n_gother, "; Target- =", n_tminus, "; All =", n_all, "\n"))
        if(n_gplus + n_gother + n_tminus != n_all) {
            cat(red("Sum not equal to cell count!!!"))
            break
        }
    }
    seurat_obj
}


# Mark cells that are positive for given guides as target_positive, all others as target_negative
mark_target_pos_neg <- function(seurat_obj, perturbed_cells_by_guide, guides, print_counts = T) {
    all_cells = Cells(seurat_obj)
    perturbed_cells = c()
    dummy_perturbed = c()
    
    for(guide in guides) {
        dummy_perturbed = unlist(perturbed_cells_by_guide[[guide]])
        perturbed_cells = union(perturbed_cells, dummy_perturbed)    
        #cat(guide," ", length(perturbed_cells), "\n")
    }

    unperturbed_cells = unlist(setdiff(all_cells, perturbed_cells))
    Idents(seurat_obj) <- "target_negative"
    seurat_obj <- SetIdent(seurat_obj, cells = perturbed_cells, value = "target_positive") 

    if(print_counts) {
        n_gplus   = length(perturbed_cells)
        n_gminus = length(unperturbed_cells)    
        cat(blue("Guide+ =",n_gplus, "; Guide- =", n_gminus, "\n"))
    }
    
    seurat_obj
}

# Mark cells that are positive for given guides as 'target_positive', unperturbed as 
# 'unperturbed', all others as 'target_negative'. The 'unperturbed' are the cells with 
# either no guides at all or only positive for non-targeting guides. I wanted use them as 
# a gold-standard set of controls, but their number is pretty low. 
mark_pert_unpert <- function(seurat_obj, perturbed_cells_by_guide, df_guides, guides, print_counts = F) {
    guide_classes = unique(filter(df_guide, guide1 %in% guides | guide2 %in% guides)$class)
    if('nt_ctrl' %in% guide_classes) {
        stop("You can't pass non-targeting guides to mark_pert_unpert()")
    }
    select_targeting_guides = df_guides$class == 'targeting'
    df_targeting  = df_guides[select_targeting_guides, ]
    targeting_guides = c(df_targeting$guide1, df_targeting$guide2)
    seurat_dummy <- mark_target_pos_neg(seurat_obj, perturbed_cells_by_guide, targeting_guides, print_counts = F)
    unperturbed_cells <- Cells(subset(seurat_dummy, idents='target_negative'))
    
    seurat_dummy <- mark_target_pos_neg(seurat_obj, perturbed_cells_by_guide, guides, print_counts = F)
    seurat_dummy <- SetIdent(seurat_dummy, cells = unperturbed_cells, value = "unperturbed")
    if(print_counts) {
        n_unpert = length(Cells(subset(seurat_dummy, idents='unperturbed')))
        n_pert = length(Cells(subset(seurat_dummy, idents='target_positive')))
        cat(blue("Guide+ =",n_pert, "; Unperturbed =", n_unpert, "\n"))
    }
    seurat_dummy
}

# Mark two different set of guides, A and B
mark_a_b <- function(seurat_obj, arg_perturbed_cells_by_guide, guides_a, guides_b) {
    all_cells = Cells(seurat_obj)
    a_cells = c()
    b_cells = c()
    dummy = c()
    
    Idents(seurat_obj) <- ""
    for(guide in guides_a) {
        dummy = unlist(arg_perturbed_cells_by_guide[[guide]])
        a_cells = union(a_cells, dummy)
    }

    dummy = c() 
    for(guide in guides_b) {
        dummy = unlist(arg_perturbed_cells_by_guide[[guide]])
        b_cells = union(b_cells, dummy)
    }

    # Discard overlapping cells
    overlapping_cells = intersect(a_cells, b_cells)
    a_cells = setdiff(a_cells, overlapping_cells)
    b_cells = setdiff(b_cells, overlapping_cells)
    
    n_a = length(a_cells)
    n_b = length(b_cells)
    cat(blue("A =", n_a, "; B =", n_b, "\n"))
    
    seurat_obj <- SetIdent(seurat_obj, cells = a_cells, value = "A") 
    seurat_obj <- SetIdent(seurat_obj, cells = b_cells, value = "B") 
    seurat_obj
}

# Find guides that belong to a class/subclass
get_guides_by_subclass <- function(df_guides, column, value) {
    # column variable can be either 'class' or 'subclass'
    select_subclass = df_guides[, column] == value
    df_subclass     = df_guides[select_subclass, ]
    subclass_guides = c(df_subclass$guide1, df_subclass$guide2)
    subclass_guides    
}

# Produce Vlnplot for given targets
vlnplot_for_targets <- function(seurat_obj, df_guide, perturbed_cells_by_guide, targets) {
    plt_list = list()
    for(i in seq_along(targets)) {
        target = targets[i]
        print(target)
        guides = get_guides_by_subclass(df_guide, 'alias', target)
        seurat_dummy <- mark_target_pos_neg(
            seurat_obj, 
            perturbed_cells_by_guide, 
            guides, 
            print_counts = T
        )

        options(repr.plot.width=5, repr.plot.height=4)
        plt <- VlnPlot(
            object = seurat_dummy,
            features =  target, 
            idents = NULL, 
            pt.size = 0., 
            sort = F, 
            ncol = 1,    
        ) + geom_boxplot(width=1, color="black", alpha=0.2) + theme(legend.position = 'none')  #+ stat_summary(fun = "mean", colour = "blue")
        plt_list[[i]] = plt
    }
    plt_list
}

# Produce RidgePlot for given targets
ridgeplot_for_targets <- function(seurat_obj, df_guide, perturbed_cells_by_guide, targets) {
    plt_list = list()
    for(i in seq_along(targets)) {
        target = targets[i]
        print(target)
        guides = get_guides_by_subclass(df_guide, 'alias', target)
        seurat_rna <- mark_target_pos_neg(
            seurat_obj, 
            perturbed_cells_by_guide, 
            guides, 
            print_counts = T
        )

        options(repr.plot.width=20, repr.plot.height=4)
        plt <- RidgePlot(
            object = seurat_rna,
            features =  targets, 
            idents = NULL,             
            sort = F, 
            ncol = length(targets),    
        ) + + theme(legend.position = 'none')
        plt_list[[i]] = plt
    }
    plt_list
}

# Get perturbed cells
get_perturbed_cells <- function(seurat_obj, donor_id) {
    seurat_obj = subset(seurat_obj, subset = donor == donor_id)
    libraries = unique(seurat_obj$library)
    seurat_libs = list()

    for(i in seq_along(libraries)){ 
        lib = libraries[i]
        seurat_libs[[i]] = subset(seurat_obj, subset = library == lib)
    }
    names(seurat_libs) <- libraries
    perturbed_cells_by_guide = list()

    for(i in 1:nrow(df_thresholds)){  
        perturbed_cells_in_all_libs = list()
        guide = df_thresholds$guide[i]
        # Loop over libraries
        for(lib in libraries){        
            seurat_lib = seurat_libs[[lib]]
            threshold = df_thresholds[i, lib]        
            #cat(blue(guide, lib, threshold, "\n"))
            cells_in_lib = Cells(seurat_lib)        
            sgrna_counts = seurat_lib[['sgRNA']]@counts
            select_perturbed = sgrna_counts[guide, cells_in_lib] >= threshold
            perturbed_cells_in_library = cells_in_lib[select_perturbed]
            #cat(length(cells_in_lib), "in", lib, guide, length(perturbed_cells_in_library), "cells >", threshold, "\n")        
            perturbed_cells_in_all_libs = append(perturbed_cells_in_all_libs, perturbed_cells_in_library)
        }
        perturbed_cells_by_guide[[i]] = perturbed_cells_in_all_libs
    }
    names(perturbed_cells_by_guide) <- df_thresholds$guide
    perturbed_cells_by_guide
}