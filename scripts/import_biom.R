# Script to handle single sample and single OTU biom (`https://github.com/joey711/phyloseq/issues/1532`)

# Adapted from biomformat
new_biom_data <- function(x, rows, columns, parallel=FALSE) {
  if(rlang::is_missing(rows)) {
    rows = 1:nrow(x)
  }
  if(rlang::is_missing(columns)) {
    columns = 1:ncol(x)
  }
  if( identical(length(rows), 0) ){
    stop("argument `rows` must have non-zero length.")
  }
  if( identical(length(columns), 0) ){
    stop("argument `columns` must have non-zero length.")
  }
  # Begin matrix section
  if( identical(x$matrix_type, "dense") ){
    # Begin dense section
    # If matrix is stored as dense, create "vanilla" R matrix, m
    m = plyr::laply(x$data[rows], function(i) i[columns], .parallel=parallel)
    if( length(rows) > 1L &
                       length(columns) > 1L &
                                       matrix_element_type(x) %in% c("int", "float")
                                             ){
      # If either dimension is length-one, don't call coerce to "Matrix"
      # Note that laply() does still work in this case.
      # If both dimension lengths > 1 & data is numeric,
      # attempt to coerce to Matrix-inherited class,
      # Mainly because it might still be sparse and this is a good way
      # to handle it in R.
      m = Matrix::Matrix(m)
    }
  } else {
    # Begin sparse section
    ## Initialize sparse matrix as either Matrix or matrix, depending on data class
    biom_shape = biom_shape(x)
    if(matrix_element_type(x) %in% c("int", "float")){
      # If data is numeric, initialize with Matrix (numeric data only)
      m = Matrix::Matrix(0, nrow=nrow(x), ncol=ncol(x), sparse=TRUE)
      # Create an assignment data.frame
      adf = plyr::ldply(x$data)
    } else {
      # Else, matrix_element_type must be "unicode" for a unicode string.
      # Use a standard R character matrix
      m = matrix(NA_character_, nrow(x), ncol(x))
      # Create an assignment data.frame.
      # Is slightly more complicated for sparse JSON w/ character values
      adf = plyr::ldply(x$data, function(x){
                data.frame(r=x[[1]], c=x[[2]], data=x[[3]], stringsAsFactors=FALSE)
            })
    }
    colnames(adf) <- c("r", "c", "data")
    # indices start at 0 in biom sparse format,
    # and are first two columns
    adf[, c("r", "c")] <- cbind(r = as.integer(adf$r) + 1L,
                                                                c = as.integer(adf$c) + 1L)
    # Subset to just indices that are in both arguments `rows` and `columns`
    adf <- adf[(adf$r %in% rows & adf$c %in% columns), ]
    # Fill in data values in matrix, m.
    # Vectorized for speed using matrix indexing.
    # See help("Extract") for details about matrix indexing. Diff than 2-vec index.
    m[as(adf[, 1:2], "matrix")] <- adf$data
    # Subset this biggest-size m to just `rows` and `columns`
    m = m[rows, columns]
  # End sparse section
  }
  # If either dimension is length-one
  if( identical(length(rows), 1L) | identical(length(columns), 1L) ){
    if( identical(length(rows), 1L) ){
      # If row is length-one (1 taxon)
      m = t(as.matrix(m))
    } else {
      # If column is length-one (1 sample)
      m = as.matrix(m)
    }
  }
  # Add row and column names
  rownames(m) <- sapply(x$rows[rows], function(i) i$id )
  colnames(m) <- sapply(x$columns[columns], function(i) i$id )
  return(m)
}

# Adapted from phyloseq
new_import_biom <- function(BIOMfilename,
        treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
        parseFunction=parse_taxonomy_default, parallel=FALSE, version=1.0, ...){

        # initialize the argument-list for phyloseq. Start empty.
        argumentlist <- list()

        # Read the data
        if(class(BIOMfilename)=="character"){
                x = read_biom(biom_file=BIOMfilename)
        } else if (class(BIOMfilename)=="biom"){
                x = BIOMfilename
        } else {
                stop("import_biom requires a 'character' string to a biom file or a 'biom-class' object")
        }

        ########################################
        # OTU table:
        ########################################
        otutab = otu_table(as(new_biom_data(x), "matrix"), taxa_are_rows=TRUE)
        argumentlist <- c(argumentlist, list(otutab))

        ########################################
        # Taxonomy Table
        ########################################
        # Need to check if taxonomy information is empty (minimal BIOM file)
        if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
                taxtab <- NULL
        } else {
                # parse once each character vector, save as a list
                taxlist = lapply(x$rows, function(i){
                        parseFunction(i$metadata$taxonomy)
                })
                names(taxlist) = sapply(x$rows, function(i){i$id})
                        taxtab = build_tax_table(taxlist)
        }
        argumentlist <- c(argumentlist, list(taxtab))

        ########################################
        # Sample Data ("columns" in QIIME/BIOM)
        ########################################
        # If there is no metadata (all NULL), then set sam_data <- NULL
        if( is.null(sample_metadata(x)) ){
                samdata <- NULL
        } else {
                samdata = sample_data(sample_metadata(x))
        }
        argumentlist <- c(argumentlist, list(samdata))

        ########################################
        # Tree data
        ########################################
        if( !is.null(treefilename) ){
                if( inherits(treefilename, "phylo") ){
                        # If argument is already a tree, don't read, just assign.
                        tree = treefilename
                } else {
                        # NULL is silently returned if tree is not read properly.
                        tree <- read_tree(treefilename, ...)
                }
                # Add to argument list or warn
                if( is.null(tree) ){
                        warning("treefilename failed import. It not included.")
                } else {
                        argumentlist <- c(argumentlist, list(tree) )
                }
        }

        ########################################
        # Reference Sequence data
        ########################################        
        if( !is.null(refseqfilename) ){
                if( inherits(refseqfilename, "XStringSet") ){
                        # If argument is already a XStringSet, don't read, just assign.
                        refseq = refseqfilename
                } else {
                        # call refseqFunction and read refseqfilename, either with or without additional args
                        if( !is.null(refseqArgs) ){
                                refseq = do.call("refseqFunction", c(list(refseqfilename), refseqArgs))
                        } else {
                                refseq = refseqFunction(refseqfilename)
                        }
                }
                argumentlist <- c(argumentlist, list(refseq) )
        }

        ########################################
        # Put together into a phyloseq object
        ########################################
        return( do.call("phyloseq", argumentlist) )

}
